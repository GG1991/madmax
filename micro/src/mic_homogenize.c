/*
   Routines for performing homogenization on RVE

   Author > Guido Giuntoli
   Date   > 18-08-2017

 */

#include "micro.h"

int mic_homogenize_taylor(MPI_Comm PROBLEM_COMM, double strain_mac[6], double strain_ave[6], double stress_ave[6])
{
  int      i, k, e, gp, ngp, npe, ierr;
  double   coor_elm[8][3], dsh[8][3], detj;
  double   *wp = NULL, stress_gp[6], DsDe[6][6];
  double   vol = -1.0, vol_tot = -1.0, stress_aux[6], strain_aux[6];
  double   wp_eff;

  for(i=0;i<nvoi;i++){
    strain_aux[i] = stress_aux[i] = 0.0;
  }
  vol = 0.0;
  for(e=0;e<nelm;e++){

    npe = eptr[e+1]-eptr[e];
    ngp = npe;
    GetElemCoord(&eind[eptr[e]], npe, coor_elm);
    GetWeight(npe, &wp);

    for(gp=0;gp<ngp;gp++){
      memset(stress_gp, 0.0, nvoi*sizeof(double));
      get_dsh(gp, npe, coor_elm, dsh, &detj);
      detj = fabs(detj);
      get_c(e, gp, strain_mac, DsDe);
      wp_eff = detj*wp[gp];
      for(i=0;i<nvoi;i++){
	for(k=0;k<nvoi;k++){
	  stress_gp[i] += DsDe[i][k]*strain_mac[k];
	}
      }
      for(i=0;i<nvoi;i++){
	stress_aux[i] += stress_gp[i]  * wp_eff;
	strain_aux[i] += strain_mac[i] * wp_eff;
      }
      vol += wp_eff;
    }
  }
  ierr = MPI_Allreduce(stress_aux, stress_ave, nvoi, MPI_DOUBLE, MPI_SUM, PROBLEM_COMM);CHKERRQ(ierr);
  ierr = MPI_Allreduce(strain_aux, strain_ave, nvoi, MPI_DOUBLE, MPI_SUM, PROBLEM_COMM);CHKERRQ(ierr);
  ierr = MPI_Allreduce(&vol, &vol_tot, 1, MPI_DOUBLE, MPI_SUM, PROBLEM_COMM);CHKERRQ(ierr);
  for(i=0;i<nvoi;i++){
    stress_ave[i] /= vol_tot;
    strain_ave[i] /= vol_tot;
  }

  return 0;
}
/****************************************************************************************************/
int mic_homogenize_unif_strains(MPI_Comm MICRO_COMM, double strain_mac[6], double strain_ave[6], double stress_ave[6])
{
  /*
     ub = E . x
  */

  int    ierr, nr_its=-1;
  int    nnods_bc, *nods_bc;
  int    i, d, k;
  int    *ix_loc, *ix_glo;
  double *ub;
  double norm=2*nr_norm_tol;
  double strain_matrix[3][3];
  double *x_arr, *b_arr;

  if(dim==2){
    strain_matrix[0][0]=strain_mac[0]; strain_matrix[0][1]=strain_mac[2];
    strain_matrix[1][0]=strain_mac[2]; strain_matrix[1][1]=strain_mac[1];
  }
  else if(dim==3){
    strain_matrix[0][0]=strain_mac[0]; strain_matrix[0][1]=strain_mac[3]; strain_matrix[0][2]=strain_mac[5];
    strain_matrix[1][0]=strain_mac[3]; strain_matrix[1][1]=strain_mac[1]; strain_matrix[1][2]=strain_mac[4];
    strain_matrix[2][0]=strain_mac[5]; strain_matrix[2][1]=strain_mac[4]; strain_matrix[2][2]=strain_mac[2];
  }

  /* get array of nods in al boundary without repetition */
  ierr = get_nods_bc( &nods_bc, &nnods_bc);
  if(ierr){
    return 1;
  }
  ix_loc = malloc(nnods_bc*dim*sizeof(int));    /* indeces to search for local coordinates */
  ix_glo = malloc(nnods_bc*dim*sizeof(int));    /* indeces to be set at boundary */
  ub     = malloc(nnods_bc*dim*sizeof(double)); /* values  to be set at boundary */
  ierr   = get_nods_index( nods_bc, nnods_bc, ix_loc, ix_glo);

  // first we fill the value ub_val in <homog_ld_lagran_t> structure
  for(i=0; i<nnods_bc; i++){
    for(d=0;d<dim;d++){
      ub[i*dim+d] = 0.0;
      for(k=0;k<dim;k++){
	ub[ i * dim + d] = ub[ i * dim + d] +  strain_matrix[d][k] * coord[ix_loc[i* dim + k]];
      }
    }
  }

  // first we fill the value ub_val in <homog_ld_lagran_t> structure
  ierr = VecGetArray(x, &x_arr);CHKERRQ(ierr);
  for(i=0; i<nnods_bc; i++){
    for(d=0;d<dim;d++){
      x_arr[ix_loc[i*dim+d]] = ub[i*dim+d]; 
    }
  }
  ierr = VecRestoreArray(x, &x_arr);CHKERRQ(ierr);

  while( nr_its < nr_max_its && norm > nr_norm_tol )
  {
    ierr = assembly_residual_sd(&x, &b);CHKERRQ(ierr);
    ierr = VecGetArray(b, &b_arr);CHKERRQ(ierr);
    for(i=0; i<nnods_bc; i++){
      for(d=0;d<dim;d++){
	b_arr[ix_loc[i*dim+d]] = 0.0;
      }
    }
    ierr = VecRestoreArray(b, &b_arr);CHKERRQ(ierr);
    ierr = VecScale(b,-1.0); CHKERRQ(ierr);

    ierr = VecNorm(b,NORM_2,&norm);CHKERRQ(ierr);
    PetscPrintf(MICRO_COMM,"|b| = %lf \n",norm);

    if( !(norm > nr_norm_tol) )break;

    ierr = assembly_jacobian_sd(&A);
    ierr = MatZeroRowsColumns(A, nnods_bc*dim, ix_glo, 1.0, NULL, NULL); CHKERRQ(ierr);
    ierr = KSPSolve(ksp,b,dx);CHKERRQ(ierr);

    ierr = VecAXPY( x, 1.0, dx); CHKERRQ(ierr);
    ierr = MatMult(A, dx, b1);CHKERRQ(ierr);
    
    if( flag_print & (1<<PRINT_PETSC) ){
      ierr = PetscViewerASCIIOpen(MICRO_COMM,"A.dat",&viewer); CHKERRQ(ierr);
      ierr = MatView(A,viewer); CHKERRQ(ierr);
      ierr = PetscViewerASCIIOpen(MICRO_COMM,"b.dat",&viewer); CHKERRQ(ierr);
      ierr = VecView(b,viewer); CHKERRQ(ierr);
      ierr = PetscViewerASCIIOpen(MICRO_COMM,"b1.dat",&viewer); CHKERRQ(ierr);
      ierr = VecView(b1,viewer); CHKERRQ(ierr);
      ierr = PetscViewerASCIIOpen(MICRO_COMM,"x.dat",&viewer); CHKERRQ(ierr);
      ierr = VecView(x,viewer); CHKERRQ(ierr);
      ierr = PetscViewerASCIIOpen(MICRO_COMM,"dx.dat",&viewer); CHKERRQ(ierr);
      ierr = VecView(x,viewer); CHKERRQ(ierr);
    }

    nr_its ++;
  }

  ierr = calc_ave_strain_stress(MICRO_COMM, &x, strain_ave, stress_ave);CHKERRQ(ierr);

  return 0;
}
/****************************************************************************************************/
int mic_homogenize_ld_lagran(MPI_Comm MICRO_COMM, double strain_mac[6], double strain_ave[6], double stress_ave[6])
{
  /*
     Here we apply the Lagrage multiplier procedure 
     to stablish the linear displacement boundary conditions. 
     The equations can be written:
     fa(u)=0
     fb(u)-ub=0
     ub-D^Te=0

         | A    ei |       | f       |     | x_val |
     J = |         |  b =  |         | x = |       |
         | ei^T 0  |       | ub-D^Te |     | l_val |

     the "f" contains fa and fb mixed depending on the node numeration
  */
//  int    ierr, nr_its=-1, nnods_bc, i, d, k, n_loc, *ixb;
//  double nr_norm_tol=1.0e-8, norm=2*nr_norm_tol, ub, strain_matrix[3][3], *x_arr, *b_arr;
//
//  if(dim==2){
//    strain_matrix[0][0]=strain_mac[0]; strain_matrix[0][1]=strain_mac[2];
//    strain_matrix[1][0]=strain_mac[2]; strain_matrix[1][1]=strain_mac[1];
//  }
//  else if(dim==3){
//    strain_matrix[0][0]=strain_mac[0]; strain_matrix[0][1]=strain_mac[3]; strain_matrix[0][2]=strain_mac[5];
//    strain_matrix[1][0]=strain_mac[3]; strain_matrix[1][1]=strain_mac[1]; strain_matrix[1][2]=strain_mac[4];
//    strain_matrix[2][0]=strain_mac[5]; strain_matrix[2][1]=strain_mac[4]; strain_matrix[2][2]=strain_mac[2];
//  }
//
//  nnods_bc = ((homog_ld_lagran_t*)homo.st)->nnods_bc;
//  ixb = ((homog_ld_lagran_t*)homo.st)->index;
//
//  // first we fill the value ub_val in <homog_ld_lagran_t> structure
//  for(i=0; i<nnods_bc; i++){
//    for(d=0;d<dim;d++){
//      n_loc = ixb[i*dim+d]/dim;
//      ub = 0.0;
//      for(k=0;k<dim;k++){
//	ub += strain_matrix[d][k]*coord[n_loc * dim + k];
//      }
//      ((homog_ld_lagran_t*)homo.st)->ub_val[i*dim+d] = ub;
//    }
//  }
//
//  while( nr_its < nr_max_its && norm > nr_norm_tol )
//  {
//
//    /* internal forces */
//    ierr = VecGetArray(x, &x_arr);CHKERRQ(ierr);
//    ierr = assembly_residual_sd(&x, &b);CHKERRQ(ierr);
//    /* fb - delta */
//    ierr = VecGetArray(b, &b_arr);CHKERRQ(ierr);
//    for(i=0; i<nnods_bc; i++){
//      for(d=0;d<dim;d++){
//	b_arr[ixb[i*dim+d]] -= x_arr[nmynods*dim+i*dim+d];
//      }
//    }
//    /* ub - D^T . e_mac */
//    for(i=0; i<nnods_bc; i++){
//      for(d=0;d<dim;d++){
//	b_arr[nmynods*dim + i*dim+d] = x_arr[ixb[i*dim+d]] - ((homog_ld_lagran_t*)homo.st)->ub_val[i*dim+d];
//      }
//    }
//    ierr = VecRestoreArray(b, &b_arr);CHKERRQ(ierr);
//    ierr = VecRestoreArray(x, &x_arr);CHKERRQ(ierr);
//
//    ierr = VecNorm(b,NORM_2,&norm);CHKERRQ(ierr);
//    PetscPrintf(MICRO_COMM,"|b| = %lf \n",norm);
//    if( !(norm > nr_norm_tol) ) break;
//    ierr = VecScale(b,-1.0); CHKERRQ(ierr);
//
//    /* Tangent matrix 
//
//         | K    1  |
//     A = |         |
//         | 1    0  |
//    */
//    ierr = assembly_jacobian_sd(&A);
//    for(i=0; i<nnods_bc; i++){
//      for(d=0;d<dim;d++){
//	MatSetValue(A,ixb[i*dim+d],nmynods*dim+i*dim+d,-1.0,INSERT_VALUES);
//      }
//    }
//    for(i=0; i<nnods_bc; i++){
//      for(d=0;d<dim;d++){
//	MatSetValue(A,nmynods*dim+i*dim+d,ixb[i*dim+d],1.0,INSERT_VALUES);
//      }
//    }
//    for(i=0; i<nnods_bc; i++){
//      for(d=0;d<dim;d++){
//	MatSetValue(A,nmynods*dim+i*dim+d,nmynods*dim+i*dim+d,0.0,INSERT_VALUES);
//      }
//    }
//    MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);
//    MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);
//  
//    ierr = KSPSolve(ksp, b, dx);CHKERRQ(ierr);
//    ierr = VecAXPY(x, 1.0, dx); CHKERRQ(ierr);
//
//    ierr = MatMult(A, dx, b1);CHKERRQ(ierr);
//    
//    if( flag_print & (1<<PRINT_PETSC) ){
//      ierr = PetscViewerASCIIOpen(MICRO_COMM,"A.dat",&viewer); CHKERRQ(ierr);
//      ierr = MatView(A,viewer); CHKERRQ(ierr);
//      ierr = PetscViewerASCIIOpen(MICRO_COMM,"b.dat",&viewer); CHKERRQ(ierr);
//      ierr = VecView(b,viewer); CHKERRQ(ierr);
//      ierr = PetscViewerASCIIOpen(MICRO_COMM,"b1.dat",&viewer); CHKERRQ(ierr);
//      ierr = VecView(b1,viewer); CHKERRQ(ierr);
//      ierr = PetscViewerASCIIOpen(MICRO_COMM,"x.dat",&viewer); CHKERRQ(ierr);
//      ierr = VecView(x,viewer); CHKERRQ(ierr);
//      ierr = PetscViewerASCIIOpen(MICRO_COMM,"dx.dat",&viewer); CHKERRQ(ierr);
//      ierr = VecView(x,viewer); CHKERRQ(ierr);
//    }
//
//    nr_its ++;
//  }
//
//  ierr = calc_ave_strain_stress(MICRO_COMM, &x, strain_ave, stress_ave);CHKERRQ(ierr);
//
  return 0;
}
/****************************************************************************************************/
int mic_homogenize(MPI_Comm MICRO_COMM, double strain_mac[6], double strain_ave[6], double stress_ave[6])
{

  /*
     Performs linear homogenization of the RVE 

     UNIF_STRAINS            > u = E . x

     HOMO_TAYLOR             > 
   */
  int ierr;

  switch(homo_type){

    case TAYLOR:
      ierr = mic_homogenize_taylor(MICRO_COMM, strain_mac, strain_ave, stress_ave);
      if(ierr){
	return 1;
      }
      break;

    case UNIF_STRAINS:
      ierr = mic_homogenize_unif_strains(MICRO_COMM, strain_mac, strain_ave, stress_ave);
      if(ierr){
	return 1;
      }
      break;

    default:
      return 1;
  }
  

  return 0;
}
/****************************************************************************************************/
int mic_calc_c_homo(MPI_Comm MICRO_COMM, double strain_mac[6], double *c_homo)
{

  /* 
  Si la micro estructura está integramente conformada por materiales
  lineales entonces este tensor será siempre el mismo para cada punto 
  de gauss en la macro escala entonces es eficiente almacenar c_homo_linear
  */

  int ierr;

  if(flag_linear_micro){

    if(first_time_c_homo_lineal_ask){
      ierr = mic_calc_c_homo_lineal(MICRO_COMM, strain_mac, c_homo_lineal);
      if(ierr){
	return 1;
      }
      first_time_c_homo_lineal_ask = 0;
    }
    c_homo = c_homo_lineal;

  }
  else{
    return 1;
  }
  return 0;
}
/****************************************************************************************************/
int mic_calc_c_homo_lineal(MPI_Comm MICRO_COMM, double strain_mac[6], double c_homo_lineal[36])
{

  int      i, j;
  int      ierr;
  double   strain[6], strain_ave[6], stress_ave[6];

  for(i=0;i<nvoi*nvoi;i++){
    c_homo_lineal[i]=0.0;
  }

  for(i=0;i<nvoi;i++){

    for(i=0;i<nvoi;i++){
      strain[i]=0.0;
    }
    strain[i]=0.005;

    ierr = mic_homogenize(MICRO_COMM, strain, strain_ave, stress_ave);
    if(ierr){
      return 1;
    }

    for(j=0;j<nvoi;j++){
      c_homo_lineal[j*nvoi+i] = stress_ave[j] / strain_ave[i];
    }

  }

  return 0;
}

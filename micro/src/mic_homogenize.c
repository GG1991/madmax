/*
   Routines for performing homogenization on RVE

   Author > Guido Giuntoli
   Date   > 18-08-2017

 */

#include "micro.h"

int micro_homogenize_taylor(MPI_Comm PROBLEM_COMM, double strain_mac[6], double strain_ave[6], double stress_ave[6])
{
  int    i, k, e, gp, ngp, npe, ierr;
  double coor_elm[8][3], dsh[8][3], detj, *wp = NULL, stress_gp[6], DsDe[6][6], disp_elm[8*3];
  double vol = -1.0, vol_tot = -1.0, stress_aux[6], strain_aux[6];
  register double wp_eff;

  for(i=0;i<6;i++){
    strain_aux[i] = stress_aux[i] = 0.0;
  }
  vol = 0.0;
  for(e=0;e<nelm;e++){

    npe = eptr[e+1]-eptr[e];
    ngp = npe;
    GetElemCoord(&eind[eptr[e]], npe, coor_elm);
    GetWeight(npe, &wp);

    // calculate <ElemResidue> by numerical integration
    for(gp=0;gp<ngp;gp++){
      memset(stress_gp,0.0,6*sizeof(double));
      get_dsh(gp, npe, coor_elm, dsh, &detj);
      GetDsDe( e, disp_elm, DsDe );
      wp_eff = detj*wp[gp];
      for(i=0;i<6;i++){
	for(k=0;k<6;k++){
	  stress_gp[i] += DsDe[i][k]*strain_mac[k];
	}
      }
      for(i=0;i<6;i++){
	stress_aux[i] += stress_gp[i] * wp_eff;
	strain_aux[i] += strain_mac[i] * wp_eff;
      }
      vol += wp_eff;
    }
  }
  ierr = MPI_Allreduce(stress_aux, stress_ave, 6, MPI_DOUBLE, MPI_SUM, PROBLEM_COMM);CHKERRQ(ierr);
  ierr = MPI_Allreduce(strain_aux, strain_ave, 6, MPI_DOUBLE, MPI_SUM, PROBLEM_COMM);CHKERRQ(ierr);
  ierr = MPI_Allreduce(&vol, &vol_tot, 1, MPI_DOUBLE, MPI_SUM, PROBLEM_COMM);CHKERRQ(ierr);
  for(i=0;i<6;i++){
    stress_ave[i] /= vol_tot;
    strain_ave[i] /= vol_tot;
  }

  return 0;
}
/****************************************************************************************************/
int micro_homogenize_linear(MPI_Comm MICRO_COMM, double strain_mac[6], double strain_ave[6], double stress_ave[6])
{

  int    ierr, nr_its=-1, max_its=6, nnods_bc, i, d, k, n_loc, *ixb;
  double norm_tol=1.0e-8, norm=2*norm_tol, ub, strain_matrix[3][3], *x_arr, *b_arr;

  if(dim==2){
    strain_matrix[0][0]=strain_mac[0]; strain_matrix[0][1]=strain_mac[2];
    strain_matrix[1][0]=strain_mac[2]; strain_matrix[1][1]=strain_mac[1];
  }
  else if(dim==3){
    strain_matrix[0][0]=strain_mac[0]; strain_matrix[0][1]=strain_mac[3]; strain_matrix[0][2]=strain_mac[5];
    strain_matrix[1][0]=strain_mac[3]; strain_matrix[1][1]=strain_mac[1]; strain_matrix[1][2]=strain_mac[4];
    strain_matrix[2][0]=strain_mac[5]; strain_matrix[2][1]=strain_mac[4]; strain_matrix[2][2]=strain_mac[2];
  }

  ierr = VecGetArray(x, &x_arr);CHKERRQ(ierr);
  nnods_bc = ((homog_ld_lagran_t*)homo.st)->nnods_bc;
  ixb = ((homog_ld_lagran_t*)homo.st)->index;

  // first we fill the value ub_val in <homog_ld_lagran_t> structure
  for(i=0; i<nnods_bc; i++){
    for(d=0;d<dim;d++){
      n_loc = ixb[i*dim+d]/dim;
      ub = 0.0;
      for(k=0;k<dim;k++){
	ub += strain_matrix[d][k]*coord[n_loc * 3 + k];
      }
      ((homog_ld_t*)homo.st)->ub_val[i*dim+d] = ub;
    }
  }
  for(i=0; i<nnods_bc; i++){
    for(d=0;d<dim;d++){
      x_arr[ixb[i*dim+d]] = ((homog_ld_t*)homo.st)->ub_val[i*dim+d]; 
    }
  }

  while( nr_its < max_its && norm > norm_tol )
  {
    ierr = assembly_residual_sd( &x, &b);CHKERRQ(ierr);
    ierr = VecGetArray(b, &b_arr);CHKERRQ(ierr);
    for(i=0; i<nnods_bc; i++){
      for(d=0;d<dim;d++){
	b_arr[ixb[i*dim+d]] = 0.0;
      }
    }
    ierr = VecRestoreArray(b, &b_arr);CHKERRQ(ierr);
    ierr = VecNorm(b,NORM_2,&norm);CHKERRQ(ierr);
    PetscPrintf(MICRO_COMM,"|b| = %lf \n",norm);
    if( !(norm > norm_tol) )break;
    ierr = VecScale(b,-1.0); CHKERRQ(ierr);
    ierr = assembly_jacobian_sd(&A);
    ierr = MatZeroRowsColumns(A, nnods_bc*dim, ixb, 1.0, NULL, NULL); CHKERRQ(ierr);
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
  int    ierr, nr_its=-1, max_its=6, nnods_bc, i, d, k, n_loc, *ixb;
  double norm_tol=1.0e-8, norm=2*norm_tol, ub, strain_matrix[3][3], *x_arr, *b_arr;

  if(dim==2){
    strain_matrix[0][0]=strain_mac[0]; strain_matrix[0][1]=strain_mac[2];
    strain_matrix[1][0]=strain_mac[2]; strain_matrix[1][1]=strain_mac[1];
  }
  else if(dim==3){
    strain_matrix[0][0]=strain_mac[0]; strain_matrix[0][1]=strain_mac[3]; strain_matrix[0][2]=strain_mac[5];
    strain_matrix[1][0]=strain_mac[3]; strain_matrix[1][1]=strain_mac[1]; strain_matrix[1][2]=strain_mac[4];
    strain_matrix[2][0]=strain_mac[5]; strain_matrix[2][1]=strain_mac[4]; strain_matrix[2][2]=strain_mac[2];
  }

  nnods_bc = ((homog_ld_lagran_t*)homo.st)->nnods_bc;
  ixb = ((homog_ld_lagran_t*)homo.st)->index;

  // first we fill the value ub_val in <homog_ld_lagran_t> structure
  for(i=0; i<nnods_bc; i++){
    for(d=0;d<dim;d++){
      n_loc = ixb[i*dim+d]/dim;
      ub = 0.0;
      for(k=0;k<dim;k++){
	ub += strain_matrix[d][k]*coord[n_loc * 3 + k];
      }
      ((homog_ld_lagran_t*)homo.st)->ub_val[i*dim+d] = ub;
    }
  }

  while( nr_its < max_its && norm > norm_tol )
  {

    /* internal forces */
    ierr = VecGetArray(x, &x_arr);CHKERRQ(ierr);
    ierr = assembly_residual_sd(&x, &b);CHKERRQ(ierr);
    /* fb - delta */
    ierr = VecGetArray(b, &b_arr);CHKERRQ(ierr);
    for(i=0; i<nnods_bc; i++){
      for(d=0;d<dim;d++){
	b_arr[ixb[i*dim+d]] -= x_arr[nmynods*dim+i*dim+d];
      }
    }
    /* ub - D^T . e_mac */
    for(i=0; i<nnods_bc; i++){
      for(d=0;d<dim;d++){
	b_arr[nmynods*dim + i*dim+d] = x_arr[ixb[i*dim+d]] - ((homog_ld_lagran_t*)homo.st)->ub_val[i*dim+d];
      }
    }
    ierr = VecRestoreArray(b, &b_arr);CHKERRQ(ierr);
    ierr = VecRestoreArray(x, &x_arr);CHKERRQ(ierr);

    ierr = VecNorm(b,NORM_2,&norm);CHKERRQ(ierr);
    PetscPrintf(MICRO_COMM,"|b| = %lf \n",norm);
    if( !(norm > norm_tol) ) break;
    ierr = VecScale(b,-1.0); CHKERRQ(ierr);

    /* Tangent matrix 

         | K    1  |
     A = |         |
         | 1    0  |
    */
    ierr = assembly_jacobian_sd(&A);
    for(i=0; i<nnods_bc; i++){
      for(d=0;d<dim;d++){
	MatSetValue(A,ixb[i*dim+d],nmynods*dim+i*dim+d,-1.0,INSERT_VALUES);
      }
    }
    for(i=0; i<nnods_bc; i++){
      for(d=0;d<dim;d++){
	MatSetValue(A,nmynods*dim+i*dim+d,ixb[i*dim+d],1.0,INSERT_VALUES);
      }
    }
    for(i=0; i<nnods_bc; i++){
      for(d=0;d<dim;d++){
	MatSetValue(A,nmynods*dim+i*dim+d,nmynods*dim+i*dim+d,0.0,INSERT_VALUES);
      }
    }
    MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);
  
    ierr = KSPSolve(ksp, b, dx);CHKERRQ(ierr);
    ierr = VecAXPY(x, 1.0, dx); CHKERRQ(ierr);

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
int micro_homogenize(MPI_Comm MICRO_COMM, double strain_mac[6], double strain_ave[6], double stress_ave[6])
{

  /*
     Performs linear homogenization of the RVE 
     Possible types are>

     HOMO_TAYLOR, HOMO_LINEAR, HOMO_PERIOD, HOMO_LD_LAGRAN_SEQ, HOMO_EXP

     LD            > Linear Displacements
     LD_LAGRAN_SEQ > Linear Displacement with Lagrangian multiplier BC

     HOMO_TAYLOR > no need of calculating a displacement field
   */
  int ierr;

  if(homo.type==HOMO_TAYLOR){
    ierr = micro_homogenize_taylor(MICRO_COMM, strain_mac, strain_ave, stress_ave);CHKERRQ(ierr);
  }
  else if(homo.type==LD){
    ierr = micro_homogenize_linear(MICRO_COMM, strain_mac, strain_ave, stress_ave);CHKERRQ(ierr);
  }
  else if(homo.type==HOMO_LINEAR_HEXA){
    ierr = micro_homogenize_linear_hexa(MICRO_COMM, strain_mac, strain_ave, stress_ave);CHKERRQ(ierr);
  }
  else if(homo.type==LD_LAGRAN_SEQ){
    ierr = mic_homogenize_ld_lagran(MICRO_COMM, strain_mac, strain_ave, stress_ave);CHKERRQ(ierr);
  }
  else{
    return 1;
  }

  return 0;
}
/****************************************************************************************************/
int mic_init_homo(void)
{
  
  if(homo.type==LD_LAGRAN_SEQ){

    /*
       a) Count how many nodes (nnods_bc) belongs 
       to the boundary search in the boundary_list
       (deleted the ones which are repeated)
       b) Allocate <index> and <ub_val>

       <index>  correspond to the global numeration where 
       the values should be assembled on the matrix and 
       vector
       <ub_val> corresponds to the product <D^T . e_mac>
    */

    homo.st = malloc(sizeof(homog_ld_lagran_t));

    int nnods_bc=0, nnods_bc_a, c, n_glo, n_loc, n_orig, i, d, *p, *aux_nod;
    node_list_t *pb,*pn;
    pb = boundary_list.head;
    while(pb){
      nnods_bc += ((boundary_t*)pb->data)->Nods.sizelist;
      pb=pb->next;
    }
    aux_nod = malloc(nnods_bc*sizeof(int));

    i=0;
    pb = boundary_list.head;
    while(pb){
      pn = ((boundary_t*)pb->data)->Nods.head;
      while(pn){
	n_orig = *(int*)pn->data;
	p = bsearch(&n_orig, mynods, nmynods, sizeof(int), cmpfunc); 
	if(!p){
	  PetscPrintf(MICRO_COMM,
	      "A boundary node (%d) seems now to not belong to this process (rank:%d)",n_orig,rank_mic);
	  return 1;
	}
	n_loc = p - mynods;       // Local numeration
	n_glo = loc2petsc[n_loc]; // PETSc numeration
	aux_nod[i]  = n_glo;
	i++;
	pn=pn->next;
      }
      pb=pb->next;
    }

    /*
       we order the vector of boundary nodes in order to delete those 
       nodes which are repeated
     */
    qsort(aux_nod, nnods_bc, sizeof(int), cmpfunc);
    for(i=0;i<nnods_bc;i++){
      if(i!=0){
	if(aux_nod[i]!=aux_nod[i-1]){
	  nnods_bc_a++;
	}
      }
      else if(i==0){
	nnods_bc_a=1;
      }
    }

    ((homog_ld_lagran_t*)homo.st)->nnods_bc = nnods_bc_a;
    ((homog_ld_lagran_t*)homo.st)->index  = malloc(nnods_bc_a*dim*sizeof(int));
    ((homog_ld_lagran_t*)homo.st)->ub_val = malloc(nnods_bc_a*dim*sizeof(double));

    c = 0;
    for(i=0;i<nnods_bc;i++){
      if(i!=0){
	if(aux_nod[i]!=aux_nod[i-1]){
	  for(d=0;d<dim;d++){
	    ((homog_ld_lagran_t*)homo.st)->index[c*dim+d] = aux_nod[i]*dim + d;
	  }
	  c++;
	}
      }
      else if(i==0){
	for(d=0;d<dim;d++){
	  ((homog_ld_lagran_t*)homo.st)->index[c*dim+d] = aux_nod[i]*dim + d;
	}
	c++;
      }
    }

  }

  if(homo.type==LD){

    /*
       a) Count how many nodes (nnods_bc) belongs 
       to the boundary search in the boundary_list
       (deleted the ones which are repeated)
       b) Allocate <index> and <ub_val>

       <index>  correspond to the global numeration where 
       the values should be assembled on the matrix and 
       vector
       <ub_val> corresponds to the product <D^T . e_mac>
    */

    homo.st = malloc(sizeof(homog_ld_t));

    int nnods_bc=0, nnods_bc_a, c, n_glo, n_loc, n_orig, i, d, *p, *aux_nod;
    node_list_t *pb,*pn;
    pb = boundary_list.head;
    while(pb){
      nnods_bc += ((boundary_t*)pb->data)->Nods.sizelist;
      pb=pb->next;
    }
    aux_nod = malloc(nnods_bc*sizeof(int));

    i=0;
    pb = boundary_list.head;
    while(pb){
      pn = ((boundary_t*)pb->data)->Nods.head;
      while(pn){
	n_orig = *(int*)pn->data;
	p = bsearch(&n_orig, mynods, nmynods, sizeof(int), cmpfunc); 
	if(!p){
	  PetscPrintf(MICRO_COMM,
	      "A boundary node (%d) seems now to not belong to this process (rank:%d)",n_orig,rank_mic);
	  return 1;
	}
	n_loc = p - mynods;       // Local numeration
	n_glo = loc2petsc[n_loc]; // PETSc numeration
	aux_nod[i]  = n_glo;
	i++;
	pn=pn->next;
      }
      pb=pb->next;
    }

    /*
       we order the vector of boundary nodes in order to delete those 
       nodes which are repeated
     */
    qsort(aux_nod, nnods_bc, sizeof(int), cmpfunc);
    for(i=0;i<nnods_bc;i++){
      if(i!=0){
	if(aux_nod[i]!=aux_nod[i-1]){
	  nnods_bc_a++;
	}
      }
      else if(i==0){
	nnods_bc_a=1;
      }
    }

    ((homog_ld_t*)homo.st)->nnods_bc = nnods_bc_a;
    ((homog_ld_t*)homo.st)->index  = malloc(nnods_bc_a*dim*sizeof(int));
    ((homog_ld_t*)homo.st)->ub_val = malloc(nnods_bc_a*dim*sizeof(double));

    c = 0;
    for(i=0;i<nnods_bc;i++){
      if(i!=0){
	if(aux_nod[i]!=aux_nod[i-1]){
	  for(d=0;d<dim;d++){
	    ((homog_ld_t*)homo.st)->index[c*dim+d] = aux_nod[i]*dim + d;
	  }
	  c++;
	}
      }
      else if(i==0){
	for(d=0;d<dim;d++){
	  ((homog_ld_t*)homo.st)->index[c*dim+d] = aux_nod[i]*dim + d;
	}
	c++;
      }
    }

  }
  return 0;
}
/****************************************************************************************************/
int micro_homogenize_linear_hexa(MPI_Comm MICRO_COMM, double strain_mac[6], double strain_ave[6], double stress_ave[6])
{

  /*
     Performs linear homogenization of the RVE cell acording to <strain_mac[6]>
   */

  int ierr, nr_its = -1; 
  //  int kspits = -1;
  double norm = -1.0, norm_tol = 1.0e-8, max_its = 3; 
  //  double ksp_norm = -1.0;

  ierr = micro_apply_bc_linear_hexa(strain_mac, &x, &A, &b, SET_DISPLACE);CHKERRQ(ierr);

  while( nr_its < max_its && norm > norm_tol )
  {
    ierr = assembly_residual_sd( &x, &b);CHKERRQ(ierr);
    ierr = micro_apply_bc_linear_hexa(strain_mac, &x, &A, &b, SET_RESIDUAL);CHKERRQ(ierr);
    ierr = VecNorm(b,NORM_2,&norm);CHKERRQ(ierr);
    if( !(norm > norm_tol) )break;
    ierr = VecScale(b,-1.0); CHKERRQ(ierr);
    ierr = assembly_jacobian_sd(&A);
    ierr = micro_apply_bc_linear_hexa(strain_mac, &x, &A, &b, SET_JACOBIAN);CHKERRQ(ierr);
    ierr = KSPSolve(ksp,b,dx);CHKERRQ(ierr);
    ierr = VecAXPY( x, 1.0, dx); CHKERRQ(ierr);
    nr_its ++;
  }

  ierr = calc_ave_strain_stress(MICRO_COMM, &x, strain_ave, stress_ave);CHKERRQ(ierr);

  return 0;
}

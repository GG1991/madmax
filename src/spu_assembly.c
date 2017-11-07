/*
   Routines for matrix assembly

   Author > Guido Giuntoli
   Date   > 03-08-2017
  
 */

#include "sputnik.h"
#include "comm.h"

int assembly_jacobian_sd(Mat *J)
{
  /* Assembly the Jacobian for Small Deformation approach */

  int         i, j, k, l, e, gp, ngp, npe;
  int         index[8*3];
  int         ierr;
  double      elm_coor[8][3];
  double      dsh[8][3];
  double      detj;
  double      Ke[8*3*8*3];
  double      B[6][3*8], strain_gp[6];
  double      c[6][6];
  double      *wp = NULL;
  double      elm_disp[8*3];
  double      wp_eff;
  double      *xvalues;
  Vec         xlocal;

  ierr = MatZeroEntries(*J);CHKERRQ(ierr);

  /* Local representation of <x> with ghost padding */
  ierr = VecGhostGetLocalForm(x,&xlocal); CHKERRQ(ierr);
  ierr = VecGetArray(xlocal, &xvalues); CHKERRQ(ierr);

  for(e=0;e<nelm;e++){

    npe = eptr[e+1]-eptr[e];
    ngp = npe;
    GetPETScIndeces( &eind[eptr[e]], npe, loc2petsc, index);
    get_elm_coor(&eind[eptr[e]], npe, elm_coor);
    GetWeight(npe, &wp);
    get_elm_disp(e, xvalues, elm_disp);

    // calculate <Ke> by numerical integration
    for(i=0;i<npe*dim*npe*dim;i++){
      Ke[i]=0.0;
    }
    for(gp=0;gp<ngp;gp++){

      ierr = get_dsh(gp, npe, elm_coor, dsh, &detj);
      detj = fabs(detj);
      ierr = GetB(npe, dsh, B);

      for(i=0;i<nvoi;i++){
	strain_gp[i] = 0.0;
	for(k=0;k<npe*dim;k++){
	  strain_gp[i] = strain_gp[i] + B[i][k] * elm_disp[k];
	}
      }
      ierr = get_c(NULL, e, gp, strain_gp, c);
      if(ierr){
	PetscPrintf(PETSC_COMM_WORLD, "%s: problem calculating constitutive tensor\n",myname);
	return 1; 
      }

      wp_eff = detj * wp[gp];
      for(i=0;i<npe*dim;i++){
	for(j=0;j<npe*dim;j++){
	  for(k=0;k<nvoi;k++){
	    for(l=0;l<nvoi;l++){
	      Ke[i*npe*dim+j] = Ke[i*npe*dim+j] + B[k][i] * c[k][l] * B[l][j] * wp_eff;
	    }
	  }
	}
      }

    }
    ierr = MatSetValues(*J, npe*dim, index, npe*dim, index, Ke, ADD_VALUES);CHKERRQ(ierr);
  }
  VecRestoreArray(xlocal,&xvalues); CHKERRQ(ierr);

  /* communication between processes */
  ierr = MatAssemblyBegin(*J, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(*J, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

  return 0;
}
/****************************************************************************************************/
int assembly_mass(Mat *M)
{
  /* Assembly the Jacobian for Small Deformation approach */

  int         i, j, d, e, gp, ngp, npe, ierr;
  int         index[8*3];
  double      elm_coor[8][3];
  double      dsh[8][3], **sh = NULL, *wp = NULL, detj, wp_eff;
  double      Me[8*3*8*3];
  double      rho_gp;

  ierr = MatZeroEntries(*M);CHKERRQ(ierr);

  for(e=0;e<nelm;e++){

    npe = eptr[e+1]-eptr[e];
    ngp = npe;
    GetPETScIndeces(&eind[eptr[e]], npe, loc2petsc, index);
    get_elm_coor(&eind[eptr[e]], npe, elm_coor);
    GetWeight(npe, &wp);

    // calculate <Ke> by numerical integration
    for(i=0;i<npe*dim*npe*dim;i++){
      Me[i]=0.0;
    }

    for(gp=0;gp<ngp;gp++){

      ierr = get_dsh(gp, npe, elm_coor, dsh, &detj);
      ierr = fem_get_sh(npe, dim, &sh);
      detj = fabs(detj);

      ierr = get_rho(NULL, e, &rho_gp);

      wp_eff = detj * wp[gp];
      for(d=0;d<dim;d++){
	for(i=0;i<npe;i++){
	  for(j=0;j<npe;j++){
	    Me[ (i*dim)*(npe*dim) + j*dim + (d*dim*npe + d) ] += rho_gp * sh[i][gp] * sh[j][gp] * wp_eff;
	  }
	}
      }

    }
    ierr = MatSetValues(*M, npe*dim, index, npe*dim, index, Me, ADD_VALUES);CHKERRQ(ierr);
  }

  /* communication between processes */
  ierr = MatAssemblyBegin(*M, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(*M, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

  return 0;
}
/****************************************************************************************************/
int assembly_residual_sd(Vec *x, Vec *b)
{
  /* Assembly the Residual for Small Deformation approach */

  int         i, k, e, gp, ngp, npe;
  int         index[8*3];
  int         ierr;
  double      elm_coor[8][3];
  double      dsh[8][3];
  double      detj;
  double      Re[8*3];
  double      B[6][3*8], stress_gp[6], strain_gp[6];
  double      DsDe[6][6];
  double      *wp = NULL;
  double      elm_disp[8*3];
  double      wp_eff;
  double      *xvalues;
  Vec         xlocal;
  material_t  *mat = NULL;

  /* Local representation of <x> with ghost padding */
  ierr = VecZeroEntries(*b); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(*x,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(*x,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostGetLocalForm(*x,&xlocal); CHKERRQ(ierr);
  ierr = VecGetArray(xlocal, &xvalues); CHKERRQ(ierr);

  for(e=0;e<nelm;e++){

    ierr = get_mat_from_elem(e, &mat);
    if(!mat){
      ierr = PetscPrintf(PETSC_COMM_WORLD, "material with physical_id %d not found\n",elm_id[e]);
      return 1;
    }

    npe = eptr[e+1]-eptr[e];
    ngp = npe;
    GetPETScIndeces(&eind[eptr[e]], npe, loc2petsc, index);
    get_elm_coor(&eind[eptr[e]], npe, elm_coor);
    GetWeight(npe, &wp);
    get_elm_disp(e, xvalues, elm_disp);

    // calculate <ElemResidue> by numerical integration

    for(i=0;i<npe*dim;i++){
      Re[i] = 0.0;
    }
    for(gp=0;gp<ngp;gp++){

      ierr = get_dsh(gp, npe, elm_coor, dsh, &detj);
      detj = fabs(detj);
      ierr = GetB(npe, dsh, B );

      for(i=0;i<nvoi;i++){
	strain_gp[i] = 0.0;
	for(k=0;k<npe*dim;k++){
	  strain_gp[i] = strain_gp[i] + B[i][k] * elm_disp[k];
	}
      }

      if(mat->type_id==MICRO){

	/* calculate stress using localization+homogenization */

	/* send instruction to micro */
	ierr = mac_send_signal(WORLD_COMM, MAC2MIC_STRAIN);
	if(ierr){
	  PetscPrintf(PETSC_COMM_WORLD, "%s: problem sending signal MAC2MIC_STRAIN to micro\n",myname);
	  return 1;
	}
	/* send strain to micro */
	ierr = mac_send_strain(WORLD_COMM, strain_gp);
	if(ierr){
	  ierr = PetscPrintf(PETSC_COMM_WORLD, "macro: problem sending strain to micro\n");
	  return 1;
	}
	/* recv stress from micro */
	ierr = mac_recv_stress(WORLD_COMM, stress_gp);
	if(ierr){
	  ierr = PetscPrintf(PETSC_COMM_WORLD, "macro: problem receiving stress from micro\n");
	  return 1;
	}

      }else{

	ierr = get_c(NULL, e, gp, strain_gp, DsDe);
	if(ierr){
	  PetscPrintf(PETSC_COMM_WORLD, "%s: problem calculating constitutive tensor\n",myname);
	  return 1; 
	}

	/* calculate stress using constitutive laws */
	for(i=0;i<nvoi;i++){
	  stress_gp[i] = 0.0;
	  for(k=0;k<nvoi;k++){
	    stress_gp[i] = stress_gp[i] + DsDe[i][k] * strain_gp[k];
	  }
	}

      }

      wp_eff = detj * wp[gp];
      for(i=0;i<npe*dim;i++){
	for(k=0;k<nvoi;k++){
	  Re[i] = Re[i] + B[k][i] * stress_gp[k] * wp_eff;
	}
      }

    }
    ierr = VecSetValues(*b, npe*dim, index, Re, ADD_VALUES);CHKERRQ(ierr);
  }
  VecRestoreArray(xlocal,&xvalues); CHKERRQ(ierr);

  /* communication between processes */
  ierr = VecAssemblyBegin(*b);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(*b);CHKERRQ(ierr);
  return 0;

}
/****************************************************************************************************/
int calc_strain_stress_energy(Vec *x, double *strain, double *stress, double *energy)
{
  /*    
	Calculate averange strain and stress tensors on each element
   */
  int    i, k, e, gp, ngp, npe;
  int    ierr;

  double elm_coor[8][3];
  double dsh[8][3];
  double detj;
  double B[6][3*8], stress_gp[6], strain_gp[6], energy_gp, stress_ave[6], strain_ave[6], energy_ave;
  double DsDe[6][6];
  double *wp = NULL;
  double elm_disp[8*3];
  double *xvalues;
  double vol = -1.0;

  Vec xlocal;

  register double wp_eff;

  /* 
     Local representation of <x> with ghost padding
  */
  ierr = VecGhostUpdateBegin(*x,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(*x,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostGetLocalForm(*x,&xlocal); CHKERRQ(ierr);
  ierr = VecGetArray(xlocal, &xvalues); CHKERRQ(ierr); CHKERRQ(ierr);

  for(e=0;e<nelm;e++){

    npe = eptr[e+1]-eptr[e];
    ngp = npe;

    get_elm_coor(&eind[eptr[e]], npe, elm_coor);
    GetWeight(npe, &wp);
    get_elm_disp( e, xvalues, elm_disp );

    vol = 0.0;
    memset(strain_ave,0.0,nvoi*sizeof(double));
    memset(stress_ave,0.0,nvoi*sizeof(double));
    energy_ave = 0.0;

    for(gp=0;gp<ngp;gp++){
      
      energy_gp = 0.0;

      get_dsh(gp, npe, elm_coor, dsh, &detj);
      detj = fabs(detj);
      GetB(npe, dsh, B );
      wp_eff = detj * wp[gp];

      for(i=0;i<nvoi;i++){
	strain_gp[i] = 0.0;
	for(k=0;k<npe*dim;k++){
	  strain_gp[i] = strain_gp[i] + B[i][k] * elm_disp[k];
	}
      }
      get_c(NULL, e, gp, strain_gp, DsDe);

      for(i=0;i<nvoi;i++){
	stress_gp[i] = 0.0;
	for(k=0;k<nvoi;k++){
	  stress_gp[i] = stress_gp[i] + DsDe[i][k] * strain_gp[k];
	}
      }
      for(k=0;k<nvoi;k++){
	energy_gp += stress_gp[k]*strain_gp[k];
      }

      for(i=0;i<nvoi;i++){
	stress_ave[i] += stress_gp[i] * wp_eff;
	strain_ave[i] += strain_gp[i] * wp_eff;
      }
      energy_ave    += energy_gp * wp_eff;

      vol += wp_eff;
    }
    for(i=0;i<nvoi;i++){
      strain[e*nvoi+i] = strain_ave[i] / vol; 
      stress[e*nvoi+i] = stress_ave[i] / vol;
    }
    energy[e] = energy_ave / vol;

  }
  VecRestoreArray(xlocal,&xvalues); CHKERRQ(ierr);

  return 0;
}
/****************************************************************************************************/
int calc_ave_strain_stress(MPI_Comm PROBLEM_COMM, Vec *x, double strain_ave[6], double stress_ave[6])
{
  /*    
	Calculate averange strain and stress tensors on each element

        Input
	<Vec x>  > distributed vector of displacements
	Ouput
	strain_ave[6] (same on every process)
	stress_ave[6] (same on every process)

   */
  int    i, k, e, gp, ngp, npe;
  int    ierr;
  int    nproc, rank;

  double elm_coor[8][3];
  double dsh[8][3];
  double detj;
  double B[6][3*8], stress_gp[6], strain_gp[6];
  double DsDe[6][6];
  double *wp = NULL;
  double elm_disp[8*3];
  double *xvalues;
  double vol = -1.0, vol_tot = -1.0;
  double stress_aux[6];
  double strain_aux[6];

  MPI_Comm_size(PROBLEM_COMM, &nproc);
  MPI_Comm_rank(PROBLEM_COMM, &rank);

  Vec xlocal;

  register double wp_eff;

  /* 
     Local representation of <x> with ghost padding
  */
  ierr = VecGhostUpdateBegin(*x,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(*x,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostGetLocalForm(*x,&xlocal); CHKERRQ(ierr);
  ierr = VecGetArray(xlocal, &xvalues); CHKERRQ(ierr); CHKERRQ(ierr);

  for(i=0;i<nvoi;i++){
    strain_aux[i] = stress_aux[i] = 0.0;
  }
  vol = 0.0;

  for(e=0;e<nelm;e++){

    npe = eptr[e+1]-eptr[e];
    ngp = npe;

    get_elm_coor(&eind[eptr[e]], npe, elm_coor);
    GetWeight(npe, &wp);
    get_elm_disp( e, xvalues, elm_disp );

    // calculate <ElemResidue> by numerical integration

    for(gp=0;gp<ngp;gp++){

      get_dsh(gp, npe, elm_coor, dsh, &detj);
      detj = fabs(detj);
      GetB( npe, dsh, B );
      wp_eff = detj*wp[gp];

      for(i=0;i<nvoi;i++){
	strain_gp[i] = 0.0;
	for(k=0;k<npe*dim;k++){
	  strain_gp[i] = strain_gp[i] + B[i][k] * elm_disp[k];
	}
      }
      get_c(NULL, e, gp, elm_disp, DsDe);

      for(i=0;i<nvoi;i++){
	stress_gp[i] = 0.0;
	for(k=0;k<nvoi;k++){
	  stress_gp[i] = stress_gp[i] + DsDe[i][k] * strain_gp[k];
	}
      }

      for(i=0;i<6;i++){
	stress_aux[i] += stress_gp[i] * wp_eff;
	strain_aux[i] += strain_gp[i] * wp_eff;
      }

      vol += wp_eff;
    }

  }
  VecRestoreArray(xlocal,&xvalues); CHKERRQ(ierr);

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
int calc_rho(MPI_Comm PROBLEM_COMM, double *rho)
{

  int    e, gp, ngp, npe, ierr;
  double elm_coor[8][3], dsh[8][3], detj, *wp=NULL, wp_eff;
  double vol = -1.0, vol_tot = -1.0;
  double rho_gp, rho_a;

  vol = 0.0;

  for(e=0;e<nelm;e++){

    npe = eptr[e+1]-eptr[e];
    ngp = npe;

    get_elm_coor(&eind[eptr[e]], npe, elm_coor);
    GetWeight(npe, &wp);

    // calculate <ElemResidue> by numerical integration

    for(gp=0;gp<ngp;gp++){

      get_dsh(gp, npe, elm_coor, dsh, &detj);
      detj = fabs(detj);
      wp_eff = detj*wp[gp];

      get_rho(NULL, e, &rho_gp);
      rho_a += rho_gp * wp_eff;

      vol += wp_eff;
    }

  }

  ierr = MPI_Allreduce(&rho_a, rho, 1, MPI_DOUBLE, MPI_SUM, PROBLEM_COMM);CHKERRQ(ierr);
  ierr = MPI_Allreduce(&vol, &vol_tot, 1, MPI_DOUBLE, MPI_SUM, PROBLEM_COMM);CHKERRQ(ierr);

  *rho = *rho / vol_tot;

  return 0;
}
/****************************************************************************************************/
int get_elm_disp(int e, double *x, double *elm_disp)
{

  int  d, n, npe;

  npe = eptr[e+1]-eptr[e];
  for(n=0;n<npe;n++){
    for(d=0;d<dim;d++){
      // para usar VecGetValues usamos la numeracion global
      elm_disp[n*dim+d] = x[eind[eptr[e]+n]*dim+d];
    }
  }
  
  return 0;
}
/****************************************************************************************************/
int get_rho(const char *name, int e, double *rho)
{
  /*  Calculates constitutive tensor */

  int         ierr;
  material_t  *mat = NULL;

  if(name!=NULL){
    ierr = get_mat_from_name(name, &mat);
    if(!mat){
      PetscPrintf(PETSC_COMM_WORLD,"material with name %s not found",name);
      return 1;
    }
  }
  else{
    ierr = get_mat_from_elem(e, &mat);
    if(!mat){
      PetscPrintf(PETSC_COMM_WORLD,"material of element %d and id %d not found", e, elm_id[e]);
      return 1;
    }
  }

  switch(mat->type_id){

    case TYPE_0:
      /* 
	 Elástico lineal 
       */
      *rho  = ((type_0*)mat->type)->rho;
      break;

    case MICRO:

      /* send instruction to micro */
      ierr = mac_send_signal(WORLD_COMM, RHO);
      if(ierr){
	ierr = PetscPrintf(PETSC_COMM_WORLD, "macro: problem sending signal RHO to micro\n");
	return 1;
      }

      /* recv rho from micro */
      ierr = mac_recv_rho(WORLD_COMM, rho);
      if(ierr){
	ierr = PetscPrintf(PETSC_COMM_WORLD, "macro: problem receiving rho from micro\n");
	return 1;
      }

      break;

    default:
      break;

  }

  return 0;
}
/****************************************************************************************************/
int get_c(const char *name, int e, int gp, double strain[6], double c[6][6])
{
  /*  Calculates constitutive tensor */

  int         i, j, ierr;
  int         macro_gp;
  material_t  *mat = NULL;
  double      la, mu, poi, you; 
  double      c_homo[36];

  if(name!=NULL){
    ierr = get_mat_from_name(name, &mat);
    if(!mat){
      PetscPrintf(PETSC_COMM_WORLD,"material with name %s not found",name);
      return 1;
    }
  }
  else{
    macro_gp = e*8+gp;
    ierr = get_mat_from_elem(e, &mat);
    if(!mat){
      PetscPrintf(PETSC_COMM_WORLD,"material of element %d and id %d not found", e, elm_id[e]);
      return 1;
    }
  }

  switch(mat->type_id){

    case TYPE_0:
      /* 
	 Elástico lineal 
       */
      la  = ((type_0*)mat->type)->lambda;
      mu  = ((type_0*)mat->type)->mu;
      poi = ((type_0*)mat->type)->poisson;
      you = ((type_0*)mat->type)->young;

      if(dim==2){
	/*
	   Plane Strain
	*/
	c[0][0]=1.0-poi; c[0][1]=poi    ; c[0][2]=0.0        ;
	c[1][0]=poi    ; c[1][1]=1.0-poi; c[1][2]=0.0        ;
	c[2][0]=0.0    ; c[2][1]=0.0    ; c[2][2]=(1-2*poi)/2;
	for(i=0;i<3;i++){
	  for(j=0;j<3;j++){
	    c[i][j] = c[i][j] * you/((1+poi)*(1-2*poi));
	  }
	}
      }
      else if(dim==3){
	c[0][0]=la+2*mu ;c[0][1]=la      ;c[0][2]=la      ;c[0][3]=0.0; c[0][4]=0.0; c[0][5]=0.0;
	c[1][0]=la      ;c[1][1]=la+2*mu ;c[1][2]=la      ;c[1][3]=0.0; c[1][4]=0.0; c[1][5]=0.0;
	c[2][0]=la      ;c[2][1]=la      ;c[2][2]=la+2*mu ;c[2][3]=0.0; c[2][4]=0.0; c[2][5]=0.0;
	c[3][0]=0.0     ;c[3][1]=0.0     ;c[3][2]=0.0     ;c[3][3]=mu ; c[3][4]=0.0; c[3][5]=0.0;
	c[4][0]=0.0     ;c[4][1]=0.0     ;c[4][2]=0.0     ;c[4][3]=0.0; c[4][4]=mu ; c[4][5]=0.0;
	c[5][0]=0.0     ;c[5][1]=0.0     ;c[5][2]=0.0     ;c[5][3]=0.0; c[5][4]=0.0; c[5][5]=mu ;
      }

      break;

    case MICRO:

      /* send instruction to micro */
      ierr = mac_send_signal(WORLD_COMM, C_HOMO);
      if(ierr){
	ierr = PetscPrintf(PETSC_COMM_WORLD, "macro: problem sending signal C_HOMO to micro\n");
	return 1;
      }

      /* send strain to micro */
      ierr = mac_send_strain(WORLD_COMM, strain);
      if(ierr){
	ierr = PetscPrintf(PETSC_COMM_WORLD, "macro: problem sending strain to micro\n");
	return 1;
      }

      /* send macro_gp to micro */
      ierr = mac_send_macro_gp(WORLD_COMM, &macro_gp);
      if(ierr){
	ierr = PetscPrintf(PETSC_COMM_WORLD, "macro: problem sending macro_gp to micro\n");
	return 1;
      }

      /* recv c_homo from micro */
      ierr = mac_recv_c_homo(WORLD_COMM, nvoi, c_homo);
      if(ierr){
	ierr = PetscPrintf(PETSC_COMM_WORLD, "macro: problem receiving c_homo from micro\n");
	return 1;
      }

      for(i=0;i<nvoi;i++){
	for(j=0;j<nvoi;j++)
	  c[i][j] = c_homo[i*nvoi+j];
      }

      break;

    default:
      break;

  }
  return 0;
}
/****************************************************************************************************/
int get_mat_from_elem(int e, material_t **mat)
{
  /* 
     -fiber_middle <radious>
   */
  int          id;
  node_list_t  *pn;

  pn = material_list.head;
  while(pn){

    if(flag_fiber_cilin)
    {
      /* fiber in the middle */
      if(is_inside_fiber_cilin(e)){
	if(!strcmp(((material_t*)pn->data)->name,"FIBER")) break;
      }
      else{
	if(!strcmp(((material_t*)pn->data)->name,"MATRIX")) break;
      }
    }
    else{
      /* normal case */
      id = elm_id[e];
      if( ((material_t*)pn->data)->GmshID == id ) break;
    }
    pn = pn->next;
  }
  if(!pn){
    PetscPrintf(PETSC_COMM_WORLD, "%s: problem finding material from element %d\n", myname, e);
    return 1;
  }

  *mat = (material_t*)pn->data;

  return 0;
}
/****************************************************************************************************/
int get_mat_from_name(const char *name, material_t **mat)
{
  /* 
   */
  node_list_t  *pn;

  pn = material_list.head;
  while(pn){
    if(!strcmp(((material_t*)pn->data)->name,name)) break;
    pn = pn->next;
  }
  if(!pn){
    PetscPrintf(PETSC_COMM_WORLD, "%s: problem finding material from name %s\n", myname, name);
    return 1;
  }

  *mat = (material_t*)pn->data;

  return 0;
}
/****************************************************************************************************/
int is_inside_fiber_cilin(int e)
{
  int    i, j, d;
  double l=-1, centroid[3];
  double deviation[2];
  get_centroid(e, centroid);

  for(i=0;i<nx_fibers;i++){
    for(j=0;j<ny_fibers;j++){
      l=0.0;
      deviation[0] = fiber_cilin_center_devi[0] -LX/2 +(LX/nx_fibers)/2 + i*(LX/nx_fibers);
      deviation[1] = fiber_cilin_center_devi[1] -LY/2 +(LY/ny_fibers)/2 + j*(LY/ny_fibers);
      for(d=0;d<2;d++){
	l = l + pow(centroid[d]-(center_domain[d]+deviation[d]),2);
      }
      l = sqrt(l);
      if (l<=fiber_cilin_r) return 1;
    }
  }
  return 0;
}
/****************************************************************************************************/
int get_centroid(int e, double centroid[3])
{
  double   elem_coor[8][3];
  int      i, d;
  int      npe = eptr[e+1]-eptr[e];

  get_elem_coor(e, elem_coor);
  for(d=0;d<dim;d++){
    centroid[d] = 0.0;
    for(i=0;i<npe;i++){
      centroid[d] = centroid[d] + elem_coor[i][d];
    }
  }
  for(d=0;d<dim;d++){
    centroid[d] = centroid[d] / npe;
  }

  return 0;
}
/****************************************************************************************************/
int get_elem_coor(int e, double elem_coor[8][3])
{
  int i, d;
  int npe = eptr[e+1]-eptr[e];

  for(i=0;i<npe;i++){
    for(d=0;d<dim;d++){
      elem_coor[i][d] = coord[ eind[eptr[e]+i]*dim + d ];
    }
  }

  return 0;
}
/****************************************************************************************************/
int get_elem_vol(int e, double *vol)
{
  int    gp, npe, ngp;
  double elem_coor[8][3];
  double *wp;
  double dsh[8][3];
  double detj;

  npe = eptr[e+1]-eptr[e];
  ngp = npe;

  get_elm_coor(&eind[eptr[e]], npe, elem_coor);
  GetWeight(npe, &wp);

  *vol = 0;
  for(gp=0;gp<ngp;gp++){
    get_dsh(gp, npe, elem_coor, dsh, &detj);
    detj = fabs(detj);
    *vol += detj*wp[gp];
  }
  return 0;
}
/****************************************************************************************************/
int GetWeight(int npe, double **wp)
{
  *wp = FemGetPointer2Weight(npe, dim);
  return 0;
}
/****************************************************************************************************/
int GetB(int npe, double dsh[8][3], double B[6][3*8] )
{
  /*
     e = Bu    e=[exx eyy ezz exy eyz exz]

  e=[
  exx exy exz
  eyx eyy eyz
  ezx ezy ezz
  ]

  */
  int i;

  if(dim==2){
    for(i=0;i<npe;i++){
      B[0][i*2+0] = dsh[i][0]; B[0][i*2+1] = 0.0      ;
      B[1][i*2+0] = 0.0      ; B[1][i*2+1] = dsh[i][1];
      B[2][i*2+0] = dsh[i][1]; B[2][i*2+1] = dsh[i][0];
    }
  }
  else if(dim==3){
    for(i=0;i<npe;i++){
      B[0][i*3+0] = dsh[i][0]; B[0][i*3+1] = 0.0       ; B[0][i*3+2] = 0.0      ; 
      B[1][i*3+0] = 0.0      ; B[1][i*3+1] = dsh[i][1] ; B[1][i*3+2] = 0.0      ; 
      B[2][i*3+0] = 0.0      ; B[2][i*3+1] = 0.0       ; B[2][i*3+2] = dsh[i][2]; 
      B[3][i*3+0] = dsh[i][1]; B[3][i*3+1] = dsh[i][0] ; B[3][i*3+2] = 0.0      ; 
      B[4][i*3+0] = 0.0      ; B[4][i*3+1] = dsh[i][2] ; B[4][i*3+2] = dsh[i][1];
      B[5][i*3+0] = dsh[i][2]; B[5][i*3+1] = 0.0       ; B[5][i*3+2] = dsh[i][0]; 

    }
  }

  return 0;
}
/****************************************************************************************************/
int get_dsh(int gp, int npe, double coor[8][3], double dsh[8][3], double *detj)
{

  double ***ShapeDerivsMaster;

  double jac[3][3];
  double ijac[3][3];

  ShapeDerivsMaster = FemGetPointer2ShapeDerivsMaster(npe, dim);
  if(ShapeDerivsMaster == NULL) return 1;

  fem_calc_jac(dim, coor, ShapeDerivsMaster, npe, gp, jac);
  fem_invjac(dim, jac, ijac, detj);
  fem_trans_dsh(dim, ijac, npe, gp, ShapeDerivsMaster, dsh);

  return 0;
}
/****************************************************************************************************/
int GetPETScIndeces(int *LocalNod, int npe, int *local2PETSc, int *PETScIndex)
{
  /* 
     Gives the indeces to Gather fields and assembly on
     PETSc matrices and vectors

     Input>
     LocalNod    >  array with local node numeration
     npe         > number of nodes of the element
     local2PETSc > renumbering vetor to transform from local to PETSc

     Output>
     PETScIndex > Array with PETSc numeration

   */
  int i, d;

  for(i=0;i<npe;i++){
    for(d=0;d<dim;d++){
      PETScIndex[i*dim+d] = local2PETSc[LocalNod[i]]*dim+d;
    }
  }
  return 0;
}
/****************************************************************************************************/
int get_elm_coor(int *LocalNod, int n, double elm_coor[8][3])
{
  /*
      Returns the coordinates of the vertex of that element
   */

  int i, d;

  for(i=0;i<n;i++){
    for(d=0;d<dim;d++){
      elm_coor[i][d] = coord[ LocalNod[i]*dim + d ];
    }
  }

  return 0;
}

/*
   Routines for matrix assembly

   Author > Guido Giuntoli
   Date   > 03-08-2017
  
 */

#include "sputnik.h"
#include "macmic.h"

int assembly_jacobian_sd(Mat *J)
{
  /* Assembly the Jacobian for Small Deformation approach */

  int         i, j, k, l, e, gp, ngp, npe;
  int         PETScIdx[8*3];
  int         ierr;
  double      ElemCoord[8][3];
  double      dsh[8][3];
  double      detj;
  double      Ke[8*3*8*3];
  double      B[6][3*8], strain_gp[6];
  double      c[6][6];
  double      *wp = NULL;
  double      ElemDispls[8*3];
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
    GetPETScIndeces( &eind[eptr[e]], npe, loc2petsc, PETScIdx);
    GetElemCoord(&eind[eptr[e]], npe, ElemCoord);
    GetWeight(npe, &wp);
    GetElemenDispls(e, xvalues, ElemDispls);

    // calculate <Ke> by numerical integration
    memset(Ke, 0.0, (npe*dim*npe*dim)*sizeof(double));
    for(gp=0;gp<ngp;gp++){

      ierr = get_dsh(gp, npe, ElemCoord, dsh, &detj);
      detj = fabs(detj);
      ierr = GetB(npe, dsh, B);

      for(i=0;i<nvoi;i++){
	strain_gp[i] = 0.0;
	for(k=0;k<npe*dim;k++){
	  strain_gp[i] = strain_gp[i] + B[i][k] * ElemDispls[k];
	}
      }
      ierr = get_c(e, strain_gp, c);
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
    ierr = MatSetValues(*J, npe*dim, PETScIdx, npe*dim, PETScIdx, Ke, ADD_VALUES);CHKERRQ(ierr);
  }
  VecRestoreArray(xlocal,&xvalues); CHKERRQ(ierr);

  /* communication between processes */
  ierr = MatAssemblyBegin(*J, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(*J, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

  return 0;
}
/****************************************************************************************************/
int assembly_residual_sd(Vec *x_old, Vec *Residue)
{
  /* Assembly the Residual for Small Deformation approach */

  int         i, k, e, gp, ngp, npe;
  int         PETScIdx[8*3];
  int         ierr;
  double      ElemCoord[8][3];
  double      dsh[8][3];
  double      detj;
  double      Re[8*3];
  double      B[6][3*8], stress_gp[6], strain_gp[6];
  double      DsDe[6][6];
  double      *wp = NULL;
  double      ElemDispls[8*3];
  double      *xvalues;
  double      wp_eff;
  Vec         xlocal;
  material_t  *material;

  /* Local representation of <x> with ghost padding */
  ierr = VecZeroEntries(*Residue); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(*x_old,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(*x_old,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostGetLocalForm(*x_old,&xlocal); CHKERRQ(ierr);
  ierr = VecGetArray(xlocal, &xvalues); CHKERRQ(ierr);

  for(e=0;e<nelm;e++){

    material = GetMaterial(e);
    if(!material){
      ierr = PetscPrintf(PETSC_COMM_WORLD, "material with physical_id %d not found\n",PhysicalID[e]);
      return 1;
    }

    npe = eptr[e+1]-eptr[e];
    ngp = npe;
    GetPETScIndeces(&eind[eptr[e]], npe, loc2petsc, PETScIdx);
    GetElemCoord(&eind[eptr[e]], npe, ElemCoord);
    GetWeight(npe, &wp);
    GetElemenDispls(e, xvalues, ElemDispls);

    // calculate <ElemResidue> by numerical integration

    for(i=0;i<npe*dim;i++){
      Re[i] = 0.0;
    }
    for(gp=0;gp<ngp;gp++){

      ierr = get_dsh(gp, npe, ElemCoord, dsh, &detj);
      detj = fabs(detj);
      ierr = GetB(npe, dsh, B );

      for(i=0;i<nvoi;i++){
	strain_gp[i] = 0.0;
	for(k=0;k<npe*dim;k++){
	  strain_gp[i] = strain_gp[i] + B[i][k] * ElemDispls[k];
	}
      }
      ierr = get_c(e, strain_gp, DsDe);
      if(ierr){
	PetscPrintf(PETSC_COMM_WORLD, "%s: problem calculating constitutive tensor\n",myname);
	return 1; 
      }

      if(material->typeID==MICRO){

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
    ierr = VecSetValues(*Residue, npe*dim, PETScIdx, Re, ADD_VALUES);CHKERRQ(ierr);
  }
  VecRestoreArray(xlocal,&xvalues); CHKERRQ(ierr);

  /* communication between processes */
  ierr = VecAssemblyBegin(*Residue);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(*Residue);CHKERRQ(ierr);
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

  double ElemCoord[8][3];
  double dsh[8][3];
  double detj;
  double B[6][3*8], stress_gp[6], strain_gp[6], energy_gp, stress_ave[6], strain_ave[6], energy_ave;
  double DsDe[6][6];
  double *wp = NULL;
  double ElemDispls[8*3];
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

    GetElemCoord(&eind[eptr[e]], npe, ElemCoord);
    GetWeight(npe, &wp);
    GetElemenDispls( e, xvalues, ElemDispls );

    vol = 0.0;
    memset(strain_ave,0.0,nvoi*sizeof(double));
    memset(stress_ave,0.0,nvoi*sizeof(double));
    energy_ave = 0.0;

    for(gp=0;gp<ngp;gp++){
      
      energy_gp = 0.0;

      get_dsh(gp, npe, ElemCoord, dsh, &detj);
      detj = fabs(detj);
      GetB( npe, dsh, B );
      wp_eff = detj * wp[gp];

      for(i=0;i<nvoi;i++){
	strain_gp[i] = 0.0;
	for(k=0;k<npe*dim;k++){
	  strain_gp[i] = strain_gp[i] + B[i][k] * ElemDispls[k];
	}
      }
      get_c(e, strain_gp, DsDe);

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

      vol += detj*wp[gp];
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

  double ElemCoord[8][3];
  double dsh[8][3];
  double detj;
  double B[6][3*8], stress_gp[6], strain_gp[6];
  double DsDe[6][6];
  double *wp = NULL;
  double ElemDispls[8*3];
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

    GetElemCoord(&eind[eptr[e]], npe, ElemCoord);
    GetWeight(npe, &wp);
    GetElemenDispls( e, xvalues, ElemDispls );

    // calculate <ElemResidue> by numerical integration

    for(gp=0;gp<ngp;gp++){

      get_dsh(gp, npe, ElemCoord, dsh, &detj);
      detj = fabs(detj);
      GetB( npe, dsh, B );
      wp_eff = detj*wp[gp];

      for(i=0;i<nvoi;i++){
	strain_gp[i] = 0.0;
	for(k=0;k<npe*dim;k++){
	  strain_gp[i] = strain_gp[i] + B[i][k] * ElemDispls[k];
	}
      }
      get_c( e, ElemDispls, DsDe);

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
int GetElemenDispls(int e, double *x, double *ElemDispls )
{

  int  d, n, npe;

  npe = eptr[e+1]-eptr[e];
  for(n=0;n<npe;n++){
    for(d=0;d<dim;d++){
      // para usar VecGetValues usamos la numeracion global
      ElemDispls[n*dim+d] = x[eind[eptr[e]+n]*dim+d];
    }
  }
  
  return 0;
}
/****************************************************************************************************/
int get_c(int e, double strain[6], double c[6][6] )
{
  /*  Calculates constitutive tensor */

  double  la, mu, poi, you; 
  double  c_homo[36];
  int     i, j, ierr;

  material_t *material = GetMaterial(e);
  if(!material){
    PetscPrintf(PETSC_COMM_WORLD,"material with physical_id %d not found",PhysicalID[e]);
    return 1;
  }

  switch(material->typeID){

    case TYPE00:
      /* 
	 Elástico lineal 
       */
      la = ((type_00*)material->type)->lambda;
      mu = ((type_00*)material->type)->mu;
      poi = ((type_00*)material->type)->poisson;
      you = ((type_00*)material->type)->young;

      if(dim==2){
	/*
	   Plane Strain
	*/
	c[0][0]=1.0; c[0][1]=poi; c[0][2]=0.0;
	c[1][0]=poi; c[1][1]=1.0; c[1][2]=0.0;
	c[2][0]=0.0; c[2][1]=0.0; c[2][2]=(1-poi)/2;
	for(i=0;i<3;i++){
	  for(j=0;j<3;j++){
	    c[i][j] = c[i][j] * you/(1-pow(poi,2));
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
      /* recv stress from micro */
      ierr = mac_recv_c_homo(WORLD_COMM, c_homo);
      if(ierr){
	ierr = PetscPrintf(PETSC_COMM_WORLD, "macro: problem receiving c_homo from micro\n");
	return 1;
      }

      for(i=0;i<nvoi;i++){
	for(j=0;j<nvoi;j++){
	  c[i][j] = c_homo[i*nvoi+j];
	}
      }

      break;

    default:
      break;

  }
  return 0;
}
/****************************************************************************************************/
material_t * GetMaterial(int e)
{
  /* 
     Search in the <material_list> if for the one
     that has <GmshIDToSearch>
     this is normally done except if one of the 
     following options it is activated
     -fiber_middle <radious>
   */
  int         id;
  node_list_t *pn;

  pn = material_list.head;
  while(pn)
  {
    if(flag_fiber_cilin){
      if(is_inside_fiber_cilin(e)){
	if(!strcmp(((material_t*)pn->data)->name,"FIBER")) break;
      }
      else{
	if(!strcmp(((material_t*)pn->data)->name,"MATRIX")) break;
      }
    }
    else{
      id = PhysicalID[e];
      if( ((material_t*)pn->data)->GmshID == id ) break;
    }
    pn = pn->next;
  }
  if(!pn) return NULL;
  return (material_t*)pn->data;
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
int GetElemCoord(int *LocalNod, int n, double ElemCoord[8][3])
{
  /*
      Returns the coordinates of the vertex of that element
   */

  int i, d;

  for(i=0;i<n;i++){
    for(d=0;d<dim;d++){
      ElemCoord[i][d] = coord[ LocalNod[i]*dim + d ];
    }
  }

  return 0;
}

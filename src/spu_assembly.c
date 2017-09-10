/*
   Routines for matrix assembly

   Author > Guido Giuntoli
   Date   > 03-08-2017
  
 */

#include "sputnik.h"
#include "macmic.h"

int assembly_jacobian_sd(Mat *J)
{

  /*    
	Assembly the Jacobian for Small Deformation
	approach.
   */

  int    i, j, k, e, gp, ngp, npe;
  int    PETScIdx[8*3];
  int    ierr;

  double ElemCoord[8][3];
  double ShapeDerivs[8][3];
  double DetJac;
  double ElemMatrix[8*3 * 8*3];
  double B[6][3*8], Baux[6][3*8];
  double DsDe[6][6];
  double *wp = NULL;
  double ElemDispls[8*3];

  register double wp_eff;

  ierr = MatZeroEntries(*J);CHKERRQ(ierr);

  for(e=0;e<nelm;e++){
    npe = eptr[e+1]-eptr[e];
    ngp = npe;
    GetPETScIndeces( &eind[eptr[e]], npe, loc2petsc, PETScIdx);
    GetElemCoord(&eind[eptr[e]], npe, ElemCoord);
    GetWeight(npe, &wp);

    // calculate <ElemMatrix> by numerical integration

    memset(ElemMatrix, 0.0, (npe*dim*npe*dim)*sizeof(double));
    for(gp=0;gp<ngp;gp++){

      get_dsh(gp, npe, ElemCoord, ShapeDerivs, &DetJac);
      GetB( npe, ShapeDerivs, B );
      GetDsDe( e, ElemDispls, DsDe );

      for(i=0;i<nvoi;i++){
	for(j=0;j<npe*dim;j++){
	  Baux[i][j]=0.0;
	  for(k=0;k<nvoi;k++){
	    Baux[i][j] += DsDe[i][k]*B[k][j];
	  }
	}
      }

      wp_eff = DetJac * wp[gp];
      for(i=0;i<npe*dim;i++){
	for(j=0;j<npe*dim;j++){
	  for(k=0;k<nvoi;k++){
	    ElemMatrix[i*npe*dim+j] += B[k][i]*Baux[k][j] * wp_eff;
	  }
	}
      }

    }
    ierr = MatSetValues(*J, npe*dim, PETScIdx, npe*dim, PETScIdx, ElemMatrix, ADD_VALUES);CHKERRQ(ierr);

  }

  /* communication between processes */
  ierr = MatAssemblyBegin(*J, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(*J, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

  return 0;
}
/****************************************************************************************************/
int assembly_residual_sd(Vec *x_old, Vec *Residue)
{
  /*
     Assembly the Residual for Small Deformation
     approach.
   */

  int    i, k, e, gp, ngp, npe;
  int    PETScIdx[8*3];
  int    ierr;

  double ElemCoord[8][3];
  double ShapeDerivs[8][3];
  double DetJac;
  double ElemResidual[8*3];
  double B[6][3*8], stress_gp[6], strain_gp[6];
  double DsDe[6][6];
  double *wp = NULL;
  double ElemDispls[8*3];
  double *xvalues;

  Vec xlocal;
  material_t *material;

  register double wp_eff;

  /* 
     Local representation of <x> with ghost padding
  */
  ierr = VecZeroEntries(*Residue); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(*x_old,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(*x_old,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostGetLocalForm(*x_old,&xlocal); CHKERRQ(ierr);
  ierr = VecGetArray(xlocal, &xvalues); CHKERRQ(ierr); CHKERRQ(ierr);

  for(e=0;e<nelm;e++){

    material = GetMaterial(PhysicalID[e]);
    if(!material) SETERRQ1(PETSC_COMM_SELF,1,"material with physical_id %d not found",PhysicalID[e]);

    npe = eptr[e+1]-eptr[e];
    ngp = npe;
    GetPETScIndeces( &eind[eptr[e]], npe, loc2petsc, PETScIdx);
    GetElemCoord(&eind[eptr[e]], npe, ElemCoord);
    GetWeight(npe, &wp);
    GetElemenDispls(e, xvalues, ElemDispls);

    // calculate <ElemResidue> by numerical integration

    memset(ElemResidual, 0.0, (npe*dim)*sizeof(double));
    for(gp=0;gp<ngp;gp++){

      memset(strain_gp, 0.0, nvoi*sizeof(double));
      memset(stress_gp, 0.0, nvoi*sizeof(double));
      get_dsh(gp, npe, ElemCoord, ShapeDerivs, &DetJac);
      GetB( npe, ShapeDerivs, B );
      GetDsDe( e, ElemDispls, DsDe );

      for(i=0;i<nvoi;i++){
	for(k=0;k<npe*dim;k++){
	  strain_gp[i] += B[i][k]*ElemDispls[k];
	}
      }

      if(material->typeID==MICRO){
	ierr = mac_send_signal(WORLD_COMM, MAC2MIC_STRAIN);CHKERRQ(ierr);
	ierr = mac_send_strain(WORLD_COMM, strain_gp);CHKERRQ(ierr);
	ierr = mac_recv_stress(WORLD_COMM, stress_gp);CHKERRQ(ierr);
      }else{
	for(i=0;i<nvoi;i++){
	  for(k=0;k<nvoi;k++){
	    stress_gp[i] += DsDe[i][k]*strain_gp[k];
	  }
	}
      }

      wp_eff = DetJac * wp[gp];
      for(i=0;i<npe*dim;i++){
	for(k=0;k<nvoi;k++){
	  ElemResidual[i] += B[k][i]*stress_gp[k] * wp_eff;
	}
      }

    }
    ierr = VecSetValues(*Residue, npe*dim, PETScIdx, ElemResidual, ADD_VALUES);CHKERRQ(ierr);
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
  double ShapeDerivs[8][3];
  double det_jac;
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
      
      memset(strain_gp,0.0,nvoi*sizeof(double));
      memset(stress_gp,0.0,nvoi*sizeof(double));
      energy_gp = 0.0;

      get_dsh(gp, npe, ElemCoord, ShapeDerivs, &det_jac);
      GetB( npe, ShapeDerivs, B );
      GetDsDe( e, ElemDispls, DsDe );
      wp_eff = det_jac*wp[gp];

      for(i=0;i<nvoi;i++){
	for(k=0;k<npe*dim;k++){
	  strain_gp[i] += B[i][k]*ElemDispls[k];
	}
      }

      for(i=0;i<nvoi;i++){
	for(k=0;k<nvoi;k++){
	  stress_gp[i] += DsDe[i][k]*strain_gp[k];
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

      vol += det_jac*wp[gp];
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
  double ShapeDerivs[8][3];
  double DetJac;
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

      memset(strain_gp, 0.0, nvoi*sizeof(double));
      memset(stress_gp  , 0.0, nvoi*sizeof(double));

      get_dsh(gp, npe, ElemCoord, ShapeDerivs, &DetJac);
      GetB( npe, ShapeDerivs, B );
      GetDsDe( e, ElemDispls, DsDe );
      wp_eff = DetJac*wp[gp];

      for(i=0;i<nvoi;i++){
	for(k=0;k<npe*dim;k++){
	  strain_gp[i] += B[i][k]*ElemDispls[k];
	}
      }

      for(i=0;i<nvoi;i++){
	for(k=0;k<nvoi;k++){
	  stress_gp[i] += DsDe[i][k]*strain_gp[k];
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
int GetDsDe( int e, double *ElemDisp, double DsDe[6][6] )
{
  /*  
      Calculates constitutive tensor
      according to the element type
   */
  double la, mu, poi, you; 
  int i, j;

  material_t *material = GetMaterial(PhysicalID[e]);
  if(!material) SETERRQ1(PETSC_COMM_SELF,1,"material with physical_id %d not found",PhysicalID[e]);

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
	DsDe[0][0]=1.0; DsDe[0][1]=poi; DsDe[0][2]=0.0;
	DsDe[1][0]=poi; DsDe[1][1]=1.0; DsDe[1][2]=0.0;
	DsDe[2][0]=0.0; DsDe[2][1]=0.0; DsDe[2][2]=(1-poi)/2;
	for(i=0;i<3;i++){
	  for(j=0;j<3;j++){
	    DsDe[i][j] = DsDe[i][j] * you/(1-pow(poi,2));
	  }
	}
      }
      else if(dim==3){
	DsDe[0][0]=la+2*mu ;DsDe[0][1]=la      ;DsDe[0][2]=la      ;DsDe[0][3]=0.0; DsDe[0][4]=0.0; DsDe[0][5]=0.0;
	DsDe[1][0]=la      ;DsDe[1][1]=la+2*mu ;DsDe[1][2]=la      ;DsDe[1][3]=0.0; DsDe[1][4]=0.0; DsDe[1][5]=0.0;
	DsDe[2][0]=la      ;DsDe[2][1]=la      ;DsDe[2][2]=la+2*mu ;DsDe[2][3]=0.0; DsDe[2][4]=0.0; DsDe[2][5]=0.0;
	DsDe[3][0]=0.0     ;DsDe[3][1]=0.0     ;DsDe[3][2]=0.0     ;DsDe[3][3]=mu ; DsDe[3][4]=0.0; DsDe[3][5]=0.0;
	DsDe[4][0]=0.0     ;DsDe[4][1]=0.0     ;DsDe[4][2]=0.0     ;DsDe[4][3]=0.0; DsDe[4][4]=mu ; DsDe[4][5]=0.0;
	DsDe[5][0]=0.0     ;DsDe[5][1]=0.0     ;DsDe[5][2]=0.0     ;DsDe[5][3]=0.0; DsDe[5][4]=0.0; DsDe[5][5]=mu ;
      }

      break;

    case MICRO:

      if(macmic.type==COUP_1){
	for(i=0;i<6;i++){
	  for(j=0;j<6;j++){
	    DsDe[i][j]=((mac_coup_1_t*)macmic.coup)->homo_cij[i*6+j];
	  }
	}
      }

      break;

    default:
      break;

  }
  return 0;
}
/****************************************************************************************************/
material_t * GetMaterial(int GmshIDToSearch)
{
  /* 
     Search in the <material_list> if for the one
     that has <GmshIDToSearch>
   */

  node_list_t *pn;
  pn = material_list.head;
  while(pn)
  {
    if( ((material_t*)pn->data)->GmshID == GmshIDToSearch ) break;
    pn = pn->next;
  }
  if(!pn) return NULL;
  return (material_t*)pn->data;
}
/****************************************************************************************************/
int GetWeight(int npe, double **wp)
{
  *wp = FemGetPointer2Weight(npe, dim);
  return 0;
}
/****************************************************************************************************/
int GetB(int npe, double ShapeDerivs[8][3], double B[6][3*8] )
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
      B[0][i*2+0] = ShapeDerivs[i][0]; B[0][i*2+1] = 0.0         ;
      B[1][i*2+0] = 0.0              ; B[1][i*2+1] = ShapeDerivs[i][1];
      B[2][i*2+0] = ShapeDerivs[i][1]; B[2][i*2+1] = ShapeDerivs[i][0];
    }
  }
  else if(dim==3){
    for(i=0;i<npe;i++){
      B[0][i*3+0] = ShapeDerivs[i][0]; B[0][i*3+1] = 0.0               ; B[0][i*3+2] = 0.0              ; 
      B[1][i*3+0] = 0.0              ; B[1][i*3+1] = ShapeDerivs[i][1] ; B[1][i*3+2] = 0.0              ; 
      B[2][i*3+0] = 0.0              ; B[2][i*3+1] = 0.0               ; B[2][i*3+2] = ShapeDerivs[i][2]; 
      B[3][i*3+0] = ShapeDerivs[i][1]; B[3][i*3+1] = ShapeDerivs[i][0] ; B[3][i*3+2] = 0.0              ; 
      B[4][i*3+0] = 0.0              ; B[4][i*3+1] = ShapeDerivs[i][2] ; B[4][i*3+2] = ShapeDerivs[i][1];
      B[5][i*3+0] = ShapeDerivs[i][2]; B[5][i*3+1] = 0.0               ; B[5][i*3+2] = ShapeDerivs[i][0]; 

    }
  }

  return 0;
}
/****************************************************************************************************/
int get_dsh(int gp, int npe, double coor[8][3], double ShapeDerivs[8][3], double *DetJac)
{

  double ***ShapeDerivsMaster;

  double jac[3][3];
  double ijac[3][3];

  ShapeDerivsMaster = FemGetPointer2ShapeDerivsMaster(npe, dim);
  if(ShapeDerivsMaster == NULL) return 1;

  fem_calc_jac(dim, coor, ShapeDerivsMaster, npe, gp, jac);
  fem_invjac(dim, jac, ijac, DetJac);
  fem_trans_dsh(dim, ijac, npe, gp, ShapeDerivsMaster, ShapeDerivs);

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

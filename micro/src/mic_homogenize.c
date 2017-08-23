/*
   Routines for performing homogenization on RVE

   Author > Guido Giuntoli
   Date   > 18-08-2017

 */

#include "micro.h"

int micro_homogenize_taylor(MPI_Comm PROBLEM_COMM, double strain_mac[6], double strain_ave[6], double stress_ave[6])
{

  int    i, k, e, gp, ngp, npe, ierr;

  double coor_elm[8][3];
  double dsh[8][3];
  double detj;
  double stress_gp[6];
  double DsDe[6][6];
  double *wp = NULL;
  double disp_elm[8*3];
  double vol = -1.0, vol_tot = -1.0;
  double stress_aux[6], strain_aux[6];
  register double wp_eff;

  for(i=0;i<6;i++) strain_aux[i] = stress_aux[i] = 0.0;
 
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
int micro_homogenize_linear(MPI_Comm MICRO_COMM, int i, double strain_bc[6], double strain_ave[6], double stress_ave[6])
{

  /*
     Performs linear homogenization of the RVE cell acording to <strain_bc[6]>
   */

  int    ierr, nr_its = -1, kspits = -1;
  double norm = -1.0, NormTol = 1.0e-8, NRMaxIts = 3, kspnorm = -1.0;

  ierr = PetscLogEventBegin(EVENT_SET_DISP_BOU,0,0,0,0);CHKERRQ(ierr);
  ierr = PetscPrintf(MICRO_COMM,"\nExperiment i=%d\n",i);CHKERRQ(ierr);
  ierr = micro_apply_bc(i, strain_bc, &x, &A, &b, SET_DISPLACE);CHKERRQ(ierr);
  if( flag_print & (1<<PRINT_PETSC) ){
    ierr = PetscViewerASCIIOpen(MICRO_COMM,"x.dat",&viewer); CHKERRQ(ierr);
    ierr = VecView(x,viewer); CHKERRQ(ierr);
  }
  ierr = PetscLogEventEnd(EVENT_SET_DISP_BOU,0,0,0,0);CHKERRQ(ierr);

  nr_its = 0; norm = 2*NormTol;
  while( nr_its < NRMaxIts && norm > NormTol )
  {
    /*
       Assemblying Residual
     */
    ierr = PetscLogEventBegin(EVENT_ASSEMBLY_RES,0,0,0,0);CHKERRQ(ierr);
    ierr = PetscPrintf(MICRO_COMM,"Assembling Residual ");CHKERRQ(ierr);
    ierr = AssemblyResidualSmallDeformation( &x, &b);CHKERRQ(ierr);
    ierr = micro_apply_bc(i, strain_bc, &x, &A, &b, SET_RESIDUAL);
    if( flag_print & (1<<PRINT_PETSC) ){
      ierr = PetscViewerASCIIOpen(MICRO_COMM,"b.dat",&viewer); CHKERRQ(ierr);
      ierr = VecView(b,viewer); CHKERRQ(ierr);
    }
    ierr = VecNorm(b,NORM_2,&norm);CHKERRQ(ierr);
    ierr = PetscPrintf(MICRO_COMM,"|b| = %e\n",norm);CHKERRQ(ierr);
    ierr = VecScale(b,-1.0); CHKERRQ(ierr);
    ierr = PetscLogEventEnd(EVENT_ASSEMBLY_RES,0,0,0,0);CHKERRQ(ierr);
    if( !(norm > NormTol) )break;
    /*
       Assemblying Jacobian
     */
    ierr = PetscLogEventBegin(EVENT_ASSEMBLY_JAC,0,0,0,0);CHKERRQ(ierr);
    ierr = PetscPrintf(MICRO_COMM,"Assembling Jacobian\n");
    ierr = AssemblyJacobianSmallDeformation(&A);
    ierr = micro_apply_bc(i, strain_bc, &x, &A, &b, SET_JACOBIAN);CHKERRQ(ierr);
    if( flag_print & (1<<PRINT_PETSC) ){
      ierr = PetscViewerASCIIOpen(MICRO_COMM,"A.dat",&viewer); CHKERRQ(ierr);
      ierr = MatView(A,viewer); CHKERRQ(ierr);
    }
    ierr = PetscLogEventEnd(EVENT_ASSEMBLY_JAC,0,0,0,0);CHKERRQ(ierr);
    /*
       Solving Problem
     */
    ierr = PetscLogEventBegin(EVENT_SOLVE_SYSTEM,0,0,0,0);CHKERRQ(ierr);
    ierr = PetscPrintf(MICRO_COMM,"Solving Linear System ");
    ierr = KSPSolve(ksp,b,dx);CHKERRQ(ierr);
    ierr = KSPGetIterationNumber(ksp,&kspits);CHKERRQ(ierr);
    ierr = KSPGetConvergedReason(ksp,&reason);CHKERRQ(ierr);
    ierr = KSPGetResidualNorm(ksp,&kspnorm);CHKERRQ(ierr);
    ierr = VecAXPY( x, 1.0, dx); CHKERRQ(ierr);
    if( flag_print & (1<<PRINT_PETSC) ){
      ierr = PetscViewerASCIIOpen(MICRO_COMM,"dx.dat",&viewer); CHKERRQ(ierr);
      ierr = VecView(dx,viewer); CHKERRQ(ierr);
      ierr = PetscViewerASCIIOpen(MICRO_COMM,"x.dat",&viewer); CHKERRQ(ierr);
      ierr = VecView(x,viewer); CHKERRQ(ierr);
    }
    ierr = PetscPrintf(MICRO_COMM,"Iterations %D Norm %e reason %d\n",kspits, kspnorm, reason);CHKERRQ(ierr);
    ierr = PetscLogEventBegin(EVENT_SOLVE_SYSTEM,0,0,0,0);CHKERRQ(ierr);

    nr_its ++;
  }
  ierr = PetscPrintf(MICRO_COMM,"NR its : %d\n",nr_its);CHKERRQ(ierr);

  /*
     Calculate the average stress, strain and constitutive tensor
   */
  ierr = SpuAveStressAndStrain(MICRO_COMM, &x, strain_ave, stress_ave);CHKERRQ(ierr);

  return 0;
}
/****************************************************************************************************/
int micro_homogenize(MPI_Comm MICRO_COMM, double strain_mac[6], double strain_ave[6], double stress_ave[6])
{

  /*
     Performs linear homogenization of the RVE 
     Possible types are>

     HOMO_TAYLOR  1
     HOMO_LINEAR  2
     HOMO_PERIOD  3
     HOMO_EXP     9

     HOMO_TAYLOR > no need of calculating a displacement field
     
   */
  int    ierr, nr_its = -1, kspits = -1;
  double norm = -1.0, NormTol = 1.0e-8, NRMaxIts = 3, kspnorm = -1.0;

  if(homo_type==HOMO_TAYLOR){
    ierr = micro_homogenize_taylor(MICRO_COMM, strain_mac, strain_ave, stress_ave);CHKERRQ(ierr);
  }
  else{

    ierr = PetscLogEventBegin(EVENT_SET_DISP_BOU,0,0,0,0);CHKERRQ(ierr);
    /*
       depending on the <kind> we apply a specific boundary condition
     */
    if(homo_type==HOMO_EXP){
      ierr = micro_apply_bc_linear(strain_mac, &x, &A, &b, SET_DISPLACE);CHKERRQ(ierr);
    }
    else if(homo_type==HOMO_PERIOD){
      SETERRQ1(MICRO_COMM,1,"<kind> %d not supported",homo_type);
    }
    else{
      SETERRQ1(MICRO_COMM,1,"<kind> %d not supported",homo_type);
    }

    if( flag_print & (1<<PRINT_PETSC) ){
      ierr = PetscViewerASCIIOpen(MICRO_COMM,"x.dat",&viewer); CHKERRQ(ierr);
      ierr = VecView(x,viewer); CHKERRQ(ierr);
    }
    ierr = PetscLogEventEnd(EVENT_SET_DISP_BOU,0,0,0,0);CHKERRQ(ierr);

    nr_its = 0; norm = 2*NormTol;
    while( nr_its < NRMaxIts && norm > NormTol )
    {
      /*
	 Assemblying Residual
       */
      ierr = PetscLogEventBegin(EVENT_ASSEMBLY_RES,0,0,0,0);CHKERRQ(ierr);
      if(flag_coupling)
	ierr = PetscPrintf(MICRO_COMM,"Assembling Residual ");CHKERRQ(ierr);
      ierr = AssemblyResidualSmallDeformation( &x, &b);CHKERRQ(ierr);
      
      if(homo_type==HOMO_EXP){
	ierr = micro_apply_bc_linear(strain_mac, &x, &A, &b, SET_RESIDUAL);CHKERRQ(ierr);
      }
      else if(homo_type==HOMO_PERIOD){
	SETERRQ1(MICRO_COMM,1,"<kind> %d not supported",homo_type);
      }
      else{
	SETERRQ1(MICRO_COMM,1,"<kind> %d not supported",homo_type);
      }

      if( flag_print & (1<<PRINT_PETSC) ){
	ierr = PetscViewerASCIIOpen(MICRO_COMM,"b.dat",&viewer); CHKERRQ(ierr);
	ierr = VecView(b,viewer); CHKERRQ(ierr);
      }
      ierr = VecNorm(b,NORM_2,&norm);CHKERRQ(ierr);
      ierr = PetscPrintf(MICRO_COMM,"|b| = %e\n",norm);CHKERRQ(ierr);
      ierr = VecScale(b,-1.0); CHKERRQ(ierr);
      ierr = PetscLogEventEnd(EVENT_ASSEMBLY_RES,0,0,0,0);CHKERRQ(ierr);
      if( !(norm > NormTol) )break;
      /*
	 Assemblying Jacobian
       */
      ierr = PetscLogEventBegin(EVENT_ASSEMBLY_JAC,0,0,0,0);CHKERRQ(ierr);
      if(flag_coupling)
	ierr = PetscPrintf(MICRO_COMM,"Assembling Jacobian\n");
      ierr = AssemblyJacobianSmallDeformation(&A);

      if(homo_type==HOMO_EXP){
	ierr = micro_apply_bc_linear(strain_mac, &x, &A, &b, SET_JACOBIAN);CHKERRQ(ierr);
      }
      else if(homo_type==HOMO_PERIOD){
	SETERRQ1(MICRO_COMM,1,"<kind> %d not supported",homo_type);
      }
      else{
	SETERRQ1(MICRO_COMM,1,"<kind> %d not supported",homo_type);
      }

      if( flag_print & (1<<PRINT_PETSC) ){
	ierr = PetscViewerASCIIOpen(MICRO_COMM,"A.dat",&viewer); CHKERRQ(ierr);
	ierr = MatView(A,viewer); CHKERRQ(ierr);
      }
      ierr = PetscLogEventEnd(EVENT_ASSEMBLY_JAC,0,0,0,0);CHKERRQ(ierr);
      /*
	 Solving Problem
       */
      ierr = PetscLogEventBegin(EVENT_SOLVE_SYSTEM,0,0,0,0);CHKERRQ(ierr);
      if(flag_coupling)
	ierr = PetscPrintf(MICRO_COMM,"Solving Linear System ");
      ierr = KSPSolve(ksp,b,dx);CHKERRQ(ierr);
      ierr = KSPGetIterationNumber(ksp,&kspits);CHKERRQ(ierr);
      ierr = KSPGetConvergedReason(ksp,&reason);CHKERRQ(ierr);
      ierr = KSPGetResidualNorm(ksp,&kspnorm);CHKERRQ(ierr);
      ierr = VecAXPY( x, 1.0, dx); CHKERRQ(ierr);
      if( flag_print & (1<<PRINT_PETSC) ){
	ierr = PetscViewerASCIIOpen(MICRO_COMM,"dx.dat",&viewer); CHKERRQ(ierr);
	ierr = VecView(dx,viewer); CHKERRQ(ierr);
	ierr = PetscViewerASCIIOpen(MICRO_COMM,"x.dat",&viewer); CHKERRQ(ierr);
	ierr = VecView(x,viewer); CHKERRQ(ierr);
      }
      if(flag_coupling)
	ierr = PetscPrintf(MICRO_COMM,"Iterations %D Norm %e reason %d\n",kspits, kspnorm, reason);CHKERRQ(ierr);
      ierr = PetscLogEventBegin(EVENT_SOLVE_SYSTEM,0,0,0,0);CHKERRQ(ierr);

      nr_its ++;
    }
    ierr = PetscPrintf(MICRO_COMM,"NR its : %d\n",nr_its);CHKERRQ(ierr);

    /*
       Calculate the average stress, strain and constitutive tensor
     */
    ierr = SpuAveStressAndStrain(MICRO_COMM, &x, strain_ave, stress_ave);CHKERRQ(ierr);
  }

  return 0;
}
/****************************************************************************************************/

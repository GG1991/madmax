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
  int    ierr, nr_its=-1, max_its=3;
  double norm_tol=1.0e-8, norm=2*norm_tol;
//  double ksp_norm=-1.0;

  ierr = micro_apply_bc_linear(strain_mac, &x, &A, &b, SET_DISPLACE);CHKERRQ(ierr);

  while( nr_its < max_its && norm > norm_tol )
  {
    ierr = assembly_residual_sd( &x, &b);CHKERRQ(ierr);
    ierr = micro_apply_bc_linear(strain_mac, &x, &A, &b, SET_RESIDUAL);CHKERRQ(ierr);
    ierr = VecNorm(b,NORM_2,&norm);CHKERRQ(ierr);
    if( !(norm > norm_tol) )break;
    ierr = VecScale(b,-1.0); CHKERRQ(ierr);
    ierr = assembly_jacobian_sd(&A);
    ierr = micro_apply_bc_linear(strain_mac, &x, &A, &b, SET_JACOBIAN);CHKERRQ(ierr);
    ierr = KSPSolve(ksp,b,dx);CHKERRQ(ierr);
    ierr = VecAXPY( x, 1.0, dx); CHKERRQ(ierr);
    nr_its ++;
  }

  ierr = calc_ave_strain_stress(MICRO_COMM, &x, strain_ave, stress_ave);CHKERRQ(ierr);

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
  int ierr;

  if(homo.type==HOMO_TAYLOR){
    ierr = micro_homogenize_taylor(MICRO_COMM, strain_mac, strain_ave, stress_ave);CHKERRQ(ierr);
  }
  else if(homo.type==HOMO_LINEAR){
    ierr = micro_homogenize_linear(MICRO_COMM, strain_mac, strain_ave, stress_ave);CHKERRQ(ierr);
  }
  else if(homo.type==HOMO_LINEAR_HEXA){
    ierr = micro_homogenize_linear_hexa(MICRO_COMM, strain_mac, strain_ave, stress_ave);CHKERRQ(ierr);
  }
  else{
    return 1;
  }

  return 0;
}
/****************************************************************************************************/
int mic_init_homo(void)
{
  
  if(homo.type==LD_LAGRAN){
    /*
       a) Count how many nodes (nnods_bc) belongs 
       to the boundary search in the boundary_list
       b) Allocate <index> and <ub_val>
    */

    homo.st = malloc(sizeof(homog_ld_lagran_t));

    int nnods_bc=0;
    node_list_t *pb;
    pb = boundary_list.head;
    while(pb){
      nnods_bc += ((boundary_t*)pb->data)->Nods.sizelist;
      pb=pb->next;
    }
    ((homog_ld_lagran_t*)homo.st)->nnods_bc = nnods_bc;


  }

  return 0;
}

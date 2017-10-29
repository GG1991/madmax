/*
   Routines for performing homogenization on RVE

   Author > Guido Giuntoli
   Date   > 18-08-2017

 */

#include "micro.h"

int get_local_index( int e, int *loc_index );
int assembly_residual_struct( void );

int mic_homogenize_taylor(MPI_Comm MICRO_COMM, double strain_mac[6], double strain_ave[6], double stress_ave[6])
{
  /* 
     This work only if inclusion names starts with FIBER and matrix name start with MATRIX
     we assume only on kind of them 
    
     This kind of homogenization do not need to perform 3 (dim=2) or 6 (dim=3) experiments
  */
  int          i, j, e, ierr;
  double       c_i[6][6];
  double       c_m[6][6];
  double       c[6][6];
  double       vol_i = 0.0, vol_ia = 0.0;  // inclusion volume
  double       vol_m = 0.0, vol_ma = 0.0;  // matrix volume
  double       vol_t = 0.0;  // total volume
  double       vol_e = 0.0;  // element volume
  material_t   *mat = NULL; 

  if(first_time_homo){
    for(e=0;e<nelm;e++){

      ierr = get_mat_from_elem(e, &mat);
      ierr = get_elem_vol(e, &vol_e);
      if(!strncmp(mat->name,"FIBER",5)){
	vol_ia += vol_e;
      }
      else if(!strncmp(mat->name,"MATRIX",6)){
	vol_ma += vol_e;
      }
      else{
	PetscPrintf(MICRO_COMM,
	"not possible to compute a Taylor homogenization with a material of type\n",mat->name);
	return 1;
      }

    }
    ierr = MPI_Allreduce(&vol_ia, &vol_i, 1, MPI_DOUBLE, MPI_SUM, MICRO_COMM);CHKERRQ(ierr);
    ierr = MPI_Allreduce(&vol_ma, &vol_m, 1, MPI_DOUBLE, MPI_SUM, MICRO_COMM);CHKERRQ(ierr);
    vol_t = vol_i + vol_m;
    vi = vol_i / vol_t; // inclusion fraction
    vm = vol_m / vol_t; // matrix fraction
  }
  get_c("FIBER" , 0, 0, strain_mac, c_i);  //returns c_i of FIBER
  get_c("MATRIX", 0, 0, strain_mac, c_m);  //returns c_m of MATRIX
 
  if(homo_type==TAYLOR_P){
    for(i=0;i<nvoi;i++){
      for(j=0;j<nvoi;j++){
	c[i][j] = vi * c_i[i][j] + vm * c_m[i][j];
      }
    }
  }
  else if(homo_type==TAYLOR_S){

    int              s;
    double           c_ia[36], c_ma[36]; // matrices in array
    double           c_ii[36], c_mi[36]; // inverted matrices
    double           c_a[36] , c_ai[36];
    gsl_matrix_view  m, mi;
    gsl_permutation  *p;

    p  = gsl_permutation_alloc(nvoi);

    for(i = 0; i < nvoi; ++i){
      for(j = 0; j < nvoi; ++j){
	c_ia[i*nvoi + j] = c_i[i][j];
	c_ma[i*nvoi + j] = c_m[i][j];
      }
    }
    
    m  = gsl_matrix_view_array(c_ia,nvoi,nvoi);
    mi = gsl_matrix_view_array(c_ii,nvoi,nvoi);
    
    gsl_linalg_LU_decomp (&m.matrix, p, &s);    
    gsl_linalg_LU_invert (&m.matrix, p, &mi.matrix);
    
    for(i = 0; i < nvoi; ++i){
      for(j = 0; j < nvoi; ++j){
	c_ii[i*nvoi+j] = gsl_matrix_get(&mi.matrix,i,j);
      }
    }

    m  = gsl_matrix_view_array(c_ma,nvoi,nvoi);
    mi = gsl_matrix_view_array(c_mi,nvoi,nvoi);
    
    gsl_linalg_LU_decomp (&m.matrix, p, &s);    
    gsl_linalg_LU_invert (&m.matrix, p, &mi.matrix);
    
    for(i = 0; i < nvoi; ++i){
      for(j = 0; j < nvoi; ++j){
	c_mi[i*nvoi+j] = gsl_matrix_get(&mi.matrix,i,j);
      }
    }
    for(i = 0; i < nvoi*nvoi; ++i){
      c_a[i] = vi * c_ii[i] + vm * c_mi[i];
    }

    m  = gsl_matrix_view_array(c_a ,nvoi,nvoi);
    mi = gsl_matrix_view_array(c_ai,nvoi,nvoi);
    
    gsl_linalg_LU_decomp (&m.matrix, p, &s);    
    gsl_linalg_LU_invert (&m.matrix, p, &mi.matrix);
    
    for(i = 0; i < nvoi; ++i){
      for(j = 0; j < nvoi; ++j){
	c[i][j] = gsl_matrix_get(&mi.matrix,i,j);
      }
    }
     
    gsl_permutation_free (p);
  }
  for(i=0;i<nvoi;i++){
    strain_ave[i] = strain_mac[i];
    stress_ave[i] = 0.0;
    for(j=0;j<nvoi;j++){
      stress_ave[i] += c[i][j] * strain_mac[j];
    }
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
    strain_matrix[0][0]=strain_mac[0]    ; strain_matrix[0][1]=strain_mac[2]*0.5;
    strain_matrix[1][0]=strain_mac[2]*0.5; strain_matrix[1][1]=strain_mac[1]    ;
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
    if(!flag_coupling)
      PetscPrintf(MICRO_COMM,"|b| = %lf \n",norm);

    if( !(norm > nr_norm_tol) )break;

    if(first_time_homo || !flag_linear_micro){
      ierr = assembly_jacobian_sd(&A);
    }
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
int mic_homog_us(MPI_Comm MICRO_COMM, double strain_mac[6], double strain_ave[6], double stress_ave[6])
{

  /* 
     homogenization routine with uniform strain BC
     for structured mesh 

     a) alloc matrix/vector if it is first time
     b) assembly residue -> evaluate norm
     c) assembly jacobian
     d) solve system 
     e) iterate to b) with NR method
     f) get average properties

     partimos el dominio en tandas por el eje y

   */
  if( flag_first_alloc == true ){

    flag_first_alloc = false;

    int rank, nproc;

    MPI_Comm_size(MICRO_COMM, &nproc);
    MPI_Comm_rank(MICRO_COMM, &rank);

    //int nyl = ny / nproc; // ny for every process
    int nyl = ny / nproc + ((ny % nproc > rank)?1:0); // ny for every process
    int nl  = ( dim == 2 ) ? nyl*nx : nyl*nx*nz;      // local nodes
    int nnz = (dim==2)? 18:81;                        // nonzeros per row

    MatCreate(MICRO_COMM,&A);
    MatSetSizes(A, nl*dim, nl*dim, nn*dim, nn*dim);
    MatSetFromOptions(A);
    MatSeqAIJSetPreallocation(A, nnz, NULL);
    MatMPIAIJSetPreallocation(A, nnz, NULL, nnz, NULL);
    MatSetOption(A, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);

    int istart, iend;
    MatGetOwnershipRange(A,&istart,&iend);

    int nghost, *ghost_index;;
    if( nproc == 1 )
      nghost = 0;
    else
      nghost = ( dim == 2 )?nx:nx*nz;

    ghost_index = malloc(nghost*dim*sizeof(int));

    int i;
    if( rank == 0 ){
      for( i = 0 ; i < nghost*dim  ; i++ )
	ghost_index[i] = istart + (nyl-1)*nghost + i;
    }
    else if( rank == (nproc - 1) ){
      for( i = 0 ; i < nghost*dim  ; i++ )
	ghost_index[i] = istart + i;
    }
    else{
      for( i = 0 ; i < nghost*dim  ; i++ ){
	ghost_index[i] = istart + i;
	ghost_index[i] = istart + (nyl-1)*nghost + i;
      }
    }

    VecCreateGhost(MICRO_COMM, nl*dim, nn*dim, nghost*dim, ghost_index, &x);
    VecDuplicate(x,&dx);
    VecDuplicate(x,&b);

    free(ghost_index);

    /* Initilize shape functions, derivatives, jacobian, b_matrix */

    int nsh = ( dim == 2 ) ? 4 : 8;
    double h[3]; h[0] = hx; h[1] = hx; h[3] = hz;

    fem_init_struct( &struct_sh, &struct_dsh, &struct_wp, h, dim);
    
    /* alloc the B matrix */
    struct_bmat = malloc( nvoi*nsh*dim* sizeof(double));
    /* alloc the local index vector for assembly */
    loc_index   = malloc( dim*npe * sizeof(int));

  } // first time for allocation

  /* 
     Begin Newton-Raphson Iterations 
   */
  int nr_its = 0, norm = nr_norm_tol*2;
  while( nr_its < nr_max_its && norm > nr_norm_tol )
  {
    assembly_residual_struct(); // assembly "b" (residue) using "x" (displacement)
  }


  return 0;
}
/****************************************************************************************************/
int assembly_residual_struct(void)
{

  VecZeroEntries(b);
  VecGhostUpdateBegin(x,INSERT_VALUES,SCATTER_FORWARD);
  VecGhostUpdateEnd(x,INSERT_VALUES,SCATTER_FORWARD);

  Vec     xlocal;
  double *xvalues;

  VecGhostGetLocalForm(x,&xlocal);
  VecGetArray(xlocal, &xvalues);

  int e;
  double *elem_disp = malloc( dim*npe * sizeof(double));

  for( e = 0 ; e < nelm ; e++ ){

    /* get the local indeces of the element vertex nodes */
    get_local_index(e, loc_index);

    /* get the elemental displacements */
    int i;
    for( i = 0 ; i < npe*dim ; i++ )
      elem_disp[i] = xvalues[loc_index[i]];

  }

  return 0;
}
/****************************************************************************************************/
int get_local_index( int e, int *loc_index )
{
  
  int d;
  if( dim == 2 ){
    for( d = 0 ; d < dim ; d++ ){
      loc_index[ 0*dim + d ] = (e + 0)      * dim + d ;
      loc_index[ 1*dim + d ] = (e + 1)      * dim + d ;
      loc_index[ 2*dim + d ] = (e + nx + 0) * dim + d ;
      loc_index[ 3*dim + d ] = (e + nx + 1) * dim + d ;
    }
  }
  return 0;
}
/****************************************************************************************************/
int mic_homogenize(MPI_Comm MICRO_COMM, double strain_mac[6], double strain_ave[6], double stress_ave[6])
{

  /*
     Performs linear homogenization of the RVE 

     UNIF_STRAINS            > u = E . x

     HOMO_TAYLOR_S           > 
   */
  int ierr;

  if(homo_type==TAYLOR_P || homo_type==TAYLOR_S){

    ierr = mic_homogenize_taylor(MICRO_COMM, strain_mac, strain_ave, stress_ave);
    if(ierr){
      return 1;
    }
  }
  else if(homo_type==UNIF_STRAINS){

    ierr = mic_homog_us(MICRO_COMM, strain_mac, strain_ave, stress_ave);
    if(ierr) return 1;

  }
  if(first_time_homo) first_time_homo = 0;

  return 0;
}
/****************************************************************************************************/
int mic_calc_c_homo(MPI_Comm MICRO_COMM, double strain_mac[6], double c_homo[36])
{

  /* 
     Si la micro estructura está integramente conformada por materiales
     lineales entonces este tensor será siempre el mismo para cada punto 
     de gauss en la macro escala entonces es eficiente almacenar c_homo_linear
  */

  int i, ierr;

  if(flag_linear_micro){

    if(first_time_homo){
      ierr = mic_calc_c_homo_lineal(MICRO_COMM, c_homo_lineal);
      if(ierr){
	return 1;
      }
    }
    for(i=0;i<nvoi*nvoi;i++){
      c_homo[i] = c_homo_lineal[i];
    }

  }
  else{
    return 1;
  }
  return 0;
}
/****************************************************************************************************/
int mic_calc_c_homo_lineal(MPI_Comm MICRO_COMM, double c_homo_lineal[36])
{

  int       i, j;
  int       ierr;
  double    strain[6];
  double    strain_ave[6];
  double    stress_ave[6];

  for(i=0;i<nvoi*nvoi;i++){
    c_homo_lineal[i]=0.0;
  }

  for(i=0;i<nvoi;i++){

    for(j=0;j<nvoi;j++){
      strain[j]=0.0;
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
/****************************************************************************************************/
int mic_calc_stress_ave(MPI_Comm MICRO_COMM, double strain_mac[6], double strain_ave[6], double stress_ave[6])
{

  /* 
     Si la micro estructura está integramente conformada por materiales
     lineales entonces este tensor será siempre el mismo para cada punto 
     de gauss en la macro escala entonces es eficiente almacenar c_homo_linear
  */

  int i, j, ierr;

  if(flag_linear_micro){

    for(i=0;i<nvoi;i++){
      strain_ave[i] = strain_mac[i];
      stress_ave[i] = 0.0;
      for(j=0;j<nvoi;j++){
	stress_ave[i] += c_homo_lineal[i*nvoi+j] * strain_mac[j];
      }
    }
    return 0;

  }
  else{

    ierr = mic_homogenize(MICRO_COMM, strain_mac, strain_ave, stress_ave);
    if(ierr){
      return 1;
    }
    return 0;

  }
  return 0;
}
/****************************************************************************************************/
int mic_check_linear_material(void)
{

  flag_linear_micro = 1;

  return 0;
}

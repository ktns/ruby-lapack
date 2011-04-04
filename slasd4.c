#include "rb_lapack.h"

extern VOID slasd4_(integer *n, integer *i, real *d, real *z, real *delta, real *rho, real *sigma, real *work, integer *info);

static VALUE
rb_slasd4(int argc, VALUE *argv, VALUE self){
  VALUE rb_i;
  integer i; 
  VALUE rb_d;
  real *d; 
  VALUE rb_z;
  real *z; 
  VALUE rb_rho;
  real rho; 
  VALUE rb_delta;
  real *delta; 
  VALUE rb_sigma;
  real sigma; 
  VALUE rb_info;
  integer info; 
  real *work;

  integer n;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  delta, sigma, info = NumRu::Lapack.slasd4( i, d, z, rho)\n    or\n  NumRu::Lapack.slasd4  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE SLASD4( N, I, D, Z, DELTA, RHO, SIGMA, WORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  This subroutine computes the square root of the I-th updated\n*  eigenvalue of a positive symmetric rank-one modification to\n*  a positive diagonal matrix whose entries are given as the squares\n*  of the corresponding entries in the array d, and that\n*\n*         0 <= D(i) < D(j)  for  i < j\n*\n*  and that RHO > 0. This is arranged by the calling routine, and is\n*  no loss in generality.  The rank-one modified system is thus\n*\n*         diag( D ) * diag( D ) +  RHO *  Z * Z_transpose.\n*\n*  where we assume the Euclidean norm of Z is 1.\n*\n*  The method consists of approximating the rational functions in the\n*  secular equation by simpler interpolating rational functions.\n*\n\n*  Arguments\n*  =========\n*\n*  N      (input) INTEGER\n*         The length of all arrays.\n*\n*  I      (input) INTEGER\n*         The index of the eigenvalue to be computed.  1 <= I <= N.\n*\n*  D      (input) REAL array, dimension ( N )\n*         The original eigenvalues.  It is assumed that they are in\n*         order, 0 <= D(I) < D(J)  for I < J.\n*\n*  Z      (input) REAL array, dimension (N)\n*         The components of the updating vector.\n*\n*  DELTA  (output) REAL array, dimension (N)\n*         If N .ne. 1, DELTA contains (D(j) - sigma_I) in its  j-th\n*         component.  If N = 1, then DELTA(1) = 1.  The vector DELTA\n*         contains the information necessary to construct the\n*         (singular) eigenvectors.\n*\n*  RHO    (input) REAL\n*         The scalar in the symmetric updating formula.\n*\n*  SIGMA  (output) REAL\n*         The computed sigma_I, the I-th updated eigenvalue.\n*\n*  WORK   (workspace) REAL array, dimension (N)\n*         If N .ne. 1, WORK contains (D(j) + sigma_I) in its  j-th\n*         component.  If N = 1, then WORK( 1 ) = 1.\n*\n*  INFO   (output) INTEGER\n*         = 0:  successful exit\n*         > 0:  if INFO = 1, the updating process failed.\n*\n*  Internal Parameters\n*  ===================\n*\n*  Logical variable ORGATI (origin-at-i?) is used for distinguishing\n*  whether D(i) or D(i+1) is treated as the origin.\n*\n*            ORGATI = .true.    origin at i\n*            ORGATI = .false.   origin at i+1\n*\n*  Logical variable SWTCH3 (switch-for-3-poles?) is for noting\n*  if we are working with THREE poles!\n*\n*  MAXIT is the maximum number of iterations allowed for each\n*  eigenvalue.\n*\n\n*  Further Details\n*  ===============\n*\n*  Based on contributions by\n*     Ren-Cang Li, Computer Science Division, University of California\n*     at Berkeley, USA\n*\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 4)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 4)", argc);
  rb_i = argv[0];
  rb_d = argv[1];
  rb_z = argv[2];
  rb_rho = argv[3];

  rho = (real)NUM2DBL(rb_rho);
  if (!NA_IsNArray(rb_z))
    rb_raise(rb_eArgError, "z (3th argument) must be NArray");
  if (NA_RANK(rb_z) != 1)
    rb_raise(rb_eArgError, "rank of z (3th argument) must be %d", 1);
  n = NA_SHAPE0(rb_z);
  if (NA_TYPE(rb_z) != NA_SFLOAT)
    rb_z = na_change_type(rb_z, NA_SFLOAT);
  z = NA_PTR_TYPE(rb_z, real*);
  if (!NA_IsNArray(rb_d))
    rb_raise(rb_eArgError, "d (2th argument) must be NArray");
  if (NA_RANK(rb_d) != 1)
    rb_raise(rb_eArgError, "rank of d (2th argument) must be %d", 1);
  if (NA_SHAPE0(rb_d) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of d must be the same as shape 0 of z");
  if (NA_TYPE(rb_d) != NA_SFLOAT)
    rb_d = na_change_type(rb_d, NA_SFLOAT);
  d = NA_PTR_TYPE(rb_d, real*);
  i = NUM2INT(rb_i);
  {
    int shape[1];
    shape[0] = n;
    rb_delta = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  delta = NA_PTR_TYPE(rb_delta, real*);
  work = ALLOC_N(real, (n));

  slasd4_(&n, &i, d, z, delta, &rho, &sigma, work, &info);

  free(work);
  rb_sigma = rb_float_new((double)sigma);
  rb_info = INT2NUM(info);
  return rb_ary_new3(3, rb_delta, rb_sigma, rb_info);
}

void
init_lapack_slasd4(VALUE mLapack){
  rb_define_module_function(mLapack, "slasd4", rb_slasd4, -1);
}

#include "rb_lapack.h"

extern VOID dlaed4_(integer *n, integer *i, doublereal *d, doublereal *z, doublereal *delta, doublereal *rho, doublereal *dlam, integer *info);

static VALUE sHelp, sUsage;

static VALUE
rblapack_dlaed4(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_i;
  integer i; 
  VALUE rblapack_d;
  doublereal *d; 
  VALUE rblapack_z;
  doublereal *z; 
  VALUE rblapack_rho;
  doublereal rho; 
  VALUE rblapack_delta;
  doublereal *delta; 
  VALUE rblapack_dlam;
  doublereal dlam; 
  VALUE rblapack_info;
  integer info; 

  integer n;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  delta, dlam, info = NumRu::Lapack.dlaed4( i, d, z, rho, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE DLAED4( N, I, D, Z, DELTA, RHO, DLAM, INFO )\n\n*  Purpose\n*  =======\n*\n*  This subroutine computes the I-th updated eigenvalue of a symmetric\n*  rank-one modification to a diagonal matrix whose elements are\n*  given in the array d, and that\n*\n*             D(i) < D(j)  for  i < j\n*\n*  and that RHO > 0.  This is arranged by the calling routine, and is\n*  no loss in generality.  The rank-one modified system is thus\n*\n*             diag( D )  +  RHO *  Z * Z_transpose.\n*\n*  where we assume the Euclidean norm of Z is 1.\n*\n*  The method consists of approximating the rational functions in the\n*  secular equation by simpler interpolating rational functions.\n*\n\n*  Arguments\n*  =========\n*\n*  N      (input) INTEGER\n*         The length of all arrays.\n*\n*  I      (input) INTEGER\n*         The index of the eigenvalue to be computed.  1 <= I <= N.\n*\n*  D      (input) DOUBLE PRECISION array, dimension (N)\n*         The original eigenvalues.  It is assumed that they are in\n*         order, D(I) < D(J)  for I < J.\n*\n*  Z      (input) DOUBLE PRECISION array, dimension (N)\n*         The components of the updating vector.\n*\n*  DELTA  (output) DOUBLE PRECISION array, dimension (N)\n*         If N .GT. 2, DELTA contains (D(j) - lambda_I) in its  j-th\n*         component.  If N = 1, then DELTA(1) = 1. If N = 2, see DLAED5\n*         for detail. The vector DELTA contains the information necessary\n*         to construct the eigenvectors by DLAED3 and DLAED9.\n*\n*  RHO    (input) DOUBLE PRECISION\n*         The scalar in the symmetric updating formula.\n*\n*  DLAM   (output) DOUBLE PRECISION\n*         The computed lambda_I, the I-th updated eigenvalue.\n*\n*  INFO   (output) INTEGER\n*         = 0:  successful exit\n*         > 0:  if INFO = 1, the updating process failed.\n*\n*  Internal Parameters\n*  ===================\n*\n*  Logical variable ORGATI (origin-at-i?) is used for distinguishing\n*  whether D(i) or D(i+1) is treated as the origin.\n*\n*            ORGATI = .true.    origin at i\n*            ORGATI = .false.   origin at i+1\n*\n*   Logical variable SWTCH3 (switch-for-3-poles?) is for noting\n*   if we are working with THREE poles!\n*\n*   MAXIT is the maximum number of iterations allowed for each\n*   eigenvalue.\n*\n\n*  Further Details\n*  ===============\n*\n*  Based on contributions by\n*     Ren-Cang Li, Computer Science Division, University of California\n*     at Berkeley, USA\n*\n*  =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  delta, dlam, info = NumRu::Lapack.dlaed4( i, d, z, rho, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 4)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 4)", argc);
  rblapack_i = argv[0];
  rblapack_d = argv[1];
  rblapack_z = argv[2];
  rblapack_rho = argv[3];
  if (rb_options != Qnil) {
  }

  rho = NUM2DBL(rblapack_rho);
  if (!NA_IsNArray(rblapack_z))
    rb_raise(rb_eArgError, "z (3th argument) must be NArray");
  if (NA_RANK(rblapack_z) != 1)
    rb_raise(rb_eArgError, "rank of z (3th argument) must be %d", 1);
  n = NA_SHAPE0(rblapack_z);
  if (NA_TYPE(rblapack_z) != NA_DFLOAT)
    rblapack_z = na_change_type(rblapack_z, NA_DFLOAT);
  z = NA_PTR_TYPE(rblapack_z, doublereal*);
  if (!NA_IsNArray(rblapack_d))
    rb_raise(rb_eArgError, "d (2th argument) must be NArray");
  if (NA_RANK(rblapack_d) != 1)
    rb_raise(rb_eArgError, "rank of d (2th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_d) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of d must be the same as shape 0 of z");
  if (NA_TYPE(rblapack_d) != NA_DFLOAT)
    rblapack_d = na_change_type(rblapack_d, NA_DFLOAT);
  d = NA_PTR_TYPE(rblapack_d, doublereal*);
  i = NUM2INT(rblapack_i);
  {
    int shape[1];
    shape[0] = n;
    rblapack_delta = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  delta = NA_PTR_TYPE(rblapack_delta, doublereal*);

  dlaed4_(&n, &i, d, z, delta, &rho, &dlam, &info);

  rblapack_dlam = rb_float_new((double)dlam);
  rblapack_info = INT2NUM(info);
  return rb_ary_new3(3, rblapack_delta, rblapack_dlam, rblapack_info);
}

void
init_lapack_dlaed4(VALUE mLapack){
  rb_define_module_function(mLapack, "dlaed4", rblapack_dlaed4, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}

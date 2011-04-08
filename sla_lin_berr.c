#include "rb_lapack.h"

extern VOID sla_lin_berr_(integer *n, integer *nz, integer *nrhs, real *res, real *ayb, real *berr);

static VALUE sHelp, sUsage;

static VALUE
rblapack_sla_lin_berr(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_nz;
  integer nz; 
  VALUE rblapack_res;
  real *res; 
  VALUE rblapack_ayb;
  real *ayb; 
  VALUE rblapack_berr;
  real *berr; 

  integer n;
  integer nrhs;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  berr = NumRu::Lapack.sla_lin_berr( nz, res, ayb, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE SLA_LIN_BERR ( N, NZ, NRHS, RES, AYB, BERR )\n\n*  Purpose\n*  =======\n*\n*     SLA_LIN_BERR computes componentwise relative backward error from\n*     the formula\n*         max(i) ( abs(R(i)) / ( abs(op(A_s))*abs(Y) + abs(B_s) )(i) )\n*     where abs(Z) is the componentwise absolute value of the matrix\n*     or vector Z.\n*\n\n*  Arguments\n*  ==========\n*\n*     N       (input) INTEGER\n*     The number of linear equations, i.e., the order of the\n*     matrix A.  N >= 0.\n*\n*     NZ      (input) INTEGER\n*     We add (NZ+1)*SLAMCH( 'Safe minimum' ) to R(i) in the numerator to\n*     guard against spuriously zero residuals. Default value is N.\n*\n*     NRHS    (input) INTEGER\n*     The number of right hand sides, i.e., the number of columns\n*     of the matrices AYB, RES, and BERR.  NRHS >= 0.\n*\n*     RES    (input) REAL array, dimension (N,NRHS)\n*     The residual matrix, i.e., the matrix R in the relative backward\n*     error formula above.\n*\n*     AYB    (input) REAL array, dimension (N, NRHS)\n*     The denominator in the relative backward error formula above, i.e.,\n*     the matrix abs(op(A_s))*abs(Y) + abs(B_s). The matrices A, Y, and B\n*     are from iterative refinement (see sla_gerfsx_extended.f).\n*     \n*     BERR   (output) REAL array, dimension (NRHS)\n*     The componentwise relative backward error from the formula above.\n*\n\n*  =====================================================================\n*\n*     .. Local Scalars ..\n      REAL               TMP\n      INTEGER            I, J\n*     ..\n*     .. Intrinsic Functions ..\n      INTRINSIC          ABS, MAX\n*     ..\n*     .. External Functions ..\n      EXTERNAL           SLAMCH\n      REAL               SLAMCH\n      REAL               SAFE1\n*     ..\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  berr = NumRu::Lapack.sla_lin_berr( nz, res, ayb, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 3)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 3)", argc);
  rblapack_nz = argv[0];
  rblapack_res = argv[1];
  rblapack_ayb = argv[2];
  if (rb_options != Qnil) {
  }

  if (!NA_IsNArray(rblapack_res))
    rb_raise(rb_eArgError, "res (2th argument) must be NArray");
  if (NA_RANK(rblapack_res) != 2)
    rb_raise(rb_eArgError, "rank of res (2th argument) must be %d", 2);
  nrhs = NA_SHAPE1(rblapack_res);
  n = NA_SHAPE0(rblapack_res);
  if (NA_TYPE(rblapack_res) != NA_SFLOAT)
    rblapack_res = na_change_type(rblapack_res, NA_SFLOAT);
  res = NA_PTR_TYPE(rblapack_res, real*);
  nz = NUM2INT(rblapack_nz);
  if (!NA_IsNArray(rblapack_ayb))
    rb_raise(rb_eArgError, "ayb (3th argument) must be NArray");
  if (NA_RANK(rblapack_ayb) != 2)
    rb_raise(rb_eArgError, "rank of ayb (3th argument) must be %d", 2);
  if (NA_SHAPE1(rblapack_ayb) != nrhs)
    rb_raise(rb_eRuntimeError, "shape 1 of ayb must be the same as shape 1 of res");
  if (NA_SHAPE0(rblapack_ayb) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of ayb must be the same as shape 0 of res");
  if (NA_TYPE(rblapack_ayb) != NA_SFLOAT)
    rblapack_ayb = na_change_type(rblapack_ayb, NA_SFLOAT);
  ayb = NA_PTR_TYPE(rblapack_ayb, real*);
  {
    int shape[1];
    shape[0] = nrhs;
    rblapack_berr = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  berr = NA_PTR_TYPE(rblapack_berr, real*);

  sla_lin_berr_(&n, &nz, &nrhs, res, ayb, berr);

  return rblapack_berr;
}

void
init_lapack_sla_lin_berr(VALUE mLapack){
  rb_define_module_function(mLapack, "sla_lin_berr", rblapack_sla_lin_berr, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}

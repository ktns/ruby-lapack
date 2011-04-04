#include "rb_lapack.h"

extern VOID sla_lin_berr_(integer *n, integer *nz, integer *nrhs, real *res, real *ayb, real *berr);

static VALUE
rb_sla_lin_berr(int argc, VALUE *argv, VALUE self){
  VALUE rb_nz;
  integer nz; 
  VALUE rb_res;
  real *res; 
  VALUE rb_ayb;
  real *ayb; 
  VALUE rb_berr;
  real *berr; 

  integer n;
  integer nrhs;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  berr = NumRu::Lapack.sla_lin_berr( nz, res, ayb)\n    or\n  NumRu::Lapack.sla_lin_berr  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE SLA_LIN_BERR ( N, NZ, NRHS, RES, AYB, BERR )\n\n*  Purpose\n*  =======\n*\n*     SLA_LIN_BERR computes componentwise relative backward error from\n*     the formula\n*         max(i) ( abs(R(i)) / ( abs(op(A_s))*abs(Y) + abs(B_s) )(i) )\n*     where abs(Z) is the componentwise absolute value of the matrix\n*     or vector Z.\n*\n\n*  Arguments\n*  ==========\n*\n*     N       (input) INTEGER\n*     The number of linear equations, i.e., the order of the\n*     matrix A.  N >= 0.\n*\n*     NZ      (input) INTEGER\n*     We add (NZ+1)*SLAMCH( 'Safe minimum' ) to R(i) in the numerator to\n*     guard against spuriously zero residuals. Default value is N.\n*\n*     NRHS    (input) INTEGER\n*     The number of right hand sides, i.e., the number of columns\n*     of the matrices AYB, RES, and BERR.  NRHS >= 0.\n*\n*     RES    (input) REAL array, dimension (N,NRHS)\n*     The residual matrix, i.e., the matrix R in the relative backward\n*     error formula above.\n*\n*     AYB    (input) REAL array, dimension (N, NRHS)\n*     The denominator in the relative backward error formula above, i.e.,\n*     the matrix abs(op(A_s))*abs(Y) + abs(B_s). The matrices A, Y, and B\n*     are from iterative refinement (see sla_gerfsx_extended.f).\n*     \n*     BERR   (output) REAL array, dimension (NRHS)\n*     The componentwise relative backward error from the formula above.\n*\n\n*  =====================================================================\n*\n*     .. Local Scalars ..\n      REAL               TMP\n      INTEGER            I, J\n*     ..\n*     .. Intrinsic Functions ..\n      INTRINSIC          ABS, MAX\n*     ..\n*     .. External Functions ..\n      EXTERNAL           SLAMCH\n      REAL               SLAMCH\n      REAL               SAFE1\n*     ..\n\n");
    return Qnil;
  }
  if (argc != 3)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 3)", argc);
  rb_nz = argv[0];
  rb_res = argv[1];
  rb_ayb = argv[2];

  if (!NA_IsNArray(rb_res))
    rb_raise(rb_eArgError, "res (2th argument) must be NArray");
  if (NA_RANK(rb_res) != 2)
    rb_raise(rb_eArgError, "rank of res (2th argument) must be %d", 2);
  nrhs = NA_SHAPE1(rb_res);
  n = NA_SHAPE0(rb_res);
  if (NA_TYPE(rb_res) != NA_SFLOAT)
    rb_res = na_change_type(rb_res, NA_SFLOAT);
  res = NA_PTR_TYPE(rb_res, real*);
  nz = NUM2INT(rb_nz);
  if (!NA_IsNArray(rb_ayb))
    rb_raise(rb_eArgError, "ayb (3th argument) must be NArray");
  if (NA_RANK(rb_ayb) != 2)
    rb_raise(rb_eArgError, "rank of ayb (3th argument) must be %d", 2);
  if (NA_SHAPE1(rb_ayb) != nrhs)
    rb_raise(rb_eRuntimeError, "shape 1 of ayb must be the same as shape 1 of res");
  if (NA_SHAPE0(rb_ayb) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of ayb must be the same as shape 0 of res");
  if (NA_TYPE(rb_ayb) != NA_SFLOAT)
    rb_ayb = na_change_type(rb_ayb, NA_SFLOAT);
  ayb = NA_PTR_TYPE(rb_ayb, real*);
  {
    int shape[1];
    shape[0] = nrhs;
    rb_berr = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  berr = NA_PTR_TYPE(rb_berr, real*);

  sla_lin_berr_(&n, &nz, &nrhs, res, ayb, berr);

  return rb_berr;
}

void
init_lapack_sla_lin_berr(VALUE mLapack){
  rb_define_module_function(mLapack, "sla_lin_berr", rb_sla_lin_berr, -1);
}

#include "rb_lapack.h"

static VALUE
rb_dla_lin_berr(int argc, VALUE *argv, VALUE self){
  VALUE rb_nz;
  integer nz; 
  VALUE rb_ayb;
  doublereal *ayb; 
  VALUE rb_res;
  doublereal *res; 
  VALUE rb_berr;
  doublereal berr; 

  integer n;
  integer nrhs;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  res, berr = NumRu::Lapack.dla_lin_berr( nz, ayb)\n    or\n  NumRu::Lapack.dla_lin_berr  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE DLA_LIN_BERR ( N, NZ, NRHS, RES, AYB, BERR )\n\n*  Purpose\n*  =======\n*\n*     DLA_LIN_BERR computes componentwise relative backward error from\n*     the formula\n*         max(i) ( abs(R(i)) / ( abs(op(A_s))*abs(Y) + abs(B_s) )(i) )\n*     where abs(Z) is the componentwise absolute value of the matrix\n*     or vector Z.\n*\n\n*  Arguments\n*  ==========\n*\n*     N       (input) INTEGER\n*     The number of linear equations, i.e., the order of the\n*     matrix A.  N >= 0.\n*\n*     NZ      (input) INTEGER\n*     We add (NZ+1)*SLAMCH( 'Safe minimum' ) to R(i) in the numerator to\n*     guard against spuriously zero residuals. Default value is N.\n*\n*     NRHS    (input) INTEGER\n*     The number of right hand sides, i.e., the number of columns\n*     of the matrices AYB, RES, and BERR.  NRHS >= 0.\n*\n*     RES    (input) DOUBLE PRECISION array, dimension (N,NRHS)\n*     The residual matrix, i.e., the matrix R in the relative backward\n*     error formula above.\n*\n*     AYB    (input) DOUBLE PRECISION array, dimension (N, NRHS)\n*     The denominator in the relative backward error formula above, i.e.,\n*     the matrix abs(op(A_s))*abs(Y) + abs(B_s). The matrices A, Y, and B\n*     are from iterative refinement (see dla_gerfsx_extended.f).\n*     \n*     RES    (output) DOUBLE PRECISION array, dimension (NRHS)\n*     The componentwise relative backward error from the formula above.\n*\n\n*  =====================================================================\n*\n*     .. Local Scalars ..\n      DOUBLE PRECISION   TMP\n      INTEGER            I, J\n*     ..\n*     .. Intrinsic Functions ..\n      INTRINSIC          ABS, MAX\n*     ..\n*     .. External Functions ..\n      EXTERNAL           DLAMCH\n      DOUBLE PRECISION   DLAMCH\n      DOUBLE PRECISION   SAFE1\n*     ..\n\n");
    return Qnil;
  }
  if (argc != 2)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 2)", argc);
  rb_nz = argv[0];
  rb_ayb = argv[1];

  nz = NUM2INT(rb_nz);
  if (!NA_IsNArray(rb_ayb))
    rb_raise(rb_eArgError, "ayb (2th argument) must be NArray");
  if (NA_RANK(rb_ayb) != 2)
    rb_raise(rb_eArgError, "rank of ayb (2th argument) must be %d", 2);
  n = NA_SHAPE0(rb_ayb);
  nrhs = NA_SHAPE1(rb_ayb);
  if (NA_TYPE(rb_ayb) != NA_DFLOAT)
    rb_ayb = na_change_type(rb_ayb, NA_DFLOAT);
  ayb = NA_PTR_TYPE(rb_ayb, doublereal*);
  {
    int shape[1];
    shape[0] = nrhs;
    rb_res = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  res = NA_PTR_TYPE(rb_res, doublereal*);

  dla_lin_berr_(&n, &nz, &nrhs, res, ayb, &berr);

  rb_berr = rb_float_new((double)berr);
  return rb_ary_new3(2, rb_res, rb_berr);
}

void
init_lapack_dla_lin_berr(VALUE mLapack){
  rb_define_module_function(mLapack, "dla_lin_berr", rb_dla_lin_berr, -1);
}

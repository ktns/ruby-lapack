#include "rb_lapack.h"

extern VOID slapll_(integer *n, real *x, integer *incx, real *y, integer *incy, real *ssmin);

static VALUE
rb_slapll(int argc, VALUE *argv, VALUE self){
  VALUE rb_n;
  integer n; 
  VALUE rb_x;
  real *x; 
  VALUE rb_incx;
  integer incx; 
  VALUE rb_y;
  real *y; 
  VALUE rb_incy;
  integer incy; 
  VALUE rb_ssmin;
  real ssmin; 
  VALUE rb_x_out__;
  real *x_out__;
  VALUE rb_y_out__;
  real *y_out__;


  if (argc == 0) {
    printf("%s\n", "USAGE:\n  ssmin, x, y = NumRu::Lapack.slapll( n, x, incx, y, incy)\n    or\n  NumRu::Lapack.slapll  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE SLAPLL( N, X, INCX, Y, INCY, SSMIN )\n\n*  Purpose\n*  =======\n*\n*  Given two column vectors X and Y, let\n*\n*                       A = ( X Y ).\n*\n*  The subroutine first computes the QR factorization of A = Q*R,\n*  and then computes the SVD of the 2-by-2 upper triangular matrix R.\n*  The smaller singular value of R is returned in SSMIN, which is used\n*  as the measurement of the linear dependency of the vectors X and Y.\n*\n\n*  Arguments\n*  =========\n*\n*  N       (input) INTEGER\n*          The length of the vectors X and Y.\n*\n*  X       (input/output) REAL array,\n*                         dimension (1+(N-1)*INCX)\n*          On entry, X contains the N-vector X.\n*          On exit, X is overwritten.\n*\n*  INCX    (input) INTEGER\n*          The increment between successive elements of X. INCX > 0.\n*\n*  Y       (input/output) REAL array,\n*                         dimension (1+(N-1)*INCY)\n*          On entry, Y contains the N-vector Y.\n*          On exit, Y is overwritten.\n*\n*  INCY    (input) INTEGER\n*          The increment between successive elements of Y. INCY > 0.\n*\n*  SSMIN   (output) REAL\n*          The smallest singular value of the N-by-2 matrix A = ( X Y ).\n*\n\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 5)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 5)", argc);
  rb_n = argv[0];
  rb_x = argv[1];
  rb_incx = argv[2];
  rb_y = argv[3];
  rb_incy = argv[4];

  incy = NUM2INT(rb_incy);
  incx = NUM2INT(rb_incx);
  n = NUM2INT(rb_n);
  if (!NA_IsNArray(rb_y))
    rb_raise(rb_eArgError, "y (4th argument) must be NArray");
  if (NA_RANK(rb_y) != 1)
    rb_raise(rb_eArgError, "rank of y (4th argument) must be %d", 1);
  if (NA_SHAPE0(rb_y) != (1+(n-1)*incy))
    rb_raise(rb_eRuntimeError, "shape 0 of y must be %d", 1+(n-1)*incy);
  if (NA_TYPE(rb_y) != NA_SFLOAT)
    rb_y = na_change_type(rb_y, NA_SFLOAT);
  y = NA_PTR_TYPE(rb_y, real*);
  if (!NA_IsNArray(rb_x))
    rb_raise(rb_eArgError, "x (2th argument) must be NArray");
  if (NA_RANK(rb_x) != 1)
    rb_raise(rb_eArgError, "rank of x (2th argument) must be %d", 1);
  if (NA_SHAPE0(rb_x) != (1+(n-1)*incx))
    rb_raise(rb_eRuntimeError, "shape 0 of x must be %d", 1+(n-1)*incx);
  if (NA_TYPE(rb_x) != NA_SFLOAT)
    rb_x = na_change_type(rb_x, NA_SFLOAT);
  x = NA_PTR_TYPE(rb_x, real*);
  {
    int shape[1];
    shape[0] = 1+(n-1)*incx;
    rb_x_out__ = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  x_out__ = NA_PTR_TYPE(rb_x_out__, real*);
  MEMCPY(x_out__, x, real, NA_TOTAL(rb_x));
  rb_x = rb_x_out__;
  x = x_out__;
  {
    int shape[1];
    shape[0] = 1+(n-1)*incy;
    rb_y_out__ = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  y_out__ = NA_PTR_TYPE(rb_y_out__, real*);
  MEMCPY(y_out__, y, real, NA_TOTAL(rb_y));
  rb_y = rb_y_out__;
  y = y_out__;

  slapll_(&n, x, &incx, y, &incy, &ssmin);

  rb_ssmin = rb_float_new((double)ssmin);
  return rb_ary_new3(3, rb_ssmin, rb_x, rb_y);
}

void
init_lapack_slapll(VALUE mLapack){
  rb_define_module_function(mLapack, "slapll", rb_slapll, -1);
}

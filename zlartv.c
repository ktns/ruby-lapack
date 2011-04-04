#include "rb_lapack.h"

extern VOID zlartv_(integer *n, doublecomplex *x, integer *incx, doublecomplex *y, integer *incy, doublereal *c, doublecomplex *s, integer *incc);

static VALUE
rb_zlartv(int argc, VALUE *argv, VALUE self){
  VALUE rb_n;
  integer n; 
  VALUE rb_x;
  doublecomplex *x; 
  VALUE rb_incx;
  integer incx; 
  VALUE rb_y;
  doublecomplex *y; 
  VALUE rb_incy;
  integer incy; 
  VALUE rb_c;
  doublereal *c; 
  VALUE rb_s;
  doublecomplex *s; 
  VALUE rb_incc;
  integer incc; 
  VALUE rb_x_out__;
  doublecomplex *x_out__;
  VALUE rb_y_out__;
  doublecomplex *y_out__;


  if (argc == 0) {
    printf("%s\n", "USAGE:\n  x, y = NumRu::Lapack.zlartv( n, x, incx, y, incy, c, s, incc)\n    or\n  NumRu::Lapack.zlartv  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE ZLARTV( N, X, INCX, Y, INCY, C, S, INCC )\n\n*  Purpose\n*  =======\n*\n*  ZLARTV applies a vector of complex plane rotations with real cosines\n*  to elements of the complex vectors x and y. For i = 1,2,...,n\n*\n*     ( x(i) ) := (        c(i)   s(i) ) ( x(i) )\n*     ( y(i) )    ( -conjg(s(i))  c(i) ) ( y(i) )\n*\n\n*  Arguments\n*  =========\n*\n*  N       (input) INTEGER\n*          The number of plane rotations to be applied.\n*\n*  X       (input/output) COMPLEX*16 array, dimension (1+(N-1)*INCX)\n*          The vector x.\n*\n*  INCX    (input) INTEGER\n*          The increment between elements of X. INCX > 0.\n*\n*  Y       (input/output) COMPLEX*16 array, dimension (1+(N-1)*INCY)\n*          The vector y.\n*\n*  INCY    (input) INTEGER\n*          The increment between elements of Y. INCY > 0.\n*\n*  C       (input) DOUBLE PRECISION array, dimension (1+(N-1)*INCC)\n*          The cosines of the plane rotations.\n*\n*  S       (input) COMPLEX*16 array, dimension (1+(N-1)*INCC)\n*          The sines of the plane rotations.\n*\n*  INCC    (input) INTEGER\n*          The increment between elements of C and S. INCC > 0.\n*\n\n*  =====================================================================\n*\n*     .. Local Scalars ..\n      INTEGER            I, IC, IX, IY\n      COMPLEX*16         XI, YI\n*     ..\n*     .. Intrinsic Functions ..\n      INTRINSIC          DCONJG\n*     ..\n\n");
    return Qnil;
  }
  if (argc != 8)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 8)", argc);
  rb_n = argv[0];
  rb_x = argv[1];
  rb_incx = argv[2];
  rb_y = argv[3];
  rb_incy = argv[4];
  rb_c = argv[5];
  rb_s = argv[6];
  rb_incc = argv[7];

  incx = NUM2INT(rb_incx);
  incy = NUM2INT(rb_incy);
  incc = NUM2INT(rb_incc);
  n = NUM2INT(rb_n);
  if (!NA_IsNArray(rb_x))
    rb_raise(rb_eArgError, "x (2th argument) must be NArray");
  if (NA_RANK(rb_x) != 1)
    rb_raise(rb_eArgError, "rank of x (2th argument) must be %d", 1);
  if (NA_SHAPE0(rb_x) != (1+(n-1)*incx))
    rb_raise(rb_eRuntimeError, "shape 0 of x must be %d", 1+(n-1)*incx);
  if (NA_TYPE(rb_x) != NA_DCOMPLEX)
    rb_x = na_change_type(rb_x, NA_DCOMPLEX);
  x = NA_PTR_TYPE(rb_x, doublecomplex*);
  if (!NA_IsNArray(rb_y))
    rb_raise(rb_eArgError, "y (4th argument) must be NArray");
  if (NA_RANK(rb_y) != 1)
    rb_raise(rb_eArgError, "rank of y (4th argument) must be %d", 1);
  if (NA_SHAPE0(rb_y) != (1+(n-1)*incy))
    rb_raise(rb_eRuntimeError, "shape 0 of y must be %d", 1+(n-1)*incy);
  if (NA_TYPE(rb_y) != NA_DCOMPLEX)
    rb_y = na_change_type(rb_y, NA_DCOMPLEX);
  y = NA_PTR_TYPE(rb_y, doublecomplex*);
  if (!NA_IsNArray(rb_c))
    rb_raise(rb_eArgError, "c (6th argument) must be NArray");
  if (NA_RANK(rb_c) != 1)
    rb_raise(rb_eArgError, "rank of c (6th argument) must be %d", 1);
  if (NA_SHAPE0(rb_c) != (1+(n-1)*incc))
    rb_raise(rb_eRuntimeError, "shape 0 of c must be %d", 1+(n-1)*incc);
  if (NA_TYPE(rb_c) != NA_DFLOAT)
    rb_c = na_change_type(rb_c, NA_DFLOAT);
  c = NA_PTR_TYPE(rb_c, doublereal*);
  if (!NA_IsNArray(rb_s))
    rb_raise(rb_eArgError, "s (7th argument) must be NArray");
  if (NA_RANK(rb_s) != 1)
    rb_raise(rb_eArgError, "rank of s (7th argument) must be %d", 1);
  if (NA_SHAPE0(rb_s) != (1+(n-1)*incc))
    rb_raise(rb_eRuntimeError, "shape 0 of s must be %d", 1+(n-1)*incc);
  if (NA_TYPE(rb_s) != NA_DCOMPLEX)
    rb_s = na_change_type(rb_s, NA_DCOMPLEX);
  s = NA_PTR_TYPE(rb_s, doublecomplex*);
  {
    int shape[1];
    shape[0] = 1+(n-1)*incx;
    rb_x_out__ = na_make_object(NA_DCOMPLEX, 1, shape, cNArray);
  }
  x_out__ = NA_PTR_TYPE(rb_x_out__, doublecomplex*);
  MEMCPY(x_out__, x, doublecomplex, NA_TOTAL(rb_x));
  rb_x = rb_x_out__;
  x = x_out__;
  {
    int shape[1];
    shape[0] = 1+(n-1)*incy;
    rb_y_out__ = na_make_object(NA_DCOMPLEX, 1, shape, cNArray);
  }
  y_out__ = NA_PTR_TYPE(rb_y_out__, doublecomplex*);
  MEMCPY(y_out__, y, doublecomplex, NA_TOTAL(rb_y));
  rb_y = rb_y_out__;
  y = y_out__;

  zlartv_(&n, x, &incx, y, &incy, c, s, &incc);

  return rb_ary_new3(2, rb_x, rb_y);
}

void
init_lapack_zlartv(VALUE mLapack){
  rb_define_module_function(mLapack, "zlartv", rb_zlartv, -1);
}

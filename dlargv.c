#include "rb_lapack.h"

extern VOID dlargv_(integer *n, doublereal *x, integer *incx, doublereal *y, integer *incy, doublereal *c, integer *incc);

static VALUE
rb_dlargv(int argc, VALUE *argv, VALUE self){
  VALUE rb_n;
  integer n; 
  VALUE rb_x;
  doublereal *x; 
  VALUE rb_incx;
  integer incx; 
  VALUE rb_y;
  doublereal *y; 
  VALUE rb_incy;
  integer incy; 
  VALUE rb_incc;
  integer incc; 
  VALUE rb_c;
  doublereal *c; 
  VALUE rb_x_out__;
  doublereal *x_out__;
  VALUE rb_y_out__;
  doublereal *y_out__;


  if (argc == 0) {
    printf("%s\n", "USAGE:\n  c, x, y = NumRu::Lapack.dlargv( n, x, incx, y, incy, incc)\n    or\n  NumRu::Lapack.dlargv  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE DLARGV( N, X, INCX, Y, INCY, C, INCC )\n\n*  Purpose\n*  =======\n*\n*  DLARGV generates a vector of real plane rotations, determined by\n*  elements of the real vectors x and y. For i = 1,2,...,n\n*\n*     (  c(i)  s(i) ) ( x(i) ) = ( a(i) )\n*     ( -s(i)  c(i) ) ( y(i) ) = (   0  )\n*\n\n*  Arguments\n*  =========\n*\n*  N       (input) INTEGER\n*          The number of plane rotations to be generated.\n*\n*  X       (input/output) DOUBLE PRECISION array,\n*                         dimension (1+(N-1)*INCX)\n*          On entry, the vector x.\n*          On exit, x(i) is overwritten by a(i), for i = 1,...,n.\n*\n*  INCX    (input) INTEGER\n*          The increment between elements of X. INCX > 0.\n*\n*  Y       (input/output) DOUBLE PRECISION array,\n*                         dimension (1+(N-1)*INCY)\n*          On entry, the vector y.\n*          On exit, the sines of the plane rotations.\n*\n*  INCY    (input) INTEGER\n*          The increment between elements of Y. INCY > 0.\n*\n*  C       (output) DOUBLE PRECISION array, dimension (1+(N-1)*INCC)\n*          The cosines of the plane rotations.\n*\n*  INCC    (input) INTEGER\n*          The increment between elements of C. INCC > 0.\n*\n\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 6)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 6)", argc);
  rb_n = argv[0];
  rb_x = argv[1];
  rb_incx = argv[2];
  rb_y = argv[3];
  rb_incy = argv[4];
  rb_incc = argv[5];

  incy = NUM2INT(rb_incy);
  incc = NUM2INT(rb_incc);
  n = NUM2INT(rb_n);
  incx = NUM2INT(rb_incx);
  if (!NA_IsNArray(rb_y))
    rb_raise(rb_eArgError, "y (4th argument) must be NArray");
  if (NA_RANK(rb_y) != 1)
    rb_raise(rb_eArgError, "rank of y (4th argument) must be %d", 1);
  if (NA_SHAPE0(rb_y) != (1+(n-1)*incy))
    rb_raise(rb_eRuntimeError, "shape 0 of y must be %d", 1+(n-1)*incy);
  if (NA_TYPE(rb_y) != NA_DFLOAT)
    rb_y = na_change_type(rb_y, NA_DFLOAT);
  y = NA_PTR_TYPE(rb_y, doublereal*);
  if (!NA_IsNArray(rb_x))
    rb_raise(rb_eArgError, "x (2th argument) must be NArray");
  if (NA_RANK(rb_x) != 1)
    rb_raise(rb_eArgError, "rank of x (2th argument) must be %d", 1);
  if (NA_SHAPE0(rb_x) != (1+(n-1)*incx))
    rb_raise(rb_eRuntimeError, "shape 0 of x must be %d", 1+(n-1)*incx);
  if (NA_TYPE(rb_x) != NA_DFLOAT)
    rb_x = na_change_type(rb_x, NA_DFLOAT);
  x = NA_PTR_TYPE(rb_x, doublereal*);
  {
    int shape[1];
    shape[0] = 1+(n-1)*incc;
    rb_c = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  c = NA_PTR_TYPE(rb_c, doublereal*);
  {
    int shape[1];
    shape[0] = 1+(n-1)*incx;
    rb_x_out__ = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  x_out__ = NA_PTR_TYPE(rb_x_out__, doublereal*);
  MEMCPY(x_out__, x, doublereal, NA_TOTAL(rb_x));
  rb_x = rb_x_out__;
  x = x_out__;
  {
    int shape[1];
    shape[0] = 1+(n-1)*incy;
    rb_y_out__ = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  y_out__ = NA_PTR_TYPE(rb_y_out__, doublereal*);
  MEMCPY(y_out__, y, doublereal, NA_TOTAL(rb_y));
  rb_y = rb_y_out__;
  y = y_out__;

  dlargv_(&n, x, &incx, y, &incy, c, &incc);

  return rb_ary_new3(3, rb_c, rb_x, rb_y);
}

void
init_lapack_dlargv(VALUE mLapack){
  rb_define_module_function(mLapack, "dlargv", rb_dlargv, -1);
}

#include "rb_lapack.h"

static VALUE
rb_zlargv(int argc, VALUE *argv, VALUE self){
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
  VALUE rb_incc;
  integer incc; 
  VALUE rb_c;
  doublereal *c; 
  VALUE rb_x_out__;
  doublecomplex *x_out__;
  VALUE rb_y_out__;
  doublecomplex *y_out__;


  if (argc == 0) {
    printf("%s\n", "USAGE:\n  c, x, y = NumRu::Lapack.zlargv( n, x, incx, y, incy, incc)\n    or\n  NumRu::Lapack.zlargv  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE ZLARGV( N, X, INCX, Y, INCY, C, INCC )\n\n*  Purpose\n*  =======\n*\n*  ZLARGV generates a vector of complex plane rotations with real\n*  cosines, determined by elements of the complex vectors x and y.\n*  For i = 1,2,...,n\n*\n*     (        c(i)   s(i) ) ( x(i) ) = ( r(i) )\n*     ( -conjg(s(i))  c(i) ) ( y(i) ) = (   0  )\n*\n*     where c(i)**2 + ABS(s(i))**2 = 1\n*\n*  The following conventions are used (these are the same as in ZLARTG,\n*  but differ from the BLAS1 routine ZROTG):\n*     If y(i)=0, then c(i)=1 and s(i)=0.\n*     If x(i)=0, then c(i)=0 and s(i) is chosen so that r(i) is real.\n*\n\n*  Arguments\n*  =========\n*\n*  N       (input) INTEGER\n*          The number of plane rotations to be generated.\n*\n*  X       (input/output) COMPLEX*16 array, dimension (1+(N-1)*INCX)\n*          On entry, the vector x.\n*          On exit, x(i) is overwritten by r(i), for i = 1,...,n.\n*\n*  INCX    (input) INTEGER\n*          The increment between elements of X. INCX > 0.\n*\n*  Y       (input/output) COMPLEX*16 array, dimension (1+(N-1)*INCY)\n*          On entry, the vector y.\n*          On exit, the sines of the plane rotations.\n*\n*  INCY    (input) INTEGER\n*          The increment between elements of Y. INCY > 0.\n*\n*  C       (output) DOUBLE PRECISION array, dimension (1+(N-1)*INCC)\n*          The cosines of the plane rotations.\n*\n*  INCC    (input) INTEGER\n*          The increment between elements of C. INCC > 0.\n*\n\n*  Further Details\n*  ======= =======\n*\n*  6-6-96 - Modified with a new algorithm by W. Kahan and J. Demmel\n*\n*  This version has a few statements commented out for thread safety\n*  (machine parameters are computed on each entry). 10 feb 03, SJH.\n*\n*  =====================================================================\n*\n\n");
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

  n = NUM2INT(rb_n);
  incx = NUM2INT(rb_incx);
  incy = NUM2INT(rb_incy);
  incc = NUM2INT(rb_incc);
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
  {
    int shape[1];
    shape[0] = 1+(n-1)*incc;
    rb_c = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  c = NA_PTR_TYPE(rb_c, doublereal*);
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

  zlargv_(&n, x, &incx, y, &incy, c, &incc);

  return rb_ary_new3(3, rb_c, rb_x, rb_y);
}

void
init_lapack_zlargv(VALUE mLapack){
  rb_define_module_function(mLapack, "zlargv", rb_zlargv, -1);
}

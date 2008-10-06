#include "rb_lapack.h"

static VALUE
rb_zlacrt(int argc, VALUE *argv, VALUE self){
  VALUE rb_cx;
  doublecomplex *cx; 
  VALUE rb_incx;
  integer incx; 
  VALUE rb_cy;
  doublecomplex *cy; 
  VALUE rb_incy;
  integer incy; 
  VALUE rb_c;
  doublecomplex c; 
  VALUE rb_s;
  doublecomplex s; 
  VALUE rb_cx_out__;
  doublecomplex *cx_out__;
  VALUE rb_cy_out__;
  doublecomplex *cy_out__;

  integer n;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  cx, cy = NumRu::Lapack.zlacrt( cx, incx, cy, incy, c, s)\n    or\n  NumRu::Lapack.zlacrt  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE ZLACRT( N, CX, INCX, CY, INCY, C, S )\n\n*  Purpose\n*  =======\n*\n*  ZLACRT performs the operation\n*\n*     (  c  s )( x )  ==> ( x )\n*     ( -s  c )( y )      ( y )\n*\n*  where c and s are complex and the vectors x and y are complex.\n*\n\n*  Arguments\n*  =========\n*\n*  N       (input) INTEGER\n*          The number of elements in the vectors CX and CY.\n*\n*  CX      (input/output) COMPLEX*16 array, dimension (N)\n*          On input, the vector x.\n*          On output, CX is overwritten with c*x + s*y.\n*\n*  INCX    (input) INTEGER\n*          The increment between successive values of CX.  INCX <> 0.\n*\n*  CY      (input/output) COMPLEX*16 array, dimension (N)\n*          On input, the vector y.\n*          On output, CY is overwritten with -s*x + c*y.\n*\n*  INCY    (input) INTEGER\n*          The increment between successive values of CY.  INCY <> 0.\n*\n*  C       (input) COMPLEX*16\n*  S       (input) COMPLEX*16\n*          C and S define the matrix\n*             [  C   S  ].\n*             [ -S   C  ]\n*\n\n* =====================================================================\n*\n*     .. Local Scalars ..\n      INTEGER            I, IX, IY\n      COMPLEX*16         CTEMP\n*     ..\n\n");
    return Qnil;
  }
  if (argc != 6)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 6)", argc);
  rb_cx = argv[0];
  rb_incx = argv[1];
  rb_cy = argv[2];
  rb_incy = argv[3];
  rb_c = argv[4];
  rb_s = argv[5];

  incx = NUM2INT(rb_incx);
  incy = NUM2INT(rb_incy);
  c.r = NUM2DBL(rb_funcall(rb_c, rb_intern("real"), 0));
  c.i = NUM2DBL(rb_funcall(rb_c, rb_intern("imag"), 0));
  s.r = NUM2DBL(rb_funcall(rb_s, rb_intern("real"), 0));
  s.i = NUM2DBL(rb_funcall(rb_s, rb_intern("imag"), 0));
  if (!NA_IsNArray(rb_cx))
    rb_raise(rb_eArgError, "cx (1th argument) must be NArray");
  if (NA_RANK(rb_cx) != 1)
    rb_raise(rb_eArgError, "rank of cx (1th argument) must be %d", 1);
  n = NA_SHAPE0(rb_cx);
  if (NA_TYPE(rb_cx) != NA_DCOMPLEX)
    rb_cx = na_change_type(rb_cx, NA_DCOMPLEX);
  cx = NA_PTR_TYPE(rb_cx, doublecomplex*);
  if (!NA_IsNArray(rb_cy))
    rb_raise(rb_eArgError, "cy (3th argument) must be NArray");
  if (NA_RANK(rb_cy) != 1)
    rb_raise(rb_eArgError, "rank of cy (3th argument) must be %d", 1);
  if (NA_SHAPE0(rb_cy) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of cy must be the same as shape 0 of cx");
  if (NA_TYPE(rb_cy) != NA_DCOMPLEX)
    rb_cy = na_change_type(rb_cy, NA_DCOMPLEX);
  cy = NA_PTR_TYPE(rb_cy, doublecomplex*);
  {
    int shape[1];
    shape[0] = n;
    rb_cx_out__ = na_make_object(NA_DCOMPLEX, 1, shape, cNArray);
  }
  cx_out__ = NA_PTR_TYPE(rb_cx_out__, doublecomplex*);
  MEMCPY(cx_out__, cx, doublecomplex, NA_TOTAL(rb_cx));
  rb_cx = rb_cx_out__;
  cx = cx_out__;
  {
    int shape[1];
    shape[0] = n;
    rb_cy_out__ = na_make_object(NA_DCOMPLEX, 1, shape, cNArray);
  }
  cy_out__ = NA_PTR_TYPE(rb_cy_out__, doublecomplex*);
  MEMCPY(cy_out__, cy, doublecomplex, NA_TOTAL(rb_cy));
  rb_cy = rb_cy_out__;
  cy = cy_out__;

  zlacrt_(&n, cx, &incx, cy, &incy, &c, &s);

  return rb_ary_new3(2, rb_cx, rb_cy);
}

void
init_lapack_zlacrt(VALUE mLapack){
  rb_define_module_function(mLapack, "zlacrt", rb_zlacrt, -1);
}

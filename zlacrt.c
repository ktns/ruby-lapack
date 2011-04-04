#include "rb_lapack.h"

extern VOID zlacrt_(integer *n, doublecomplex *cx, integer *incx, doublecomplex *cy, integer *incy, doublecomplex *c, doublecomplex *s);

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
    printf("%s\n", "USAGE:\n  cx, cy = NumRu::Lapack.zlacrt( cx, incx, cy, incy, c, s)\n    or\n  NumRu::Lapack.zlacrt  # print help\n\n\nFORTRAN MANUAL\n\n");
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

  if (!NA_IsNArray(rb_cy))
    rb_raise(rb_eArgError, "cy (3th argument) must be NArray");
  if (NA_RANK(rb_cy) != 1)
    rb_raise(rb_eArgError, "rank of cy (3th argument) must be %d", 1);
  n = NA_SHAPE0(rb_cy);
  if (NA_TYPE(rb_cy) != NA_DCOMPLEX)
    rb_cy = na_change_type(rb_cy, NA_DCOMPLEX);
  cy = NA_PTR_TYPE(rb_cy, doublecomplex*);
  c.r = NUM2DBL(rb_funcall(rb_c, rb_intern("real"), 0));
  c.i = NUM2DBL(rb_funcall(rb_c, rb_intern("imag"), 0));
  incx = NUM2INT(rb_incx);
  incy = NUM2INT(rb_incy);
  s.r = NUM2DBL(rb_funcall(rb_s, rb_intern("real"), 0));
  s.i = NUM2DBL(rb_funcall(rb_s, rb_intern("imag"), 0));
  if (!NA_IsNArray(rb_cx))
    rb_raise(rb_eArgError, "cx (1th argument) must be NArray");
  if (NA_RANK(rb_cx) != 1)
    rb_raise(rb_eArgError, "rank of cx (1th argument) must be %d", 1);
  if (NA_SHAPE0(rb_cx) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of cx must be the same as shape 0 of cy");
  if (NA_TYPE(rb_cx) != NA_DCOMPLEX)
    rb_cx = na_change_type(rb_cx, NA_DCOMPLEX);
  cx = NA_PTR_TYPE(rb_cx, doublecomplex*);
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

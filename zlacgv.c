#include "rb_lapack.h"

extern VOID zlacgv_(integer *n, doublecomplex *x, integer *incx);

static VALUE
rb_zlacgv(int argc, VALUE *argv, VALUE self){
  VALUE rb_n;
  integer n; 
  VALUE rb_x;
  doublecomplex *x; 
  VALUE rb_incx;
  integer incx; 
  VALUE rb_x_out__;
  doublecomplex *x_out__;


  if (argc == 0) {
    printf("%s\n", "USAGE:\n  x = NumRu::Lapack.zlacgv( n, x, incx)\n    or\n  NumRu::Lapack.zlacgv  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 3)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 3)", argc);
  rb_n = argv[0];
  rb_x = argv[1];
  rb_incx = argv[2];

  incx = NUM2INT(rb_incx);
  n = NUM2INT(rb_n);
  if (!NA_IsNArray(rb_x))
    rb_raise(rb_eArgError, "x (2th argument) must be NArray");
  if (NA_RANK(rb_x) != 1)
    rb_raise(rb_eArgError, "rank of x (2th argument) must be %d", 1);
  if (NA_SHAPE0(rb_x) != (1+(n-1)*abs(incx)))
    rb_raise(rb_eRuntimeError, "shape 0 of x must be %d", 1+(n-1)*abs(incx));
  if (NA_TYPE(rb_x) != NA_DCOMPLEX)
    rb_x = na_change_type(rb_x, NA_DCOMPLEX);
  x = NA_PTR_TYPE(rb_x, doublecomplex*);
  {
    int shape[1];
    shape[0] = 1+(n-1)*abs(incx);
    rb_x_out__ = na_make_object(NA_DCOMPLEX, 1, shape, cNArray);
  }
  x_out__ = NA_PTR_TYPE(rb_x_out__, doublecomplex*);
  MEMCPY(x_out__, x, doublecomplex, NA_TOTAL(rb_x));
  rb_x = rb_x_out__;
  x = x_out__;

  zlacgv_(&n, x, &incx);

  return rb_x;
}

void
init_lapack_zlacgv(VALUE mLapack){
  rb_define_module_function(mLapack, "zlacgv", rb_zlacgv, -1);
}

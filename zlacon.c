#include "rb_lapack.h"

extern VOID zlacon_(integer *n, doublecomplex *v, doublecomplex *x, doublereal *est, integer *kase);

static VALUE
rb_zlacon(int argc, VALUE *argv, VALUE self){
  VALUE rb_x;
  doublecomplex *x; 
  VALUE rb_est;
  doublereal est; 
  VALUE rb_kase;
  integer kase; 
  VALUE rb_x_out__;
  doublecomplex *x_out__;
  doublecomplex *v;

  integer n;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  x, est, kase = NumRu::Lapack.zlacon( x, est, kase)\n    or\n  NumRu::Lapack.zlacon  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 3)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 3)", argc);
  rb_x = argv[0];
  rb_est = argv[1];
  rb_kase = argv[2];

  if (!NA_IsNArray(rb_x))
    rb_raise(rb_eArgError, "x (1th argument) must be NArray");
  if (NA_RANK(rb_x) != 1)
    rb_raise(rb_eArgError, "rank of x (1th argument) must be %d", 1);
  n = NA_SHAPE0(rb_x);
  if (NA_TYPE(rb_x) != NA_DCOMPLEX)
    rb_x = na_change_type(rb_x, NA_DCOMPLEX);
  x = NA_PTR_TYPE(rb_x, doublecomplex*);
  est = NUM2DBL(rb_est);
  kase = NUM2INT(rb_kase);
  {
    int shape[1];
    shape[0] = n;
    rb_x_out__ = na_make_object(NA_DCOMPLEX, 1, shape, cNArray);
  }
  x_out__ = NA_PTR_TYPE(rb_x_out__, doublecomplex*);
  MEMCPY(x_out__, x, doublecomplex, NA_TOTAL(rb_x));
  rb_x = rb_x_out__;
  x = x_out__;
  v = ALLOC_N(doublecomplex, (n));

  zlacon_(&n, v, x, &est, &kase);

  free(v);
  rb_est = rb_float_new((double)est);
  rb_kase = INT2NUM(kase);
  return rb_ary_new3(3, rb_x, rb_est, rb_kase);
}

void
init_lapack_zlacon(VALUE mLapack){
  rb_define_module_function(mLapack, "zlacon", rb_zlacon, -1);
}

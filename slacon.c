#include "rb_lapack.h"

extern VOID slacon_(integer *n, real *v, real *x, integer *isgn, real *est, integer *kase);

static VALUE
rb_slacon(int argc, VALUE *argv, VALUE self){
  VALUE rb_x;
  real *x; 
  VALUE rb_est;
  real est; 
  VALUE rb_kase;
  integer kase; 
  VALUE rb_x_out__;
  real *x_out__;
  real *v;
  integer *isgn;

  integer n;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  x, est, kase = NumRu::Lapack.slacon( x, est, kase)\n    or\n  NumRu::Lapack.slacon  # print help\n\n\nFORTRAN MANUAL\n\n");
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
  if (NA_TYPE(rb_x) != NA_SFLOAT)
    rb_x = na_change_type(rb_x, NA_SFLOAT);
  x = NA_PTR_TYPE(rb_x, real*);
  est = (real)NUM2DBL(rb_est);
  kase = NUM2INT(rb_kase);
  {
    int shape[1];
    shape[0] = n;
    rb_x_out__ = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  x_out__ = NA_PTR_TYPE(rb_x_out__, real*);
  MEMCPY(x_out__, x, real, NA_TOTAL(rb_x));
  rb_x = rb_x_out__;
  x = x_out__;
  v = ALLOC_N(real, (n));
  isgn = ALLOC_N(integer, (n));

  slacon_(&n, v, x, isgn, &est, &kase);

  free(v);
  free(isgn);
  rb_est = rb_float_new((double)est);
  rb_kase = INT2NUM(kase);
  return rb_ary_new3(3, rb_x, rb_est, rb_kase);
}

void
init_lapack_slacon(VALUE mLapack){
  rb_define_module_function(mLapack, "slacon", rb_slacon, -1);
}

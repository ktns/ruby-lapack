#include "rb_lapack.h"

extern VOID slapmr_(logical *forwrd, integer *m, integer *n, real *x, integer *ldx, integer *k);

static VALUE
rb_slapmr(int argc, VALUE *argv, VALUE self){
  VALUE rb_forwrd;
  logical forwrd; 
  VALUE rb_x;
  real *x; 
  VALUE rb_k;
  integer *k; 
  VALUE rb_x_out__;
  real *x_out__;
  VALUE rb_k_out__;
  integer *k_out__;

  integer ldx;
  integer n;
  integer m;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  x, k = NumRu::Lapack.slapmr( forwrd, x, k)\n    or\n  NumRu::Lapack.slapmr  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 3)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 3)", argc);
  rb_forwrd = argv[0];
  rb_x = argv[1];
  rb_k = argv[2];

  if (!NA_IsNArray(rb_k))
    rb_raise(rb_eArgError, "k (3th argument) must be NArray");
  if (NA_RANK(rb_k) != 1)
    rb_raise(rb_eArgError, "rank of k (3th argument) must be %d", 1);
  m = NA_SHAPE0(rb_k);
  if (NA_TYPE(rb_k) != NA_LINT)
    rb_k = na_change_type(rb_k, NA_LINT);
  k = NA_PTR_TYPE(rb_k, integer*);
  forwrd = (rb_forwrd == Qtrue);
  if (!NA_IsNArray(rb_x))
    rb_raise(rb_eArgError, "x (2th argument) must be NArray");
  if (NA_RANK(rb_x) != 2)
    rb_raise(rb_eArgError, "rank of x (2th argument) must be %d", 2);
  n = NA_SHAPE1(rb_x);
  ldx = NA_SHAPE0(rb_x);
  if (NA_TYPE(rb_x) != NA_SFLOAT)
    rb_x = na_change_type(rb_x, NA_SFLOAT);
  x = NA_PTR_TYPE(rb_x, real*);
  {
    int shape[2];
    shape[0] = ldx;
    shape[1] = n;
    rb_x_out__ = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  x_out__ = NA_PTR_TYPE(rb_x_out__, real*);
  MEMCPY(x_out__, x, real, NA_TOTAL(rb_x));
  rb_x = rb_x_out__;
  x = x_out__;
  {
    int shape[1];
    shape[0] = m;
    rb_k_out__ = na_make_object(NA_LINT, 1, shape, cNArray);
  }
  k_out__ = NA_PTR_TYPE(rb_k_out__, integer*);
  MEMCPY(k_out__, k, integer, NA_TOTAL(rb_k));
  rb_k = rb_k_out__;
  k = k_out__;

  slapmr_(&forwrd, &m, &n, x, &ldx, k);

  return rb_ary_new3(2, rb_x, rb_k);
}

void
init_lapack_slapmr(VALUE mLapack){
  rb_define_module_function(mLapack, "slapmr", rb_slapmr, -1);
}

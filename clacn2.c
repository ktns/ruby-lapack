#include "rb_lapack.h"

extern VOID clacn2_(integer *n, complex *v, complex *x, real *est, integer *kase, integer *isave);

static VALUE
rb_clacn2(int argc, VALUE *argv, VALUE self){
  VALUE rb_x;
  complex *x; 
  VALUE rb_est;
  real est; 
  VALUE rb_kase;
  integer kase; 
  VALUE rb_isave;
  integer *isave; 
  VALUE rb_x_out__;
  complex *x_out__;
  VALUE rb_isave_out__;
  integer *isave_out__;
  complex *v;

  integer n;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  x, est, kase, isave = NumRu::Lapack.clacn2( x, est, kase, isave)\n    or\n  NumRu::Lapack.clacn2  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 4)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 4)", argc);
  rb_x = argv[0];
  rb_est = argv[1];
  rb_kase = argv[2];
  rb_isave = argv[3];

  if (!NA_IsNArray(rb_x))
    rb_raise(rb_eArgError, "x (1th argument) must be NArray");
  if (NA_RANK(rb_x) != 1)
    rb_raise(rb_eArgError, "rank of x (1th argument) must be %d", 1);
  n = NA_SHAPE0(rb_x);
  if (NA_TYPE(rb_x) != NA_SCOMPLEX)
    rb_x = na_change_type(rb_x, NA_SCOMPLEX);
  x = NA_PTR_TYPE(rb_x, complex*);
  est = (real)NUM2DBL(rb_est);
  if (!NA_IsNArray(rb_isave))
    rb_raise(rb_eArgError, "isave (4th argument) must be NArray");
  if (NA_RANK(rb_isave) != 1)
    rb_raise(rb_eArgError, "rank of isave (4th argument) must be %d", 1);
  if (NA_SHAPE0(rb_isave) != (3))
    rb_raise(rb_eRuntimeError, "shape 0 of isave must be %d", 3);
  if (NA_TYPE(rb_isave) != NA_LINT)
    rb_isave = na_change_type(rb_isave, NA_LINT);
  isave = NA_PTR_TYPE(rb_isave, integer*);
  kase = NUM2INT(rb_kase);
  {
    int shape[1];
    shape[0] = n;
    rb_x_out__ = na_make_object(NA_SCOMPLEX, 1, shape, cNArray);
  }
  x_out__ = NA_PTR_TYPE(rb_x_out__, complex*);
  MEMCPY(x_out__, x, complex, NA_TOTAL(rb_x));
  rb_x = rb_x_out__;
  x = x_out__;
  {
    int shape[1];
    shape[0] = 3;
    rb_isave_out__ = na_make_object(NA_LINT, 1, shape, cNArray);
  }
  isave_out__ = NA_PTR_TYPE(rb_isave_out__, integer*);
  MEMCPY(isave_out__, isave, integer, NA_TOTAL(rb_isave));
  rb_isave = rb_isave_out__;
  isave = isave_out__;
  v = ALLOC_N(complex, (n));

  clacn2_(&n, v, x, &est, &kase, isave);

  free(v);
  rb_est = rb_float_new((double)est);
  rb_kase = INT2NUM(kase);
  return rb_ary_new3(4, rb_x, rb_est, rb_kase, rb_isave);
}

void
init_lapack_clacn2(VALUE mLapack){
  rb_define_module_function(mLapack, "clacn2", rb_clacn2, -1);
}

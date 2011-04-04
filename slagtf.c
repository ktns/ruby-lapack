#include "rb_lapack.h"

extern VOID slagtf_(integer *n, real *a, real *lambda, real *b, real *c, real *tol, real *d, integer *in, integer *info);

static VALUE
rb_slagtf(int argc, VALUE *argv, VALUE self){
  VALUE rb_a;
  real *a; 
  VALUE rb_lambda;
  real lambda; 
  VALUE rb_b;
  real *b; 
  VALUE rb_c;
  real *c; 
  VALUE rb_tol;
  real tol; 
  VALUE rb_d;
  real *d; 
  VALUE rb_in;
  integer *in; 
  VALUE rb_info;
  integer info; 
  VALUE rb_a_out__;
  real *a_out__;
  VALUE rb_b_out__;
  real *b_out__;
  VALUE rb_c_out__;
  real *c_out__;

  integer n;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  d, in, info, a, b, c = NumRu::Lapack.slagtf( a, lambda, b, c, tol)\n    or\n  NumRu::Lapack.slagtf  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 5)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 5)", argc);
  rb_a = argv[0];
  rb_lambda = argv[1];
  rb_b = argv[2];
  rb_c = argv[3];
  rb_tol = argv[4];

  if (!NA_IsNArray(rb_a))
    rb_raise(rb_eArgError, "a (1th argument) must be NArray");
  if (NA_RANK(rb_a) != 1)
    rb_raise(rb_eArgError, "rank of a (1th argument) must be %d", 1);
  n = NA_SHAPE0(rb_a);
  if (NA_TYPE(rb_a) != NA_SFLOAT)
    rb_a = na_change_type(rb_a, NA_SFLOAT);
  a = NA_PTR_TYPE(rb_a, real*);
  tol = (real)NUM2DBL(rb_tol);
  lambda = (real)NUM2DBL(rb_lambda);
  if (!NA_IsNArray(rb_b))
    rb_raise(rb_eArgError, "b (3th argument) must be NArray");
  if (NA_RANK(rb_b) != 1)
    rb_raise(rb_eArgError, "rank of b (3th argument) must be %d", 1);
  if (NA_SHAPE0(rb_b) != (n-1))
    rb_raise(rb_eRuntimeError, "shape 0 of b must be %d", n-1);
  if (NA_TYPE(rb_b) != NA_SFLOAT)
    rb_b = na_change_type(rb_b, NA_SFLOAT);
  b = NA_PTR_TYPE(rb_b, real*);
  if (!NA_IsNArray(rb_c))
    rb_raise(rb_eArgError, "c (4th argument) must be NArray");
  if (NA_RANK(rb_c) != 1)
    rb_raise(rb_eArgError, "rank of c (4th argument) must be %d", 1);
  if (NA_SHAPE0(rb_c) != (n-1))
    rb_raise(rb_eRuntimeError, "shape 0 of c must be %d", n-1);
  if (NA_TYPE(rb_c) != NA_SFLOAT)
    rb_c = na_change_type(rb_c, NA_SFLOAT);
  c = NA_PTR_TYPE(rb_c, real*);
  {
    int shape[1];
    shape[0] = n-2;
    rb_d = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  d = NA_PTR_TYPE(rb_d, real*);
  {
    int shape[1];
    shape[0] = n;
    rb_in = na_make_object(NA_LINT, 1, shape, cNArray);
  }
  in = NA_PTR_TYPE(rb_in, integer*);
  {
    int shape[1];
    shape[0] = n;
    rb_a_out__ = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  a_out__ = NA_PTR_TYPE(rb_a_out__, real*);
  MEMCPY(a_out__, a, real, NA_TOTAL(rb_a));
  rb_a = rb_a_out__;
  a = a_out__;
  {
    int shape[1];
    shape[0] = n-1;
    rb_b_out__ = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  b_out__ = NA_PTR_TYPE(rb_b_out__, real*);
  MEMCPY(b_out__, b, real, NA_TOTAL(rb_b));
  rb_b = rb_b_out__;
  b = b_out__;
  {
    int shape[1];
    shape[0] = n-1;
    rb_c_out__ = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  c_out__ = NA_PTR_TYPE(rb_c_out__, real*);
  MEMCPY(c_out__, c, real, NA_TOTAL(rb_c));
  rb_c = rb_c_out__;
  c = c_out__;

  slagtf_(&n, a, &lambda, b, c, &tol, d, in, &info);

  rb_info = INT2NUM(info);
  return rb_ary_new3(6, rb_d, rb_in, rb_info, rb_a, rb_b, rb_c);
}

void
init_lapack_slagtf(VALUE mLapack){
  rb_define_module_function(mLapack, "slagtf", rb_slagtf, -1);
}

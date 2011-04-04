#include "rb_lapack.h"

extern VOID slagts_(integer *job, integer *n, real *a, real *b, real *c, real *d, integer *in, real *y, real *tol, integer *info);

static VALUE
rb_slagts(int argc, VALUE *argv, VALUE self){
  VALUE rb_job;
  integer job; 
  VALUE rb_a;
  real *a; 
  VALUE rb_b;
  real *b; 
  VALUE rb_c;
  real *c; 
  VALUE rb_d;
  real *d; 
  VALUE rb_in;
  integer *in; 
  VALUE rb_y;
  real *y; 
  VALUE rb_tol;
  real tol; 
  VALUE rb_info;
  integer info; 
  VALUE rb_y_out__;
  real *y_out__;

  integer n;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  info, y, tol = NumRu::Lapack.slagts( job, a, b, c, d, in, y, tol)\n    or\n  NumRu::Lapack.slagts  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 8)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 8)", argc);
  rb_job = argv[0];
  rb_a = argv[1];
  rb_b = argv[2];
  rb_c = argv[3];
  rb_d = argv[4];
  rb_in = argv[5];
  rb_y = argv[6];
  rb_tol = argv[7];

  if (!NA_IsNArray(rb_a))
    rb_raise(rb_eArgError, "a (2th argument) must be NArray");
  if (NA_RANK(rb_a) != 1)
    rb_raise(rb_eArgError, "rank of a (2th argument) must be %d", 1);
  n = NA_SHAPE0(rb_a);
  if (NA_TYPE(rb_a) != NA_SFLOAT)
    rb_a = na_change_type(rb_a, NA_SFLOAT);
  a = NA_PTR_TYPE(rb_a, real*);
  tol = (real)NUM2DBL(rb_tol);
  if (!NA_IsNArray(rb_y))
    rb_raise(rb_eArgError, "y (7th argument) must be NArray");
  if (NA_RANK(rb_y) != 1)
    rb_raise(rb_eArgError, "rank of y (7th argument) must be %d", 1);
  if (NA_SHAPE0(rb_y) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of y must be the same as shape 0 of a");
  if (NA_TYPE(rb_y) != NA_SFLOAT)
    rb_y = na_change_type(rb_y, NA_SFLOAT);
  y = NA_PTR_TYPE(rb_y, real*);
  if (!NA_IsNArray(rb_in))
    rb_raise(rb_eArgError, "in (6th argument) must be NArray");
  if (NA_RANK(rb_in) != 1)
    rb_raise(rb_eArgError, "rank of in (6th argument) must be %d", 1);
  if (NA_SHAPE0(rb_in) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of in must be the same as shape 0 of a");
  if (NA_TYPE(rb_in) != NA_LINT)
    rb_in = na_change_type(rb_in, NA_LINT);
  in = NA_PTR_TYPE(rb_in, integer*);
  job = NUM2INT(rb_job);
  if (!NA_IsNArray(rb_c))
    rb_raise(rb_eArgError, "c (4th argument) must be NArray");
  if (NA_RANK(rb_c) != 1)
    rb_raise(rb_eArgError, "rank of c (4th argument) must be %d", 1);
  if (NA_SHAPE0(rb_c) != (n-1))
    rb_raise(rb_eRuntimeError, "shape 0 of c must be %d", n-1);
  if (NA_TYPE(rb_c) != NA_SFLOAT)
    rb_c = na_change_type(rb_c, NA_SFLOAT);
  c = NA_PTR_TYPE(rb_c, real*);
  if (!NA_IsNArray(rb_b))
    rb_raise(rb_eArgError, "b (3th argument) must be NArray");
  if (NA_RANK(rb_b) != 1)
    rb_raise(rb_eArgError, "rank of b (3th argument) must be %d", 1);
  if (NA_SHAPE0(rb_b) != (n-1))
    rb_raise(rb_eRuntimeError, "shape 0 of b must be %d", n-1);
  if (NA_TYPE(rb_b) != NA_SFLOAT)
    rb_b = na_change_type(rb_b, NA_SFLOAT);
  b = NA_PTR_TYPE(rb_b, real*);
  if (!NA_IsNArray(rb_d))
    rb_raise(rb_eArgError, "d (5th argument) must be NArray");
  if (NA_RANK(rb_d) != 1)
    rb_raise(rb_eArgError, "rank of d (5th argument) must be %d", 1);
  if (NA_SHAPE0(rb_d) != (n-2))
    rb_raise(rb_eRuntimeError, "shape 0 of d must be %d", n-2);
  if (NA_TYPE(rb_d) != NA_SFLOAT)
    rb_d = na_change_type(rb_d, NA_SFLOAT);
  d = NA_PTR_TYPE(rb_d, real*);
  {
    int shape[1];
    shape[0] = n;
    rb_y_out__ = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  y_out__ = NA_PTR_TYPE(rb_y_out__, real*);
  MEMCPY(y_out__, y, real, NA_TOTAL(rb_y));
  rb_y = rb_y_out__;
  y = y_out__;

  slagts_(&job, &n, a, b, c, d, in, y, &tol, &info);

  rb_info = INT2NUM(info);
  rb_tol = rb_float_new((double)tol);
  return rb_ary_new3(3, rb_info, rb_y, rb_tol);
}

void
init_lapack_slagts(VALUE mLapack){
  rb_define_module_function(mLapack, "slagts", rb_slagts, -1);
}

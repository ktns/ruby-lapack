#include "rb_lapack.h"

extern VOID dlagts_(integer *job, integer *n, doublereal *a, doublereal *b, doublereal *c, doublereal *d, integer *in, doublereal *y, doublereal *tol, integer *info);

static VALUE
rb_dlagts(int argc, VALUE *argv, VALUE self){
  VALUE rb_job;
  integer job; 
  VALUE rb_a;
  doublereal *a; 
  VALUE rb_b;
  doublereal *b; 
  VALUE rb_c;
  doublereal *c; 
  VALUE rb_d;
  doublereal *d; 
  VALUE rb_in;
  integer *in; 
  VALUE rb_y;
  doublereal *y; 
  VALUE rb_tol;
  doublereal tol; 
  VALUE rb_info;
  integer info; 
  VALUE rb_y_out__;
  doublereal *y_out__;

  integer n;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  info, y, tol = NumRu::Lapack.dlagts( job, a, b, c, d, in, y, tol)\n    or\n  NumRu::Lapack.dlagts  # print help\n\n\nFORTRAN MANUAL\n\n");
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
  if (NA_TYPE(rb_a) != NA_DFLOAT)
    rb_a = na_change_type(rb_a, NA_DFLOAT);
  a = NA_PTR_TYPE(rb_a, doublereal*);
  tol = NUM2DBL(rb_tol);
  if (!NA_IsNArray(rb_y))
    rb_raise(rb_eArgError, "y (7th argument) must be NArray");
  if (NA_RANK(rb_y) != 1)
    rb_raise(rb_eArgError, "rank of y (7th argument) must be %d", 1);
  if (NA_SHAPE0(rb_y) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of y must be the same as shape 0 of a");
  if (NA_TYPE(rb_y) != NA_DFLOAT)
    rb_y = na_change_type(rb_y, NA_DFLOAT);
  y = NA_PTR_TYPE(rb_y, doublereal*);
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
  if (NA_TYPE(rb_c) != NA_DFLOAT)
    rb_c = na_change_type(rb_c, NA_DFLOAT);
  c = NA_PTR_TYPE(rb_c, doublereal*);
  if (!NA_IsNArray(rb_b))
    rb_raise(rb_eArgError, "b (3th argument) must be NArray");
  if (NA_RANK(rb_b) != 1)
    rb_raise(rb_eArgError, "rank of b (3th argument) must be %d", 1);
  if (NA_SHAPE0(rb_b) != (n-1))
    rb_raise(rb_eRuntimeError, "shape 0 of b must be %d", n-1);
  if (NA_TYPE(rb_b) != NA_DFLOAT)
    rb_b = na_change_type(rb_b, NA_DFLOAT);
  b = NA_PTR_TYPE(rb_b, doublereal*);
  if (!NA_IsNArray(rb_d))
    rb_raise(rb_eArgError, "d (5th argument) must be NArray");
  if (NA_RANK(rb_d) != 1)
    rb_raise(rb_eArgError, "rank of d (5th argument) must be %d", 1);
  if (NA_SHAPE0(rb_d) != (n-2))
    rb_raise(rb_eRuntimeError, "shape 0 of d must be %d", n-2);
  if (NA_TYPE(rb_d) != NA_DFLOAT)
    rb_d = na_change_type(rb_d, NA_DFLOAT);
  d = NA_PTR_TYPE(rb_d, doublereal*);
  {
    int shape[1];
    shape[0] = n;
    rb_y_out__ = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  y_out__ = NA_PTR_TYPE(rb_y_out__, doublereal*);
  MEMCPY(y_out__, y, doublereal, NA_TOTAL(rb_y));
  rb_y = rb_y_out__;
  y = y_out__;

  dlagts_(&job, &n, a, b, c, d, in, y, &tol, &info);

  rb_info = INT2NUM(info);
  rb_tol = rb_float_new((double)tol);
  return rb_ary_new3(3, rb_info, rb_y, rb_tol);
}

void
init_lapack_dlagts(VALUE mLapack){
  rb_define_module_function(mLapack, "dlagts", rb_dlagts, -1);
}

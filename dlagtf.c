#include "rb_lapack.h"

extern VOID dlagtf_(integer *n, doublereal *a, doublereal *lambda, doublereal *b, doublereal *c, doublereal *tol, doublereal *d, integer *in, integer *info);

static VALUE
rb_dlagtf(int argc, VALUE *argv, VALUE self){
  VALUE rb_a;
  doublereal *a; 
  VALUE rb_lambda;
  doublereal lambda; 
  VALUE rb_b;
  doublereal *b; 
  VALUE rb_c;
  doublereal *c; 
  VALUE rb_tol;
  doublereal tol; 
  VALUE rb_d;
  doublereal *d; 
  VALUE rb_in;
  integer *in; 
  VALUE rb_info;
  integer info; 
  VALUE rb_a_out__;
  doublereal *a_out__;
  VALUE rb_b_out__;
  doublereal *b_out__;
  VALUE rb_c_out__;
  doublereal *c_out__;

  integer n;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  d, in, info, a, b, c = NumRu::Lapack.dlagtf( a, lambda, b, c, tol)\n    or\n  NumRu::Lapack.dlagtf  # print help\n\n\nFORTRAN MANUAL\n\n");
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
  if (NA_TYPE(rb_a) != NA_DFLOAT)
    rb_a = na_change_type(rb_a, NA_DFLOAT);
  a = NA_PTR_TYPE(rb_a, doublereal*);
  tol = NUM2DBL(rb_tol);
  lambda = NUM2DBL(rb_lambda);
  if (!NA_IsNArray(rb_b))
    rb_raise(rb_eArgError, "b (3th argument) must be NArray");
  if (NA_RANK(rb_b) != 1)
    rb_raise(rb_eArgError, "rank of b (3th argument) must be %d", 1);
  if (NA_SHAPE0(rb_b) != (n-1))
    rb_raise(rb_eRuntimeError, "shape 0 of b must be %d", n-1);
  if (NA_TYPE(rb_b) != NA_DFLOAT)
    rb_b = na_change_type(rb_b, NA_DFLOAT);
  b = NA_PTR_TYPE(rb_b, doublereal*);
  if (!NA_IsNArray(rb_c))
    rb_raise(rb_eArgError, "c (4th argument) must be NArray");
  if (NA_RANK(rb_c) != 1)
    rb_raise(rb_eArgError, "rank of c (4th argument) must be %d", 1);
  if (NA_SHAPE0(rb_c) != (n-1))
    rb_raise(rb_eRuntimeError, "shape 0 of c must be %d", n-1);
  if (NA_TYPE(rb_c) != NA_DFLOAT)
    rb_c = na_change_type(rb_c, NA_DFLOAT);
  c = NA_PTR_TYPE(rb_c, doublereal*);
  {
    int shape[1];
    shape[0] = n-2;
    rb_d = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  d = NA_PTR_TYPE(rb_d, doublereal*);
  {
    int shape[1];
    shape[0] = n;
    rb_in = na_make_object(NA_LINT, 1, shape, cNArray);
  }
  in = NA_PTR_TYPE(rb_in, integer*);
  {
    int shape[1];
    shape[0] = n;
    rb_a_out__ = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  a_out__ = NA_PTR_TYPE(rb_a_out__, doublereal*);
  MEMCPY(a_out__, a, doublereal, NA_TOTAL(rb_a));
  rb_a = rb_a_out__;
  a = a_out__;
  {
    int shape[1];
    shape[0] = n-1;
    rb_b_out__ = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  b_out__ = NA_PTR_TYPE(rb_b_out__, doublereal*);
  MEMCPY(b_out__, b, doublereal, NA_TOTAL(rb_b));
  rb_b = rb_b_out__;
  b = b_out__;
  {
    int shape[1];
    shape[0] = n-1;
    rb_c_out__ = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  c_out__ = NA_PTR_TYPE(rb_c_out__, doublereal*);
  MEMCPY(c_out__, c, doublereal, NA_TOTAL(rb_c));
  rb_c = rb_c_out__;
  c = c_out__;

  dlagtf_(&n, a, &lambda, b, c, &tol, d, in, &info);

  rb_info = INT2NUM(info);
  return rb_ary_new3(6, rb_d, rb_in, rb_info, rb_a, rb_b, rb_c);
}

void
init_lapack_dlagtf(VALUE mLapack){
  rb_define_module_function(mLapack, "dlagtf", rb_dlagtf, -1);
}

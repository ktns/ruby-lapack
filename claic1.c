#include "rb_lapack.h"

extern VOID claic1_(integer *job, integer *j, complex *x, real *sest, complex *w, complex *gamma, real *sestpr, complex *s, complex *c);

static VALUE
rb_claic1(int argc, VALUE *argv, VALUE self){
  VALUE rb_job;
  integer job; 
  VALUE rb_x;
  complex *x; 
  VALUE rb_sest;
  real sest; 
  VALUE rb_w;
  complex *w; 
  VALUE rb_gamma;
  complex gamma; 
  VALUE rb_sestpr;
  real sestpr; 
  VALUE rb_s;
  complex s; 
  VALUE rb_c;
  complex c; 

  integer j;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  sestpr, s, c = NumRu::Lapack.claic1( job, x, sest, w, gamma)\n    or\n  NumRu::Lapack.claic1  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 5)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 5)", argc);
  rb_job = argv[0];
  rb_x = argv[1];
  rb_sest = argv[2];
  rb_w = argv[3];
  rb_gamma = argv[4];

  if (!NA_IsNArray(rb_w))
    rb_raise(rb_eArgError, "w (4th argument) must be NArray");
  if (NA_RANK(rb_w) != 1)
    rb_raise(rb_eArgError, "rank of w (4th argument) must be %d", 1);
  j = NA_SHAPE0(rb_w);
  if (NA_TYPE(rb_w) != NA_SCOMPLEX)
    rb_w = na_change_type(rb_w, NA_SCOMPLEX);
  w = NA_PTR_TYPE(rb_w, complex*);
  if (!NA_IsNArray(rb_x))
    rb_raise(rb_eArgError, "x (2th argument) must be NArray");
  if (NA_RANK(rb_x) != 1)
    rb_raise(rb_eArgError, "rank of x (2th argument) must be %d", 1);
  if (NA_SHAPE0(rb_x) != j)
    rb_raise(rb_eRuntimeError, "shape 0 of x must be the same as shape 0 of w");
  if (NA_TYPE(rb_x) != NA_SCOMPLEX)
    rb_x = na_change_type(rb_x, NA_SCOMPLEX);
  x = NA_PTR_TYPE(rb_x, complex*);
  gamma.r = (real)NUM2DBL(rb_funcall(rb_gamma, rb_intern("real"), 0));
  gamma.i = (real)NUM2DBL(rb_funcall(rb_gamma, rb_intern("imag"), 0));
  job = NUM2INT(rb_job);
  sest = (real)NUM2DBL(rb_sest);

  claic1_(&job, &j, x, &sest, w, &gamma, &sestpr, &s, &c);

  rb_sestpr = rb_float_new((double)sestpr);
  rb_s = rb_funcall(rb_gv_get("Complex"), rb_intern("new"), 2, rb_float_new((double)(s.r)), rb_float_new((double)(s.i)));
  rb_c = rb_funcall(rb_gv_get("Complex"), rb_intern("new"), 2, rb_float_new((double)(c.r)), rb_float_new((double)(c.i)));
  return rb_ary_new3(3, rb_sestpr, rb_s, rb_c);
}

void
init_lapack_claic1(VALUE mLapack){
  rb_define_module_function(mLapack, "claic1", rb_claic1, -1);
}

#include "rb_lapack.h"

extern VOID dlaic1_(integer *job, integer *j, doublereal *x, doublereal *sest, doublereal *w, doublereal *gamma, doublereal *sestpr, doublereal *s, doublereal *c);

static VALUE
rb_dlaic1(int argc, VALUE *argv, VALUE self){
  VALUE rb_job;
  integer job; 
  VALUE rb_x;
  doublereal *x; 
  VALUE rb_sest;
  doublereal sest; 
  VALUE rb_w;
  doublereal *w; 
  VALUE rb_gamma;
  doublereal gamma; 
  VALUE rb_sestpr;
  doublereal sestpr; 
  VALUE rb_s;
  doublereal s; 
  VALUE rb_c;
  doublereal c; 

  integer j;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  sestpr, s, c = NumRu::Lapack.dlaic1( job, x, sest, w, gamma)\n    or\n  NumRu::Lapack.dlaic1  # print help\n\n\nFORTRAN MANUAL\n\n");
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
  if (NA_TYPE(rb_w) != NA_DFLOAT)
    rb_w = na_change_type(rb_w, NA_DFLOAT);
  w = NA_PTR_TYPE(rb_w, doublereal*);
  if (!NA_IsNArray(rb_x))
    rb_raise(rb_eArgError, "x (2th argument) must be NArray");
  if (NA_RANK(rb_x) != 1)
    rb_raise(rb_eArgError, "rank of x (2th argument) must be %d", 1);
  if (NA_SHAPE0(rb_x) != j)
    rb_raise(rb_eRuntimeError, "shape 0 of x must be the same as shape 0 of w");
  if (NA_TYPE(rb_x) != NA_DFLOAT)
    rb_x = na_change_type(rb_x, NA_DFLOAT);
  x = NA_PTR_TYPE(rb_x, doublereal*);
  gamma = NUM2DBL(rb_gamma);
  job = NUM2INT(rb_job);
  sest = NUM2DBL(rb_sest);

  dlaic1_(&job, &j, x, &sest, w, &gamma, &sestpr, &s, &c);

  rb_sestpr = rb_float_new((double)sestpr);
  rb_s = rb_float_new((double)s);
  rb_c = rb_float_new((double)c);
  return rb_ary_new3(3, rb_sestpr, rb_s, rb_c);
}

void
init_lapack_dlaic1(VALUE mLapack){
  rb_define_module_function(mLapack, "dlaic1", rb_dlaic1, -1);
}

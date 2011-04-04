#include "rb_lapack.h"

extern VOID dlag2_(doublereal *a, integer *lda, doublereal *b, integer *ldb, doublereal *safmin, doublereal *scale1, doublereal *scale2, doublereal *wr1, doublereal *wr2, doublereal *wi);

static VALUE
rb_dlag2(int argc, VALUE *argv, VALUE self){
  VALUE rb_a;
  doublereal *a; 
  VALUE rb_b;
  doublereal *b; 
  VALUE rb_safmin;
  doublereal safmin; 
  VALUE rb_scale1;
  doublereal scale1; 
  VALUE rb_scale2;
  doublereal scale2; 
  VALUE rb_wr1;
  doublereal wr1; 
  VALUE rb_wr2;
  doublereal wr2; 
  VALUE rb_wi;
  doublereal wi; 

  integer lda;
  integer ldb;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  scale1, scale2, wr1, wr2, wi = NumRu::Lapack.dlag2( a, b, safmin)\n    or\n  NumRu::Lapack.dlag2  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 3)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 3)", argc);
  rb_a = argv[0];
  rb_b = argv[1];
  rb_safmin = argv[2];

  if (!NA_IsNArray(rb_a))
    rb_raise(rb_eArgError, "a (1th argument) must be NArray");
  if (NA_RANK(rb_a) != 2)
    rb_raise(rb_eArgError, "rank of a (1th argument) must be %d", 2);
  if (NA_SHAPE1(rb_a) != (2))
    rb_raise(rb_eRuntimeError, "shape 1 of a must be %d", 2);
  lda = NA_SHAPE0(rb_a);
  if (NA_TYPE(rb_a) != NA_DFLOAT)
    rb_a = na_change_type(rb_a, NA_DFLOAT);
  a = NA_PTR_TYPE(rb_a, doublereal*);
  if (!NA_IsNArray(rb_b))
    rb_raise(rb_eArgError, "b (2th argument) must be NArray");
  if (NA_RANK(rb_b) != 2)
    rb_raise(rb_eArgError, "rank of b (2th argument) must be %d", 2);
  if (NA_SHAPE1(rb_b) != (2))
    rb_raise(rb_eRuntimeError, "shape 1 of b must be %d", 2);
  ldb = NA_SHAPE0(rb_b);
  if (NA_TYPE(rb_b) != NA_DFLOAT)
    rb_b = na_change_type(rb_b, NA_DFLOAT);
  b = NA_PTR_TYPE(rb_b, doublereal*);
  safmin = NUM2DBL(rb_safmin);

  dlag2_(a, &lda, b, &ldb, &safmin, &scale1, &scale2, &wr1, &wr2, &wi);

  rb_scale1 = rb_float_new((double)scale1);
  rb_scale2 = rb_float_new((double)scale2);
  rb_wr1 = rb_float_new((double)wr1);
  rb_wr2 = rb_float_new((double)wr2);
  rb_wi = rb_float_new((double)wi);
  return rb_ary_new3(5, rb_scale1, rb_scale2, rb_wr1, rb_wr2, rb_wi);
}

void
init_lapack_dlag2(VALUE mLapack){
  rb_define_module_function(mLapack, "dlag2", rb_dlag2, -1);
}

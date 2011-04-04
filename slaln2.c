#include "rb_lapack.h"

extern VOID slaln2_(logical *ltrans, integer *na, integer *nw, real *smin, real *ca, real *a, integer *lda, real *d1, real *d2, real *b, integer *ldb, real *wr, real *wi, real *x, integer *ldx, real *scale, real *xnorm, integer *info);

static VALUE
rb_slaln2(int argc, VALUE *argv, VALUE self){
  VALUE rb_ltrans;
  logical ltrans; 
  VALUE rb_smin;
  real smin; 
  VALUE rb_ca;
  real ca; 
  VALUE rb_a;
  real *a; 
  VALUE rb_d1;
  real d1; 
  VALUE rb_d2;
  real d2; 
  VALUE rb_b;
  real *b; 
  VALUE rb_wr;
  real wr; 
  VALUE rb_wi;
  real wi; 
  VALUE rb_x;
  real *x; 
  VALUE rb_scale;
  real scale; 
  VALUE rb_xnorm;
  real xnorm; 
  VALUE rb_info;
  integer info; 

  integer lda;
  integer na;
  integer ldb;
  integer nw;
  integer ldx;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  x, scale, xnorm, info = NumRu::Lapack.slaln2( ltrans, smin, ca, a, d1, d2, b, wr, wi)\n    or\n  NumRu::Lapack.slaln2  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 9)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 9)", argc);
  rb_ltrans = argv[0];
  rb_smin = argv[1];
  rb_ca = argv[2];
  rb_a = argv[3];
  rb_d1 = argv[4];
  rb_d2 = argv[5];
  rb_b = argv[6];
  rb_wr = argv[7];
  rb_wi = argv[8];

  smin = (real)NUM2DBL(rb_smin);
  if (!NA_IsNArray(rb_a))
    rb_raise(rb_eArgError, "a (4th argument) must be NArray");
  if (NA_RANK(rb_a) != 2)
    rb_raise(rb_eArgError, "rank of a (4th argument) must be %d", 2);
  na = NA_SHAPE1(rb_a);
  lda = NA_SHAPE0(rb_a);
  if (NA_TYPE(rb_a) != NA_SFLOAT)
    rb_a = na_change_type(rb_a, NA_SFLOAT);
  a = NA_PTR_TYPE(rb_a, real*);
  if (!NA_IsNArray(rb_b))
    rb_raise(rb_eArgError, "b (7th argument) must be NArray");
  if (NA_RANK(rb_b) != 2)
    rb_raise(rb_eArgError, "rank of b (7th argument) must be %d", 2);
  nw = NA_SHAPE1(rb_b);
  ldb = NA_SHAPE0(rb_b);
  if (NA_TYPE(rb_b) != NA_SFLOAT)
    rb_b = na_change_type(rb_b, NA_SFLOAT);
  b = NA_PTR_TYPE(rb_b, real*);
  d1 = (real)NUM2DBL(rb_d1);
  d2 = (real)NUM2DBL(rb_d2);
  ca = (real)NUM2DBL(rb_ca);
  ltrans = (rb_ltrans == Qtrue);
  wi = (real)NUM2DBL(rb_wi);
  wr = (real)NUM2DBL(rb_wr);
  ldx = na;
  {
    int shape[2];
    shape[0] = ldx;
    shape[1] = nw;
    rb_x = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  x = NA_PTR_TYPE(rb_x, real*);

  slaln2_(&ltrans, &na, &nw, &smin, &ca, a, &lda, &d1, &d2, b, &ldb, &wr, &wi, x, &ldx, &scale, &xnorm, &info);

  rb_scale = rb_float_new((double)scale);
  rb_xnorm = rb_float_new((double)xnorm);
  rb_info = INT2NUM(info);
  return rb_ary_new3(4, rb_x, rb_scale, rb_xnorm, rb_info);
}

void
init_lapack_slaln2(VALUE mLapack){
  rb_define_module_function(mLapack, "slaln2", rb_slaln2, -1);
}

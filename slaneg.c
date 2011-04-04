#include "rb_lapack.h"

extern integer slaneg_(integer *n, real *d, real *lld, real *sigma, real *pivmin, integer *r);

static VALUE
rb_slaneg(int argc, VALUE *argv, VALUE self){
  VALUE rb_d;
  real *d; 
  VALUE rb_lld;
  real *lld; 
  VALUE rb_sigma;
  real sigma; 
  VALUE rb_pivmin;
  real pivmin; 
  VALUE rb_r;
  integer r; 
  VALUE rb___out__;
  integer __out__; 

  integer n;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  __out__ = NumRu::Lapack.slaneg( d, lld, sigma, pivmin, r)\n    or\n  NumRu::Lapack.slaneg  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 5)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 5)", argc);
  rb_d = argv[0];
  rb_lld = argv[1];
  rb_sigma = argv[2];
  rb_pivmin = argv[3];
  rb_r = argv[4];

  pivmin = (real)NUM2DBL(rb_pivmin);
  if (!NA_IsNArray(rb_d))
    rb_raise(rb_eArgError, "d (1th argument) must be NArray");
  if (NA_RANK(rb_d) != 1)
    rb_raise(rb_eArgError, "rank of d (1th argument) must be %d", 1);
  n = NA_SHAPE0(rb_d);
  if (NA_TYPE(rb_d) != NA_SFLOAT)
    rb_d = na_change_type(rb_d, NA_SFLOAT);
  d = NA_PTR_TYPE(rb_d, real*);
  sigma = (real)NUM2DBL(rb_sigma);
  r = NUM2INT(rb_r);
  if (!NA_IsNArray(rb_lld))
    rb_raise(rb_eArgError, "lld (2th argument) must be NArray");
  if (NA_RANK(rb_lld) != 1)
    rb_raise(rb_eArgError, "rank of lld (2th argument) must be %d", 1);
  if (NA_SHAPE0(rb_lld) != (n-1))
    rb_raise(rb_eRuntimeError, "shape 0 of lld must be %d", n-1);
  if (NA_TYPE(rb_lld) != NA_SFLOAT)
    rb_lld = na_change_type(rb_lld, NA_SFLOAT);
  lld = NA_PTR_TYPE(rb_lld, real*);

  __out__ = slaneg_(&n, d, lld, &sigma, &pivmin, &r);

  rb___out__ = INT2NUM(__out__);
  return rb___out__;
}

void
init_lapack_slaneg(VALUE mLapack){
  rb_define_module_function(mLapack, "slaneg", rb_slaneg, -1);
}

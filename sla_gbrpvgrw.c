#include "rb_lapack.h"

extern real sla_gbrpvgrw_(integer *n, integer *kl, integer *ku, integer *ncols, real *ab, integer *ldab, real *afb, integer *ldafb);

static VALUE
rb_sla_gbrpvgrw(int argc, VALUE *argv, VALUE self){
  VALUE rb_kl;
  integer kl; 
  VALUE rb_ku;
  integer ku; 
  VALUE rb_ncols;
  integer ncols; 
  VALUE rb_ab;
  real *ab; 
  VALUE rb_afb;
  real *afb; 
  VALUE rb___out__;
  real __out__; 

  integer ldab;
  integer n;
  integer ldafb;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  __out__ = NumRu::Lapack.sla_gbrpvgrw( kl, ku, ncols, ab, afb)\n    or\n  NumRu::Lapack.sla_gbrpvgrw  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 5)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 5)", argc);
  rb_kl = argv[0];
  rb_ku = argv[1];
  rb_ncols = argv[2];
  rb_ab = argv[3];
  rb_afb = argv[4];

  if (!NA_IsNArray(rb_ab))
    rb_raise(rb_eArgError, "ab (4th argument) must be NArray");
  if (NA_RANK(rb_ab) != 2)
    rb_raise(rb_eArgError, "rank of ab (4th argument) must be %d", 2);
  n = NA_SHAPE1(rb_ab);
  ldab = NA_SHAPE0(rb_ab);
  if (NA_TYPE(rb_ab) != NA_SFLOAT)
    rb_ab = na_change_type(rb_ab, NA_SFLOAT);
  ab = NA_PTR_TYPE(rb_ab, real*);
  kl = NUM2INT(rb_kl);
  ncols = NUM2INT(rb_ncols);
  if (!NA_IsNArray(rb_afb))
    rb_raise(rb_eArgError, "afb (5th argument) must be NArray");
  if (NA_RANK(rb_afb) != 2)
    rb_raise(rb_eArgError, "rank of afb (5th argument) must be %d", 2);
  if (NA_SHAPE1(rb_afb) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of afb must be the same as shape 1 of ab");
  ldafb = NA_SHAPE0(rb_afb);
  if (NA_TYPE(rb_afb) != NA_SFLOAT)
    rb_afb = na_change_type(rb_afb, NA_SFLOAT);
  afb = NA_PTR_TYPE(rb_afb, real*);
  ku = NUM2INT(rb_ku);

  __out__ = sla_gbrpvgrw_(&n, &kl, &ku, &ncols, ab, &ldab, afb, &ldafb);

  rb___out__ = rb_float_new((double)__out__);
  return rb___out__;
}

void
init_lapack_sla_gbrpvgrw(VALUE mLapack){
  rb_define_module_function(mLapack, "sla_gbrpvgrw", rb_sla_gbrpvgrw, -1);
}

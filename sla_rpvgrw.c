#include "rb_lapack.h"

extern real sla_rpvgrw_(integer *n, integer *ncols, real *a, integer *lda, real *af, integer *ldaf);

static VALUE
rb_sla_rpvgrw(int argc, VALUE *argv, VALUE self){
  VALUE rb_ncols;
  integer ncols; 
  VALUE rb_a;
  real *a; 
  VALUE rb_af;
  real *af; 
  VALUE rb___out__;
  real __out__; 

  integer lda;
  integer n;
  integer ldaf;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  __out__ = NumRu::Lapack.sla_rpvgrw( ncols, a, af)\n    or\n  NumRu::Lapack.sla_rpvgrw  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 3)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 3)", argc);
  rb_ncols = argv[0];
  rb_a = argv[1];
  rb_af = argv[2];

  if (!NA_IsNArray(rb_a))
    rb_raise(rb_eArgError, "a (2th argument) must be NArray");
  if (NA_RANK(rb_a) != 2)
    rb_raise(rb_eArgError, "rank of a (2th argument) must be %d", 2);
  n = NA_SHAPE1(rb_a);
  lda = NA_SHAPE0(rb_a);
  if (NA_TYPE(rb_a) != NA_SFLOAT)
    rb_a = na_change_type(rb_a, NA_SFLOAT);
  a = NA_PTR_TYPE(rb_a, real*);
  ncols = NUM2INT(rb_ncols);
  if (!NA_IsNArray(rb_af))
    rb_raise(rb_eArgError, "af (3th argument) must be NArray");
  if (NA_RANK(rb_af) != 2)
    rb_raise(rb_eArgError, "rank of af (3th argument) must be %d", 2);
  if (NA_SHAPE1(rb_af) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of af must be the same as shape 1 of a");
  ldaf = NA_SHAPE0(rb_af);
  if (NA_TYPE(rb_af) != NA_SFLOAT)
    rb_af = na_change_type(rb_af, NA_SFLOAT);
  af = NA_PTR_TYPE(rb_af, real*);

  __out__ = sla_rpvgrw_(&n, &ncols, a, &lda, af, &ldaf);

  rb___out__ = rb_float_new((double)__out__);
  return rb___out__;
}

void
init_lapack_sla_rpvgrw(VALUE mLapack){
  rb_define_module_function(mLapack, "sla_rpvgrw", rb_sla_rpvgrw, -1);
}

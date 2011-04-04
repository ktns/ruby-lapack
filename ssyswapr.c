#include "rb_lapack.h"

extern VOID ssyswapr_(char *uplo, integer *n, real *a, integer *i1, integer *i2);

static VALUE
rb_ssyswapr(int argc, VALUE *argv, VALUE self){
  VALUE rb_uplo;
  char uplo; 
  VALUE rb_a;
  real *a; 
  VALUE rb_i1;
  integer i1; 
  VALUE rb_i2;
  integer i2; 
  VALUE rb_a_out__;
  real *a_out__;

  integer lda;
  integer n;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  a = NumRu::Lapack.ssyswapr( uplo, a, i1, i2)\n    or\n  NumRu::Lapack.ssyswapr  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 4)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 4)", argc);
  rb_uplo = argv[0];
  rb_a = argv[1];
  rb_i1 = argv[2];
  rb_i2 = argv[3];

  if (!NA_IsNArray(rb_a))
    rb_raise(rb_eArgError, "a (2th argument) must be NArray");
  if (NA_RANK(rb_a) != 2)
    rb_raise(rb_eArgError, "rank of a (2th argument) must be %d", 2);
  n = NA_SHAPE1(rb_a);
  lda = NA_SHAPE0(rb_a);
  if (NA_TYPE(rb_a) != NA_SFLOAT)
    rb_a = na_change_type(rb_a, NA_SFLOAT);
  a = NA_PTR_TYPE(rb_a, real*);
  i1 = NUM2INT(rb_i1);
  i2 = NUM2INT(rb_i2);
  uplo = StringValueCStr(rb_uplo)[0];
  {
    int shape[2];
    shape[0] = lda;
    shape[1] = n;
    rb_a_out__ = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  a_out__ = NA_PTR_TYPE(rb_a_out__, real*);
  MEMCPY(a_out__, a, real, NA_TOTAL(rb_a));
  rb_a = rb_a_out__;
  a = a_out__;

  ssyswapr_(&uplo, &n, a, &i1, &i2);

  return rb_a;
}

void
init_lapack_ssyswapr(VALUE mLapack){
  rb_define_module_function(mLapack, "ssyswapr", rb_ssyswapr, -1);
}

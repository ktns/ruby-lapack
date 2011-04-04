#include "rb_lapack.h"

extern VOID strttf_(char *transr, char *uplo, integer *n, real *a, integer *lda, real *arf, integer *info);

static VALUE
rb_strttf(int argc, VALUE *argv, VALUE self){
  VALUE rb_transr;
  char transr; 
  VALUE rb_uplo;
  char uplo; 
  VALUE rb_a;
  real *a; 
  VALUE rb_arf;
  real *arf; 
  VALUE rb_info;
  integer info; 

  integer lda;
  integer n;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  arf, info = NumRu::Lapack.strttf( transr, uplo, a)\n    or\n  NumRu::Lapack.strttf  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 3)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 3)", argc);
  rb_transr = argv[0];
  rb_uplo = argv[1];
  rb_a = argv[2];

  if (!NA_IsNArray(rb_a))
    rb_raise(rb_eArgError, "a (3th argument) must be NArray");
  if (NA_RANK(rb_a) != 2)
    rb_raise(rb_eArgError, "rank of a (3th argument) must be %d", 2);
  n = NA_SHAPE1(rb_a);
  lda = NA_SHAPE0(rb_a);
  if (NA_TYPE(rb_a) != NA_SFLOAT)
    rb_a = na_change_type(rb_a, NA_SFLOAT);
  a = NA_PTR_TYPE(rb_a, real*);
  uplo = StringValueCStr(rb_uplo)[0];
  transr = StringValueCStr(rb_transr)[0];
  {
    int shape[1];
    shape[0] = n*(n+1)/2;
    rb_arf = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  arf = NA_PTR_TYPE(rb_arf, real*);

  strttf_(&transr, &uplo, &n, a, &lda, arf, &info);

  rb_info = INT2NUM(info);
  return rb_ary_new3(2, rb_arf, rb_info);
}

void
init_lapack_strttf(VALUE mLapack){
  rb_define_module_function(mLapack, "strttf", rb_strttf, -1);
}

#include "rb_lapack.h"

extern VOID strttp_(char *uplo, integer *n, real *a, integer *lda, real *ap, integer *info);

static VALUE
rb_strttp(int argc, VALUE *argv, VALUE self){
  VALUE rb_uplo;
  char uplo; 
  VALUE rb_a;
  real *a; 
  VALUE rb_ap;
  real *ap; 
  VALUE rb_info;
  integer info; 

  integer lda;
  integer n;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  ap, info = NumRu::Lapack.strttp( uplo, a)\n    or\n  NumRu::Lapack.strttp  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 2)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 2)", argc);
  rb_uplo = argv[0];
  rb_a = argv[1];

  if (!NA_IsNArray(rb_a))
    rb_raise(rb_eArgError, "a (2th argument) must be NArray");
  if (NA_RANK(rb_a) != 2)
    rb_raise(rb_eArgError, "rank of a (2th argument) must be %d", 2);
  n = NA_SHAPE1(rb_a);
  lda = NA_SHAPE0(rb_a);
  if (NA_TYPE(rb_a) != NA_SFLOAT)
    rb_a = na_change_type(rb_a, NA_SFLOAT);
  a = NA_PTR_TYPE(rb_a, real*);
  uplo = StringValueCStr(rb_uplo)[0];
  {
    int shape[1];
    shape[0] = n*(n+1)/2;
    rb_ap = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  ap = NA_PTR_TYPE(rb_ap, real*);

  strttp_(&uplo, &n, a, &lda, ap, &info);

  rb_info = INT2NUM(info);
  return rb_ary_new3(2, rb_ap, rb_info);
}

void
init_lapack_strttp(VALUE mLapack){
  rb_define_module_function(mLapack, "strttp", rb_strttp, -1);
}
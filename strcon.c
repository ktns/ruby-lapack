#include "rb_lapack.h"

extern VOID strcon_(char *norm, char *uplo, char *diag, integer *n, real *a, integer *lda, real *rcond, real *work, integer *iwork, integer *info);

static VALUE
rb_strcon(int argc, VALUE *argv, VALUE self){
  VALUE rb_norm;
  char norm; 
  VALUE rb_uplo;
  char uplo; 
  VALUE rb_diag;
  char diag; 
  VALUE rb_a;
  real *a; 
  VALUE rb_rcond;
  real rcond; 
  VALUE rb_info;
  integer info; 
  real *work;
  integer *iwork;

  integer lda;
  integer n;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  rcond, info = NumRu::Lapack.strcon( norm, uplo, diag, a)\n    or\n  NumRu::Lapack.strcon  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 4)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 4)", argc);
  rb_norm = argv[0];
  rb_uplo = argv[1];
  rb_diag = argv[2];
  rb_a = argv[3];

  if (!NA_IsNArray(rb_a))
    rb_raise(rb_eArgError, "a (4th argument) must be NArray");
  if (NA_RANK(rb_a) != 2)
    rb_raise(rb_eArgError, "rank of a (4th argument) must be %d", 2);
  n = NA_SHAPE1(rb_a);
  lda = NA_SHAPE0(rb_a);
  if (NA_TYPE(rb_a) != NA_SFLOAT)
    rb_a = na_change_type(rb_a, NA_SFLOAT);
  a = NA_PTR_TYPE(rb_a, real*);
  diag = StringValueCStr(rb_diag)[0];
  norm = StringValueCStr(rb_norm)[0];
  uplo = StringValueCStr(rb_uplo)[0];
  work = ALLOC_N(real, (3*n));
  iwork = ALLOC_N(integer, (n));

  strcon_(&norm, &uplo, &diag, &n, a, &lda, &rcond, work, iwork, &info);

  free(work);
  free(iwork);
  rb_rcond = rb_float_new((double)rcond);
  rb_info = INT2NUM(info);
  return rb_ary_new3(2, rb_rcond, rb_info);
}

void
init_lapack_strcon(VALUE mLapack){
  rb_define_module_function(mLapack, "strcon", rb_strcon, -1);
}

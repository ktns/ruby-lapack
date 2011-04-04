#include "rb_lapack.h"

extern VOID sgecon_(char *norm, integer *n, real *a, integer *lda, real *anorm, real *rcond, real *work, integer *iwork, integer *info);

static VALUE
rb_sgecon(int argc, VALUE *argv, VALUE self){
  VALUE rb_norm;
  char norm; 
  VALUE rb_a;
  real *a; 
  VALUE rb_anorm;
  real anorm; 
  VALUE rb_rcond;
  real rcond; 
  VALUE rb_info;
  integer info; 
  real *work;
  integer *iwork;

  integer lda;
  integer n;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  rcond, info = NumRu::Lapack.sgecon( norm, a, anorm)\n    or\n  NumRu::Lapack.sgecon  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 3)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 3)", argc);
  rb_norm = argv[0];
  rb_a = argv[1];
  rb_anorm = argv[2];

  anorm = (real)NUM2DBL(rb_anorm);
  if (!NA_IsNArray(rb_a))
    rb_raise(rb_eArgError, "a (2th argument) must be NArray");
  if (NA_RANK(rb_a) != 2)
    rb_raise(rb_eArgError, "rank of a (2th argument) must be %d", 2);
  n = NA_SHAPE1(rb_a);
  lda = NA_SHAPE0(rb_a);
  if (NA_TYPE(rb_a) != NA_SFLOAT)
    rb_a = na_change_type(rb_a, NA_SFLOAT);
  a = NA_PTR_TYPE(rb_a, real*);
  norm = StringValueCStr(rb_norm)[0];
  work = ALLOC_N(real, (4*n));
  iwork = ALLOC_N(integer, (n));

  sgecon_(&norm, &n, a, &lda, &anorm, &rcond, work, iwork, &info);

  free(work);
  free(iwork);
  rb_rcond = rb_float_new((double)rcond);
  rb_info = INT2NUM(info);
  return rb_ary_new3(2, rb_rcond, rb_info);
}

void
init_lapack_sgecon(VALUE mLapack){
  rb_define_module_function(mLapack, "sgecon", rb_sgecon, -1);
}

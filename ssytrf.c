#include "rb_lapack.h"

extern VOID ssytrf_(char *uplo, integer *n, real *a, integer *lda, integer *ipiv, real *work, integer *lwork, integer *info);

static VALUE
rb_ssytrf(int argc, VALUE *argv, VALUE self){
  VALUE rb_uplo;
  char uplo; 
  VALUE rb_a;
  real *a; 
  VALUE rb_lwork;
  integer lwork; 
  VALUE rb_ipiv;
  integer *ipiv; 
  VALUE rb_work;
  real *work; 
  VALUE rb_info;
  integer info; 
  VALUE rb_a_out__;
  real *a_out__;

  integer lda;
  integer n;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  ipiv, work, info, a = NumRu::Lapack.ssytrf( uplo, a, lwork)\n    or\n  NumRu::Lapack.ssytrf  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 3)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 3)", argc);
  rb_uplo = argv[0];
  rb_a = argv[1];
  rb_lwork = argv[2];

  if (!NA_IsNArray(rb_a))
    rb_raise(rb_eArgError, "a (2th argument) must be NArray");
  if (NA_RANK(rb_a) != 2)
    rb_raise(rb_eArgError, "rank of a (2th argument) must be %d", 2);
  n = NA_SHAPE1(rb_a);
  lda = NA_SHAPE0(rb_a);
  if (NA_TYPE(rb_a) != NA_SFLOAT)
    rb_a = na_change_type(rb_a, NA_SFLOAT);
  a = NA_PTR_TYPE(rb_a, real*);
  lwork = NUM2INT(rb_lwork);
  uplo = StringValueCStr(rb_uplo)[0];
  {
    int shape[1];
    shape[0] = n;
    rb_ipiv = na_make_object(NA_LINT, 1, shape, cNArray);
  }
  ipiv = NA_PTR_TYPE(rb_ipiv, integer*);
  {
    int shape[1];
    shape[0] = MAX(1,lwork);
    rb_work = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  work = NA_PTR_TYPE(rb_work, real*);
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

  ssytrf_(&uplo, &n, a, &lda, ipiv, work, &lwork, &info);

  rb_info = INT2NUM(info);
  return rb_ary_new3(4, rb_ipiv, rb_work, rb_info, rb_a);
}

void
init_lapack_ssytrf(VALUE mLapack){
  rb_define_module_function(mLapack, "ssytrf", rb_ssytrf, -1);
}

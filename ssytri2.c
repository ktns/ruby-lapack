#include "rb_lapack.h"

extern VOID ssytri2_(char *uplo, integer *n, real *a, integer *lda, integer *ipiv, real *work, integer *lwork, integer *info);

static VALUE
rb_ssytri2(int argc, VALUE *argv, VALUE self){
  VALUE rb_uplo;
  char uplo; 
  VALUE rb_a;
  real *a; 
  VALUE rb_ipiv;
  integer *ipiv; 
  VALUE rb_work;
  real *work; 
  VALUE rb_info;
  integer info; 
  VALUE rb_a_out__;
  real *a_out__;
  VALUE rb_work_out__;
  real *work_out__;
  integer c__1;
  integer nb;
  integer c__m1;

  integer lda;
  integer n;
  integer lwork;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  info, a, work = NumRu::Lapack.ssytri2( uplo, a, ipiv, work)\n    or\n  NumRu::Lapack.ssytri2  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 4)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 4)", argc);
  rb_uplo = argv[0];
  rb_a = argv[1];
  rb_ipiv = argv[2];
  rb_work = argv[3];

  if (!NA_IsNArray(rb_ipiv))
    rb_raise(rb_eArgError, "ipiv (3th argument) must be NArray");
  if (NA_RANK(rb_ipiv) != 1)
    rb_raise(rb_eArgError, "rank of ipiv (3th argument) must be %d", 1);
  n = NA_SHAPE0(rb_ipiv);
  if (NA_TYPE(rb_ipiv) != NA_LINT)
    rb_ipiv = na_change_type(rb_ipiv, NA_LINT);
  ipiv = NA_PTR_TYPE(rb_ipiv, integer*);
  if (!NA_IsNArray(rb_a))
    rb_raise(rb_eArgError, "a (2th argument) must be NArray");
  if (NA_RANK(rb_a) != 2)
    rb_raise(rb_eArgError, "rank of a (2th argument) must be %d", 2);
  if (NA_SHAPE1(rb_a) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of a must be the same as shape 0 of ipiv");
  lda = NA_SHAPE0(rb_a);
  if (NA_TYPE(rb_a) != NA_SFLOAT)
    rb_a = na_change_type(rb_a, NA_SFLOAT);
  a = NA_PTR_TYPE(rb_a, real*);
  c__1 = 1;
  if (!NA_IsNArray(rb_work))
    rb_raise(rb_eArgError, "work (4th argument) must be NArray");
  if (NA_RANK(rb_work) != 1)
    rb_raise(rb_eArgError, "rank of work (4th argument) must be %d", 1);
  lwork = NA_SHAPE0(rb_work);
  if (NA_TYPE(rb_work) != NA_SFLOAT)
    rb_work = na_change_type(rb_work, NA_SFLOAT);
  work = NA_PTR_TYPE(rb_work, real*);
  uplo = StringValueCStr(rb_uplo)[0];
  c__m1 = -1;
  nb = ilaenv_(&c__1, "SSYTRF", &uplo, &n, &c__m1, &c__m1, &c__m1);
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
  {
    int shape[1];
    shape[0] = lwork;
    rb_work_out__ = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  work_out__ = NA_PTR_TYPE(rb_work_out__, real*);
  MEMCPY(work_out__, work, real, NA_TOTAL(rb_work));
  rb_work = rb_work_out__;
  work = work_out__;

  ssytri2_(&uplo, &n, a, &lda, ipiv, work, &lwork, &info);

  rb_info = INT2NUM(info);
  return rb_ary_new3(3, rb_info, rb_a, rb_work);
}

void
init_lapack_ssytri2(VALUE mLapack){
  rb_define_module_function(mLapack, "ssytri2", rb_ssytri2, -1);
}

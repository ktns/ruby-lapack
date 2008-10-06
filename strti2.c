#include "rb_lapack.h"

static VALUE
rb_strti2(int argc, VALUE *argv, VALUE self){
  VALUE rb_uplo;
  char uplo; 
  VALUE rb_diag;
  char diag; 
  VALUE rb_a;
  real *a; 
  VALUE rb_info;
  integer info; 
  VALUE rb_a_out__;
  real *a_out__;

  integer lda;
  integer n;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  info, a = NumRu::Lapack.strti2( uplo, diag, a)\n    or\n  NumRu::Lapack.strti2  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE STRTI2( UPLO, DIAG, N, A, LDA, INFO )\n\n*  Purpose\n*  =======\n*\n*  STRTI2 computes the inverse of a real upper or lower triangular\n*  matrix.\n*\n*  This is the Level 2 BLAS version of the algorithm.\n*\n\n*  Arguments\n*  =========\n*\n*  UPLO    (input) CHARACTER*1\n*          Specifies whether the matrix A is upper or lower triangular.\n*          = 'U':  Upper triangular\n*          = 'L':  Lower triangular\n*\n*  DIAG    (input) CHARACTER*1\n*          Specifies whether or not the matrix A is unit triangular.\n*          = 'N':  Non-unit triangular\n*          = 'U':  Unit triangular\n*\n*  N       (input) INTEGER\n*          The order of the matrix A.  N >= 0.\n*\n*  A       (input/output) REAL array, dimension (LDA,N)\n*          On entry, the triangular matrix A.  If UPLO = 'U', the\n*          leading n by n upper triangular part of the array A contains\n*          the upper triangular matrix, and the strictly lower\n*          triangular part of A is not referenced.  If UPLO = 'L', the\n*          leading n by n lower triangular part of the array A contains\n*          the lower triangular matrix, and the strictly upper\n*          triangular part of A is not referenced.  If DIAG = 'U', the\n*          diagonal elements of A are also not referenced and are\n*          assumed to be 1.\n*\n*          On exit, the (triangular) inverse of the original matrix, in\n*          the same storage format.\n*\n*  LDA     (input) INTEGER\n*          The leading dimension of the array A.  LDA >= max(1,N).\n*\n*  INFO    (output) INTEGER\n*          = 0: successful exit\n*          < 0: if INFO = -k, the k-th argument had an illegal value\n*\n\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 3)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 3)", argc);
  rb_uplo = argv[0];
  rb_diag = argv[1];
  rb_a = argv[2];

  uplo = StringValueCStr(rb_uplo)[0];
  diag = StringValueCStr(rb_diag)[0];
  if (!NA_IsNArray(rb_a))
    rb_raise(rb_eArgError, "a (3th argument) must be NArray");
  if (NA_RANK(rb_a) != 2)
    rb_raise(rb_eArgError, "rank of a (3th argument) must be %d", 2);
  lda = NA_SHAPE0(rb_a);
  n = NA_SHAPE1(rb_a);
  if (NA_TYPE(rb_a) != NA_SFLOAT)
    rb_a = na_change_type(rb_a, NA_SFLOAT);
  a = NA_PTR_TYPE(rb_a, real*);
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

  strti2_(&uplo, &diag, &n, a, &lda, &info);

  rb_info = INT2NUM(info);
  return rb_ary_new3(2, rb_info, rb_a);
}

void
init_lapack_strti2(VALUE mLapack){
  rb_define_module_function(mLapack, "strti2", rb_strti2, -1);
}

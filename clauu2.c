#include "rb_lapack.h"

extern VOID clauu2_(char *uplo, integer *n, complex *a, integer *lda, integer *info);

static VALUE
rb_clauu2(int argc, VALUE *argv, VALUE self){
  VALUE rb_uplo;
  char uplo; 
  VALUE rb_a;
  complex *a; 
  VALUE rb_info;
  integer info; 
  VALUE rb_a_out__;
  complex *a_out__;

  integer lda;
  integer n;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  info, a = NumRu::Lapack.clauu2( uplo, a)\n    or\n  NumRu::Lapack.clauu2  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE CLAUU2( UPLO, N, A, LDA, INFO )\n\n*  Purpose\n*  =======\n*\n*  CLAUU2 computes the product U * U' or L' * L, where the triangular\n*  factor U or L is stored in the upper or lower triangular part of\n*  the array A.\n*\n*  If UPLO = 'U' or 'u' then the upper triangle of the result is stored,\n*  overwriting the factor U in A.\n*  If UPLO = 'L' or 'l' then the lower triangle of the result is stored,\n*  overwriting the factor L in A.\n*\n*  This is the unblocked form of the algorithm, calling Level 2 BLAS.\n*\n\n*  Arguments\n*  =========\n*\n*  UPLO    (input) CHARACTER*1\n*          Specifies whether the triangular factor stored in the array A\n*          is upper or lower triangular:\n*          = 'U':  Upper triangular\n*          = 'L':  Lower triangular\n*\n*  N       (input) INTEGER\n*          The order of the triangular factor U or L.  N >= 0.\n*\n*  A       (input/output) COMPLEX array, dimension (LDA,N)\n*          On entry, the triangular factor U or L.\n*          On exit, if UPLO = 'U', the upper triangle of A is\n*          overwritten with the upper triangle of the product U * U';\n*          if UPLO = 'L', the lower triangle of A is overwritten with\n*          the lower triangle of the product L' * L.\n*\n*  LDA     (input) INTEGER\n*          The leading dimension of the array A.  LDA >= max(1,N).\n*\n*  INFO    (output) INTEGER\n*          = 0: successful exit\n*          < 0: if INFO = -k, the k-th argument had an illegal value\n*\n\n*  =====================================================================\n*\n\n");
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
  if (NA_TYPE(rb_a) != NA_SCOMPLEX)
    rb_a = na_change_type(rb_a, NA_SCOMPLEX);
  a = NA_PTR_TYPE(rb_a, complex*);
  uplo = StringValueCStr(rb_uplo)[0];
  {
    int shape[2];
    shape[0] = lda;
    shape[1] = n;
    rb_a_out__ = na_make_object(NA_SCOMPLEX, 2, shape, cNArray);
  }
  a_out__ = NA_PTR_TYPE(rb_a_out__, complex*);
  MEMCPY(a_out__, a, complex, NA_TOTAL(rb_a));
  rb_a = rb_a_out__;
  a = a_out__;

  clauu2_(&uplo, &n, a, &lda, &info);

  rb_info = INT2NUM(info);
  return rb_ary_new3(2, rb_info, rb_a);
}

void
init_lapack_clauu2(VALUE mLapack){
  rb_define_module_function(mLapack, "clauu2", rb_clauu2, -1);
}

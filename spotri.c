#include "rb_lapack.h"

static VALUE
rb_spotri(int argc, VALUE *argv, VALUE self){
  VALUE rb_uplo;
  char uplo; 
  VALUE rb_a;
  real *a; 
  VALUE rb_info;
  integer info; 
  VALUE rb_a_out__;
  real *a_out__;

  integer lda;
  integer n;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  info, a = NumRu::Lapack.spotri( uplo, a)\n    or\n  NumRu::Lapack.spotri  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE SPOTRI( UPLO, N, A, LDA, INFO )\n\n*  Purpose\n*  =======\n*\n*  SPOTRI computes the inverse of a real symmetric positive definite\n*  matrix A using the Cholesky factorization A = U**T*U or A = L*L**T\n*  computed by SPOTRF.\n*\n\n*  Arguments\n*  =========\n*\n*  UPLO    (input) CHARACTER*1\n*          = 'U':  Upper triangle of A is stored;\n*          = 'L':  Lower triangle of A is stored.\n*\n*  N       (input) INTEGER\n*          The order of the matrix A.  N >= 0.\n*\n*  A       (input/output) REAL array, dimension (LDA,N)\n*          On entry, the triangular factor U or L from the Cholesky\n*          factorization A = U**T*U or A = L*L**T, as computed by\n*          SPOTRF.\n*          On exit, the upper or lower triangle of the (symmetric)\n*          inverse of A, overwriting the input factor U or L.\n*\n*  LDA     (input) INTEGER\n*          The leading dimension of the array A.  LDA >= max(1,N).\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit\n*          < 0:  if INFO = -i, the i-th argument had an illegal value\n*          > 0:  if INFO = i, the (i,i) element of the factor U or L is\n*                zero, and the inverse could not be computed.\n*\n\n*  =====================================================================\n*\n*     .. External Functions ..\n      LOGICAL            LSAME\n      EXTERNAL           LSAME\n*     ..\n*     .. External Subroutines ..\n      EXTERNAL           SLAUUM, STRTRI, XERBLA\n*     ..\n*     .. Intrinsic Functions ..\n      INTRINSIC          MAX\n*     ..\n\n");
    return Qnil;
  }
  if (argc != 2)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 2)", argc);
  rb_uplo = argv[0];
  rb_a = argv[1];

  uplo = StringValueCStr(rb_uplo)[0];
  if (!NA_IsNArray(rb_a))
    rb_raise(rb_eArgError, "a (2th argument) must be NArray");
  if (NA_RANK(rb_a) != 2)
    rb_raise(rb_eArgError, "rank of a (2th argument) must be %d", 2);
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

  spotri_(&uplo, &n, a, &lda, &info);

  rb_info = INT2NUM(info);
  return rb_ary_new3(2, rb_info, rb_a);
}

void
init_lapack_spotri(VALUE mLapack){
  rb_define_module_function(mLapack, "spotri", rb_spotri, -1);
}

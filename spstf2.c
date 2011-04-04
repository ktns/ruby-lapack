#include "rb_lapack.h"

extern VOID spstf2_(char *uplo, integer *n, real *a, integer *lda, integer *piv, integer *rank, real *tol, real *work, integer *info);

static VALUE
rb_spstf2(int argc, VALUE *argv, VALUE self){
  VALUE rb_uplo;
  char uplo; 
  VALUE rb_a;
  real *a; 
  VALUE rb_tol;
  real tol; 
  VALUE rb_piv;
  integer *piv; 
  VALUE rb_rank;
  integer rank; 
  VALUE rb_info;
  integer info; 
  VALUE rb_a_out__;
  real *a_out__;
  real *work;

  integer lda;
  integer n;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  piv, rank, info, a = NumRu::Lapack.spstf2( uplo, a, tol)\n    or\n  NumRu::Lapack.spstf2  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE SPSTF2( UPLO, N, A, LDA, PIV, RANK, TOL, WORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  SPSTF2 computes the Cholesky factorization with complete\n*  pivoting of a real symmetric positive semidefinite matrix A.\n*\n*  The factorization has the form\n*     P' * A * P = U' * U ,  if UPLO = 'U',\n*     P' * A * P = L  * L',  if UPLO = 'L',\n*  where U is an upper triangular matrix and L is lower triangular, and\n*  P is stored as vector PIV.\n*\n*  This algorithm does not attempt to check that A is positive\n*  semidefinite. This version of the algorithm calls level 2 BLAS.\n*\n\n*  Arguments\n*  =========\n*\n*  UPLO    (input) CHARACTER*1\n*          Specifies whether the upper or lower triangular part of the\n*          symmetric matrix A is stored.\n*          = 'U':  Upper triangular\n*          = 'L':  Lower triangular\n*\n*  N       (input) INTEGER\n*          The order of the matrix A.  N >= 0.\n*\n*  A       (input/output) REAL array, dimension (LDA,N)\n*          On entry, the symmetric matrix A.  If UPLO = 'U', the leading\n*          n by n upper triangular part of A contains the upper\n*          triangular part of the matrix A, and the strictly lower\n*          triangular part of A is not referenced.  If UPLO = 'L', the\n*          leading n by n lower triangular part of A contains the lower\n*          triangular part of the matrix A, and the strictly upper\n*          triangular part of A is not referenced.\n*\n*          On exit, if INFO = 0, the factor U or L from the Cholesky\n*          factorization as above.\n*\n*  PIV     (output) INTEGER array, dimension (N)\n*          PIV is such that the nonzero entries are P( PIV(K), K ) = 1.\n*\n*  RANK    (output) INTEGER\n*          The rank of A given by the number of steps the algorithm\n*          completed.\n*\n*  TOL     (input) REAL\n*          User defined tolerance. If TOL < 0, then N*U*MAX( A( K,K ) )\n*          will be used. The algorithm terminates at the (K-1)st step\n*          if the pivot <= TOL.\n*\n*  LDA     (input) INTEGER\n*          The leading dimension of the array A.  LDA >= max(1,N).\n*\n*  WORK    (workspace) REAL array, dimension (2*N)\n*          Work space.\n*\n*  INFO    (output) INTEGER\n*          < 0: If INFO = -K, the K-th argument had an illegal value,\n*          = 0: algorithm completed successfully, and\n*          > 0: the matrix A is either rank deficient with computed rank\n*               as returned in RANK, or is indefinite.  See Section 7 of\n*               LAPACK Working Note #161 for further information.\n*\n\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 3)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 3)", argc);
  rb_uplo = argv[0];
  rb_a = argv[1];
  rb_tol = argv[2];

  if (!NA_IsNArray(rb_a))
    rb_raise(rb_eArgError, "a (2th argument) must be NArray");
  if (NA_RANK(rb_a) != 2)
    rb_raise(rb_eArgError, "rank of a (2th argument) must be %d", 2);
  n = NA_SHAPE1(rb_a);
  lda = NA_SHAPE0(rb_a);
  if (NA_TYPE(rb_a) != NA_SFLOAT)
    rb_a = na_change_type(rb_a, NA_SFLOAT);
  a = NA_PTR_TYPE(rb_a, real*);
  tol = (real)NUM2DBL(rb_tol);
  uplo = StringValueCStr(rb_uplo)[0];
  {
    int shape[1];
    shape[0] = n;
    rb_piv = na_make_object(NA_LINT, 1, shape, cNArray);
  }
  piv = NA_PTR_TYPE(rb_piv, integer*);
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
  work = ALLOC_N(real, (2*n));

  spstf2_(&uplo, &n, a, &lda, piv, &rank, &tol, work, &info);

  free(work);
  rb_rank = INT2NUM(rank);
  rb_info = INT2NUM(info);
  return rb_ary_new3(4, rb_piv, rb_rank, rb_info, rb_a);
}

void
init_lapack_spstf2(VALUE mLapack){
  rb_define_module_function(mLapack, "spstf2", rb_spstf2, -1);
}

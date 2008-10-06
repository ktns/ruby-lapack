#include "rb_lapack.h"

static VALUE
rb_csytf2(int argc, VALUE *argv, VALUE self){
  VALUE rb_uplo;
  char uplo; 
  VALUE rb_a;
  complex *a; 
  VALUE rb_ipiv;
  integer *ipiv; 
  VALUE rb_info;
  integer info; 
  VALUE rb_a_out__;
  complex *a_out__;

  integer lda;
  integer n;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  ipiv, info, a = NumRu::Lapack.csytf2( uplo, a)\n    or\n  NumRu::Lapack.csytf2  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE CSYTF2( UPLO, N, A, LDA, IPIV, INFO )\n\n*  Purpose\n*  =======\n*\n*  CSYTF2 computes the factorization of a complex symmetric matrix A\n*  using the Bunch-Kaufman diagonal pivoting method:\n*\n*     A = U*D*U'  or  A = L*D*L'\n*\n*  where U (or L) is a product of permutation and unit upper (lower)\n*  triangular matrices, U' is the transpose of U, and D is symmetric and\n*  block diagonal with 1-by-1 and 2-by-2 diagonal blocks.\n*\n*  This is the unblocked version of the algorithm, calling Level 2 BLAS.\n*\n\n*  Arguments\n*  =========\n*\n*  UPLO    (input) CHARACTER*1\n*          Specifies whether the upper or lower triangular part of the\n*          symmetric matrix A is stored:\n*          = 'U':  Upper triangular\n*          = 'L':  Lower triangular\n*\n*  N       (input) INTEGER\n*          The order of the matrix A.  N >= 0.\n*\n*  A       (input/output) COMPLEX array, dimension (LDA,N)\n*          On entry, the symmetric matrix A.  If UPLO = 'U', the leading\n*          n-by-n upper triangular part of A contains the upper\n*          triangular part of the matrix A, and the strictly lower\n*          triangular part of A is not referenced.  If UPLO = 'L', the\n*          leading n-by-n lower triangular part of A contains the lower\n*          triangular part of the matrix A, and the strictly upper\n*          triangular part of A is not referenced.\n*\n*          On exit, the block diagonal matrix D and the multipliers used\n*          to obtain the factor U or L (see below for further details).\n*\n*  LDA     (input) INTEGER\n*          The leading dimension of the array A.  LDA >= max(1,N).\n*\n*  IPIV    (output) INTEGER array, dimension (N)\n*          Details of the interchanges and the block structure of D.\n*          If IPIV(k) > 0, then rows and columns k and IPIV(k) were\n*          interchanged and D(k,k) is a 1-by-1 diagonal block.\n*          If UPLO = 'U' and IPIV(k) = IPIV(k-1) < 0, then rows and\n*          columns k-1 and -IPIV(k) were interchanged and D(k-1:k,k-1:k)\n*          is a 2-by-2 diagonal block.  If UPLO = 'L' and IPIV(k) =\n*          IPIV(k+1) < 0, then rows and columns k+1 and -IPIV(k) were\n*          interchanged and D(k:k+1,k:k+1) is a 2-by-2 diagonal block.\n*\n*  INFO    (output) INTEGER\n*          = 0: successful exit\n*          < 0: if INFO = -k, the k-th argument had an illegal value\n*          > 0: if INFO = k, D(k,k) is exactly zero.  The factorization\n*               has been completed, but the block diagonal matrix D is\n*               exactly singular, and division by zero will occur if it\n*               is used to solve a system of equations.\n*\n\n*  Further Details\n*  ===============\n*\n*  09-29-06 - patch from\n*    Bobby Cheng, MathWorks\n*\n*    Replace l.209 and l.377\n*         IF( MAX( ABSAKK, COLMAX ).EQ.ZERO ) THEN\n*    by\n*         IF( (MAX( ABSAKK, COLMAX ).EQ.ZERO) .OR. SISNAN(ABSAKK) ) THEN\n*\n*  1-96 - Based on modifications by J. Lewis, Boeing Computer Services\n*         Company\n*\n*  If UPLO = 'U', then A = U*D*U', where\n*     U = P(n)*U(n)* ... *P(k)U(k)* ...,\n*  i.e., U is a product of terms P(k)*U(k), where k decreases from n to\n*  1 in steps of 1 or 2, and D is a block diagonal matrix with 1-by-1\n*  and 2-by-2 diagonal blocks D(k).  P(k) is a permutation matrix as\n*  defined by IPIV(k), and U(k) is a unit upper triangular matrix, such\n*  that if the diagonal block D(k) is of order s (s = 1 or 2), then\n*\n*             (   I    v    0   )   k-s\n*     U(k) =  (   0    I    0   )   s\n*             (   0    0    I   )   n-k\n*                k-s   s   n-k\n*\n*  If s = 1, D(k) overwrites A(k,k), and v overwrites A(1:k-1,k).\n*  If s = 2, the upper triangle of D(k) overwrites A(k-1,k-1), A(k-1,k),\n*  and A(k,k), and v overwrites A(1:k-2,k-1:k).\n*\n*  If UPLO = 'L', then A = L*D*L', where\n*     L = P(1)*L(1)* ... *P(k)*L(k)* ...,\n*  i.e., L is a product of terms P(k)*L(k), where k increases from 1 to\n*  n in steps of 1 or 2, and D is a block diagonal matrix with 1-by-1\n*  and 2-by-2 diagonal blocks D(k).  P(k) is a permutation matrix as\n*  defined by IPIV(k), and L(k) is a unit lower triangular matrix, such\n*  that if the diagonal block D(k) is of order s (s = 1 or 2), then\n*\n*             (   I    0     0   )  k-1\n*     L(k) =  (   0    I     0   )  s\n*             (   0    v     I   )  n-k-s+1\n*                k-1   s  n-k-s+1\n*\n*  If s = 1, D(k) overwrites A(k,k), and v overwrites A(k+1:n,k).\n*  If s = 2, the lower triangle of D(k) overwrites A(k,k), A(k+1,k),\n*  and A(k+1,k+1), and v overwrites A(k+2:n,k:k+1).\n*\n*  =====================================================================\n*\n\n");
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
  if (NA_TYPE(rb_a) != NA_SCOMPLEX)
    rb_a = na_change_type(rb_a, NA_SCOMPLEX);
  a = NA_PTR_TYPE(rb_a, complex*);
  {
    int shape[1];
    shape[0] = n;
    rb_ipiv = na_make_object(NA_LINT, 1, shape, cNArray);
  }
  ipiv = NA_PTR_TYPE(rb_ipiv, integer*);
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

  csytf2_(&uplo, &n, a, &lda, ipiv, &info);

  rb_info = INT2NUM(info);
  return rb_ary_new3(3, rb_ipiv, rb_info, rb_a);
}

void
init_lapack_csytf2(VALUE mLapack){
  rb_define_module_function(mLapack, "csytf2", rb_csytf2, -1);
}

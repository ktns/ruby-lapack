#include "rb_lapack.h"

static VALUE
rb_csysv(int argc, VALUE *argv, VALUE self){
  VALUE rb_uplo;
  char uplo; 
  VALUE rb_a;
  complex *a; 
  VALUE rb_b;
  complex *b; 
  VALUE rb_lwork;
  integer lwork; 
  VALUE rb_ipiv;
  integer *ipiv; 
  VALUE rb_work;
  complex *work; 
  VALUE rb_info;
  integer info; 
  VALUE rb_a_out__;
  complex *a_out__;
  VALUE rb_b_out__;
  complex *b_out__;

  integer lda;
  integer n;
  integer ldb;
  integer nrhs;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  ipiv, work, info, a, b = NumRu::Lapack.csysv( uplo, a, b, lwork)\n    or\n  NumRu::Lapack.csysv  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE CSYSV( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, WORK, LWORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  CSYSV computes the solution to a complex system of linear equations\n*     A * X = B,\n*  where A is an N-by-N symmetric matrix and X and B are N-by-NRHS\n*  matrices.\n*\n*  The diagonal pivoting method is used to factor A as\n*     A = U * D * U**T,  if UPLO = 'U', or\n*     A = L * D * L**T,  if UPLO = 'L',\n*  where U (or L) is a product of permutation and unit upper (lower)\n*  triangular matrices, and D is symmetric and block diagonal with \n*  1-by-1 and 2-by-2 diagonal blocks.  The factored form of A is then\n*  used to solve the system of equations A * X = B.\n*\n\n*  Arguments\n*  =========\n*\n*  UPLO    (input) CHARACTER*1\n*          = 'U':  Upper triangle of A is stored;\n*          = 'L':  Lower triangle of A is stored.\n*\n*  N       (input) INTEGER\n*          The number of linear equations, i.e., the order of the\n*          matrix A.  N >= 0.\n*\n*  NRHS    (input) INTEGER\n*          The number of right hand sides, i.e., the number of columns\n*          of the matrix B.  NRHS >= 0.\n*\n*  A       (input/output) COMPLEX array, dimension (LDA,N)\n*          On entry, the symmetric matrix A.  If UPLO = 'U', the leading\n*          N-by-N upper triangular part of A contains the upper\n*          triangular part of the matrix A, and the strictly lower\n*          triangular part of A is not referenced.  If UPLO = 'L', the\n*          leading N-by-N lower triangular part of A contains the lower\n*          triangular part of the matrix A, and the strictly upper\n*          triangular part of A is not referenced.\n*\n*          On exit, if INFO = 0, the block diagonal matrix D and the\n*          multipliers used to obtain the factor U or L from the\n*          factorization A = U*D*U**T or A = L*D*L**T as computed by\n*          CSYTRF.\n*\n*  LDA     (input) INTEGER\n*          The leading dimension of the array A.  LDA >= max(1,N).\n*\n*  IPIV    (output) INTEGER array, dimension (N)\n*          Details of the interchanges and the block structure of D, as\n*          determined by CSYTRF.  If IPIV(k) > 0, then rows and columns\n*          k and IPIV(k) were interchanged, and D(k,k) is a 1-by-1\n*          diagonal block.  If UPLO = 'U' and IPIV(k) = IPIV(k-1) < 0,\n*          then rows and columns k-1 and -IPIV(k) were interchanged and\n*          D(k-1:k,k-1:k) is a 2-by-2 diagonal block.  If UPLO = 'L' and\n*          IPIV(k) = IPIV(k+1) < 0, then rows and columns k+1 and\n*          -IPIV(k) were interchanged and D(k:k+1,k:k+1) is a 2-by-2\n*          diagonal block.\n*\n*  B       (input/output) COMPLEX array, dimension (LDB,NRHS)\n*          On entry, the N-by-NRHS right hand side matrix B.\n*          On exit, if INFO = 0, the N-by-NRHS solution matrix X.\n*\n*  LDB     (input) INTEGER\n*          The leading dimension of the array B.  LDB >= max(1,N).\n*\n*  WORK    (workspace/output) COMPLEX array, dimension (MAX(1,LWORK))\n*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.\n*\n*  LWORK   (input) INTEGER\n*          The length of WORK.  LWORK >= 1, and for best performance\n*          LWORK >= max(1,N*NB), where NB is the optimal blocksize for\n*          CSYTRF.\n*\n*          If LWORK = -1, then a workspace query is assumed; the routine\n*          only calculates the optimal size of the WORK array, returns\n*          this value as the first entry of the WORK array, and no error\n*          message related to LWORK is issued by XERBLA.\n*\n*  INFO    (output) INTEGER\n*          = 0: successful exit\n*          < 0: if INFO = -i, the i-th argument had an illegal value\n*          > 0: if INFO = i, D(i,i) is exactly zero.  The factorization\n*               has been completed, but the block diagonal matrix D is\n*               exactly singular, so the solution could not be computed.\n*\n\n*  =====================================================================\n*\n*     .. Local Scalars ..\n      LOGICAL            LQUERY\n      INTEGER            LWKOPT, NB\n*     ..\n*     .. External Functions ..\n      LOGICAL            LSAME\n      INTEGER            ILAENV\n      EXTERNAL           ILAENV, LSAME\n*     ..\n*     .. External Subroutines ..\n      EXTERNAL           CSYTRF, CSYTRS, XERBLA\n*     ..\n*     .. Intrinsic Functions ..\n      INTRINSIC          MAX\n*     ..\n\n");
    return Qnil;
  }
  if (argc != 4)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 4)", argc);
  rb_uplo = argv[0];
  rb_a = argv[1];
  rb_b = argv[2];
  rb_lwork = argv[3];

  uplo = StringValueCStr(rb_uplo)[0];
  lwork = NUM2INT(rb_lwork);
  if (!NA_IsNArray(rb_a))
    rb_raise(rb_eArgError, "a (2th argument) must be NArray");
  if (NA_RANK(rb_a) != 2)
    rb_raise(rb_eArgError, "rank of a (2th argument) must be %d", 2);
  lda = NA_SHAPE0(rb_a);
  n = NA_SHAPE1(rb_a);
  if (NA_TYPE(rb_a) != NA_SCOMPLEX)
    rb_a = na_change_type(rb_a, NA_SCOMPLEX);
  a = NA_PTR_TYPE(rb_a, complex*);
  if (!NA_IsNArray(rb_b))
    rb_raise(rb_eArgError, "b (3th argument) must be NArray");
  if (NA_RANK(rb_b) != 2)
    rb_raise(rb_eArgError, "rank of b (3th argument) must be %d", 2);
  ldb = NA_SHAPE0(rb_b);
  nrhs = NA_SHAPE1(rb_b);
  if (NA_TYPE(rb_b) != NA_SCOMPLEX)
    rb_b = na_change_type(rb_b, NA_SCOMPLEX);
  b = NA_PTR_TYPE(rb_b, complex*);
  {
    int shape[1];
    shape[0] = n;
    rb_ipiv = na_make_object(NA_LINT, 1, shape, cNArray);
  }
  ipiv = NA_PTR_TYPE(rb_ipiv, integer*);
  {
    int shape[1];
    shape[0] = MAX(1,lwork);
    rb_work = na_make_object(NA_SCOMPLEX, 1, shape, cNArray);
  }
  work = NA_PTR_TYPE(rb_work, complex*);
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
  {
    int shape[2];
    shape[0] = ldb;
    shape[1] = nrhs;
    rb_b_out__ = na_make_object(NA_SCOMPLEX, 2, shape, cNArray);
  }
  b_out__ = NA_PTR_TYPE(rb_b_out__, complex*);
  MEMCPY(b_out__, b, complex, NA_TOTAL(rb_b));
  rb_b = rb_b_out__;
  b = b_out__;

  csysv_(&uplo, &n, &nrhs, a, &lda, ipiv, b, &ldb, work, &lwork, &info);

  rb_info = INT2NUM(info);
  return rb_ary_new3(5, rb_ipiv, rb_work, rb_info, rb_a, rb_b);
}

void
init_lapack_csysv(VALUE mLapack){
  rb_define_module_function(mLapack, "csysv", rb_csysv, -1);
}

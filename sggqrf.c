#include "rb_lapack.h"

extern VOID sggqrf_(integer *n, integer *m, integer *p, real *a, integer *lda, real *taua, real *b, integer *ldb, real *taub, real *work, integer *lwork, integer *info);

static VALUE sHelp, sUsage;

static VALUE
rblapack_sggqrf(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_n;
  integer n; 
  VALUE rblapack_a;
  real *a; 
  VALUE rblapack_b;
  real *b; 
  VALUE rblapack_lwork;
  integer lwork; 
  VALUE rblapack_taua;
  real *taua; 
  VALUE rblapack_taub;
  real *taub; 
  VALUE rblapack_work;
  real *work; 
  VALUE rblapack_info;
  integer info; 
  VALUE rblapack_a_out__;
  real *a_out__;
  VALUE rblapack_b_out__;
  real *b_out__;

  integer lda;
  integer m;
  integer ldb;
  integer p;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  taua, taub, work, info, a, b = NumRu::Lapack.sggqrf( n, a, b, lwork, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE SGGQRF( N, M, P, A, LDA, TAUA, B, LDB, TAUB, WORK, LWORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  SGGQRF computes a generalized QR factorization of an N-by-M matrix A\n*  and an N-by-P matrix B:\n*\n*              A = Q*R,        B = Q*T*Z,\n*\n*  where Q is an N-by-N orthogonal matrix, Z is a P-by-P orthogonal\n*  matrix, and R and T assume one of the forms:\n*\n*  if N >= M,  R = ( R11 ) M  ,   or if N < M,  R = ( R11  R12 ) N,\n*                  (  0  ) N-M                         N   M-N\n*                     M\n*\n*  where R11 is upper triangular, and\n*\n*  if N <= P,  T = ( 0  T12 ) N,   or if N > P,  T = ( T11 ) N-P,\n*                   P-N  N                           ( T21 ) P\n*                                                       P\n*\n*  where T12 or T21 is upper triangular.\n*\n*  In particular, if B is square and nonsingular, the GQR factorization\n*  of A and B implicitly gives the QR factorization of inv(B)*A:\n*\n*               inv(B)*A = Z'*(inv(T)*R)\n*\n*  where inv(B) denotes the inverse of the matrix B, and Z' denotes the\n*  transpose of the matrix Z.\n*\n\n*  Arguments\n*  =========\n*\n*  N       (input) INTEGER\n*          The number of rows of the matrices A and B. N >= 0.\n*\n*  M       (input) INTEGER\n*          The number of columns of the matrix A.  M >= 0.\n*\n*  P       (input) INTEGER\n*          The number of columns of the matrix B.  P >= 0.\n*\n*  A       (input/output) REAL array, dimension (LDA,M)\n*          On entry, the N-by-M matrix A.\n*          On exit, the elements on and above the diagonal of the array\n*          contain the min(N,M)-by-M upper trapezoidal matrix R (R is\n*          upper triangular if N >= M); the elements below the diagonal,\n*          with the array TAUA, represent the orthogonal matrix Q as a\n*          product of min(N,M) elementary reflectors (see Further\n*          Details).\n*\n*  LDA     (input) INTEGER\n*          The leading dimension of the array A. LDA >= max(1,N).\n*\n*  TAUA    (output) REAL array, dimension (min(N,M))\n*          The scalar factors of the elementary reflectors which\n*          represent the orthogonal matrix Q (see Further Details).\n*\n*  B       (input/output) REAL array, dimension (LDB,P)\n*          On entry, the N-by-P matrix B.\n*          On exit, if N <= P, the upper triangle of the subarray\n*          B(1:N,P-N+1:P) contains the N-by-N upper triangular matrix T;\n*          if N > P, the elements on and above the (N-P)-th subdiagonal\n*          contain the N-by-P upper trapezoidal matrix T; the remaining\n*          elements, with the array TAUB, represent the orthogonal\n*          matrix Z as a product of elementary reflectors (see Further\n*          Details).\n*\n*  LDB     (input) INTEGER\n*          The leading dimension of the array B. LDB >= max(1,N).\n*\n*  TAUB    (output) REAL array, dimension (min(N,P))\n*          The scalar factors of the elementary reflectors which\n*          represent the orthogonal matrix Z (see Further Details).\n*\n*  WORK    (workspace/output) REAL array, dimension (MAX(1,LWORK))\n*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.\n*\n*  LWORK   (input) INTEGER\n*          The dimension of the array WORK. LWORK >= max(1,N,M,P).\n*          For optimum performance LWORK >= max(N,M,P)*max(NB1,NB2,NB3),\n*          where NB1 is the optimal blocksize for the QR factorization\n*          of an N-by-M matrix, NB2 is the optimal blocksize for the\n*          RQ factorization of an N-by-P matrix, and NB3 is the optimal\n*          blocksize for a call of SORMQR.\n*\n*          If LWORK = -1, then a workspace query is assumed; the routine\n*          only calculates the optimal size of the WORK array, returns\n*          this value as the first entry of the WORK array, and no error\n*          message related to LWORK is issued by XERBLA.\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit\n*          < 0:  if INFO = -i, the i-th argument had an illegal value.\n*\n\n*  Further Details\n*  ===============\n*\n*  The matrix Q is represented as a product of elementary reflectors\n*\n*     Q = H(1) H(2) . . . H(k), where k = min(n,m).\n*\n*  Each H(i) has the form\n*\n*     H(i) = I - taua * v * v'\n*\n*  where taua is a real scalar, and v is a real vector with\n*  v(1:i-1) = 0 and v(i) = 1; v(i+1:n) is stored on exit in A(i+1:n,i),\n*  and taua in TAUA(i).\n*  To form Q explicitly, use LAPACK subroutine SORGQR.\n*  To use Q to update another matrix, use LAPACK subroutine SORMQR.\n*\n*  The matrix Z is represented as a product of elementary reflectors\n*\n*     Z = H(1) H(2) . . . H(k), where k = min(n,p).\n*\n*  Each H(i) has the form\n*\n*     H(i) = I - taub * v * v'\n*\n*  where taub is a real scalar, and v is a real vector with\n*  v(p-k+i+1:p) = 0 and v(p-k+i) = 1; v(1:p-k+i-1) is stored on exit in\n*  B(n-k+i,1:p-k+i-1), and taub in TAUB(i).\n*  To form Z explicitly, use LAPACK subroutine SORGRQ.\n*  To use Z to update another matrix, use LAPACK subroutine SORMRQ.\n*\n*  =====================================================================\n*\n*     .. Local Scalars ..\n      LOGICAL            LQUERY\n      INTEGER            LOPT, LWKOPT, NB, NB1, NB2, NB3\n*     ..\n*     .. External Subroutines ..\n      EXTERNAL           SGEQRF, SGERQF, SORMQR, XERBLA\n*     ..\n*     .. External Functions ..\n      INTEGER            ILAENV\n      EXTERNAL           ILAENV \n*     ..\n*     .. Intrinsic Functions ..\n      INTRINSIC          INT, MAX, MIN\n*     ..\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  taua, taub, work, info, a, b = NumRu::Lapack.sggqrf( n, a, b, lwork, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 4)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 4)", argc);
  rblapack_n = argv[0];
  rblapack_a = argv[1];
  rblapack_b = argv[2];
  rblapack_lwork = argv[3];
  if (rb_options != Qnil) {
  }

  if (!NA_IsNArray(rblapack_a))
    rb_raise(rb_eArgError, "a (2th argument) must be NArray");
  if (NA_RANK(rblapack_a) != 2)
    rb_raise(rb_eArgError, "rank of a (2th argument) must be %d", 2);
  m = NA_SHAPE1(rblapack_a);
  lda = NA_SHAPE0(rblapack_a);
  if (NA_TYPE(rblapack_a) != NA_SFLOAT)
    rblapack_a = na_change_type(rblapack_a, NA_SFLOAT);
  a = NA_PTR_TYPE(rblapack_a, real*);
  if (!NA_IsNArray(rblapack_b))
    rb_raise(rb_eArgError, "b (3th argument) must be NArray");
  if (NA_RANK(rblapack_b) != 2)
    rb_raise(rb_eArgError, "rank of b (3th argument) must be %d", 2);
  p = NA_SHAPE1(rblapack_b);
  ldb = NA_SHAPE0(rblapack_b);
  if (NA_TYPE(rblapack_b) != NA_SFLOAT)
    rblapack_b = na_change_type(rblapack_b, NA_SFLOAT);
  b = NA_PTR_TYPE(rblapack_b, real*);
  n = NUM2INT(rblapack_n);
  lwork = NUM2INT(rblapack_lwork);
  {
    int shape[1];
    shape[0] = MIN(n,m);
    rblapack_taua = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  taua = NA_PTR_TYPE(rblapack_taua, real*);
  {
    int shape[1];
    shape[0] = MIN(n,p);
    rblapack_taub = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  taub = NA_PTR_TYPE(rblapack_taub, real*);
  {
    int shape[1];
    shape[0] = MAX(1,lwork);
    rblapack_work = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  work = NA_PTR_TYPE(rblapack_work, real*);
  {
    int shape[2];
    shape[0] = lda;
    shape[1] = m;
    rblapack_a_out__ = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  a_out__ = NA_PTR_TYPE(rblapack_a_out__, real*);
  MEMCPY(a_out__, a, real, NA_TOTAL(rblapack_a));
  rblapack_a = rblapack_a_out__;
  a = a_out__;
  {
    int shape[2];
    shape[0] = ldb;
    shape[1] = p;
    rblapack_b_out__ = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  b_out__ = NA_PTR_TYPE(rblapack_b_out__, real*);
  MEMCPY(b_out__, b, real, NA_TOTAL(rblapack_b));
  rblapack_b = rblapack_b_out__;
  b = b_out__;

  sggqrf_(&n, &m, &p, a, &lda, taua, b, &ldb, taub, work, &lwork, &info);

  rblapack_info = INT2NUM(info);
  return rb_ary_new3(6, rblapack_taua, rblapack_taub, rblapack_work, rblapack_info, rblapack_a, rblapack_b);
}

void
init_lapack_sggqrf(VALUE mLapack){
  rb_define_module_function(mLapack, "sggqrf", rblapack_sggqrf, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}

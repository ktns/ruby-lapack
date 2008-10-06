#include "rb_lapack.h"

static VALUE
rb_sggrqf(int argc, VALUE *argv, VALUE self){
  VALUE rb_m;
  integer m; 
  VALUE rb_p;
  integer p; 
  VALUE rb_a;
  real *a; 
  VALUE rb_b;
  real *b; 
  VALUE rb_lwork;
  integer lwork; 
  VALUE rb_taua;
  real *taua; 
  VALUE rb_taub;
  real *taub; 
  VALUE rb_work;
  real *work; 
  VALUE rb_info;
  integer info; 
  VALUE rb_a_out__;
  real *a_out__;
  VALUE rb_b_out__;
  real *b_out__;

  integer lda;
  integer n;
  integer ldb;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  taua, taub, work, info, a, b = NumRu::Lapack.sggrqf( m, p, a, b, lwork)\n    or\n  NumRu::Lapack.sggrqf  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE SGGRQF( M, P, N, A, LDA, TAUA, B, LDB, TAUB, WORK, LWORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  SGGRQF computes a generalized RQ factorization of an M-by-N matrix A\n*  and a P-by-N matrix B:\n*\n*              A = R*Q,        B = Z*T*Q,\n*\n*  where Q is an N-by-N orthogonal matrix, Z is a P-by-P orthogonal\n*  matrix, and R and T assume one of the forms:\n*\n*  if M <= N,  R = ( 0  R12 ) M,   or if M > N,  R = ( R11 ) M-N,\n*                   N-M  M                           ( R21 ) N\n*                                                       N\n*\n*  where R12 or R21 is upper triangular, and\n*\n*  if P >= N,  T = ( T11 ) N  ,   or if P < N,  T = ( T11  T12 ) P,\n*                  (  0  ) P-N                         P   N-P\n*                     N\n*\n*  where T11 is upper triangular.\n*\n*  In particular, if B is square and nonsingular, the GRQ factorization\n*  of A and B implicitly gives the RQ factorization of A*inv(B):\n*\n*               A*inv(B) = (R*inv(T))*Z'\n*\n*  where inv(B) denotes the inverse of the matrix B, and Z' denotes the\n*  transpose of the matrix Z.\n*\n\n*  Arguments\n*  =========\n*\n*  M       (input) INTEGER\n*          The number of rows of the matrix A.  M >= 0.\n*\n*  P       (input) INTEGER\n*          The number of rows of the matrix B.  P >= 0.\n*\n*  N       (input) INTEGER\n*          The number of columns of the matrices A and B. N >= 0.\n*\n*  A       (input/output) REAL array, dimension (LDA,N)\n*          On entry, the M-by-N matrix A.\n*          On exit, if M <= N, the upper triangle of the subarray\n*          A(1:M,N-M+1:N) contains the M-by-M upper triangular matrix R;\n*          if M > N, the elements on and above the (M-N)-th subdiagonal\n*          contain the M-by-N upper trapezoidal matrix R; the remaining\n*          elements, with the array TAUA, represent the orthogonal\n*          matrix Q as a product of elementary reflectors (see Further\n*          Details).\n*\n*  LDA     (input) INTEGER\n*          The leading dimension of the array A. LDA >= max(1,M).\n*\n*  TAUA    (output) REAL array, dimension (min(M,N))\n*          The scalar factors of the elementary reflectors which\n*          represent the orthogonal matrix Q (see Further Details).\n*\n*  B       (input/output) REAL array, dimension (LDB,N)\n*          On entry, the P-by-N matrix B.\n*          On exit, the elements on and above the diagonal of the array\n*          contain the min(P,N)-by-N upper trapezoidal matrix T (T is\n*          upper triangular if P >= N); the elements below the diagonal,\n*          with the array TAUB, represent the orthogonal matrix Z as a\n*          product of elementary reflectors (see Further Details).\n*\n*  LDB     (input) INTEGER\n*          The leading dimension of the array B. LDB >= max(1,P).\n*\n*  TAUB    (output) REAL array, dimension (min(P,N))\n*          The scalar factors of the elementary reflectors which\n*          represent the orthogonal matrix Z (see Further Details).\n*\n*  WORK    (workspace/output) REAL array, dimension (MAX(1,LWORK))\n*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.\n*\n*  LWORK   (input) INTEGER\n*          The dimension of the array WORK. LWORK >= max(1,N,M,P).\n*          For optimum performance LWORK >= max(N,M,P)*max(NB1,NB2,NB3),\n*          where NB1 is the optimal blocksize for the RQ factorization\n*          of an M-by-N matrix, NB2 is the optimal blocksize for the\n*          QR factorization of a P-by-N matrix, and NB3 is the optimal\n*          blocksize for a call of SORMRQ.\n*\n*          If LWORK = -1, then a workspace query is assumed; the routine\n*          only calculates the optimal size of the WORK array, returns\n*          this value as the first entry of the WORK array, and no error\n*          message related to LWORK is issued by XERBLA.\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit\n*          < 0:  if INF0= -i, the i-th argument had an illegal value.\n*\n\n*  Further Details\n*  ===============\n*\n*  The matrix Q is represented as a product of elementary reflectors\n*\n*     Q = H(1) H(2) . . . H(k), where k = min(m,n).\n*\n*  Each H(i) has the form\n*\n*     H(i) = I - taua * v * v'\n*\n*  where taua is a real scalar, and v is a real vector with\n*  v(n-k+i+1:n) = 0 and v(n-k+i) = 1; v(1:n-k+i-1) is stored on exit in\n*  A(m-k+i,1:n-k+i-1), and taua in TAUA(i).\n*  To form Q explicitly, use LAPACK subroutine SORGRQ.\n*  To use Q to update another matrix, use LAPACK subroutine SORMRQ.\n*\n*  The matrix Z is represented as a product of elementary reflectors\n*\n*     Z = H(1) H(2) . . . H(k), where k = min(p,n).\n*\n*  Each H(i) has the form\n*\n*     H(i) = I - taub * v * v'\n*\n*  where taub is a real scalar, and v is a real vector with\n*  v(1:i-1) = 0 and v(i) = 1; v(i+1:p) is stored on exit in B(i+1:p,i),\n*  and taub in TAUB(i).\n*  To form Z explicitly, use LAPACK subroutine SORGQR.\n*  To use Z to update another matrix, use LAPACK subroutine SORMQR.\n*\n*  =====================================================================\n*\n*     .. Local Scalars ..\n      LOGICAL            LQUERY\n      INTEGER            LOPT, LWKOPT, NB, NB1, NB2, NB3\n*     ..\n*     .. External Subroutines ..\n      EXTERNAL           SGEQRF, SGERQF, SORMRQ, XERBLA\n*     ..\n*     .. External Functions ..\n      INTEGER            ILAENV \n      EXTERNAL           ILAENV \n*     ..\n*     .. Intrinsic Functions ..\n      INTRINSIC          INT, MAX, MIN\n*     ..\n\n");
    return Qnil;
  }
  if (argc != 5)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 5)", argc);
  rb_m = argv[0];
  rb_p = argv[1];
  rb_a = argv[2];
  rb_b = argv[3];
  rb_lwork = argv[4];

  m = NUM2INT(rb_m);
  p = NUM2INT(rb_p);
  lwork = NUM2INT(rb_lwork);
  if (!NA_IsNArray(rb_a))
    rb_raise(rb_eArgError, "a (3th argument) must be NArray");
  if (NA_RANK(rb_a) != 2)
    rb_raise(rb_eArgError, "rank of a (3th argument) must be %d", 2);
  lda = NA_SHAPE0(rb_a);
  n = NA_SHAPE1(rb_a);
  if (NA_TYPE(rb_a) != NA_SFLOAT)
    rb_a = na_change_type(rb_a, NA_SFLOAT);
  a = NA_PTR_TYPE(rb_a, real*);
  if (!NA_IsNArray(rb_b))
    rb_raise(rb_eArgError, "b (4th argument) must be NArray");
  if (NA_RANK(rb_b) != 2)
    rb_raise(rb_eArgError, "rank of b (4th argument) must be %d", 2);
  ldb = NA_SHAPE0(rb_b);
  if (NA_SHAPE1(rb_b) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of b must be the same as shape 1 of a");
  if (NA_TYPE(rb_b) != NA_SFLOAT)
    rb_b = na_change_type(rb_b, NA_SFLOAT);
  b = NA_PTR_TYPE(rb_b, real*);
  {
    int shape[1];
    shape[0] = MIN(m,n);
    rb_taua = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  taua = NA_PTR_TYPE(rb_taua, real*);
  {
    int shape[1];
    shape[0] = MIN(p,n);
    rb_taub = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  taub = NA_PTR_TYPE(rb_taub, real*);
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
  {
    int shape[2];
    shape[0] = ldb;
    shape[1] = n;
    rb_b_out__ = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  b_out__ = NA_PTR_TYPE(rb_b_out__, real*);
  MEMCPY(b_out__, b, real, NA_TOTAL(rb_b));
  rb_b = rb_b_out__;
  b = b_out__;

  sggrqf_(&m, &p, &n, a, &lda, taua, b, &ldb, taub, work, &lwork, &info);

  rb_info = INT2NUM(info);
  return rb_ary_new3(6, rb_taua, rb_taub, rb_work, rb_info, rb_a, rb_b);
}

void
init_lapack_sggrqf(VALUE mLapack){
  rb_define_module_function(mLapack, "sggrqf", rb_sggrqf, -1);
}

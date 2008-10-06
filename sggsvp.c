#include "rb_lapack.h"

static VALUE
rb_sggsvp(int argc, VALUE *argv, VALUE self){
  VALUE rb_jobu;
  char jobu; 
  VALUE rb_jobv;
  char jobv; 
  VALUE rb_jobq;
  char jobq; 
  VALUE rb_p;
  integer p; 
  VALUE rb_a;
  real *a; 
  VALUE rb_b;
  real *b; 
  VALUE rb_tola;
  real tola; 
  VALUE rb_tolb;
  real tolb; 
  VALUE rb_k;
  integer k; 
  VALUE rb_l;
  integer l; 
  VALUE rb_u;
  real *u; 
  VALUE rb_v;
  real *v; 
  VALUE rb_q;
  real *q; 
  VALUE rb_info;
  integer info; 
  VALUE rb_a_out__;
  real *a_out__;
  VALUE rb_b_out__;
  real *b_out__;
  integer *iwork;
  real *tau;
  real *work;

  integer lda;
  integer n;
  integer ldb;
  integer ldu;
  integer m;
  integer ldv;
  integer ldq;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  k, l, u, v, q, info, a, b = NumRu::Lapack.sggsvp( jobu, jobv, jobq, p, a, b, tola, tolb)\n    or\n  NumRu::Lapack.sggsvp  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE SGGSVP( JOBU, JOBV, JOBQ, M, P, N, A, LDA, B, LDB, TOLA, TOLB, K, L, U, LDU, V, LDV, Q, LDQ, IWORK, TAU, WORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  SGGSVP computes orthogonal matrices U, V and Q such that\n*\n*                   N-K-L  K    L\n*   U'*A*Q =     K ( 0    A12  A13 )  if M-K-L >= 0;\n*                L ( 0     0   A23 )\n*            M-K-L ( 0     0    0  )\n*\n*                   N-K-L  K    L\n*          =     K ( 0    A12  A13 )  if M-K-L < 0;\n*              M-K ( 0     0   A23 )\n*\n*                 N-K-L  K    L\n*   V'*B*Q =   L ( 0     0   B13 )\n*            P-L ( 0     0    0  )\n*\n*  where the K-by-K matrix A12 and L-by-L matrix B13 are nonsingular\n*  upper triangular; A23 is L-by-L upper triangular if M-K-L >= 0,\n*  otherwise A23 is (M-K)-by-L upper trapezoidal.  K+L = the effective\n*  numerical rank of the (M+P)-by-N matrix (A',B')'.  Z' denotes the\n*  transpose of Z.\n*\n*  This decomposition is the preprocessing step for computing the\n*  Generalized Singular Value Decomposition (GSVD), see subroutine\n*  SGGSVD.\n*\n\n*  Arguments\n*  =========\n*\n*  JOBU    (input) CHARACTER*1\n*          = 'U':  Orthogonal matrix U is computed;\n*          = 'N':  U is not computed.\n*\n*  JOBV    (input) CHARACTER*1\n*          = 'V':  Orthogonal matrix V is computed;\n*          = 'N':  V is not computed.\n*\n*  JOBQ    (input) CHARACTER*1\n*          = 'Q':  Orthogonal matrix Q is computed;\n*          = 'N':  Q is not computed.\n*\n*  M       (input) INTEGER\n*          The number of rows of the matrix A.  M >= 0.\n*\n*  P       (input) INTEGER\n*          The number of rows of the matrix B.  P >= 0.\n*\n*  N       (input) INTEGER\n*          The number of columns of the matrices A and B.  N >= 0.\n*\n*  A       (input/output) REAL array, dimension (LDA,N)\n*          On entry, the M-by-N matrix A.\n*          On exit, A contains the triangular (or trapezoidal) matrix\n*          described in the Purpose section.\n*\n*  LDA     (input) INTEGER\n*          The leading dimension of the array A. LDA >= max(1,M).\n*\n*  B       (input/output) REAL array, dimension (LDB,N)\n*          On entry, the P-by-N matrix B.\n*          On exit, B contains the triangular matrix described in\n*          the Purpose section.\n*\n*  LDB     (input) INTEGER\n*          The leading dimension of the array B. LDB >= max(1,P).\n*\n*  TOLA    (input) REAL\n*  TOLB    (input) REAL\n*          TOLA and TOLB are the thresholds to determine the effective\n*          numerical rank of matrix B and a subblock of A. Generally,\n*          they are set to\n*             TOLA = MAX(M,N)*norm(A)*MACHEPS,\n*             TOLB = MAX(P,N)*norm(B)*MACHEPS.\n*          The size of TOLA and TOLB may affect the size of backward\n*          errors of the decomposition.\n*\n*  K       (output) INTEGER\n*  L       (output) INTEGER\n*          On exit, K and L specify the dimension of the subblocks\n*          described in Purpose.\n*          K + L = effective numerical rank of (A',B')'.\n*\n*  U       (output) REAL array, dimension (LDU,M)\n*          If JOBU = 'U', U contains the orthogonal matrix U.\n*          If JOBU = 'N', U is not referenced.\n*\n*  LDU     (input) INTEGER\n*          The leading dimension of the array U. LDU >= max(1,M) if\n*          JOBU = 'U'; LDU >= 1 otherwise.\n*\n*  V       (output) REAL array, dimension (LDV,M)\n*          If JOBV = 'V', V contains the orthogonal matrix V.\n*          If JOBV = 'N', V is not referenced.\n*\n*  LDV     (input) INTEGER\n*          The leading dimension of the array V. LDV >= max(1,P) if\n*          JOBV = 'V'; LDV >= 1 otherwise.\n*\n*  Q       (output) REAL array, dimension (LDQ,N)\n*          If JOBQ = 'Q', Q contains the orthogonal matrix Q.\n*          If JOBQ = 'N', Q is not referenced.\n*\n*  LDQ     (input) INTEGER\n*          The leading dimension of the array Q. LDQ >= max(1,N) if\n*          JOBQ = 'Q'; LDQ >= 1 otherwise.\n*\n*  IWORK   (workspace) INTEGER array, dimension (N)\n*\n*  TAU     (workspace) REAL array, dimension (N)\n*\n*  WORK    (workspace) REAL array, dimension (max(3*N,M,P))\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit\n*          < 0:  if INFO = -i, the i-th argument had an illegal value.\n*\n*\n\n*  Further Details\n*  ===============\n*\n*  The subroutine uses LAPACK subroutine SGEQPF for the QR factorization\n*  with column pivoting to detect the effective numerical rank of the\n*  a matrix. It may be replaced by a better rank determination strategy.\n*\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 8)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 8)", argc);
  rb_jobu = argv[0];
  rb_jobv = argv[1];
  rb_jobq = argv[2];
  rb_p = argv[3];
  rb_a = argv[4];
  rb_b = argv[5];
  rb_tola = argv[6];
  rb_tolb = argv[7];

  jobu = StringValueCStr(rb_jobu)[0];
  jobv = StringValueCStr(rb_jobv)[0];
  jobq = StringValueCStr(rb_jobq)[0];
  p = NUM2INT(rb_p);
  tola = (real)NUM2DBL(rb_tola);
  tolb = (real)NUM2DBL(rb_tolb);
  if (!NA_IsNArray(rb_a))
    rb_raise(rb_eArgError, "a (5th argument) must be NArray");
  if (NA_RANK(rb_a) != 2)
    rb_raise(rb_eArgError, "rank of a (5th argument) must be %d", 2);
  lda = NA_SHAPE0(rb_a);
  n = NA_SHAPE1(rb_a);
  if (NA_TYPE(rb_a) != NA_SFLOAT)
    rb_a = na_change_type(rb_a, NA_SFLOAT);
  a = NA_PTR_TYPE(rb_a, real*);
  if (!NA_IsNArray(rb_b))
    rb_raise(rb_eArgError, "b (6th argument) must be NArray");
  if (NA_RANK(rb_b) != 2)
    rb_raise(rb_eArgError, "rank of b (6th argument) must be %d", 2);
  ldb = NA_SHAPE0(rb_b);
  if (NA_SHAPE1(rb_b) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of b must be the same as shape 1 of a");
  if (NA_TYPE(rb_b) != NA_SFLOAT)
    rb_b = na_change_type(rb_b, NA_SFLOAT);
  b = NA_PTR_TYPE(rb_b, real*);
  m = lda;
  ldu = lsame_(&jobu,"U") ? MAX(1,m) : 1;
  {
    int shape[2];
    shape[0] = ldu;
    shape[1] = m;
    rb_u = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  u = NA_PTR_TYPE(rb_u, real*);
  ldv = lsame_(&jobv,"V") ? MAX(1,p) : 1;
  {
    int shape[2];
    shape[0] = ldv;
    shape[1] = m;
    rb_v = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  v = NA_PTR_TYPE(rb_v, real*);
  ldq = lsame_(&jobq,"Q") ? MAX(1,n) : 1;
  {
    int shape[2];
    shape[0] = ldq;
    shape[1] = n;
    rb_q = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  q = NA_PTR_TYPE(rb_q, real*);
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
  iwork = ALLOC_N(integer, (n));
  tau = ALLOC_N(real, (n));
  work = ALLOC_N(real, (MAX(3*n,m)*(p)));

  sggsvp_(&jobu, &jobv, &jobq, &m, &p, &n, a, &lda, b, &ldb, &tola, &tolb, &k, &l, u, &ldu, v, &ldv, q, &ldq, iwork, tau, work, &info);

  free(iwork);
  free(tau);
  free(work);
  rb_k = INT2NUM(k);
  rb_l = INT2NUM(l);
  rb_info = INT2NUM(info);
  return rb_ary_new3(8, rb_k, rb_l, rb_u, rb_v, rb_q, rb_info, rb_a, rb_b);
}

void
init_lapack_sggsvp(VALUE mLapack){
  rb_define_module_function(mLapack, "sggsvp", rb_sggsvp, -1);
}

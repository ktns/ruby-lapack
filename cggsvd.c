#include "rb_lapack.h"

extern VOID cggsvd_(char *jobu, char *jobv, char *jobq, integer *m, integer *n, integer *p, integer *k, integer *l, complex *a, integer *lda, complex *b, integer *ldb, real *alpha, real *beta, complex *u, integer *ldu, complex *v, integer *ldv, complex *q, integer *ldq, complex *work, real *rwork, integer *iwork, integer *info);

static VALUE sHelp, sUsage;

static VALUE
rblapack_cggsvd(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_jobu;
  char jobu; 
  VALUE rblapack_jobv;
  char jobv; 
  VALUE rblapack_jobq;
  char jobq; 
  VALUE rblapack_a;
  complex *a; 
  VALUE rblapack_b;
  complex *b; 
  VALUE rblapack_k;
  integer k; 
  VALUE rblapack_l;
  integer l; 
  VALUE rblapack_alpha;
  real *alpha; 
  VALUE rblapack_beta;
  real *beta; 
  VALUE rblapack_u;
  complex *u; 
  VALUE rblapack_v;
  complex *v; 
  VALUE rblapack_q;
  complex *q; 
  VALUE rblapack_iwork;
  integer *iwork; 
  VALUE rblapack_info;
  integer info; 
  VALUE rblapack_a_out__;
  complex *a_out__;
  VALUE rblapack_b_out__;
  complex *b_out__;
  complex *work;
  real *rwork;

  integer lda;
  integer n;
  integer ldb;
  integer ldu;
  integer m;
  integer ldv;
  integer p;
  integer ldq;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  k, l, alpha, beta, u, v, q, iwork, info, a, b = NumRu::Lapack.cggsvd( jobu, jobv, jobq, a, b, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE CGGSVD( JOBU, JOBV, JOBQ, M, N, P, K, L, A, LDA, B, LDB, ALPHA, BETA, U, LDU, V, LDV, Q, LDQ, WORK, RWORK, IWORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  CGGSVD computes the generalized singular value decomposition (GSVD)\n*  of an M-by-N complex matrix A and P-by-N complex matrix B:\n*\n*        U'*A*Q = D1*( 0 R ),    V'*B*Q = D2*( 0 R )\n*\n*  where U, V and Q are unitary matrices, and Z' means the conjugate\n*  transpose of Z.  Let K+L = the effective numerical rank of the\n*  matrix (A',B')', then R is a (K+L)-by-(K+L) nonsingular upper\n*  triangular matrix, D1 and D2 are M-by-(K+L) and P-by-(K+L) \"diagonal\"\n*  matrices and of the following structures, respectively:\n*\n*  If M-K-L >= 0,\n*\n*                      K  L\n*         D1 =     K ( I  0 )\n*                  L ( 0  C )\n*              M-K-L ( 0  0 )\n*\n*                    K  L\n*         D2 =   L ( 0  S )\n*              P-L ( 0  0 )\n*\n*                  N-K-L  K    L\n*    ( 0 R ) = K (  0   R11  R12 )\n*              L (  0    0   R22 )\n*  where\n*\n*    C = diag( ALPHA(K+1), ... , ALPHA(K+L) ),\n*    S = diag( BETA(K+1),  ... , BETA(K+L) ),\n*    C**2 + S**2 = I.\n*\n*    R is stored in A(1:K+L,N-K-L+1:N) on exit.\n*\n*  If M-K-L < 0,\n*\n*                    K M-K K+L-M\n*         D1 =   K ( I  0    0   )\n*              M-K ( 0  C    0   )\n*\n*                      K M-K K+L-M\n*         D2 =   M-K ( 0  S    0  )\n*              K+L-M ( 0  0    I  )\n*                P-L ( 0  0    0  )\n*\n*                     N-K-L  K   M-K  K+L-M\n*    ( 0 R ) =     K ( 0    R11  R12  R13  )\n*                M-K ( 0     0   R22  R23  )\n*              K+L-M ( 0     0    0   R33  )\n*\n*  where\n*\n*    C = diag( ALPHA(K+1), ... , ALPHA(M) ),\n*    S = diag( BETA(K+1),  ... , BETA(M) ),\n*    C**2 + S**2 = I.\n*\n*    (R11 R12 R13 ) is stored in A(1:M, N-K-L+1:N), and R33 is stored\n*    ( 0  R22 R23 )\n*    in B(M-K+1:L,N+M-K-L+1:N) on exit.\n*\n*  The routine computes C, S, R, and optionally the unitary\n*  transformation matrices U, V and Q.\n*\n*  In particular, if B is an N-by-N nonsingular matrix, then the GSVD of\n*  A and B implicitly gives the SVD of A*inv(B):\n*                       A*inv(B) = U*(D1*inv(D2))*V'.\n*  If ( A',B')' has orthnormal columns, then the GSVD of A and B is also\n*  equal to the CS decomposition of A and B. Furthermore, the GSVD can\n*  be used to derive the solution of the eigenvalue problem:\n*                       A'*A x = lambda* B'*B x.\n*  In some literature, the GSVD of A and B is presented in the form\n*                   U'*A*X = ( 0 D1 ),   V'*B*X = ( 0 D2 )\n*  where U and V are orthogonal and X is nonsingular, and D1 and D2 are\n*  ``diagonal''.  The former GSVD form can be converted to the latter\n*  form by taking the nonsingular matrix X as\n*\n*                        X = Q*(  I   0    )\n*                              (  0 inv(R) )\n*\n\n*  Arguments\n*  =========\n*\n*  JOBU    (input) CHARACTER*1\n*          = 'U':  Unitary matrix U is computed;\n*          = 'N':  U is not computed.\n*\n*  JOBV    (input) CHARACTER*1\n*          = 'V':  Unitary matrix V is computed;\n*          = 'N':  V is not computed.\n*\n*  JOBQ    (input) CHARACTER*1\n*          = 'Q':  Unitary matrix Q is computed;\n*          = 'N':  Q is not computed.\n*\n*  M       (input) INTEGER\n*          The number of rows of the matrix A.  M >= 0.\n*\n*  N       (input) INTEGER\n*          The number of columns of the matrices A and B.  N >= 0.\n*\n*  P       (input) INTEGER\n*          The number of rows of the matrix B.  P >= 0.\n*\n*  K       (output) INTEGER\n*  L       (output) INTEGER\n*          On exit, K and L specify the dimension of the subblocks\n*          described in Purpose.\n*          K + L = effective numerical rank of (A',B')'.\n*\n*  A       (input/output) COMPLEX array, dimension (LDA,N)\n*          On entry, the M-by-N matrix A.\n*          On exit, A contains the triangular matrix R, or part of R.\n*          See Purpose for details.\n*\n*  LDA     (input) INTEGER\n*          The leading dimension of the array A. LDA >= max(1,M).\n*\n*  B       (input/output) COMPLEX array, dimension (LDB,N)\n*          On entry, the P-by-N matrix B.\n*          On exit, B contains part of the triangular matrix R if\n*          M-K-L < 0.  See Purpose for details.\n*\n*  LDB     (input) INTEGER\n*          The leading dimension of the array B. LDB >= max(1,P).\n*\n*  ALPHA   (output) REAL array, dimension (N)\n*  BETA    (output) REAL array, dimension (N)\n*          On exit, ALPHA and BETA contain the generalized singular\n*          value pairs of A and B;\n*            ALPHA(1:K) = 1,\n*            BETA(1:K)  = 0,\n*          and if M-K-L >= 0,\n*            ALPHA(K+1:K+L) = C,\n*            BETA(K+1:K+L)  = S,\n*          or if M-K-L < 0,\n*            ALPHA(K+1:M)= C, ALPHA(M+1:K+L)= 0\n*            BETA(K+1:M) = S, BETA(M+1:K+L) = 1\n*          and\n*            ALPHA(K+L+1:N) = 0\n*            BETA(K+L+1:N)  = 0\n*\n*  U       (output) COMPLEX array, dimension (LDU,M)\n*          If JOBU = 'U', U contains the M-by-M unitary matrix U.\n*          If JOBU = 'N', U is not referenced.\n*\n*  LDU     (input) INTEGER\n*          The leading dimension of the array U. LDU >= max(1,M) if\n*          JOBU = 'U'; LDU >= 1 otherwise.\n*\n*  V       (output) COMPLEX array, dimension (LDV,P)\n*          If JOBV = 'V', V contains the P-by-P unitary matrix V.\n*          If JOBV = 'N', V is not referenced.\n*\n*  LDV     (input) INTEGER\n*          The leading dimension of the array V. LDV >= max(1,P) if\n*          JOBV = 'V'; LDV >= 1 otherwise.\n*\n*  Q       (output) COMPLEX array, dimension (LDQ,N)\n*          If JOBQ = 'Q', Q contains the N-by-N unitary matrix Q.\n*          If JOBQ = 'N', Q is not referenced.\n*\n*  LDQ     (input) INTEGER\n*          The leading dimension of the array Q. LDQ >= max(1,N) if\n*          JOBQ = 'Q'; LDQ >= 1 otherwise.\n*\n*  WORK    (workspace) COMPLEX array, dimension (max(3*N,M,P)+N)\n*\n*  RWORK   (workspace) REAL array, dimension (2*N)\n*\n*  IWORK   (workspace/output) INTEGER array, dimension (N)\n*          On exit, IWORK stores the sorting information. More\n*          precisely, the following loop will sort ALPHA\n*             for I = K+1, min(M,K+L)\n*                 swap ALPHA(I) and ALPHA(IWORK(I))\n*             endfor\n*          such that ALPHA(1) >= ALPHA(2) >= ... >= ALPHA(N).\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit.\n*          < 0:  if INFO = -i, the i-th argument had an illegal value.\n*          > 0:  if INFO = 1, the Jacobi-type procedure failed to\n*                converge.  For further details, see subroutine CTGSJA.\n*\n*  Internal Parameters\n*  ===================\n*\n*  TOLA    REAL\n*  TOLB    REAL\n*          TOLA and TOLB are the thresholds to determine the effective\n*          rank of (A',B')'. Generally, they are set to\n*                   TOLA = MAX(M,N)*norm(A)*MACHEPS,\n*                   TOLB = MAX(P,N)*norm(B)*MACHEPS.\n*          The size of TOLA and TOLB may affect the size of backward\n*          errors of the decomposition.\n*\n\n*  Further Details\n*  ===============\n*\n*  2-96 Based on modifications by\n*     Ming Gu and Huan Ren, Computer Science Division, University of\n*     California at Berkeley, USA\n*\n*  =====================================================================\n*\n*     .. Local Scalars ..\n      LOGICAL            WANTQ, WANTU, WANTV\n      INTEGER            I, IBND, ISUB, J, NCYCLE\n      REAL               ANORM, BNORM, SMAX, TEMP, TOLA, TOLB, ULP, UNFL\n*     ..\n*     .. External Functions ..\n      LOGICAL            LSAME\n      REAL               CLANGE, SLAMCH\n      EXTERNAL           LSAME, CLANGE, SLAMCH\n*     ..\n*     .. External Subroutines ..\n      EXTERNAL           CGGSVP, CTGSJA, SCOPY, XERBLA\n*     ..\n*     .. Intrinsic Functions ..\n      INTRINSIC          MAX, MIN\n*     ..\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  k, l, alpha, beta, u, v, q, iwork, info, a, b = NumRu::Lapack.cggsvd( jobu, jobv, jobq, a, b, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 5)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 5)", argc);
  rblapack_jobu = argv[0];
  rblapack_jobv = argv[1];
  rblapack_jobq = argv[2];
  rblapack_a = argv[3];
  rblapack_b = argv[4];
  if (rb_options != Qnil) {
  }

  jobq = StringValueCStr(rblapack_jobq)[0];
  if (!NA_IsNArray(rblapack_a))
    rb_raise(rb_eArgError, "a (4th argument) must be NArray");
  if (NA_RANK(rblapack_a) != 2)
    rb_raise(rb_eArgError, "rank of a (4th argument) must be %d", 2);
  n = NA_SHAPE1(rblapack_a);
  lda = NA_SHAPE0(rblapack_a);
  if (lda != (m))
    rb_raise(rb_eRuntimeError, "shape 0 of a must be %d", m);
  m = lda;
  if (NA_TYPE(rblapack_a) != NA_SCOMPLEX)
    rblapack_a = na_change_type(rblapack_a, NA_SCOMPLEX);
  a = NA_PTR_TYPE(rblapack_a, complex*);
  if (!NA_IsNArray(rblapack_b))
    rb_raise(rb_eArgError, "b (5th argument) must be NArray");
  if (NA_RANK(rblapack_b) != 2)
    rb_raise(rb_eArgError, "rank of b (5th argument) must be %d", 2);
  if (NA_SHAPE1(rblapack_b) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of b must be the same as shape 1 of a");
  ldb = NA_SHAPE0(rblapack_b);
  if (ldb != (p))
    rb_raise(rb_eRuntimeError, "shape 0 of b must be %d", p);
  p = ldb;
  if (NA_TYPE(rblapack_b) != NA_SCOMPLEX)
    rblapack_b = na_change_type(rblapack_b, NA_SCOMPLEX);
  b = NA_PTR_TYPE(rblapack_b, complex*);
  jobu = StringValueCStr(rblapack_jobu)[0];
  jobv = StringValueCStr(rblapack_jobv)[0];
  ldq = lsame_(&jobq,"Q") ? MAX(1,n) : 1;
  ldu = lsame_(&jobu,"U") ? MAX(1,m) : 1;
  ldv = lsame_(&jobv,"V") ? MAX(1,p) : 1;
  lda = m;
  ldb = p;
  {
    int shape[1];
    shape[0] = n;
    rblapack_alpha = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  alpha = NA_PTR_TYPE(rblapack_alpha, real*);
  {
    int shape[1];
    shape[0] = n;
    rblapack_beta = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  beta = NA_PTR_TYPE(rblapack_beta, real*);
  {
    int shape[2];
    shape[0] = ldu;
    shape[1] = m;
    rblapack_u = na_make_object(NA_SCOMPLEX, 2, shape, cNArray);
  }
  u = NA_PTR_TYPE(rblapack_u, complex*);
  {
    int shape[2];
    shape[0] = ldv;
    shape[1] = p;
    rblapack_v = na_make_object(NA_SCOMPLEX, 2, shape, cNArray);
  }
  v = NA_PTR_TYPE(rblapack_v, complex*);
  {
    int shape[2];
    shape[0] = ldq;
    shape[1] = n;
    rblapack_q = na_make_object(NA_SCOMPLEX, 2, shape, cNArray);
  }
  q = NA_PTR_TYPE(rblapack_q, complex*);
  {
    int shape[1];
    shape[0] = n;
    rblapack_iwork = na_make_object(NA_LINT, 1, shape, cNArray);
  }
  iwork = NA_PTR_TYPE(rblapack_iwork, integer*);
  {
    int shape[2];
    shape[0] = lda;
    shape[1] = n;
    rblapack_a_out__ = na_make_object(NA_SCOMPLEX, 2, shape, cNArray);
  }
  a_out__ = NA_PTR_TYPE(rblapack_a_out__, complex*);
  MEMCPY(a_out__, a, complex, NA_TOTAL(rblapack_a));
  rblapack_a = rblapack_a_out__;
  a = a_out__;
  {
    int shape[2];
    shape[0] = ldb;
    shape[1] = n;
    rblapack_b_out__ = na_make_object(NA_SCOMPLEX, 2, shape, cNArray);
  }
  b_out__ = NA_PTR_TYPE(rblapack_b_out__, complex*);
  MEMCPY(b_out__, b, complex, NA_TOTAL(rblapack_b));
  rblapack_b = rblapack_b_out__;
  b = b_out__;
  work = ALLOC_N(complex, (MAX(3*n,m)*(p)+n));
  rwork = ALLOC_N(real, (2*n));

  cggsvd_(&jobu, &jobv, &jobq, &m, &n, &p, &k, &l, a, &lda, b, &ldb, alpha, beta, u, &ldu, v, &ldv, q, &ldq, work, rwork, iwork, &info);

  free(work);
  free(rwork);
  rblapack_k = INT2NUM(k);
  rblapack_l = INT2NUM(l);
  rblapack_info = INT2NUM(info);
  return rb_ary_new3(11, rblapack_k, rblapack_l, rblapack_alpha, rblapack_beta, rblapack_u, rblapack_v, rblapack_q, rblapack_iwork, rblapack_info, rblapack_a, rblapack_b);
}

void
init_lapack_cggsvd(VALUE mLapack){
  rb_define_module_function(mLapack, "cggsvd", rblapack_cggsvd, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}

#include "rb_lapack.h"

extern VOID dggsvp_(char *jobu, char *jobv, char *jobq, integer *m, integer *p, integer *n, doublereal *a, integer *lda, doublereal *b, integer *ldb, doublereal *tola, doublereal *tolb, integer *k, integer *l, doublereal *u, integer *ldu, doublereal *v, integer *ldv, doublereal *q, integer *ldq, integer *iwork, doublereal *tau, doublereal *work, integer *info);

static VALUE sHelp, sUsage;

static VALUE
rblapack_dggsvp(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_jobu;
  char jobu; 
  VALUE rblapack_jobv;
  char jobv; 
  VALUE rblapack_jobq;
  char jobq; 
  VALUE rblapack_a;
  doublereal *a; 
  VALUE rblapack_b;
  doublereal *b; 
  VALUE rblapack_tola;
  doublereal tola; 
  VALUE rblapack_tolb;
  doublereal tolb; 
  VALUE rblapack_k;
  integer k; 
  VALUE rblapack_l;
  integer l; 
  VALUE rblapack_u;
  doublereal *u; 
  VALUE rblapack_v;
  doublereal *v; 
  VALUE rblapack_q;
  doublereal *q; 
  VALUE rblapack_info;
  integer info; 
  VALUE rblapack_a_out__;
  doublereal *a_out__;
  VALUE rblapack_b_out__;
  doublereal *b_out__;
  integer *iwork;
  doublereal *tau;
  doublereal *work;

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
      printf("%s\n", "USAGE:\n  k, l, u, v, q, info, a, b = NumRu::Lapack.dggsvp( jobu, jobv, jobq, a, b, tola, tolb, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE DGGSVP( JOBU, JOBV, JOBQ, M, P, N, A, LDA, B, LDB, TOLA, TOLB, K, L, U, LDU, V, LDV, Q, LDQ, IWORK, TAU, WORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  DGGSVP computes orthogonal matrices U, V and Q such that\n*\n*                   N-K-L  K    L\n*   U'*A*Q =     K ( 0    A12  A13 )  if M-K-L >= 0;\n*                L ( 0     0   A23 )\n*            M-K-L ( 0     0    0  )\n*\n*                   N-K-L  K    L\n*          =     K ( 0    A12  A13 )  if M-K-L < 0;\n*              M-K ( 0     0   A23 )\n*\n*                 N-K-L  K    L\n*   V'*B*Q =   L ( 0     0   B13 )\n*            P-L ( 0     0    0  )\n*\n*  where the K-by-K matrix A12 and L-by-L matrix B13 are nonsingular\n*  upper triangular; A23 is L-by-L upper triangular if M-K-L >= 0,\n*  otherwise A23 is (M-K)-by-L upper trapezoidal.  K+L = the effective\n*  numerical rank of the (M+P)-by-N matrix (A',B')'.  Z' denotes the\n*  transpose of Z.\n*\n*  This decomposition is the preprocessing step for computing the\n*  Generalized Singular Value Decomposition (GSVD), see subroutine\n*  DGGSVD.\n*\n\n*  Arguments\n*  =========\n*\n*  JOBU    (input) CHARACTER*1\n*          = 'U':  Orthogonal matrix U is computed;\n*          = 'N':  U is not computed.\n*\n*  JOBV    (input) CHARACTER*1\n*          = 'V':  Orthogonal matrix V is computed;\n*          = 'N':  V is not computed.\n*\n*  JOBQ    (input) CHARACTER*1\n*          = 'Q':  Orthogonal matrix Q is computed;\n*          = 'N':  Q is not computed.\n*\n*  M       (input) INTEGER\n*          The number of rows of the matrix A.  M >= 0.\n*\n*  P       (input) INTEGER\n*          The number of rows of the matrix B.  P >= 0.\n*\n*  N       (input) INTEGER\n*          The number of columns of the matrices A and B.  N >= 0.\n*\n*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)\n*          On entry, the M-by-N matrix A.\n*          On exit, A contains the triangular (or trapezoidal) matrix\n*          described in the Purpose section.\n*\n*  LDA     (input) INTEGER\n*          The leading dimension of the array A. LDA >= max(1,M).\n*\n*  B       (input/output) DOUBLE PRECISION array, dimension (LDB,N)\n*          On entry, the P-by-N matrix B.\n*          On exit, B contains the triangular matrix described in\n*          the Purpose section.\n*\n*  LDB     (input) INTEGER\n*          The leading dimension of the array B. LDB >= max(1,P).\n*\n*  TOLA    (input) DOUBLE PRECISION\n*  TOLB    (input) DOUBLE PRECISION\n*          TOLA and TOLB are the thresholds to determine the effective\n*          numerical rank of matrix B and a subblock of A. Generally,\n*          they are set to\n*             TOLA = MAX(M,N)*norm(A)*MAZHEPS,\n*             TOLB = MAX(P,N)*norm(B)*MAZHEPS.\n*          The size of TOLA and TOLB may affect the size of backward\n*          errors of the decomposition.\n*\n*  K       (output) INTEGER\n*  L       (output) INTEGER\n*          On exit, K and L specify the dimension of the subblocks\n*          described in Purpose.\n*          K + L = effective numerical rank of (A',B')'.\n*\n*  U       (output) DOUBLE PRECISION array, dimension (LDU,M)\n*          If JOBU = 'U', U contains the orthogonal matrix U.\n*          If JOBU = 'N', U is not referenced.\n*\n*  LDU     (input) INTEGER\n*          The leading dimension of the array U. LDU >= max(1,M) if\n*          JOBU = 'U'; LDU >= 1 otherwise.\n*\n*  V       (output) DOUBLE PRECISION array, dimension (LDV,P)\n*          If JOBV = 'V', V contains the orthogonal matrix V.\n*          If JOBV = 'N', V is not referenced.\n*\n*  LDV     (input) INTEGER\n*          The leading dimension of the array V. LDV >= max(1,P) if\n*          JOBV = 'V'; LDV >= 1 otherwise.\n*\n*  Q       (output) DOUBLE PRECISION array, dimension (LDQ,N)\n*          If JOBQ = 'Q', Q contains the orthogonal matrix Q.\n*          If JOBQ = 'N', Q is not referenced.\n*\n*  LDQ     (input) INTEGER\n*          The leading dimension of the array Q. LDQ >= max(1,N) if\n*          JOBQ = 'Q'; LDQ >= 1 otherwise.\n*\n*  IWORK   (workspace) INTEGER array, dimension (N)\n*\n*  TAU     (workspace) DOUBLE PRECISION array, dimension (N)\n*\n*  WORK    (workspace) DOUBLE PRECISION array, dimension (max(3*N,M,P))\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit\n*          < 0:  if INFO = -i, the i-th argument had an illegal value.\n*\n*\n\n*  Further Details\n*  ===============\n*\n*  The subroutine uses LAPACK subroutine DGEQPF for the QR factorization\n*  with column pivoting to detect the effective numerical rank of the\n*  a matrix. It may be replaced by a better rank determination strategy.\n*\n*  =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  k, l, u, v, q, info, a, b = NumRu::Lapack.dggsvp( jobu, jobv, jobq, a, b, tola, tolb, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 7)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 7)", argc);
  rblapack_jobu = argv[0];
  rblapack_jobv = argv[1];
  rblapack_jobq = argv[2];
  rblapack_a = argv[3];
  rblapack_b = argv[4];
  rblapack_tola = argv[5];
  rblapack_tolb = argv[6];
  if (rb_options != Qnil) {
  }

  jobq = StringValueCStr(rblapack_jobq)[0];
  if (!NA_IsNArray(rblapack_a))
    rb_raise(rb_eArgError, "a (4th argument) must be NArray");
  if (NA_RANK(rblapack_a) != 2)
    rb_raise(rb_eArgError, "rank of a (4th argument) must be %d", 2);
  n = NA_SHAPE1(rblapack_a);
  lda = NA_SHAPE0(rblapack_a);
  if (NA_TYPE(rblapack_a) != NA_DFLOAT)
    rblapack_a = na_change_type(rblapack_a, NA_DFLOAT);
  a = NA_PTR_TYPE(rblapack_a, doublereal*);
  tolb = NUM2DBL(rblapack_tolb);
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
  if (NA_TYPE(rblapack_b) != NA_DFLOAT)
    rblapack_b = na_change_type(rblapack_b, NA_DFLOAT);
  b = NA_PTR_TYPE(rblapack_b, doublereal*);
  tola = NUM2DBL(rblapack_tola);
  jobu = StringValueCStr(rblapack_jobu)[0];
  jobv = StringValueCStr(rblapack_jobv)[0];
  ldu = lsame_(&jobu,"U") ? MAX(1,m) : 1;
  ldq = lsame_(&jobq,"Q") ? MAX(1,n) : 1;
  ldv = lsame_(&jobv,"V") ? MAX(1,p) : 1;
  m = lda;
  ldb = p;
  {
    int shape[2];
    shape[0] = ldu;
    shape[1] = m;
    rblapack_u = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  u = NA_PTR_TYPE(rblapack_u, doublereal*);
  {
    int shape[2];
    shape[0] = ldv;
    shape[1] = p;
    rblapack_v = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  v = NA_PTR_TYPE(rblapack_v, doublereal*);
  {
    int shape[2];
    shape[0] = ldq;
    shape[1] = n;
    rblapack_q = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  q = NA_PTR_TYPE(rblapack_q, doublereal*);
  {
    int shape[2];
    shape[0] = lda;
    shape[1] = n;
    rblapack_a_out__ = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  a_out__ = NA_PTR_TYPE(rblapack_a_out__, doublereal*);
  MEMCPY(a_out__, a, doublereal, NA_TOTAL(rblapack_a));
  rblapack_a = rblapack_a_out__;
  a = a_out__;
  {
    int shape[2];
    shape[0] = ldb;
    shape[1] = n;
    rblapack_b_out__ = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  b_out__ = NA_PTR_TYPE(rblapack_b_out__, doublereal*);
  MEMCPY(b_out__, b, doublereal, NA_TOTAL(rblapack_b));
  rblapack_b = rblapack_b_out__;
  b = b_out__;
  iwork = ALLOC_N(integer, (n));
  tau = ALLOC_N(doublereal, (n));
  work = ALLOC_N(doublereal, (MAX(MAX(3*n,m),p)));

  dggsvp_(&jobu, &jobv, &jobq, &m, &p, &n, a, &lda, b, &ldb, &tola, &tolb, &k, &l, u, &ldu, v, &ldv, q, &ldq, iwork, tau, work, &info);

  free(iwork);
  free(tau);
  free(work);
  rblapack_k = INT2NUM(k);
  rblapack_l = INT2NUM(l);
  rblapack_info = INT2NUM(info);
  return rb_ary_new3(8, rblapack_k, rblapack_l, rblapack_u, rblapack_v, rblapack_q, rblapack_info, rblapack_a, rblapack_b);
}

void
init_lapack_dggsvp(VALUE mLapack){
  rb_define_module_function(mLapack, "dggsvp", rblapack_dggsvp, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}

#include "rb_lapack.h"

extern VOID dtgsja_(char *jobu, char *jobv, char *jobq, integer *m, integer *p, integer *n, integer *k, integer *l, doublereal *a, integer *lda, doublereal *b, integer *ldb, doublereal *tola, doublereal *tolb, doublereal *alpha, doublereal *beta, doublereal *u, integer *ldu, doublereal *v, integer *ldv, doublereal *q, integer *ldq, doublereal *work, integer *ncycle, integer *info);

static VALUE sHelp, sUsage;

static VALUE
rblapack_dtgsja(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_jobu;
  char jobu; 
  VALUE rblapack_jobv;
  char jobv; 
  VALUE rblapack_jobq;
  char jobq; 
  VALUE rblapack_k;
  integer k; 
  VALUE rblapack_l;
  integer l; 
  VALUE rblapack_a;
  doublereal *a; 
  VALUE rblapack_b;
  doublereal *b; 
  VALUE rblapack_tola;
  doublereal tola; 
  VALUE rblapack_tolb;
  doublereal tolb; 
  VALUE rblapack_u;
  doublereal *u; 
  VALUE rblapack_v;
  doublereal *v; 
  VALUE rblapack_q;
  doublereal *q; 
  VALUE rblapack_alpha;
  doublereal *alpha; 
  VALUE rblapack_beta;
  doublereal *beta; 
  VALUE rblapack_ncycle;
  integer ncycle; 
  VALUE rblapack_info;
  integer info; 
  VALUE rblapack_a_out__;
  doublereal *a_out__;
  VALUE rblapack_b_out__;
  doublereal *b_out__;
  VALUE rblapack_u_out__;
  doublereal *u_out__;
  VALUE rblapack_v_out__;
  doublereal *v_out__;
  VALUE rblapack_q_out__;
  doublereal *q_out__;
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
      printf("%s\n", "USAGE:\n  alpha, beta, ncycle, info, a, b, u, v, q = NumRu::Lapack.dtgsja( jobu, jobv, jobq, k, l, a, b, tola, tolb, u, v, q, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE DTGSJA( JOBU, JOBV, JOBQ, M, P, N, K, L, A, LDA, B, LDB, TOLA, TOLB, ALPHA, BETA, U, LDU, V, LDV, Q, LDQ, WORK, NCYCLE, INFO )\n\n*  Purpose\n*  =======\n*\n*  DTGSJA computes the generalized singular value decomposition (GSVD)\n*  of two real upper triangular (or trapezoidal) matrices A and B.\n*\n*  On entry, it is assumed that matrices A and B have the following\n*  forms, which may be obtained by the preprocessing subroutine DGGSVP\n*  from a general M-by-N matrix A and P-by-N matrix B:\n*\n*               N-K-L  K    L\n*     A =    K ( 0    A12  A13 ) if M-K-L >= 0;\n*            L ( 0     0   A23 )\n*        M-K-L ( 0     0    0  )\n*\n*             N-K-L  K    L\n*     A =  K ( 0    A12  A13 ) if M-K-L < 0;\n*        M-K ( 0     0   A23 )\n*\n*             N-K-L  K    L\n*     B =  L ( 0     0   B13 )\n*        P-L ( 0     0    0  )\n*\n*  where the K-by-K matrix A12 and L-by-L matrix B13 are nonsingular\n*  upper triangular; A23 is L-by-L upper triangular if M-K-L >= 0,\n*  otherwise A23 is (M-K)-by-L upper trapezoidal.\n*\n*  On exit,\n*\n*              U'*A*Q = D1*( 0 R ),    V'*B*Q = D2*( 0 R ),\n*\n*  where U, V and Q are orthogonal matrices, Z' denotes the transpose\n*  of Z, R is a nonsingular upper triangular matrix, and D1 and D2 are\n*  ``diagonal'' matrices, which are of the following structures:\n*\n*  If M-K-L >= 0,\n*\n*                      K  L\n*         D1 =     K ( I  0 )\n*                  L ( 0  C )\n*              M-K-L ( 0  0 )\n*\n*                    K  L\n*         D2 = L   ( 0  S )\n*              P-L ( 0  0 )\n*\n*                 N-K-L  K    L\n*    ( 0 R ) = K (  0   R11  R12 ) K\n*              L (  0    0   R22 ) L\n*\n*  where\n*\n*    C = diag( ALPHA(K+1), ... , ALPHA(K+L) ),\n*    S = diag( BETA(K+1),  ... , BETA(K+L) ),\n*    C**2 + S**2 = I.\n*\n*    R is stored in A(1:K+L,N-K-L+1:N) on exit.\n*\n*  If M-K-L < 0,\n*\n*                 K M-K K+L-M\n*      D1 =   K ( I  0    0   )\n*           M-K ( 0  C    0   )\n*\n*                   K M-K K+L-M\n*      D2 =   M-K ( 0  S    0   )\n*           K+L-M ( 0  0    I   )\n*             P-L ( 0  0    0   )\n*\n*                 N-K-L  K   M-K  K+L-M\n* ( 0 R ) =    K ( 0    R11  R12  R13  )\n*            M-K ( 0     0   R22  R23  )\n*          K+L-M ( 0     0    0   R33  )\n*\n*  where\n*  C = diag( ALPHA(K+1), ... , ALPHA(M) ),\n*  S = diag( BETA(K+1),  ... , BETA(M) ),\n*  C**2 + S**2 = I.\n*\n*  R = ( R11 R12 R13 ) is stored in A(1:M, N-K-L+1:N) and R33 is stored\n*      (  0  R22 R23 )\n*  in B(M-K+1:L,N+M-K-L+1:N) on exit.\n*\n*  The computation of the orthogonal transformation matrices U, V or Q\n*  is optional.  These matrices may either be formed explicitly, or they\n*  may be postmultiplied into input matrices U1, V1, or Q1.\n*\n\n*  Arguments\n*  =========\n*\n*  JOBU    (input) CHARACTER*1\n*          = 'U':  U must contain an orthogonal matrix U1 on entry, and\n*                  the product U1*U is returned;\n*          = 'I':  U is initialized to the unit matrix, and the\n*                  orthogonal matrix U is returned;\n*          = 'N':  U is not computed.\n*\n*  JOBV    (input) CHARACTER*1\n*          = 'V':  V must contain an orthogonal matrix V1 on entry, and\n*                  the product V1*V is returned;\n*          = 'I':  V is initialized to the unit matrix, and the\n*                  orthogonal matrix V is returned;\n*          = 'N':  V is not computed.\n*\n*  JOBQ    (input) CHARACTER*1\n*          = 'Q':  Q must contain an orthogonal matrix Q1 on entry, and\n*                  the product Q1*Q is returned;\n*          = 'I':  Q is initialized to the unit matrix, and the\n*                  orthogonal matrix Q is returned;\n*          = 'N':  Q is not computed.\n*\n*  M       (input) INTEGER\n*          The number of rows of the matrix A.  M >= 0.\n*\n*  P       (input) INTEGER\n*          The number of rows of the matrix B.  P >= 0.\n*\n*  N       (input) INTEGER\n*          The number of columns of the matrices A and B.  N >= 0.\n*\n*  K       (input) INTEGER\n*  L       (input) INTEGER\n*          K and L specify the subblocks in the input matrices A and B:\n*          A23 = A(K+1:MIN(K+L,M),N-L+1:N) and B13 = B(1:L,N-L+1:N)\n*          of A and B, whose GSVD is going to be computed by DTGSJA.\n*          See Further Details.\n*\n*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)\n*          On entry, the M-by-N matrix A.\n*          On exit, A(N-K+1:N,1:MIN(K+L,M) ) contains the triangular\n*          matrix R or part of R.  See Purpose for details.\n*\n*  LDA     (input) INTEGER\n*          The leading dimension of the array A. LDA >= max(1,M).\n*\n*  B       (input/output) DOUBLE PRECISION array, dimension (LDB,N)\n*          On entry, the P-by-N matrix B.\n*          On exit, if necessary, B(M-K+1:L,N+M-K-L+1:N) contains\n*          a part of R.  See Purpose for details.\n*\n*  LDB     (input) INTEGER\n*          The leading dimension of the array B. LDB >= max(1,P).\n*\n*  TOLA    (input) DOUBLE PRECISION\n*  TOLB    (input) DOUBLE PRECISION\n*          TOLA and TOLB are the convergence criteria for the Jacobi-\n*          Kogbetliantz iteration procedure. Generally, they are the\n*          same as used in the preprocessing step, say\n*              TOLA = max(M,N)*norm(A)*MAZHEPS,\n*              TOLB = max(P,N)*norm(B)*MAZHEPS.\n*\n*  ALPHA   (output) DOUBLE PRECISION array, dimension (N)\n*  BETA    (output) DOUBLE PRECISION array, dimension (N)\n*          On exit, ALPHA and BETA contain the generalized singular\n*          value pairs of A and B;\n*            ALPHA(1:K) = 1,\n*            BETA(1:K)  = 0,\n*          and if M-K-L >= 0,\n*            ALPHA(K+1:K+L) = diag(C),\n*            BETA(K+1:K+L)  = diag(S),\n*          or if M-K-L < 0,\n*            ALPHA(K+1:M)= C, ALPHA(M+1:K+L)= 0\n*            BETA(K+1:M) = S, BETA(M+1:K+L) = 1.\n*          Furthermore, if K+L < N,\n*            ALPHA(K+L+1:N) = 0 and\n*            BETA(K+L+1:N)  = 0.\n*\n*  U       (input/output) DOUBLE PRECISION array, dimension (LDU,M)\n*          On entry, if JOBU = 'U', U must contain a matrix U1 (usually\n*          the orthogonal matrix returned by DGGSVP).\n*          On exit,\n*          if JOBU = 'I', U contains the orthogonal matrix U;\n*          if JOBU = 'U', U contains the product U1*U.\n*          If JOBU = 'N', U is not referenced.\n*\n*  LDU     (input) INTEGER\n*          The leading dimension of the array U. LDU >= max(1,M) if\n*          JOBU = 'U'; LDU >= 1 otherwise.\n*\n*  V       (input/output) DOUBLE PRECISION array, dimension (LDV,P)\n*          On entry, if JOBV = 'V', V must contain a matrix V1 (usually\n*          the orthogonal matrix returned by DGGSVP).\n*          On exit,\n*          if JOBV = 'I', V contains the orthogonal matrix V;\n*          if JOBV = 'V', V contains the product V1*V.\n*          If JOBV = 'N', V is not referenced.\n*\n*  LDV     (input) INTEGER\n*          The leading dimension of the array V. LDV >= max(1,P) if\n*          JOBV = 'V'; LDV >= 1 otherwise.\n*\n*  Q       (input/output) DOUBLE PRECISION array, dimension (LDQ,N)\n*          On entry, if JOBQ = 'Q', Q must contain a matrix Q1 (usually\n*          the orthogonal matrix returned by DGGSVP).\n*          On exit,\n*          if JOBQ = 'I', Q contains the orthogonal matrix Q;\n*          if JOBQ = 'Q', Q contains the product Q1*Q.\n*          If JOBQ = 'N', Q is not referenced.\n*\n*  LDQ     (input) INTEGER\n*          The leading dimension of the array Q. LDQ >= max(1,N) if\n*          JOBQ = 'Q'; LDQ >= 1 otherwise.\n*\n*  WORK    (workspace) DOUBLE PRECISION array, dimension (2*N)\n*\n*  NCYCLE  (output) INTEGER\n*          The number of cycles required for convergence.\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit\n*          < 0:  if INFO = -i, the i-th argument had an illegal value.\n*          = 1:  the procedure does not converge after MAXIT cycles.\n*\n*  Internal Parameters\n*  ===================\n*\n*  MAXIT   INTEGER\n*          MAXIT specifies the total loops that the iterative procedure\n*          may take. If after MAXIT cycles, the routine fails to\n*          converge, we return INFO = 1.\n*\n\n*  Further Details\n*  ===============\n*\n*  DTGSJA essentially uses a variant of Kogbetliantz algorithm to reduce\n*  min(L,M-K)-by-L triangular (or trapezoidal) matrix A23 and L-by-L\n*  matrix B13 to the form:\n*\n*           U1'*A13*Q1 = C1*R1; V1'*B13*Q1 = S1*R1,\n*\n*  where U1, V1 and Q1 are orthogonal matrix, and Z' is the transpose\n*  of Z.  C1 and S1 are diagonal matrices satisfying\n*\n*                C1**2 + S1**2 = I,\n*\n*  and R1 is an L-by-L nonsingular upper triangular matrix.\n*\n*  =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  alpha, beta, ncycle, info, a, b, u, v, q = NumRu::Lapack.dtgsja( jobu, jobv, jobq, k, l, a, b, tola, tolb, u, v, q, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 12)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 12)", argc);
  rblapack_jobu = argv[0];
  rblapack_jobv = argv[1];
  rblapack_jobq = argv[2];
  rblapack_k = argv[3];
  rblapack_l = argv[4];
  rblapack_a = argv[5];
  rblapack_b = argv[6];
  rblapack_tola = argv[7];
  rblapack_tolb = argv[8];
  rblapack_u = argv[9];
  rblapack_v = argv[10];
  rblapack_q = argv[11];
  if (rb_options != Qnil) {
  }

  if (!NA_IsNArray(rblapack_v))
    rb_raise(rb_eArgError, "v (11th argument) must be NArray");
  if (NA_RANK(rblapack_v) != 2)
    rb_raise(rb_eArgError, "rank of v (11th argument) must be %d", 2);
  p = NA_SHAPE1(rblapack_v);
  ldv = NA_SHAPE0(rblapack_v);
  if (NA_TYPE(rblapack_v) != NA_DFLOAT)
    rblapack_v = na_change_type(rblapack_v, NA_DFLOAT);
  v = NA_PTR_TYPE(rblapack_v, doublereal*);
  k = NUM2INT(rblapack_k);
  jobq = StringValueCStr(rblapack_jobq)[0];
  if (!NA_IsNArray(rblapack_a))
    rb_raise(rb_eArgError, "a (6th argument) must be NArray");
  if (NA_RANK(rblapack_a) != 2)
    rb_raise(rb_eArgError, "rank of a (6th argument) must be %d", 2);
  n = NA_SHAPE1(rblapack_a);
  lda = NA_SHAPE0(rblapack_a);
  if (NA_TYPE(rblapack_a) != NA_DFLOAT)
    rblapack_a = na_change_type(rblapack_a, NA_DFLOAT);
  a = NA_PTR_TYPE(rblapack_a, doublereal*);
  l = NUM2INT(rblapack_l);
  if (!NA_IsNArray(rblapack_b))
    rb_raise(rb_eArgError, "b (7th argument) must be NArray");
  if (NA_RANK(rblapack_b) != 2)
    rb_raise(rb_eArgError, "rank of b (7th argument) must be %d", 2);
  if (NA_SHAPE1(rblapack_b) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of b must be the same as shape 1 of a");
  ldb = NA_SHAPE0(rblapack_b);
  if (NA_TYPE(rblapack_b) != NA_DFLOAT)
    rblapack_b = na_change_type(rblapack_b, NA_DFLOAT);
  b = NA_PTR_TYPE(rblapack_b, doublereal*);
  jobu = StringValueCStr(rblapack_jobu)[0];
  jobv = StringValueCStr(rblapack_jobv)[0];
  if (!NA_IsNArray(rblapack_q))
    rb_raise(rb_eArgError, "q (12th argument) must be NArray");
  if (NA_RANK(rblapack_q) != 2)
    rb_raise(rb_eArgError, "rank of q (12th argument) must be %d", 2);
  if (NA_SHAPE1(rblapack_q) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of q must be the same as shape 1 of a");
  ldq = NA_SHAPE0(rblapack_q);
  if (NA_TYPE(rblapack_q) != NA_DFLOAT)
    rblapack_q = na_change_type(rblapack_q, NA_DFLOAT);
  q = NA_PTR_TYPE(rblapack_q, doublereal*);
  tola = NUM2DBL(rblapack_tola);
  tolb = NUM2DBL(rblapack_tolb);
  if (!NA_IsNArray(rblapack_u))
    rb_raise(rb_eArgError, "u (10th argument) must be NArray");
  if (NA_RANK(rblapack_u) != 2)
    rb_raise(rb_eArgError, "rank of u (10th argument) must be %d", 2);
  m = NA_SHAPE1(rblapack_u);
  ldu = NA_SHAPE0(rblapack_u);
  if (NA_TYPE(rblapack_u) != NA_DFLOAT)
    rblapack_u = na_change_type(rblapack_u, NA_DFLOAT);
  u = NA_PTR_TYPE(rblapack_u, doublereal*);
  {
    int shape[1];
    shape[0] = n;
    rblapack_alpha = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  alpha = NA_PTR_TYPE(rblapack_alpha, doublereal*);
  {
    int shape[1];
    shape[0] = n;
    rblapack_beta = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  beta = NA_PTR_TYPE(rblapack_beta, doublereal*);
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
  {
    int shape[2];
    shape[0] = ldu;
    shape[1] = m;
    rblapack_u_out__ = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  u_out__ = NA_PTR_TYPE(rblapack_u_out__, doublereal*);
  MEMCPY(u_out__, u, doublereal, NA_TOTAL(rblapack_u));
  rblapack_u = rblapack_u_out__;
  u = u_out__;
  {
    int shape[2];
    shape[0] = ldv;
    shape[1] = p;
    rblapack_v_out__ = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  v_out__ = NA_PTR_TYPE(rblapack_v_out__, doublereal*);
  MEMCPY(v_out__, v, doublereal, NA_TOTAL(rblapack_v));
  rblapack_v = rblapack_v_out__;
  v = v_out__;
  {
    int shape[2];
    shape[0] = ldq;
    shape[1] = n;
    rblapack_q_out__ = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  q_out__ = NA_PTR_TYPE(rblapack_q_out__, doublereal*);
  MEMCPY(q_out__, q, doublereal, NA_TOTAL(rblapack_q));
  rblapack_q = rblapack_q_out__;
  q = q_out__;
  work = ALLOC_N(doublereal, (2*n));

  dtgsja_(&jobu, &jobv, &jobq, &m, &p, &n, &k, &l, a, &lda, b, &ldb, &tola, &tolb, alpha, beta, u, &ldu, v, &ldv, q, &ldq, work, &ncycle, &info);

  free(work);
  rblapack_ncycle = INT2NUM(ncycle);
  rblapack_info = INT2NUM(info);
  return rb_ary_new3(9, rblapack_alpha, rblapack_beta, rblapack_ncycle, rblapack_info, rblapack_a, rblapack_b, rblapack_u, rblapack_v, rblapack_q);
}

void
init_lapack_dtgsja(VALUE mLapack){
  rb_define_module_function(mLapack, "dtgsja", rblapack_dtgsja, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}

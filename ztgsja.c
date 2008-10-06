#include "rb_lapack.h"

static VALUE
rb_ztgsja(int argc, VALUE *argv, VALUE self){
  VALUE rb_jobu;
  char jobu; 
  VALUE rb_jobv;
  char jobv; 
  VALUE rb_jobq;
  char jobq; 
  VALUE rb_k;
  integer k; 
  VALUE rb_l;
  integer l; 
  VALUE rb_a;
  doublecomplex *a; 
  VALUE rb_b;
  doublecomplex *b; 
  VALUE rb_tola;
  doublereal tola; 
  VALUE rb_tolb;
  doublereal tolb; 
  VALUE rb_u;
  doublecomplex *u; 
  VALUE rb_v;
  doublecomplex *v; 
  VALUE rb_q;
  doublecomplex *q; 
  VALUE rb_alpha;
  doublereal *alpha; 
  VALUE rb_beta;
  doublereal *beta; 
  VALUE rb_ncycle;
  integer ncycle; 
  VALUE rb_info;
  integer info; 
  VALUE rb_a_out__;
  doublecomplex *a_out__;
  VALUE rb_b_out__;
  doublecomplex *b_out__;
  VALUE rb_u_out__;
  doublecomplex *u_out__;
  VALUE rb_v_out__;
  doublecomplex *v_out__;
  VALUE rb_q_out__;
  doublecomplex *q_out__;
  doublecomplex *work;

  integer lda;
  integer n;
  integer ldb;
  integer ldu;
  integer m;
  integer ldv;
  integer p;
  integer ldq;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  alpha, beta, ncycle, info, a, b, u, v, q = NumRu::Lapack.ztgsja( jobu, jobv, jobq, k, l, a, b, tola, tolb, u, v, q)\n    or\n  NumRu::Lapack.ztgsja  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE ZTGSJA( JOBU, JOBV, JOBQ, M, P, N, K, L, A, LDA, B, LDB, TOLA, TOLB, ALPHA, BETA, U, LDU, V, LDV, Q, LDQ, WORK, NCYCLE, INFO )\n\n*  Purpose\n*  =======\n*\n*  ZTGSJA computes the generalized singular value decomposition (GSVD)\n*  of two complex upper triangular (or trapezoidal) matrices A and B.\n*\n*  On entry, it is assumed that matrices A and B have the following\n*  forms, which may be obtained by the preprocessing subroutine ZGGSVP\n*  from a general M-by-N matrix A and P-by-N matrix B:\n*\n*               N-K-L  K    L\n*     A =    K ( 0    A12  A13 ) if M-K-L >= 0;\n*            L ( 0     0   A23 )\n*        M-K-L ( 0     0    0  )\n*\n*             N-K-L  K    L\n*     A =  K ( 0    A12  A13 ) if M-K-L < 0;\n*        M-K ( 0     0   A23 )\n*\n*             N-K-L  K    L\n*     B =  L ( 0     0   B13 )\n*        P-L ( 0     0    0  )\n*\n*  where the K-by-K matrix A12 and L-by-L matrix B13 are nonsingular\n*  upper triangular; A23 is L-by-L upper triangular if M-K-L >= 0,\n*  otherwise A23 is (M-K)-by-L upper trapezoidal.\n*\n*  On exit,\n*\n*         U'*A*Q = D1*( 0 R ),    V'*B*Q = D2*( 0 R ),\n*\n*  where U, V and Q are unitary matrices, Z' denotes the conjugate\n*  transpose of Z, R is a nonsingular upper triangular matrix, and D1\n*  and D2 are ``diagonal'' matrices, which are of the following\n*  structures:\n*\n*  If M-K-L >= 0,\n*\n*                      K  L\n*         D1 =     K ( I  0 )\n*                  L ( 0  C )\n*              M-K-L ( 0  0 )\n*\n*                     K  L\n*         D2 = L   ( 0  S )\n*              P-L ( 0  0 )\n*\n*                 N-K-L  K    L\n*    ( 0 R ) = K (  0   R11  R12 ) K\n*              L (  0    0   R22 ) L\n*\n*  where\n*\n*    C = diag( ALPHA(K+1), ... , ALPHA(K+L) ),\n*    S = diag( BETA(K+1),  ... , BETA(K+L) ),\n*    C**2 + S**2 = I.\n*\n*    R is stored in A(1:K+L,N-K-L+1:N) on exit.\n*\n*  If M-K-L < 0,\n*\n*                 K M-K K+L-M\n*      D1 =   K ( I  0    0   )\n*           M-K ( 0  C    0   )\n*\n*                   K M-K K+L-M\n*      D2 =   M-K ( 0  S    0   )\n*           K+L-M ( 0  0    I   )\n*             P-L ( 0  0    0   )\n*\n*                 N-K-L  K   M-K  K+L-M\n* ( 0 R ) =    K ( 0    R11  R12  R13  )\n*            M-K ( 0     0   R22  R23  )\n*          K+L-M ( 0     0    0   R33  )\n*\n*  where\n*  C = diag( ALPHA(K+1), ... , ALPHA(M) ),\n*  S = diag( BETA(K+1),  ... , BETA(M) ),\n*  C**2 + S**2 = I.\n*\n*  R = ( R11 R12 R13 ) is stored in A(1:M, N-K-L+1:N) and R33 is stored\n*      (  0  R22 R23 )\n*  in B(M-K+1:L,N+M-K-L+1:N) on exit.\n*\n*  The computation of the unitary transformation matrices U, V or Q\n*  is optional.  These matrices may either be formed explicitly, or they\n*  may be postmultiplied into input matrices U1, V1, or Q1.\n*\n\n*  Arguments\n*  =========\n*\n*  JOBU    (input) CHARACTER*1\n*          = 'U':  U must contain a unitary matrix U1 on entry, and\n*                  the product U1*U is returned;\n*          = 'I':  U is initialized to the unit matrix, and the\n*                  unitary matrix U is returned;\n*          = 'N':  U is not computed.\n*\n*  JOBV    (input) CHARACTER*1\n*          = 'V':  V must contain a unitary matrix V1 on entry, and\n*                  the product V1*V is returned;\n*          = 'I':  V is initialized to the unit matrix, and the\n*                  unitary matrix V is returned;\n*          = 'N':  V is not computed.\n*\n*  JOBQ    (input) CHARACTER*1\n*          = 'Q':  Q must contain a unitary matrix Q1 on entry, and\n*                  the product Q1*Q is returned;\n*          = 'I':  Q is initialized to the unit matrix, and the\n*                  unitary matrix Q is returned;\n*          = 'N':  Q is not computed.\n*\n*  M       (input) INTEGER\n*          The number of rows of the matrix A.  M >= 0.\n*\n*  P       (input) INTEGER\n*          The number of rows of the matrix B.  P >= 0.\n*\n*  N       (input) INTEGER\n*          The number of columns of the matrices A and B.  N >= 0.\n*\n*  K       (input) INTEGER\n*  L       (input) INTEGER\n*          K and L specify the subblocks in the input matrices A and B:\n*          A23 = A(K+1:MIN(K+L,M),N-L+1:N) and B13 = B(1:L,,N-L+1:N)\n*          of A and B, whose GSVD is going to be computed by ZTGSJA.\n*          See Further details.\n*\n*  A       (input/output) COMPLEX*16 array, dimension (LDA,N)\n*          On entry, the M-by-N matrix A.\n*          On exit, A(N-K+1:N,1:MIN(K+L,M) ) contains the triangular\n*          matrix R or part of R.  See Purpose for details.\n*\n*  LDA     (input) INTEGER\n*          The leading dimension of the array A. LDA >= max(1,M).\n*\n*  B       (input/output) COMPLEX*16 array, dimension (LDB,N)\n*          On entry, the P-by-N matrix B.\n*          On exit, if necessary, B(M-K+1:L,N+M-K-L+1:N) contains\n*          a part of R.  See Purpose for details.\n*\n*  LDB     (input) INTEGER\n*          The leading dimension of the array B. LDB >= max(1,P).\n*\n*  TOLA    (input) DOUBLE PRECISION\n*  TOLB    (input) DOUBLE PRECISION\n*          TOLA and TOLB are the convergence criteria for the Jacobi-\n*          Kogbetliantz iteration procedure. Generally, they are the\n*          same as used in the preprocessing step, say\n*              TOLA = MAX(M,N)*norm(A)*MAZHEPS,\n*              TOLB = MAX(P,N)*norm(B)*MAZHEPS.\n*\n*  ALPHA   (output) DOUBLE PRECISION array, dimension (N)\n*  BETA    (output) DOUBLE PRECISION array, dimension (N)\n*          On exit, ALPHA and BETA contain the generalized singular\n*          value pairs of A and B;\n*            ALPHA(1:K) = 1,\n*            BETA(1:K)  = 0,\n*          and if M-K-L >= 0,\n*            ALPHA(K+1:K+L) = diag(C),\n*            BETA(K+1:K+L)  = diag(S),\n*          or if M-K-L < 0,\n*            ALPHA(K+1:M)= C, ALPHA(M+1:K+L)= 0\n*            BETA(K+1:M) = S, BETA(M+1:K+L) = 1.\n*          Furthermore, if K+L < N,\n*            ALPHA(K+L+1:N) = 0\n*            BETA(K+L+1:N)  = 0.\n*\n*  U       (input/output) COMPLEX*16 array, dimension (LDU,M)\n*          On entry, if JOBU = 'U', U must contain a matrix U1 (usually\n*          the unitary matrix returned by ZGGSVP).\n*          On exit,\n*          if JOBU = 'I', U contains the unitary matrix U;\n*          if JOBU = 'U', U contains the product U1*U.\n*          If JOBU = 'N', U is not referenced.\n*\n*  LDU     (input) INTEGER\n*          The leading dimension of the array U. LDU >= max(1,M) if\n*          JOBU = 'U'; LDU >= 1 otherwise.\n*\n*  V       (input/output) COMPLEX*16 array, dimension (LDV,P)\n*          On entry, if JOBV = 'V', V must contain a matrix V1 (usually\n*          the unitary matrix returned by ZGGSVP).\n*          On exit,\n*          if JOBV = 'I', V contains the unitary matrix V;\n*          if JOBV = 'V', V contains the product V1*V.\n*          If JOBV = 'N', V is not referenced.\n*\n*  LDV     (input) INTEGER\n*          The leading dimension of the array V. LDV >= max(1,P) if\n*          JOBV = 'V'; LDV >= 1 otherwise.\n*\n*  Q       (input/output) COMPLEX*16 array, dimension (LDQ,N)\n*          On entry, if JOBQ = 'Q', Q must contain a matrix Q1 (usually\n*          the unitary matrix returned by ZGGSVP).\n*          On exit,\n*          if JOBQ = 'I', Q contains the unitary matrix Q;\n*          if JOBQ = 'Q', Q contains the product Q1*Q.\n*          If JOBQ = 'N', Q is not referenced.\n*\n*  LDQ     (input) INTEGER\n*          The leading dimension of the array Q. LDQ >= max(1,N) if\n*          JOBQ = 'Q'; LDQ >= 1 otherwise.\n*\n*  WORK    (workspace) COMPLEX*16 array, dimension (2*N)\n*\n*  NCYCLE  (output) INTEGER\n*          The number of cycles required for convergence.\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit\n*          < 0:  if INFO = -i, the i-th argument had an illegal value.\n*          = 1:  the procedure does not converge after MAXIT cycles.\n*\n*  Internal Parameters\n*  ===================\n*\n*  MAXIT   INTEGER\n*          MAXIT specifies the total loops that the iterative procedure\n*          may take. If after MAXIT cycles, the routine fails to\n*          converge, we return INFO = 1.\n*\n\n*  Further Details\n*  ===============\n*\n*  ZTGSJA essentially uses a variant of Kogbetliantz algorithm to reduce\n*  min(L,M-K)-by-L triangular (or trapezoidal) matrix A23 and L-by-L\n*  matrix B13 to the form:\n*\n*           U1'*A13*Q1 = C1*R1; V1'*B13*Q1 = S1*R1,\n*\n*  where U1, V1 and Q1 are unitary matrix, and Z' is the conjugate\n*  transpose of Z.  C1 and S1 are diagonal matrices satisfying\n*\n*                C1**2 + S1**2 = I,\n*\n*  and R1 is an L-by-L nonsingular upper triangular matrix.\n*\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 12)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 12)", argc);
  rb_jobu = argv[0];
  rb_jobv = argv[1];
  rb_jobq = argv[2];
  rb_k = argv[3];
  rb_l = argv[4];
  rb_a = argv[5];
  rb_b = argv[6];
  rb_tola = argv[7];
  rb_tolb = argv[8];
  rb_u = argv[9];
  rb_v = argv[10];
  rb_q = argv[11];

  jobu = StringValueCStr(rb_jobu)[0];
  jobv = StringValueCStr(rb_jobv)[0];
  jobq = StringValueCStr(rb_jobq)[0];
  k = NUM2INT(rb_k);
  l = NUM2INT(rb_l);
  tola = NUM2DBL(rb_tola);
  tolb = NUM2DBL(rb_tolb);
  if (!NA_IsNArray(rb_a))
    rb_raise(rb_eArgError, "a (6th argument) must be NArray");
  if (NA_RANK(rb_a) != 2)
    rb_raise(rb_eArgError, "rank of a (6th argument) must be %d", 2);
  lda = NA_SHAPE0(rb_a);
  n = NA_SHAPE1(rb_a);
  if (NA_TYPE(rb_a) != NA_DCOMPLEX)
    rb_a = na_change_type(rb_a, NA_DCOMPLEX);
  a = NA_PTR_TYPE(rb_a, doublecomplex*);
  if (!NA_IsNArray(rb_b))
    rb_raise(rb_eArgError, "b (7th argument) must be NArray");
  if (NA_RANK(rb_b) != 2)
    rb_raise(rb_eArgError, "rank of b (7th argument) must be %d", 2);
  ldb = NA_SHAPE0(rb_b);
  if (NA_SHAPE1(rb_b) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of b must be the same as shape 1 of a");
  if (NA_TYPE(rb_b) != NA_DCOMPLEX)
    rb_b = na_change_type(rb_b, NA_DCOMPLEX);
  b = NA_PTR_TYPE(rb_b, doublecomplex*);
  if (!NA_IsNArray(rb_u))
    rb_raise(rb_eArgError, "u (10th argument) must be NArray");
  if (NA_RANK(rb_u) != 2)
    rb_raise(rb_eArgError, "rank of u (10th argument) must be %d", 2);
  ldu = NA_SHAPE0(rb_u);
  m = NA_SHAPE1(rb_u);
  if (NA_TYPE(rb_u) != NA_DCOMPLEX)
    rb_u = na_change_type(rb_u, NA_DCOMPLEX);
  u = NA_PTR_TYPE(rb_u, doublecomplex*);
  if (!NA_IsNArray(rb_v))
    rb_raise(rb_eArgError, "v (11th argument) must be NArray");
  if (NA_RANK(rb_v) != 2)
    rb_raise(rb_eArgError, "rank of v (11th argument) must be %d", 2);
  ldv = NA_SHAPE0(rb_v);
  p = NA_SHAPE1(rb_v);
  if (NA_TYPE(rb_v) != NA_DCOMPLEX)
    rb_v = na_change_type(rb_v, NA_DCOMPLEX);
  v = NA_PTR_TYPE(rb_v, doublecomplex*);
  if (!NA_IsNArray(rb_q))
    rb_raise(rb_eArgError, "q (12th argument) must be NArray");
  if (NA_RANK(rb_q) != 2)
    rb_raise(rb_eArgError, "rank of q (12th argument) must be %d", 2);
  ldq = NA_SHAPE0(rb_q);
  if (NA_SHAPE1(rb_q) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of q must be the same as shape 1 of a");
  if (NA_TYPE(rb_q) != NA_DCOMPLEX)
    rb_q = na_change_type(rb_q, NA_DCOMPLEX);
  q = NA_PTR_TYPE(rb_q, doublecomplex*);
  {
    int shape[1];
    shape[0] = n;
    rb_alpha = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  alpha = NA_PTR_TYPE(rb_alpha, doublereal*);
  {
    int shape[1];
    shape[0] = n;
    rb_beta = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  beta = NA_PTR_TYPE(rb_beta, doublereal*);
  {
    int shape[2];
    shape[0] = lda;
    shape[1] = n;
    rb_a_out__ = na_make_object(NA_DCOMPLEX, 2, shape, cNArray);
  }
  a_out__ = NA_PTR_TYPE(rb_a_out__, doublecomplex*);
  MEMCPY(a_out__, a, doublecomplex, NA_TOTAL(rb_a));
  rb_a = rb_a_out__;
  a = a_out__;
  {
    int shape[2];
    shape[0] = ldb;
    shape[1] = n;
    rb_b_out__ = na_make_object(NA_DCOMPLEX, 2, shape, cNArray);
  }
  b_out__ = NA_PTR_TYPE(rb_b_out__, doublecomplex*);
  MEMCPY(b_out__, b, doublecomplex, NA_TOTAL(rb_b));
  rb_b = rb_b_out__;
  b = b_out__;
  {
    int shape[2];
    shape[0] = ldu;
    shape[1] = m;
    rb_u_out__ = na_make_object(NA_DCOMPLEX, 2, shape, cNArray);
  }
  u_out__ = NA_PTR_TYPE(rb_u_out__, doublecomplex*);
  MEMCPY(u_out__, u, doublecomplex, NA_TOTAL(rb_u));
  rb_u = rb_u_out__;
  u = u_out__;
  {
    int shape[2];
    shape[0] = ldv;
    shape[1] = p;
    rb_v_out__ = na_make_object(NA_DCOMPLEX, 2, shape, cNArray);
  }
  v_out__ = NA_PTR_TYPE(rb_v_out__, doublecomplex*);
  MEMCPY(v_out__, v, doublecomplex, NA_TOTAL(rb_v));
  rb_v = rb_v_out__;
  v = v_out__;
  {
    int shape[2];
    shape[0] = ldq;
    shape[1] = n;
    rb_q_out__ = na_make_object(NA_DCOMPLEX, 2, shape, cNArray);
  }
  q_out__ = NA_PTR_TYPE(rb_q_out__, doublecomplex*);
  MEMCPY(q_out__, q, doublecomplex, NA_TOTAL(rb_q));
  rb_q = rb_q_out__;
  q = q_out__;
  work = ALLOC_N(doublecomplex, (2*n));

  ztgsja_(&jobu, &jobv, &jobq, &m, &p, &n, &k, &l, a, &lda, b, &ldb, &tola, &tolb, alpha, beta, u, &ldu, v, &ldv, q, &ldq, work, &ncycle, &info);

  free(work);
  rb_ncycle = INT2NUM(ncycle);
  rb_info = INT2NUM(info);
  return rb_ary_new3(9, rb_alpha, rb_beta, rb_ncycle, rb_info, rb_a, rb_b, rb_u, rb_v, rb_q);
}

void
init_lapack_ztgsja(VALUE mLapack){
  rb_define_module_function(mLapack, "ztgsja", rb_ztgsja, -1);
}

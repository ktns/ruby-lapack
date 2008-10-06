#include "rb_lapack.h"

static VALUE
rb_dhgeqz(int argc, VALUE *argv, VALUE self){
  VALUE rb_job;
  char job; 
  VALUE rb_compq;
  char compq; 
  VALUE rb_compz;
  char compz; 
  VALUE rb_ilo;
  integer ilo; 
  VALUE rb_ihi;
  integer ihi; 
  VALUE rb_h;
  doublereal *h; 
  VALUE rb_t;
  doublereal *t; 
  VALUE rb_q;
  doublereal *q; 
  VALUE rb_z;
  doublereal *z; 
  VALUE rb_lwork;
  integer lwork; 
  VALUE rb_alphar;
  doublereal *alphar; 
  VALUE rb_alphai;
  doublereal *alphai; 
  VALUE rb_beta;
  doublereal *beta; 
  VALUE rb_work;
  doublereal *work; 
  VALUE rb_info;
  integer info; 
  VALUE rb_h_out__;
  doublereal *h_out__;
  VALUE rb_t_out__;
  doublereal *t_out__;
  VALUE rb_q_out__;
  doublereal *q_out__;
  VALUE rb_z_out__;
  doublereal *z_out__;

  integer ldh;
  integer n;
  integer ldt;
  integer ldq;
  integer ldz;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  alphar, alphai, beta, work, info, h, t, q, z = NumRu::Lapack.dhgeqz( job, compq, compz, ilo, ihi, h, t, q, z, lwork)\n    or\n  NumRu::Lapack.dhgeqz  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE DHGEQZ( JOB, COMPQ, COMPZ, N, ILO, IHI, H, LDH, T, LDT, ALPHAR, ALPHAI, BETA, Q, LDQ, Z, LDZ, WORK, LWORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  DHGEQZ computes the eigenvalues of a real matrix pair (H,T),\n*  where H is an upper Hessenberg matrix and T is upper triangular,\n*  using the double-shift QZ method.\n*  Matrix pairs of this type are produced by the reduction to\n*  generalized upper Hessenberg form of a real matrix pair (A,B):\n*\n*     A = Q1*H*Z1**T,  B = Q1*T*Z1**T,\n*\n*  as computed by DGGHRD.\n*\n*  If JOB='S', then the Hessenberg-triangular pair (H,T) is\n*  also reduced to generalized Schur form,\n*  \n*     H = Q*S*Z**T,  T = Q*P*Z**T,\n*  \n*  where Q and Z are orthogonal matrices, P is an upper triangular\n*  matrix, and S is a quasi-triangular matrix with 1-by-1 and 2-by-2\n*  diagonal blocks.\n*\n*  The 1-by-1 blocks correspond to real eigenvalues of the matrix pair\n*  (H,T) and the 2-by-2 blocks correspond to complex conjugate pairs of\n*  eigenvalues.\n*\n*  Additionally, the 2-by-2 upper triangular diagonal blocks of P\n*  corresponding to 2-by-2 blocks of S are reduced to positive diagonal\n*  form, i.e., if S(j+1,j) is non-zero, then P(j+1,j) = P(j,j+1) = 0,\n*  P(j,j) > 0, and P(j+1,j+1) > 0.\n*\n*  Optionally, the orthogonal matrix Q from the generalized Schur\n*  factorization may be postmultiplied into an input matrix Q1, and the\n*  orthogonal matrix Z may be postmultiplied into an input matrix Z1.\n*  If Q1 and Z1 are the orthogonal matrices from DGGHRD that reduced\n*  the matrix pair (A,B) to generalized upper Hessenberg form, then the\n*  output matrices Q1*Q and Z1*Z are the orthogonal factors from the\n*  generalized Schur factorization of (A,B):\n*\n*     A = (Q1*Q)*S*(Z1*Z)**T,  B = (Q1*Q)*P*(Z1*Z)**T.\n*  \n*  To avoid overflow, eigenvalues of the matrix pair (H,T) (equivalently,\n*  of (A,B)) are computed as a pair of values (alpha,beta), where alpha is\n*  complex and beta real.\n*  If beta is nonzero, lambda = alpha / beta is an eigenvalue of the\n*  generalized nonsymmetric eigenvalue problem (GNEP)\n*     A*x = lambda*B*x\n*  and if alpha is nonzero, mu = beta / alpha is an eigenvalue of the\n*  alternate form of the GNEP\n*     mu*A*y = B*y.\n*  Real eigenvalues can be read directly from the generalized Schur\n*  form: \n*    alpha = S(i,i), beta = P(i,i).\n*\n*  Ref: C.B. Moler & G.W. Stewart, \"An Algorithm for Generalized Matrix\n*       Eigenvalue Problems\", SIAM J. Numer. Anal., 10(1973),\n*       pp. 241--256.\n*\n\n*  Arguments\n*  =========\n*\n*  JOB     (input) CHARACTER*1\n*          = 'E': Compute eigenvalues only;\n*          = 'S': Compute eigenvalues and the Schur form. \n*\n*  COMPQ   (input) CHARACTER*1\n*          = 'N': Left Schur vectors (Q) are not computed;\n*          = 'I': Q is initialized to the unit matrix and the matrix Q\n*                 of left Schur vectors of (H,T) is returned;\n*          = 'V': Q must contain an orthogonal matrix Q1 on entry and\n*                 the product Q1*Q is returned.\n*\n*  COMPZ   (input) CHARACTER*1\n*          = 'N': Right Schur vectors (Z) are not computed;\n*          = 'I': Z is initialized to the unit matrix and the matrix Z\n*                 of right Schur vectors of (H,T) is returned;\n*          = 'V': Z must contain an orthogonal matrix Z1 on entry and\n*                 the product Z1*Z is returned.\n*\n*  N       (input) INTEGER\n*          The order of the matrices H, T, Q, and Z.  N >= 0.\n*\n*  ILO     (input) INTEGER\n*  IHI     (input) INTEGER\n*          ILO and IHI mark the rows and columns of H which are in\n*          Hessenberg form.  It is assumed that A is already upper\n*          triangular in rows and columns 1:ILO-1 and IHI+1:N.\n*          If N > 0, 1 <= ILO <= IHI <= N; if N = 0, ILO=1 and IHI=0.\n*\n*  H       (input/output) DOUBLE PRECISION array, dimension (LDH, N)\n*          On entry, the N-by-N upper Hessenberg matrix H.\n*          On exit, if JOB = 'S', H contains the upper quasi-triangular\n*          matrix S from the generalized Schur factorization;\n*          2-by-2 diagonal blocks (corresponding to complex conjugate\n*          pairs of eigenvalues) are returned in standard form, with\n*          H(i,i) = H(i+1,i+1) and H(i+1,i)*H(i,i+1) < 0.\n*          If JOB = 'E', the diagonal blocks of H match those of S, but\n*          the rest of H is unspecified.\n*\n*  LDH     (input) INTEGER\n*          The leading dimension of the array H.  LDH >= max( 1, N ).\n*\n*  T       (input/output) DOUBLE PRECISION array, dimension (LDT, N)\n*          On entry, the N-by-N upper triangular matrix T.\n*          On exit, if JOB = 'S', T contains the upper triangular\n*          matrix P from the generalized Schur factorization;\n*          2-by-2 diagonal blocks of P corresponding to 2-by-2 blocks of S\n*          are reduced to positive diagonal form, i.e., if H(j+1,j) is\n*          non-zero, then T(j+1,j) = T(j,j+1) = 0, T(j,j) > 0, and\n*          T(j+1,j+1) > 0.\n*          If JOB = 'E', the diagonal blocks of T match those of P, but\n*          the rest of T is unspecified.\n*\n*  LDT     (input) INTEGER\n*          The leading dimension of the array T.  LDT >= max( 1, N ).\n*\n*  ALPHAR  (output) DOUBLE PRECISION array, dimension (N)\n*          The real parts of each scalar alpha defining an eigenvalue\n*          of GNEP.\n*\n*  ALPHAI  (output) DOUBLE PRECISION array, dimension (N)\n*          The imaginary parts of each scalar alpha defining an\n*          eigenvalue of GNEP.\n*          If ALPHAI(j) is zero, then the j-th eigenvalue is real; if\n*          positive, then the j-th and (j+1)-st eigenvalues are a\n*          complex conjugate pair, with ALPHAI(j+1) = -ALPHAI(j).\n*\n*  BETA    (output) DOUBLE PRECISION array, dimension (N)\n*          The scalars beta that define the eigenvalues of GNEP.\n*          Together, the quantities alpha = (ALPHAR(j),ALPHAI(j)) and\n*          beta = BETA(j) represent the j-th eigenvalue of the matrix\n*          pair (A,B), in one of the forms lambda = alpha/beta or\n*          mu = beta/alpha.  Since either lambda or mu may overflow,\n*          they should not, in general, be computed.\n*\n*  Q       (input/output) DOUBLE PRECISION array, dimension (LDQ, N)\n*          On entry, if COMPZ = 'V', the orthogonal matrix Q1 used in\n*          the reduction of (A,B) to generalized Hessenberg form.\n*          On exit, if COMPZ = 'I', the orthogonal matrix of left Schur\n*          vectors of (H,T), and if COMPZ = 'V', the orthogonal matrix\n*          of left Schur vectors of (A,B).\n*          Not referenced if COMPZ = 'N'.\n*\n*  LDQ     (input) INTEGER\n*          The leading dimension of the array Q.  LDQ >= 1.\n*          If COMPQ='V' or 'I', then LDQ >= N.\n*\n*  Z       (input/output) DOUBLE PRECISION array, dimension (LDZ, N)\n*          On entry, if COMPZ = 'V', the orthogonal matrix Z1 used in\n*          the reduction of (A,B) to generalized Hessenberg form.\n*          On exit, if COMPZ = 'I', the orthogonal matrix of\n*          right Schur vectors of (H,T), and if COMPZ = 'V', the\n*          orthogonal matrix of right Schur vectors of (A,B).\n*          Not referenced if COMPZ = 'N'.\n*\n*  LDZ     (input) INTEGER\n*          The leading dimension of the array Z.  LDZ >= 1.\n*          If COMPZ='V' or 'I', then LDZ >= N.\n*\n*  WORK    (workspace/output) DOUBLE PRECISION array, dimension (MAX(1,LWORK))\n*          On exit, if INFO >= 0, WORK(1) returns the optimal LWORK.\n*\n*  LWORK   (input) INTEGER\n*          The dimension of the array WORK.  LWORK >= max(1,N).\n*\n*          If LWORK = -1, then a workspace query is assumed; the routine\n*          only calculates the optimal size of the WORK array, returns\n*          this value as the first entry of the WORK array, and no error\n*          message related to LWORK is issued by XERBLA.\n*\n*  INFO    (output) INTEGER\n*          = 0: successful exit\n*          < 0: if INFO = -i, the i-th argument had an illegal value\n*          = 1,...,N: the QZ iteration did not converge.  (H,T) is not\n*                     in Schur form, but ALPHAR(i), ALPHAI(i), and\n*                     BETA(i), i=INFO+1,...,N should be correct.\n*          = N+1,...,2*N: the shift calculation failed.  (H,T) is not\n*                     in Schur form, but ALPHAR(i), ALPHAI(i), and\n*                     BETA(i), i=INFO-N+1,...,N should be correct.\n*\n\n*  Further Details\n*  ===============\n*\n*  Iteration counters:\n*\n*  JITER  -- counts iterations.\n*  IITER  -- counts iterations run since ILAST was last\n*            changed.  This is therefore reset only when a 1-by-1 or\n*            2-by-2 block deflates off the bottom.\n*\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 10)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 10)", argc);
  rb_job = argv[0];
  rb_compq = argv[1];
  rb_compz = argv[2];
  rb_ilo = argv[3];
  rb_ihi = argv[4];
  rb_h = argv[5];
  rb_t = argv[6];
  rb_q = argv[7];
  rb_z = argv[8];
  rb_lwork = argv[9];

  job = StringValueCStr(rb_job)[0];
  compq = StringValueCStr(rb_compq)[0];
  compz = StringValueCStr(rb_compz)[0];
  ilo = NUM2INT(rb_ilo);
  ihi = NUM2INT(rb_ihi);
  lwork = NUM2INT(rb_lwork);
  if (!NA_IsNArray(rb_h))
    rb_raise(rb_eArgError, "h (6th argument) must be NArray");
  if (NA_RANK(rb_h) != 2)
    rb_raise(rb_eArgError, "rank of h (6th argument) must be %d", 2);
  ldh = NA_SHAPE0(rb_h);
  n = NA_SHAPE1(rb_h);
  if (NA_TYPE(rb_h) != NA_DFLOAT)
    rb_h = na_change_type(rb_h, NA_DFLOAT);
  h = NA_PTR_TYPE(rb_h, doublereal*);
  if (!NA_IsNArray(rb_t))
    rb_raise(rb_eArgError, "t (7th argument) must be NArray");
  if (NA_RANK(rb_t) != 2)
    rb_raise(rb_eArgError, "rank of t (7th argument) must be %d", 2);
  ldt = NA_SHAPE0(rb_t);
  if (NA_SHAPE1(rb_t) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of t must be the same as shape 1 of h");
  if (NA_TYPE(rb_t) != NA_DFLOAT)
    rb_t = na_change_type(rb_t, NA_DFLOAT);
  t = NA_PTR_TYPE(rb_t, doublereal*);
  if (!NA_IsNArray(rb_q))
    rb_raise(rb_eArgError, "q (8th argument) must be NArray");
  if (NA_RANK(rb_q) != 2)
    rb_raise(rb_eArgError, "rank of q (8th argument) must be %d", 2);
  ldq = NA_SHAPE0(rb_q);
  if (NA_SHAPE1(rb_q) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of q must be the same as shape 1 of h");
  if (NA_TYPE(rb_q) != NA_DFLOAT)
    rb_q = na_change_type(rb_q, NA_DFLOAT);
  q = NA_PTR_TYPE(rb_q, doublereal*);
  if (!NA_IsNArray(rb_z))
    rb_raise(rb_eArgError, "z (9th argument) must be NArray");
  if (NA_RANK(rb_z) != 2)
    rb_raise(rb_eArgError, "rank of z (9th argument) must be %d", 2);
  ldz = NA_SHAPE0(rb_z);
  if (NA_SHAPE1(rb_z) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of z must be the same as shape 1 of h");
  if (NA_TYPE(rb_z) != NA_DFLOAT)
    rb_z = na_change_type(rb_z, NA_DFLOAT);
  z = NA_PTR_TYPE(rb_z, doublereal*);
  {
    int shape[1];
    shape[0] = n;
    rb_alphar = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  alphar = NA_PTR_TYPE(rb_alphar, doublereal*);
  {
    int shape[1];
    shape[0] = n;
    rb_alphai = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  alphai = NA_PTR_TYPE(rb_alphai, doublereal*);
  {
    int shape[1];
    shape[0] = n;
    rb_beta = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  beta = NA_PTR_TYPE(rb_beta, doublereal*);
  {
    int shape[1];
    shape[0] = MAX(1,lwork);
    rb_work = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  work = NA_PTR_TYPE(rb_work, doublereal*);
  {
    int shape[2];
    shape[0] = ldh;
    shape[1] = n;
    rb_h_out__ = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  h_out__ = NA_PTR_TYPE(rb_h_out__, doublereal*);
  MEMCPY(h_out__, h, doublereal, NA_TOTAL(rb_h));
  rb_h = rb_h_out__;
  h = h_out__;
  {
    int shape[2];
    shape[0] = ldt;
    shape[1] = n;
    rb_t_out__ = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  t_out__ = NA_PTR_TYPE(rb_t_out__, doublereal*);
  MEMCPY(t_out__, t, doublereal, NA_TOTAL(rb_t));
  rb_t = rb_t_out__;
  t = t_out__;
  {
    int shape[2];
    shape[0] = ldq;
    shape[1] = n;
    rb_q_out__ = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  q_out__ = NA_PTR_TYPE(rb_q_out__, doublereal*);
  MEMCPY(q_out__, q, doublereal, NA_TOTAL(rb_q));
  rb_q = rb_q_out__;
  q = q_out__;
  {
    int shape[2];
    shape[0] = ldz;
    shape[1] = n;
    rb_z_out__ = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  z_out__ = NA_PTR_TYPE(rb_z_out__, doublereal*);
  MEMCPY(z_out__, z, doublereal, NA_TOTAL(rb_z));
  rb_z = rb_z_out__;
  z = z_out__;

  dhgeqz_(&job, &compq, &compz, &n, &ilo, &ihi, h, &ldh, t, &ldt, alphar, alphai, beta, q, &ldq, z, &ldz, work, &lwork, &info);

  rb_info = INT2NUM(info);
  return rb_ary_new3(9, rb_alphar, rb_alphai, rb_beta, rb_work, rb_info, rb_h, rb_t, rb_q, rb_z);
}

void
init_lapack_dhgeqz(VALUE mLapack){
  rb_define_module_function(mLapack, "dhgeqz", rb_dhgeqz, -1);
}

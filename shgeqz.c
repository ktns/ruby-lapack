#include "rb_lapack.h"

extern VOID shgeqz_(char *job, char *compq, char *compz, integer *n, integer *ilo, integer *ihi, real *h, integer *ldh, real *t, integer *ldt, real *alphar, real *alphai, real *beta, real *q, integer *ldq, real *z, integer *ldz, real *work, integer *lwork, integer *info);

static VALUE sHelp, sUsage;

static VALUE
rblapack_shgeqz(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_job;
  char job; 
  VALUE rblapack_compq;
  char compq; 
  VALUE rblapack_compz;
  char compz; 
  VALUE rblapack_ilo;
  integer ilo; 
  VALUE rblapack_ihi;
  integer ihi; 
  VALUE rblapack_h;
  real *h; 
  VALUE rblapack_t;
  real *t; 
  VALUE rblapack_q;
  real *q; 
  VALUE rblapack_z;
  real *z; 
  VALUE rblapack_lwork;
  integer lwork; 
  VALUE rblapack_alphar;
  real *alphar; 
  VALUE rblapack_alphai;
  real *alphai; 
  VALUE rblapack_beta;
  real *beta; 
  VALUE rblapack_work;
  real *work; 
  VALUE rblapack_info;
  integer info; 
  VALUE rblapack_h_out__;
  real *h_out__;
  VALUE rblapack_t_out__;
  real *t_out__;
  VALUE rblapack_q_out__;
  real *q_out__;
  VALUE rblapack_z_out__;
  real *z_out__;

  integer ldh;
  integer n;
  integer ldt;
  integer ldq;
  integer ldz;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  alphar, alphai, beta, work, info, h, t, q, z = NumRu::Lapack.shgeqz( job, compq, compz, ilo, ihi, h, t, q, z, lwork, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE SHGEQZ( JOB, COMPQ, COMPZ, N, ILO, IHI, H, LDH, T, LDT, ALPHAR, ALPHAI, BETA, Q, LDQ, Z, LDZ, WORK, LWORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  SHGEQZ computes the eigenvalues of a real matrix pair (H,T),\n*  where H is an upper Hessenberg matrix and T is upper triangular,\n*  using the double-shift QZ method.\n*  Matrix pairs of this type are produced by the reduction to\n*  generalized upper Hessenberg form of a real matrix pair (A,B):\n*\n*     A = Q1*H*Z1**T,  B = Q1*T*Z1**T,\n*\n*  as computed by SGGHRD.\n*\n*  If JOB='S', then the Hessenberg-triangular pair (H,T) is\n*  also reduced to generalized Schur form,\n*  \n*     H = Q*S*Z**T,  T = Q*P*Z**T,\n*  \n*  where Q and Z are orthogonal matrices, P is an upper triangular\n*  matrix, and S is a quasi-triangular matrix with 1-by-1 and 2-by-2\n*  diagonal blocks.\n*\n*  The 1-by-1 blocks correspond to real eigenvalues of the matrix pair\n*  (H,T) and the 2-by-2 blocks correspond to complex conjugate pairs of\n*  eigenvalues.\n*\n*  Additionally, the 2-by-2 upper triangular diagonal blocks of P\n*  corresponding to 2-by-2 blocks of S are reduced to positive diagonal\n*  form, i.e., if S(j+1,j) is non-zero, then P(j+1,j) = P(j,j+1) = 0,\n*  P(j,j) > 0, and P(j+1,j+1) > 0.\n*\n*  Optionally, the orthogonal matrix Q from the generalized Schur\n*  factorization may be postmultiplied into an input matrix Q1, and the\n*  orthogonal matrix Z may be postmultiplied into an input matrix Z1.\n*  If Q1 and Z1 are the orthogonal matrices from SGGHRD that reduced\n*  the matrix pair (A,B) to generalized upper Hessenberg form, then the\n*  output matrices Q1*Q and Z1*Z are the orthogonal factors from the\n*  generalized Schur factorization of (A,B):\n*\n*     A = (Q1*Q)*S*(Z1*Z)**T,  B = (Q1*Q)*P*(Z1*Z)**T.\n*  \n*  To avoid overflow, eigenvalues of the matrix pair (H,T) (equivalently,\n*  of (A,B)) are computed as a pair of values (alpha,beta), where alpha is\n*  complex and beta real.\n*  If beta is nonzero, lambda = alpha / beta is an eigenvalue of the\n*  generalized nonsymmetric eigenvalue problem (GNEP)\n*     A*x = lambda*B*x\n*  and if alpha is nonzero, mu = beta / alpha is an eigenvalue of the\n*  alternate form of the GNEP\n*     mu*A*y = B*y.\n*  Real eigenvalues can be read directly from the generalized Schur\n*  form: \n*    alpha = S(i,i), beta = P(i,i).\n*\n*  Ref: C.B. Moler & G.W. Stewart, \"An Algorithm for Generalized Matrix\n*       Eigenvalue Problems\", SIAM J. Numer. Anal., 10(1973),\n*       pp. 241--256.\n*\n\n*  Arguments\n*  =========\n*\n*  JOB     (input) CHARACTER*1\n*          = 'E': Compute eigenvalues only;\n*          = 'S': Compute eigenvalues and the Schur form. \n*\n*  COMPQ   (input) CHARACTER*1\n*          = 'N': Left Schur vectors (Q) are not computed;\n*          = 'I': Q is initialized to the unit matrix and the matrix Q\n*                 of left Schur vectors of (H,T) is returned;\n*          = 'V': Q must contain an orthogonal matrix Q1 on entry and\n*                 the product Q1*Q is returned.\n*\n*  COMPZ   (input) CHARACTER*1\n*          = 'N': Right Schur vectors (Z) are not computed;\n*          = 'I': Z is initialized to the unit matrix and the matrix Z\n*                 of right Schur vectors of (H,T) is returned;\n*          = 'V': Z must contain an orthogonal matrix Z1 on entry and\n*                 the product Z1*Z is returned.\n*\n*  N       (input) INTEGER\n*          The order of the matrices H, T, Q, and Z.  N >= 0.\n*\n*  ILO     (input) INTEGER\n*  IHI     (input) INTEGER\n*          ILO and IHI mark the rows and columns of H which are in\n*          Hessenberg form.  It is assumed that A is already upper\n*          triangular in rows and columns 1:ILO-1 and IHI+1:N.\n*          If N > 0, 1 <= ILO <= IHI <= N; if N = 0, ILO=1 and IHI=0.\n*\n*  H       (input/output) REAL array, dimension (LDH, N)\n*          On entry, the N-by-N upper Hessenberg matrix H.\n*          On exit, if JOB = 'S', H contains the upper quasi-triangular\n*          matrix S from the generalized Schur factorization;\n*          2-by-2 diagonal blocks (corresponding to complex conjugate\n*          pairs of eigenvalues) are returned in standard form, with\n*          H(i,i) = H(i+1,i+1) and H(i+1,i)*H(i,i+1) < 0.\n*          If JOB = 'E', the diagonal blocks of H match those of S, but\n*          the rest of H is unspecified.\n*\n*  LDH     (input) INTEGER\n*          The leading dimension of the array H.  LDH >= max( 1, N ).\n*\n*  T       (input/output) REAL array, dimension (LDT, N)\n*          On entry, the N-by-N upper triangular matrix T.\n*          On exit, if JOB = 'S', T contains the upper triangular\n*          matrix P from the generalized Schur factorization;\n*          2-by-2 diagonal blocks of P corresponding to 2-by-2 blocks of S\n*          are reduced to positive diagonal form, i.e., if H(j+1,j) is\n*          non-zero, then T(j+1,j) = T(j,j+1) = 0, T(j,j) > 0, and\n*          T(j+1,j+1) > 0.\n*          If JOB = 'E', the diagonal blocks of T match those of P, but\n*          the rest of T is unspecified.\n*\n*  LDT     (input) INTEGER\n*          The leading dimension of the array T.  LDT >= max( 1, N ).\n*\n*  ALPHAR  (output) REAL array, dimension (N)\n*          The real parts of each scalar alpha defining an eigenvalue\n*          of GNEP.\n*\n*  ALPHAI  (output) REAL array, dimension (N)\n*          The imaginary parts of each scalar alpha defining an\n*          eigenvalue of GNEP.\n*          If ALPHAI(j) is zero, then the j-th eigenvalue is real; if\n*          positive, then the j-th and (j+1)-st eigenvalues are a\n*          complex conjugate pair, with ALPHAI(j+1) = -ALPHAI(j).\n*\n*  BETA    (output) REAL array, dimension (N)\n*          The scalars beta that define the eigenvalues of GNEP.\n*          Together, the quantities alpha = (ALPHAR(j),ALPHAI(j)) and\n*          beta = BETA(j) represent the j-th eigenvalue of the matrix\n*          pair (A,B), in one of the forms lambda = alpha/beta or\n*          mu = beta/alpha.  Since either lambda or mu may overflow,\n*          they should not, in general, be computed.\n*\n*  Q       (input/output) REAL array, dimension (LDQ, N)\n*          On entry, if COMPZ = 'V', the orthogonal matrix Q1 used in\n*          the reduction of (A,B) to generalized Hessenberg form.\n*          On exit, if COMPZ = 'I', the orthogonal matrix of left Schur\n*          vectors of (H,T), and if COMPZ = 'V', the orthogonal matrix\n*          of left Schur vectors of (A,B).\n*          Not referenced if COMPZ = 'N'.\n*\n*  LDQ     (input) INTEGER\n*          The leading dimension of the array Q.  LDQ >= 1.\n*          If COMPQ='V' or 'I', then LDQ >= N.\n*\n*  Z       (input/output) REAL array, dimension (LDZ, N)\n*          On entry, if COMPZ = 'V', the orthogonal matrix Z1 used in\n*          the reduction of (A,B) to generalized Hessenberg form.\n*          On exit, if COMPZ = 'I', the orthogonal matrix of\n*          right Schur vectors of (H,T), and if COMPZ = 'V', the\n*          orthogonal matrix of right Schur vectors of (A,B).\n*          Not referenced if COMPZ = 'N'.\n*\n*  LDZ     (input) INTEGER\n*          The leading dimension of the array Z.  LDZ >= 1.\n*          If COMPZ='V' or 'I', then LDZ >= N.\n*\n*  WORK    (workspace/output) REAL array, dimension (MAX(1,LWORK))\n*          On exit, if INFO >= 0, WORK(1) returns the optimal LWORK.\n*\n*  LWORK   (input) INTEGER\n*          The dimension of the array WORK.  LWORK >= max(1,N).\n*\n*          If LWORK = -1, then a workspace query is assumed; the routine\n*          only calculates the optimal size of the WORK array, returns\n*          this value as the first entry of the WORK array, and no error\n*          message related to LWORK is issued by XERBLA.\n*\n*  INFO    (output) INTEGER\n*          = 0: successful exit\n*          < 0: if INFO = -i, the i-th argument had an illegal value\n*          = 1,...,N: the QZ iteration did not converge.  (H,T) is not\n*                     in Schur form, but ALPHAR(i), ALPHAI(i), and\n*                     BETA(i), i=INFO+1,...,N should be correct.\n*          = N+1,...,2*N: the shift calculation failed.  (H,T) is not\n*                     in Schur form, but ALPHAR(i), ALPHAI(i), and\n*                     BETA(i), i=INFO-N+1,...,N should be correct.\n*\n\n*  Further Details\n*  ===============\n*\n*  Iteration counters:\n*\n*  JITER  -- counts iterations.\n*  IITER  -- counts iterations run since ILAST was last\n*            changed.  This is therefore reset only when a 1-by-1 or\n*            2-by-2 block deflates off the bottom.\n*\n*  =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  alphar, alphai, beta, work, info, h, t, q, z = NumRu::Lapack.shgeqz( job, compq, compz, ilo, ihi, h, t, q, z, lwork, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 10)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 10)", argc);
  rblapack_job = argv[0];
  rblapack_compq = argv[1];
  rblapack_compz = argv[2];
  rblapack_ilo = argv[3];
  rblapack_ihi = argv[4];
  rblapack_h = argv[5];
  rblapack_t = argv[6];
  rblapack_q = argv[7];
  rblapack_z = argv[8];
  rblapack_lwork = argv[9];
  if (rb_options != Qnil) {
  }

  ilo = NUM2INT(rblapack_ilo);
  compz = StringValueCStr(rblapack_compz)[0];
  if (!NA_IsNArray(rblapack_z))
    rb_raise(rb_eArgError, "z (9th argument) must be NArray");
  if (NA_RANK(rblapack_z) != 2)
    rb_raise(rb_eArgError, "rank of z (9th argument) must be %d", 2);
  n = NA_SHAPE1(rblapack_z);
  ldz = NA_SHAPE0(rblapack_z);
  if (NA_TYPE(rblapack_z) != NA_SFLOAT)
    rblapack_z = na_change_type(rblapack_z, NA_SFLOAT);
  z = NA_PTR_TYPE(rblapack_z, real*);
  compq = StringValueCStr(rblapack_compq)[0];
  lwork = NUM2INT(rblapack_lwork);
  if (!NA_IsNArray(rblapack_q))
    rb_raise(rb_eArgError, "q (8th argument) must be NArray");
  if (NA_RANK(rblapack_q) != 2)
    rb_raise(rb_eArgError, "rank of q (8th argument) must be %d", 2);
  if (NA_SHAPE1(rblapack_q) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of q must be the same as shape 1 of z");
  ldq = NA_SHAPE0(rblapack_q);
  if (NA_TYPE(rblapack_q) != NA_SFLOAT)
    rblapack_q = na_change_type(rblapack_q, NA_SFLOAT);
  q = NA_PTR_TYPE(rblapack_q, real*);
  job = StringValueCStr(rblapack_job)[0];
  if (!NA_IsNArray(rblapack_h))
    rb_raise(rb_eArgError, "h (6th argument) must be NArray");
  if (NA_RANK(rblapack_h) != 2)
    rb_raise(rb_eArgError, "rank of h (6th argument) must be %d", 2);
  if (NA_SHAPE1(rblapack_h) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of h must be the same as shape 1 of z");
  ldh = NA_SHAPE0(rblapack_h);
  if (NA_TYPE(rblapack_h) != NA_SFLOAT)
    rblapack_h = na_change_type(rblapack_h, NA_SFLOAT);
  h = NA_PTR_TYPE(rblapack_h, real*);
  if (!NA_IsNArray(rblapack_t))
    rb_raise(rb_eArgError, "t (7th argument) must be NArray");
  if (NA_RANK(rblapack_t) != 2)
    rb_raise(rb_eArgError, "rank of t (7th argument) must be %d", 2);
  if (NA_SHAPE1(rblapack_t) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of t must be the same as shape 1 of z");
  ldt = NA_SHAPE0(rblapack_t);
  if (NA_TYPE(rblapack_t) != NA_SFLOAT)
    rblapack_t = na_change_type(rblapack_t, NA_SFLOAT);
  t = NA_PTR_TYPE(rblapack_t, real*);
  ihi = NUM2INT(rblapack_ihi);
  {
    int shape[1];
    shape[0] = n;
    rblapack_alphar = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  alphar = NA_PTR_TYPE(rblapack_alphar, real*);
  {
    int shape[1];
    shape[0] = n;
    rblapack_alphai = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  alphai = NA_PTR_TYPE(rblapack_alphai, real*);
  {
    int shape[1];
    shape[0] = n;
    rblapack_beta = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  beta = NA_PTR_TYPE(rblapack_beta, real*);
  {
    int shape[1];
    shape[0] = MAX(1,lwork);
    rblapack_work = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  work = NA_PTR_TYPE(rblapack_work, real*);
  {
    int shape[2];
    shape[0] = ldh;
    shape[1] = n;
    rblapack_h_out__ = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  h_out__ = NA_PTR_TYPE(rblapack_h_out__, real*);
  MEMCPY(h_out__, h, real, NA_TOTAL(rblapack_h));
  rblapack_h = rblapack_h_out__;
  h = h_out__;
  {
    int shape[2];
    shape[0] = ldt;
    shape[1] = n;
    rblapack_t_out__ = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  t_out__ = NA_PTR_TYPE(rblapack_t_out__, real*);
  MEMCPY(t_out__, t, real, NA_TOTAL(rblapack_t));
  rblapack_t = rblapack_t_out__;
  t = t_out__;
  {
    int shape[2];
    shape[0] = ldq;
    shape[1] = n;
    rblapack_q_out__ = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  q_out__ = NA_PTR_TYPE(rblapack_q_out__, real*);
  MEMCPY(q_out__, q, real, NA_TOTAL(rblapack_q));
  rblapack_q = rblapack_q_out__;
  q = q_out__;
  {
    int shape[2];
    shape[0] = ldz;
    shape[1] = n;
    rblapack_z_out__ = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  z_out__ = NA_PTR_TYPE(rblapack_z_out__, real*);
  MEMCPY(z_out__, z, real, NA_TOTAL(rblapack_z));
  rblapack_z = rblapack_z_out__;
  z = z_out__;

  shgeqz_(&job, &compq, &compz, &n, &ilo, &ihi, h, &ldh, t, &ldt, alphar, alphai, beta, q, &ldq, z, &ldz, work, &lwork, &info);

  rblapack_info = INT2NUM(info);
  return rb_ary_new3(9, rblapack_alphar, rblapack_alphai, rblapack_beta, rblapack_work, rblapack_info, rblapack_h, rblapack_t, rblapack_q, rblapack_z);
}

void
init_lapack_shgeqz(VALUE mLapack){
  rb_define_module_function(mLapack, "shgeqz", rblapack_shgeqz, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}

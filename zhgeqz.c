#include "rb_lapack.h"

extern VOID zhgeqz_(char *job, char *compq, char *compz, integer *n, integer *ilo, integer *ihi, doublecomplex *h, integer *ldh, doublecomplex *t, integer *ldt, doublecomplex *alpha, doublecomplex *beta, doublecomplex *q, integer *ldq, doublecomplex *z, integer *ldz, doublecomplex *work, integer *lwork, doublereal *rwork, integer *info);

static VALUE sHelp, sUsage;

static VALUE
rblapack_zhgeqz(int argc, VALUE *argv, VALUE self){
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
  doublecomplex *h; 
  VALUE rblapack_t;
  doublecomplex *t; 
  VALUE rblapack_q;
  doublecomplex *q; 
  VALUE rblapack_z;
  doublecomplex *z; 
  VALUE rblapack_lwork;
  integer lwork; 
  VALUE rblapack_alpha;
  doublecomplex *alpha; 
  VALUE rblapack_beta;
  doublecomplex *beta; 
  VALUE rblapack_work;
  doublecomplex *work; 
  VALUE rblapack_info;
  integer info; 
  VALUE rblapack_h_out__;
  doublecomplex *h_out__;
  VALUE rblapack_t_out__;
  doublecomplex *t_out__;
  VALUE rblapack_q_out__;
  doublecomplex *q_out__;
  VALUE rblapack_z_out__;
  doublecomplex *z_out__;
  doublereal *rwork;

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
      printf("%s\n", "USAGE:\n  alpha, beta, work, info, h, t, q, z = NumRu::Lapack.zhgeqz( job, compq, compz, ilo, ihi, h, t, q, z, lwork, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE ZHGEQZ( JOB, COMPQ, COMPZ, N, ILO, IHI, H, LDH, T, LDT, ALPHA, BETA, Q, LDQ, Z, LDZ, WORK, LWORK, RWORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  ZHGEQZ computes the eigenvalues of a complex matrix pair (H,T),\n*  where H is an upper Hessenberg matrix and T is upper triangular,\n*  using the single-shift QZ method.\n*  Matrix pairs of this type are produced by the reduction to\n*  generalized upper Hessenberg form of a complex matrix pair (A,B):\n*  \n*     A = Q1*H*Z1**H,  B = Q1*T*Z1**H,\n*  \n*  as computed by ZGGHRD.\n*  \n*  If JOB='S', then the Hessenberg-triangular pair (H,T) is\n*  also reduced to generalized Schur form,\n*  \n*     H = Q*S*Z**H,  T = Q*P*Z**H,\n*  \n*  where Q and Z are unitary matrices and S and P are upper triangular.\n*  \n*  Optionally, the unitary matrix Q from the generalized Schur\n*  factorization may be postmultiplied into an input matrix Q1, and the\n*  unitary matrix Z may be postmultiplied into an input matrix Z1.\n*  If Q1 and Z1 are the unitary matrices from ZGGHRD that reduced\n*  the matrix pair (A,B) to generalized Hessenberg form, then the output\n*  matrices Q1*Q and Z1*Z are the unitary factors from the generalized\n*  Schur factorization of (A,B):\n*  \n*     A = (Q1*Q)*S*(Z1*Z)**H,  B = (Q1*Q)*P*(Z1*Z)**H.\n*  \n*  To avoid overflow, eigenvalues of the matrix pair (H,T)\n*  (equivalently, of (A,B)) are computed as a pair of complex values\n*  (alpha,beta).  If beta is nonzero, lambda = alpha / beta is an\n*  eigenvalue of the generalized nonsymmetric eigenvalue problem (GNEP)\n*     A*x = lambda*B*x\n*  and if alpha is nonzero, mu = beta / alpha is an eigenvalue of the\n*  alternate form of the GNEP\n*     mu*A*y = B*y.\n*  The values of alpha and beta for the i-th eigenvalue can be read\n*  directly from the generalized Schur form:  alpha = S(i,i),\n*  beta = P(i,i).\n*\n*  Ref: C.B. Moler & G.W. Stewart, \"An Algorithm for Generalized Matrix\n*       Eigenvalue Problems\", SIAM J. Numer. Anal., 10(1973),\n*       pp. 241--256.\n*\n\n*  Arguments\n*  =========\n*\n*  JOB     (input) CHARACTER*1\n*          = 'E': Compute eigenvalues only;\n*          = 'S': Computer eigenvalues and the Schur form.\n*\n*  COMPQ   (input) CHARACTER*1\n*          = 'N': Left Schur vectors (Q) are not computed;\n*          = 'I': Q is initialized to the unit matrix and the matrix Q\n*                 of left Schur vectors of (H,T) is returned;\n*          = 'V': Q must contain a unitary matrix Q1 on entry and\n*                 the product Q1*Q is returned.\n*\n*  COMPZ   (input) CHARACTER*1\n*          = 'N': Right Schur vectors (Z) are not computed;\n*          = 'I': Q is initialized to the unit matrix and the matrix Z\n*                 of right Schur vectors of (H,T) is returned;\n*          = 'V': Z must contain a unitary matrix Z1 on entry and\n*                 the product Z1*Z is returned.\n*\n*  N       (input) INTEGER\n*          The order of the matrices H, T, Q, and Z.  N >= 0.\n*\n*  ILO     (input) INTEGER\n*  IHI     (input) INTEGER\n*          ILO and IHI mark the rows and columns of H which are in\n*          Hessenberg form.  It is assumed that A is already upper\n*          triangular in rows and columns 1:ILO-1 and IHI+1:N.\n*          If N > 0, 1 <= ILO <= IHI <= N; if N = 0, ILO=1 and IHI=0.\n*\n*  H       (input/output) COMPLEX*16 array, dimension (LDH, N)\n*          On entry, the N-by-N upper Hessenberg matrix H.\n*          On exit, if JOB = 'S', H contains the upper triangular\n*          matrix S from the generalized Schur factorization.\n*          If JOB = 'E', the diagonal of H matches that of S, but\n*          the rest of H is unspecified.\n*\n*  LDH     (input) INTEGER\n*          The leading dimension of the array H.  LDH >= max( 1, N ).\n*\n*  T       (input/output) COMPLEX*16 array, dimension (LDT, N)\n*          On entry, the N-by-N upper triangular matrix T.\n*          On exit, if JOB = 'S', T contains the upper triangular\n*          matrix P from the generalized Schur factorization.\n*          If JOB = 'E', the diagonal of T matches that of P, but\n*          the rest of T is unspecified.\n*\n*  LDT     (input) INTEGER\n*          The leading dimension of the array T.  LDT >= max( 1, N ).\n*\n*  ALPHA   (output) COMPLEX*16 array, dimension (N)\n*          The complex scalars alpha that define the eigenvalues of\n*          GNEP.  ALPHA(i) = S(i,i) in the generalized Schur\n*          factorization.\n*\n*  BETA    (output) COMPLEX*16 array, dimension (N)\n*          The real non-negative scalars beta that define the\n*          eigenvalues of GNEP.  BETA(i) = P(i,i) in the generalized\n*          Schur factorization.\n*\n*          Together, the quantities alpha = ALPHA(j) and beta = BETA(j)\n*          represent the j-th eigenvalue of the matrix pair (A,B), in\n*          one of the forms lambda = alpha/beta or mu = beta/alpha.\n*          Since either lambda or mu may overflow, they should not,\n*          in general, be computed.\n*\n*  Q       (input/output) COMPLEX*16 array, dimension (LDQ, N)\n*          On entry, if COMPZ = 'V', the unitary matrix Q1 used in the\n*          reduction of (A,B) to generalized Hessenberg form.\n*          On exit, if COMPZ = 'I', the unitary matrix of left Schur\n*          vectors of (H,T), and if COMPZ = 'V', the unitary matrix of\n*          left Schur vectors of (A,B).\n*          Not referenced if COMPZ = 'N'.\n*\n*  LDQ     (input) INTEGER\n*          The leading dimension of the array Q.  LDQ >= 1.\n*          If COMPQ='V' or 'I', then LDQ >= N.\n*\n*  Z       (input/output) COMPLEX*16 array, dimension (LDZ, N)\n*          On entry, if COMPZ = 'V', the unitary matrix Z1 used in the\n*          reduction of (A,B) to generalized Hessenberg form.\n*          On exit, if COMPZ = 'I', the unitary matrix of right Schur\n*          vectors of (H,T), and if COMPZ = 'V', the unitary matrix of\n*          right Schur vectors of (A,B).\n*          Not referenced if COMPZ = 'N'.\n*\n*  LDZ     (input) INTEGER\n*          The leading dimension of the array Z.  LDZ >= 1.\n*          If COMPZ='V' or 'I', then LDZ >= N.\n*\n*  WORK    (workspace/output) COMPLEX*16 array, dimension (MAX(1,LWORK))\n*          On exit, if INFO >= 0, WORK(1) returns the optimal LWORK.\n*\n*  LWORK   (input) INTEGER\n*          The dimension of the array WORK.  LWORK >= max(1,N).\n*\n*          If LWORK = -1, then a workspace query is assumed; the routine\n*          only calculates the optimal size of the WORK array, returns\n*          this value as the first entry of the WORK array, and no error\n*          message related to LWORK is issued by XERBLA.\n*\n*  RWORK   (workspace) DOUBLE PRECISION array, dimension (N)\n*\n*  INFO    (output) INTEGER\n*          = 0: successful exit\n*          < 0: if INFO = -i, the i-th argument had an illegal value\n*          = 1,...,N: the QZ iteration did not converge.  (H,T) is not\n*                     in Schur form, but ALPHA(i) and BETA(i),\n*                     i=INFO+1,...,N should be correct.\n*          = N+1,...,2*N: the shift calculation failed.  (H,T) is not\n*                     in Schur form, but ALPHA(i) and BETA(i),\n*                     i=INFO-N+1,...,N should be correct.\n*\n\n*  Further Details\n*  ===============\n*\n*  We assume that complex ABS works as long as its value is less than\n*  overflow.\n*\n*  =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  alpha, beta, work, info, h, t, q, z = NumRu::Lapack.zhgeqz( job, compq, compz, ilo, ihi, h, t, q, z, lwork, [:usage => usage, :help => help])\n");
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
  if (NA_TYPE(rblapack_z) != NA_DCOMPLEX)
    rblapack_z = na_change_type(rblapack_z, NA_DCOMPLEX);
  z = NA_PTR_TYPE(rblapack_z, doublecomplex*);
  compq = StringValueCStr(rblapack_compq)[0];
  lwork = NUM2INT(rblapack_lwork);
  if (!NA_IsNArray(rblapack_q))
    rb_raise(rb_eArgError, "q (8th argument) must be NArray");
  if (NA_RANK(rblapack_q) != 2)
    rb_raise(rb_eArgError, "rank of q (8th argument) must be %d", 2);
  if (NA_SHAPE1(rblapack_q) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of q must be the same as shape 1 of z");
  ldq = NA_SHAPE0(rblapack_q);
  if (NA_TYPE(rblapack_q) != NA_DCOMPLEX)
    rblapack_q = na_change_type(rblapack_q, NA_DCOMPLEX);
  q = NA_PTR_TYPE(rblapack_q, doublecomplex*);
  job = StringValueCStr(rblapack_job)[0];
  if (!NA_IsNArray(rblapack_h))
    rb_raise(rb_eArgError, "h (6th argument) must be NArray");
  if (NA_RANK(rblapack_h) != 2)
    rb_raise(rb_eArgError, "rank of h (6th argument) must be %d", 2);
  if (NA_SHAPE1(rblapack_h) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of h must be the same as shape 1 of z");
  ldh = NA_SHAPE0(rblapack_h);
  if (NA_TYPE(rblapack_h) != NA_DCOMPLEX)
    rblapack_h = na_change_type(rblapack_h, NA_DCOMPLEX);
  h = NA_PTR_TYPE(rblapack_h, doublecomplex*);
  if (!NA_IsNArray(rblapack_t))
    rb_raise(rb_eArgError, "t (7th argument) must be NArray");
  if (NA_RANK(rblapack_t) != 2)
    rb_raise(rb_eArgError, "rank of t (7th argument) must be %d", 2);
  if (NA_SHAPE1(rblapack_t) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of t must be the same as shape 1 of z");
  ldt = NA_SHAPE0(rblapack_t);
  if (NA_TYPE(rblapack_t) != NA_DCOMPLEX)
    rblapack_t = na_change_type(rblapack_t, NA_DCOMPLEX);
  t = NA_PTR_TYPE(rblapack_t, doublecomplex*);
  ihi = NUM2INT(rblapack_ihi);
  {
    int shape[1];
    shape[0] = n;
    rblapack_alpha = na_make_object(NA_DCOMPLEX, 1, shape, cNArray);
  }
  alpha = NA_PTR_TYPE(rblapack_alpha, doublecomplex*);
  {
    int shape[1];
    shape[0] = n;
    rblapack_beta = na_make_object(NA_DCOMPLEX, 1, shape, cNArray);
  }
  beta = NA_PTR_TYPE(rblapack_beta, doublecomplex*);
  {
    int shape[1];
    shape[0] = MAX(1,lwork);
    rblapack_work = na_make_object(NA_DCOMPLEX, 1, shape, cNArray);
  }
  work = NA_PTR_TYPE(rblapack_work, doublecomplex*);
  {
    int shape[2];
    shape[0] = ldh;
    shape[1] = n;
    rblapack_h_out__ = na_make_object(NA_DCOMPLEX, 2, shape, cNArray);
  }
  h_out__ = NA_PTR_TYPE(rblapack_h_out__, doublecomplex*);
  MEMCPY(h_out__, h, doublecomplex, NA_TOTAL(rblapack_h));
  rblapack_h = rblapack_h_out__;
  h = h_out__;
  {
    int shape[2];
    shape[0] = ldt;
    shape[1] = n;
    rblapack_t_out__ = na_make_object(NA_DCOMPLEX, 2, shape, cNArray);
  }
  t_out__ = NA_PTR_TYPE(rblapack_t_out__, doublecomplex*);
  MEMCPY(t_out__, t, doublecomplex, NA_TOTAL(rblapack_t));
  rblapack_t = rblapack_t_out__;
  t = t_out__;
  {
    int shape[2];
    shape[0] = ldq;
    shape[1] = n;
    rblapack_q_out__ = na_make_object(NA_DCOMPLEX, 2, shape, cNArray);
  }
  q_out__ = NA_PTR_TYPE(rblapack_q_out__, doublecomplex*);
  MEMCPY(q_out__, q, doublecomplex, NA_TOTAL(rblapack_q));
  rblapack_q = rblapack_q_out__;
  q = q_out__;
  {
    int shape[2];
    shape[0] = ldz;
    shape[1] = n;
    rblapack_z_out__ = na_make_object(NA_DCOMPLEX, 2, shape, cNArray);
  }
  z_out__ = NA_PTR_TYPE(rblapack_z_out__, doublecomplex*);
  MEMCPY(z_out__, z, doublecomplex, NA_TOTAL(rblapack_z));
  rblapack_z = rblapack_z_out__;
  z = z_out__;
  rwork = ALLOC_N(doublereal, (n));

  zhgeqz_(&job, &compq, &compz, &n, &ilo, &ihi, h, &ldh, t, &ldt, alpha, beta, q, &ldq, z, &ldz, work, &lwork, rwork, &info);

  free(rwork);
  rblapack_info = INT2NUM(info);
  return rb_ary_new3(8, rblapack_alpha, rblapack_beta, rblapack_work, rblapack_info, rblapack_h, rblapack_t, rblapack_q, rblapack_z);
}

void
init_lapack_zhgeqz(VALUE mLapack){
  rb_define_module_function(mLapack, "zhgeqz", rblapack_zhgeqz, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}

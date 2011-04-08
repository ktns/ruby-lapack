#include "rb_lapack.h"

extern VOID strsen_(char *job, char *compq, logical *select, integer *n, real *t, integer *ldt, real *q, integer *ldq, real *wr, real *wi, integer *m, real *s, real *sep, real *work, integer *lwork, integer *iwork, integer *liwork, integer *info);

static VALUE sHelp, sUsage;

static VALUE
rblapack_strsen(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_job;
  char job; 
  VALUE rblapack_compq;
  char compq; 
  VALUE rblapack_select;
  logical *select; 
  VALUE rblapack_t;
  real *t; 
  VALUE rblapack_q;
  real *q; 
  VALUE rblapack_lwork;
  integer lwork; 
  VALUE rblapack_liwork;
  integer liwork; 
  VALUE rblapack_wr;
  real *wr; 
  VALUE rblapack_wi;
  real *wi; 
  VALUE rblapack_m;
  integer m; 
  VALUE rblapack_s;
  real s; 
  VALUE rblapack_sep;
  real sep; 
  VALUE rblapack_work;
  real *work; 
  VALUE rblapack_info;
  integer info; 
  VALUE rblapack_t_out__;
  real *t_out__;
  VALUE rblapack_q_out__;
  real *q_out__;
  integer *iwork;

  integer n;
  integer ldt;
  integer ldq;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  wr, wi, m, s, sep, work, info, t, q = NumRu::Lapack.strsen( job, compq, select, t, q, lwork, liwork, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE STRSEN( JOB, COMPQ, SELECT, N, T, LDT, Q, LDQ, WR, WI, M, S, SEP, WORK, LWORK, IWORK, LIWORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  STRSEN reorders the real Schur factorization of a real matrix\n*  A = Q*T*Q**T, so that a selected cluster of eigenvalues appears in\n*  the leading diagonal blocks of the upper quasi-triangular matrix T,\n*  and the leading columns of Q form an orthonormal basis of the\n*  corresponding right invariant subspace.\n*\n*  Optionally the routine computes the reciprocal condition numbers of\n*  the cluster of eigenvalues and/or the invariant subspace.\n*\n*  T must be in Schur canonical form (as returned by SHSEQR), that is,\n*  block upper triangular with 1-by-1 and 2-by-2 diagonal blocks; each\n*  2-by-2 diagonal block has its diagonal elemnts equal and its\n*  off-diagonal elements of opposite sign.\n*\n\n*  Arguments\n*  =========\n*\n*  JOB     (input) CHARACTER*1\n*          Specifies whether condition numbers are required for the\n*          cluster of eigenvalues (S) or the invariant subspace (SEP):\n*          = 'N': none;\n*          = 'E': for eigenvalues only (S);\n*          = 'V': for invariant subspace only (SEP);\n*          = 'B': for both eigenvalues and invariant subspace (S and\n*                 SEP).\n*\n*  COMPQ   (input) CHARACTER*1\n*          = 'V': update the matrix Q of Schur vectors;\n*          = 'N': do not update Q.\n*\n*  SELECT  (input) LOGICAL array, dimension (N)\n*          SELECT specifies the eigenvalues in the selected cluster. To\n*          select a real eigenvalue w(j), SELECT(j) must be set to\n*          .TRUE.. To select a complex conjugate pair of eigenvalues\n*          w(j) and w(j+1), corresponding to a 2-by-2 diagonal block,\n*          either SELECT(j) or SELECT(j+1) or both must be set to\n*          .TRUE.; a complex conjugate pair of eigenvalues must be\n*          either both included in the cluster or both excluded.\n*\n*  N       (input) INTEGER\n*          The order of the matrix T. N >= 0.\n*\n*  T       (input/output) REAL array, dimension (LDT,N)\n*          On entry, the upper quasi-triangular matrix T, in Schur\n*          canonical form.\n*          On exit, T is overwritten by the reordered matrix T, again in\n*          Schur canonical form, with the selected eigenvalues in the\n*          leading diagonal blocks.\n*\n*  LDT     (input) INTEGER\n*          The leading dimension of the array T. LDT >= max(1,N).\n*\n*  Q       (input/output) REAL array, dimension (LDQ,N)\n*          On entry, if COMPQ = 'V', the matrix Q of Schur vectors.\n*          On exit, if COMPQ = 'V', Q has been postmultiplied by the\n*          orthogonal transformation matrix which reorders T; the\n*          leading M columns of Q form an orthonormal basis for the\n*          specified invariant subspace.\n*          If COMPQ = 'N', Q is not referenced.\n*\n*  LDQ     (input) INTEGER\n*          The leading dimension of the array Q.\n*          LDQ >= 1; and if COMPQ = 'V', LDQ >= N.\n*\n*  WR      (output) REAL array, dimension (N)\n*  WI      (output) REAL array, dimension (N)\n*          The real and imaginary parts, respectively, of the reordered\n*          eigenvalues of T. The eigenvalues are stored in the same\n*          order as on the diagonal of T, with WR(i) = T(i,i) and, if\n*          T(i:i+1,i:i+1) is a 2-by-2 diagonal block, WI(i) > 0 and\n*          WI(i+1) = -WI(i). Note that if a complex eigenvalue is\n*          sufficiently ill-conditioned, then its value may differ\n*          significantly from its value before reordering.\n*\n*  M       (output) INTEGER\n*          The dimension of the specified invariant subspace.\n*          0 < = M <= N.\n*\n*  S       (output) REAL\n*          If JOB = 'E' or 'B', S is a lower bound on the reciprocal\n*          condition number for the selected cluster of eigenvalues.\n*          S cannot underestimate the true reciprocal condition number\n*          by more than a factor of sqrt(N). If M = 0 or N, S = 1.\n*          If JOB = 'N' or 'V', S is not referenced.\n*\n*  SEP     (output) REAL\n*          If JOB = 'V' or 'B', SEP is the estimated reciprocal\n*          condition number of the specified invariant subspace. If\n*          M = 0 or N, SEP = norm(T).\n*          If JOB = 'N' or 'E', SEP is not referenced.\n*\n*  WORK    (workspace/output) REAL array, dimension (MAX(1,LWORK))\n*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.\n*\n*  LWORK   (input) INTEGER\n*          The dimension of the array WORK.\n*          If JOB = 'N', LWORK >= max(1,N);\n*          if JOB = 'E', LWORK >= max(1,M*(N-M));\n*          if JOB = 'V' or 'B', LWORK >= max(1,2*M*(N-M)).\n*\n*          If LWORK = -1, then a workspace query is assumed; the routine\n*          only calculates the optimal size of the WORK array, returns\n*          this value as the first entry of the WORK array, and no error\n*          message related to LWORK is issued by XERBLA.\n*\n*  IWORK   (workspace) INTEGER array, dimension (MAX(1,LIWORK))\n*          On exit, if INFO = 0, IWORK(1) returns the optimal LIWORK.\n*\n*  LIWORK  (input) INTEGER\n*          The dimension of the array IWORK.\n*          If JOB = 'N' or 'E', LIWORK >= 1;\n*          if JOB = 'V' or 'B', LIWORK >= max(1,M*(N-M)).\n*\n*          If LIWORK = -1, then a workspace query is assumed; the\n*          routine only calculates the optimal size of the IWORK array,\n*          returns this value as the first entry of the IWORK array, and\n*          no error message related to LIWORK is issued by XERBLA.\n*\n*  INFO    (output) INTEGER\n*          = 0: successful exit\n*          < 0: if INFO = -i, the i-th argument had an illegal value\n*          = 1: reordering of T failed because some eigenvalues are too\n*               close to separate (the problem is very ill-conditioned);\n*               T may have been partially reordered, and WR and WI\n*               contain the eigenvalues in the same order as in T; S and\n*               SEP (if requested) are set to zero.\n*\n\n*  Further Details\n*  ===============\n*\n*  STRSEN first collects the selected eigenvalues by computing an\n*  orthogonal transformation Z to move them to the top left corner of T.\n*  In other words, the selected eigenvalues are the eigenvalues of T11\n*  in:\n*\n*                Z'*T*Z = ( T11 T12 ) n1\n*                         (  0  T22 ) n2\n*                            n1  n2\n*\n*  where N = n1+n2 and Z' means the transpose of Z. The first n1 columns\n*  of Z span the specified invariant subspace of T.\n*\n*  If T has been obtained from the real Schur factorization of a matrix\n*  A = Q*T*Q', then the reordered real Schur factorization of A is given\n*  by A = (Q*Z)*(Z'*T*Z)*(Q*Z)', and the first n1 columns of Q*Z span\n*  the corresponding invariant subspace of A.\n*\n*  The reciprocal condition number of the average of the eigenvalues of\n*  T11 may be returned in S. S lies between 0 (very badly conditioned)\n*  and 1 (very well conditioned). It is computed as follows. First we\n*  compute R so that\n*\n*                         P = ( I  R ) n1\n*                             ( 0  0 ) n2\n*                               n1 n2\n*\n*  is the projector on the invariant subspace associated with T11.\n*  R is the solution of the Sylvester equation:\n*\n*                        T11*R - R*T22 = T12.\n*\n*  Let F-norm(M) denote the Frobenius-norm of M and 2-norm(M) denote\n*  the two-norm of M. Then S is computed as the lower bound\n*\n*                      (1 + F-norm(R)**2)**(-1/2)\n*\n*  on the reciprocal of 2-norm(P), the true reciprocal condition number.\n*  S cannot underestimate 1 / 2-norm(P) by more than a factor of\n*  sqrt(N).\n*\n*  An approximate error bound for the computed average of the\n*  eigenvalues of T11 is\n*\n*                         EPS * norm(T) / S\n*\n*  where EPS is the machine precision.\n*\n*  The reciprocal condition number of the right invariant subspace\n*  spanned by the first n1 columns of Z (or of Q*Z) is returned in SEP.\n*  SEP is defined as the separation of T11 and T22:\n*\n*                     sep( T11, T22 ) = sigma-min( C )\n*\n*  where sigma-min(C) is the smallest singular value of the\n*  n1*n2-by-n1*n2 matrix\n*\n*     C  = kprod( I(n2), T11 ) - kprod( transpose(T22), I(n1) )\n*\n*  I(m) is an m by m identity matrix, and kprod denotes the Kronecker\n*  product. We estimate sigma-min(C) by the reciprocal of an estimate of\n*  the 1-norm of inverse(C). The true reciprocal 1-norm of inverse(C)\n*  cannot differ from sigma-min(C) by more than a factor of sqrt(n1*n2).\n*\n*  When SEP is small, small changes in T can cause large changes in\n*  the invariant subspace. An approximate bound on the maximum angular\n*  error in the computed right invariant subspace is\n*\n*                      EPS * norm(T) / SEP\n*\n*  =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  wr, wi, m, s, sep, work, info, t, q = NumRu::Lapack.strsen( job, compq, select, t, q, lwork, liwork, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 7)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 7)", argc);
  rblapack_job = argv[0];
  rblapack_compq = argv[1];
  rblapack_select = argv[2];
  rblapack_t = argv[3];
  rblapack_q = argv[4];
  rblapack_lwork = argv[5];
  rblapack_liwork = argv[6];
  if (rb_options != Qnil) {
  }

  liwork = NUM2INT(rblapack_liwork);
  compq = StringValueCStr(rblapack_compq)[0];
  lwork = NUM2INT(rblapack_lwork);
  if (!NA_IsNArray(rblapack_q))
    rb_raise(rb_eArgError, "q (5th argument) must be NArray");
  if (NA_RANK(rblapack_q) != 2)
    rb_raise(rb_eArgError, "rank of q (5th argument) must be %d", 2);
  n = NA_SHAPE1(rblapack_q);
  ldq = NA_SHAPE0(rblapack_q);
  if (NA_TYPE(rblapack_q) != NA_SFLOAT)
    rblapack_q = na_change_type(rblapack_q, NA_SFLOAT);
  q = NA_PTR_TYPE(rblapack_q, real*);
  job = StringValueCStr(rblapack_job)[0];
  if (!NA_IsNArray(rblapack_select))
    rb_raise(rb_eArgError, "select (3th argument) must be NArray");
  if (NA_RANK(rblapack_select) != 1)
    rb_raise(rb_eArgError, "rank of select (3th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_select) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of select must be the same as shape 1 of q");
  if (NA_TYPE(rblapack_select) != NA_LINT)
    rblapack_select = na_change_type(rblapack_select, NA_LINT);
  select = NA_PTR_TYPE(rblapack_select, logical*);
  if (!NA_IsNArray(rblapack_t))
    rb_raise(rb_eArgError, "t (4th argument) must be NArray");
  if (NA_RANK(rblapack_t) != 2)
    rb_raise(rb_eArgError, "rank of t (4th argument) must be %d", 2);
  if (NA_SHAPE1(rblapack_t) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of t must be the same as shape 1 of q");
  ldt = NA_SHAPE0(rblapack_t);
  if (NA_TYPE(rblapack_t) != NA_SFLOAT)
    rblapack_t = na_change_type(rblapack_t, NA_SFLOAT);
  t = NA_PTR_TYPE(rblapack_t, real*);
  {
    int shape[1];
    shape[0] = n;
    rblapack_wr = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  wr = NA_PTR_TYPE(rblapack_wr, real*);
  {
    int shape[1];
    shape[0] = n;
    rblapack_wi = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  wi = NA_PTR_TYPE(rblapack_wi, real*);
  {
    int shape[1];
    shape[0] = MAX(1,lwork);
    rblapack_work = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  work = NA_PTR_TYPE(rblapack_work, real*);
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
  iwork = ALLOC_N(integer, (MAX(1,liwork)));

  strsen_(&job, &compq, select, &n, t, &ldt, q, &ldq, wr, wi, &m, &s, &sep, work, &lwork, iwork, &liwork, &info);

  free(iwork);
  rblapack_m = INT2NUM(m);
  rblapack_s = rb_float_new((double)s);
  rblapack_sep = rb_float_new((double)sep);
  rblapack_info = INT2NUM(info);
  return rb_ary_new3(9, rblapack_wr, rblapack_wi, rblapack_m, rblapack_s, rblapack_sep, rblapack_work, rblapack_info, rblapack_t, rblapack_q);
}

void
init_lapack_strsen(VALUE mLapack){
  rb_define_module_function(mLapack, "strsen", rblapack_strsen, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}

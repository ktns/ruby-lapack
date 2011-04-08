#include "rb_lapack.h"

extern VOID cgegs_(char *jobvsl, char *jobvsr, integer *n, complex *a, integer *lda, complex *b, integer *ldb, complex *alpha, complex *beta, complex *vsl, integer *ldvsl, complex *vsr, integer *ldvsr, complex *work, integer *lwork, real *rwork, integer *info);

static VALUE sHelp, sUsage;

static VALUE
rblapack_cgegs(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_jobvsl;
  char jobvsl; 
  VALUE rblapack_jobvsr;
  char jobvsr; 
  VALUE rblapack_a;
  complex *a; 
  VALUE rblapack_b;
  complex *b; 
  VALUE rblapack_lwork;
  integer lwork; 
  VALUE rblapack_alpha;
  complex *alpha; 
  VALUE rblapack_beta;
  complex *beta; 
  VALUE rblapack_vsl;
  complex *vsl; 
  VALUE rblapack_vsr;
  complex *vsr; 
  VALUE rblapack_work;
  complex *work; 
  VALUE rblapack_info;
  integer info; 
  VALUE rblapack_a_out__;
  complex *a_out__;
  VALUE rblapack_b_out__;
  complex *b_out__;
  real *rwork;

  integer lda;
  integer n;
  integer ldb;
  integer ldvsl;
  integer ldvsr;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  alpha, beta, vsl, vsr, work, info, a, b = NumRu::Lapack.cgegs( jobvsl, jobvsr, a, b, lwork, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE CGEGS( JOBVSL, JOBVSR, N, A, LDA, B, LDB, ALPHA, BETA, VSL, LDVSL, VSR, LDVSR, WORK, LWORK, RWORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  This routine is deprecated and has been replaced by routine CGGES.\n*\n*  CGEGS computes the eigenvalues, Schur form, and, optionally, the\n*  left and or/right Schur vectors of a complex matrix pair (A,B).\n*  Given two square matrices A and B, the generalized Schur\n*  factorization has the form\n*  \n*     A = Q*S*Z**H,  B = Q*T*Z**H\n*  \n*  where Q and Z are unitary matrices and S and T are upper triangular.\n*  The columns of Q are the left Schur vectors\n*  and the columns of Z are the right Schur vectors.\n*  \n*  If only the eigenvalues of (A,B) are needed, the driver routine\n*  CGEGV should be used instead.  See CGEGV for a description of the\n*  eigenvalues of the generalized nonsymmetric eigenvalue problem\n*  (GNEP).\n*\n\n*  Arguments\n*  =========\n*\n*  JOBVSL   (input) CHARACTER*1\n*          = 'N':  do not compute the left Schur vectors;\n*          = 'V':  compute the left Schur vectors (returned in VSL).\n*\n*  JOBVSR   (input) CHARACTER*1\n*          = 'N':  do not compute the right Schur vectors;\n*          = 'V':  compute the right Schur vectors (returned in VSR).\n*\n*  N       (input) INTEGER\n*          The order of the matrices A, B, VSL, and VSR.  N >= 0.\n*\n*  A       (input/output) COMPLEX array, dimension (LDA, N)\n*          On entry, the matrix A.\n*          On exit, the upper triangular matrix S from the generalized\n*          Schur factorization.\n*\n*  LDA     (input) INTEGER\n*          The leading dimension of A.  LDA >= max(1,N).\n*\n*  B       (input/output) COMPLEX array, dimension (LDB, N)\n*          On entry, the matrix B.\n*          On exit, the upper triangular matrix T from the generalized\n*          Schur factorization.\n*\n*  LDB     (input) INTEGER\n*          The leading dimension of B.  LDB >= max(1,N).\n*\n*  ALPHA   (output) COMPLEX array, dimension (N)\n*          The complex scalars alpha that define the eigenvalues of\n*          GNEP.  ALPHA(j) = S(j,j), the diagonal element of the Schur\n*          form of A.\n*\n*  BETA    (output) COMPLEX array, dimension (N)\n*          The non-negative real scalars beta that define the\n*          eigenvalues of GNEP.  BETA(j) = T(j,j), the diagonal element\n*          of the triangular factor T.\n*\n*          Together, the quantities alpha = ALPHA(j) and beta = BETA(j)\n*          represent the j-th eigenvalue of the matrix pair (A,B), in\n*          one of the forms lambda = alpha/beta or mu = beta/alpha.\n*          Since either lambda or mu may overflow, they should not,\n*          in general, be computed.\n*\n*  VSL     (output) COMPLEX array, dimension (LDVSL,N)\n*          If JOBVSL = 'V', the matrix of left Schur vectors Q.\n*          Not referenced if JOBVSL = 'N'.\n*\n*  LDVSL   (input) INTEGER\n*          The leading dimension of the matrix VSL. LDVSL >= 1, and\n*          if JOBVSL = 'V', LDVSL >= N.\n*\n*  VSR     (output) COMPLEX array, dimension (LDVSR,N)\n*          If JOBVSR = 'V', the matrix of right Schur vectors Z.\n*          Not referenced if JOBVSR = 'N'.\n*\n*  LDVSR   (input) INTEGER\n*          The leading dimension of the matrix VSR. LDVSR >= 1, and\n*          if JOBVSR = 'V', LDVSR >= N.\n*\n*  WORK    (workspace/output) COMPLEX array, dimension (MAX(1,LWORK))\n*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.\n*\n*  LWORK   (input) INTEGER\n*          The dimension of the array WORK.  LWORK >= max(1,2*N).\n*          For good performance, LWORK must generally be larger.\n*          To compute the optimal value of LWORK, call ILAENV to get\n*          blocksizes (for CGEQRF, CUNMQR, and CUNGQR.)  Then compute:\n*          NB  -- MAX of the blocksizes for CGEQRF, CUNMQR, and CUNGQR;\n*          the optimal LWORK is N*(NB+1).\n*\n*          If LWORK = -1, then a workspace query is assumed; the routine\n*          only calculates the optimal size of the WORK array, returns\n*          this value as the first entry of the WORK array, and no error\n*          message related to LWORK is issued by XERBLA.\n*\n*  RWORK   (workspace) REAL array, dimension (3*N)\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit\n*          < 0:  if INFO = -i, the i-th argument had an illegal value.\n*          =1,...,N:\n*                The QZ iteration failed.  (A,B) are not in Schur\n*                form, but ALPHA(j) and BETA(j) should be correct for\n*                j=INFO+1,...,N.\n*          > N:  errors that usually indicate LAPACK problems:\n*                =N+1: error return from CGGBAL\n*                =N+2: error return from CGEQRF\n*                =N+3: error return from CUNMQR\n*                =N+4: error return from CUNGQR\n*                =N+5: error return from CGGHRD\n*                =N+6: error return from CHGEQZ (other than failed\n*                                               iteration)\n*                =N+7: error return from CGGBAK (computing VSL)\n*                =N+8: error return from CGGBAK (computing VSR)\n*                =N+9: error return from CLASCL (various places)\n*\n\n*  =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  alpha, beta, vsl, vsr, work, info, a, b = NumRu::Lapack.cgegs( jobvsl, jobvsr, a, b, lwork, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 5)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 5)", argc);
  rblapack_jobvsl = argv[0];
  rblapack_jobvsr = argv[1];
  rblapack_a = argv[2];
  rblapack_b = argv[3];
  rblapack_lwork = argv[4];
  if (rb_options != Qnil) {
  }

  jobvsl = StringValueCStr(rblapack_jobvsl)[0];
  if (!NA_IsNArray(rblapack_a))
    rb_raise(rb_eArgError, "a (3th argument) must be NArray");
  if (NA_RANK(rblapack_a) != 2)
    rb_raise(rb_eArgError, "rank of a (3th argument) must be %d", 2);
  n = NA_SHAPE1(rblapack_a);
  lda = NA_SHAPE0(rblapack_a);
  if (NA_TYPE(rblapack_a) != NA_SCOMPLEX)
    rblapack_a = na_change_type(rblapack_a, NA_SCOMPLEX);
  a = NA_PTR_TYPE(rblapack_a, complex*);
  if (!NA_IsNArray(rblapack_b))
    rb_raise(rb_eArgError, "b (4th argument) must be NArray");
  if (NA_RANK(rblapack_b) != 2)
    rb_raise(rb_eArgError, "rank of b (4th argument) must be %d", 2);
  if (NA_SHAPE1(rblapack_b) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of b must be the same as shape 1 of a");
  ldb = NA_SHAPE0(rblapack_b);
  if (NA_TYPE(rblapack_b) != NA_SCOMPLEX)
    rblapack_b = na_change_type(rblapack_b, NA_SCOMPLEX);
  b = NA_PTR_TYPE(rblapack_b, complex*);
  lwork = NUM2INT(rblapack_lwork);
  jobvsr = StringValueCStr(rblapack_jobvsr)[0];
  ldvsr = lsame_(&jobvsr,"V") ? n : 1;
  ldvsl = lsame_(&jobvsl,"V") ? n : 1;
  {
    int shape[1];
    shape[0] = n;
    rblapack_alpha = na_make_object(NA_SCOMPLEX, 1, shape, cNArray);
  }
  alpha = NA_PTR_TYPE(rblapack_alpha, complex*);
  {
    int shape[1];
    shape[0] = n;
    rblapack_beta = na_make_object(NA_SCOMPLEX, 1, shape, cNArray);
  }
  beta = NA_PTR_TYPE(rblapack_beta, complex*);
  {
    int shape[2];
    shape[0] = ldvsl;
    shape[1] = n;
    rblapack_vsl = na_make_object(NA_SCOMPLEX, 2, shape, cNArray);
  }
  vsl = NA_PTR_TYPE(rblapack_vsl, complex*);
  {
    int shape[2];
    shape[0] = ldvsr;
    shape[1] = n;
    rblapack_vsr = na_make_object(NA_SCOMPLEX, 2, shape, cNArray);
  }
  vsr = NA_PTR_TYPE(rblapack_vsr, complex*);
  {
    int shape[1];
    shape[0] = MAX(1,lwork);
    rblapack_work = na_make_object(NA_SCOMPLEX, 1, shape, cNArray);
  }
  work = NA_PTR_TYPE(rblapack_work, complex*);
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
  rwork = ALLOC_N(real, (3*n));

  cgegs_(&jobvsl, &jobvsr, &n, a, &lda, b, &ldb, alpha, beta, vsl, &ldvsl, vsr, &ldvsr, work, &lwork, rwork, &info);

  free(rwork);
  rblapack_info = INT2NUM(info);
  return rb_ary_new3(8, rblapack_alpha, rblapack_beta, rblapack_vsl, rblapack_vsr, rblapack_work, rblapack_info, rblapack_a, rblapack_b);
}

void
init_lapack_cgegs(VALUE mLapack){
  rb_define_module_function(mLapack, "cgegs", rblapack_cgegs, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}

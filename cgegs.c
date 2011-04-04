#include "rb_lapack.h"

extern VOID cgegs_(char *jobvsl, char *jobvsr, integer *n, complex *a, integer *lda, complex *b, integer *ldb, complex *alpha, complex *beta, complex *vsl, integer *ldvsl, complex *vsr, integer *ldvsr, complex *work, integer *lwork, real *rwork, integer *info);

static VALUE
rb_cgegs(int argc, VALUE *argv, VALUE self){
  VALUE rb_jobvsl;
  char jobvsl; 
  VALUE rb_jobvsr;
  char jobvsr; 
  VALUE rb_a;
  complex *a; 
  VALUE rb_b;
  complex *b; 
  VALUE rb_lwork;
  integer lwork; 
  VALUE rb_alpha;
  complex *alpha; 
  VALUE rb_beta;
  complex *beta; 
  VALUE rb_vsl;
  complex *vsl; 
  VALUE rb_vsr;
  complex *vsr; 
  VALUE rb_work;
  complex *work; 
  VALUE rb_info;
  integer info; 
  VALUE rb_a_out__;
  complex *a_out__;
  VALUE rb_b_out__;
  complex *b_out__;
  real *rwork;

  integer lda;
  integer n;
  integer ldb;
  integer ldvsl;
  integer ldvsr;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  alpha, beta, vsl, vsr, work, info, a, b = NumRu::Lapack.cgegs( jobvsl, jobvsr, a, b, lwork)\n    or\n  NumRu::Lapack.cgegs  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE CGEGS( JOBVSL, JOBVSR, N, A, LDA, B, LDB, ALPHA, BETA, VSL, LDVSL, VSR, LDVSR, WORK, LWORK, RWORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  This routine is deprecated and has been replaced by routine CGGES.\n*\n*  CGEGS computes the eigenvalues, Schur form, and, optionally, the\n*  left and or/right Schur vectors of a complex matrix pair (A,B).\n*  Given two square matrices A and B, the generalized Schur\n*  factorization has the form\n*  \n*     A = Q*S*Z**H,  B = Q*T*Z**H\n*  \n*  where Q and Z are unitary matrices and S and T are upper triangular.\n*  The columns of Q are the left Schur vectors\n*  and the columns of Z are the right Schur vectors.\n*  \n*  If only the eigenvalues of (A,B) are needed, the driver routine\n*  CGEGV should be used instead.  See CGEGV for a description of the\n*  eigenvalues of the generalized nonsymmetric eigenvalue problem\n*  (GNEP).\n*\n\n*  Arguments\n*  =========\n*\n*  JOBVSL   (input) CHARACTER*1\n*          = 'N':  do not compute the left Schur vectors;\n*          = 'V':  compute the left Schur vectors (returned in VSL).\n*\n*  JOBVSR   (input) CHARACTER*1\n*          = 'N':  do not compute the right Schur vectors;\n*          = 'V':  compute the right Schur vectors (returned in VSR).\n*\n*  N       (input) INTEGER\n*          The order of the matrices A, B, VSL, and VSR.  N >= 0.\n*\n*  A       (input/output) COMPLEX array, dimension (LDA, N)\n*          On entry, the matrix A.\n*          On exit, the upper triangular matrix S from the generalized\n*          Schur factorization.\n*\n*  LDA     (input) INTEGER\n*          The leading dimension of A.  LDA >= max(1,N).\n*\n*  B       (input/output) COMPLEX array, dimension (LDB, N)\n*          On entry, the matrix B.\n*          On exit, the upper triangular matrix T from the generalized\n*          Schur factorization.\n*\n*  LDB     (input) INTEGER\n*          The leading dimension of B.  LDB >= max(1,N).\n*\n*  ALPHA   (output) COMPLEX array, dimension (N)\n*          The complex scalars alpha that define the eigenvalues of\n*          GNEP.  ALPHA(j) = S(j,j), the diagonal element of the Schur\n*          form of A.\n*\n*  BETA    (output) COMPLEX array, dimension (N)\n*          The non-negative real scalars beta that define the\n*          eigenvalues of GNEP.  BETA(j) = T(j,j), the diagonal element\n*          of the triangular factor T.\n*\n*          Together, the quantities alpha = ALPHA(j) and beta = BETA(j)\n*          represent the j-th eigenvalue of the matrix pair (A,B), in\n*          one of the forms lambda = alpha/beta or mu = beta/alpha.\n*          Since either lambda or mu may overflow, they should not,\n*          in general, be computed.\n*\n*  VSL     (output) COMPLEX array, dimension (LDVSL,N)\n*          If JOBVSL = 'V', the matrix of left Schur vectors Q.\n*          Not referenced if JOBVSL = 'N'.\n*\n*  LDVSL   (input) INTEGER\n*          The leading dimension of the matrix VSL. LDVSL >= 1, and\n*          if JOBVSL = 'V', LDVSL >= N.\n*\n*  VSR     (output) COMPLEX array, dimension (LDVSR,N)\n*          If JOBVSR = 'V', the matrix of right Schur vectors Z.\n*          Not referenced if JOBVSR = 'N'.\n*\n*  LDVSR   (input) INTEGER\n*          The leading dimension of the matrix VSR. LDVSR >= 1, and\n*          if JOBVSR = 'V', LDVSR >= N.\n*\n*  WORK    (workspace/output) COMPLEX array, dimension (MAX(1,LWORK))\n*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.\n*\n*  LWORK   (input) INTEGER\n*          The dimension of the array WORK.  LWORK >= max(1,2*N).\n*          For good performance, LWORK must generally be larger.\n*          To compute the optimal value of LWORK, call ILAENV to get\n*          blocksizes (for CGEQRF, CUNMQR, and CUNGQR.)  Then compute:\n*          NB  -- MAX of the blocksizes for CGEQRF, CUNMQR, and CUNGQR;\n*          the optimal LWORK is N*(NB+1).\n*\n*          If LWORK = -1, then a workspace query is assumed; the routine\n*          only calculates the optimal size of the WORK array, returns\n*          this value as the first entry of the WORK array, and no error\n*          message related to LWORK is issued by XERBLA.\n*\n*  RWORK   (workspace) REAL array, dimension (3*N)\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit\n*          < 0:  if INFO = -i, the i-th argument had an illegal value.\n*          =1,...,N:\n*                The QZ iteration failed.  (A,B) are not in Schur\n*                form, but ALPHA(j) and BETA(j) should be correct for\n*                j=INFO+1,...,N.\n*          > N:  errors that usually indicate LAPACK problems:\n*                =N+1: error return from CGGBAL\n*                =N+2: error return from CGEQRF\n*                =N+3: error return from CUNMQR\n*                =N+4: error return from CUNGQR\n*                =N+5: error return from CGGHRD\n*                =N+6: error return from CHGEQZ (other than failed\n*                                               iteration)\n*                =N+7: error return from CGGBAK (computing VSL)\n*                =N+8: error return from CGGBAK (computing VSR)\n*                =N+9: error return from CLASCL (various places)\n*\n\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 5)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 5)", argc);
  rb_jobvsl = argv[0];
  rb_jobvsr = argv[1];
  rb_a = argv[2];
  rb_b = argv[3];
  rb_lwork = argv[4];

  jobvsl = StringValueCStr(rb_jobvsl)[0];
  if (!NA_IsNArray(rb_a))
    rb_raise(rb_eArgError, "a (3th argument) must be NArray");
  if (NA_RANK(rb_a) != 2)
    rb_raise(rb_eArgError, "rank of a (3th argument) must be %d", 2);
  n = NA_SHAPE1(rb_a);
  lda = NA_SHAPE0(rb_a);
  if (NA_TYPE(rb_a) != NA_SCOMPLEX)
    rb_a = na_change_type(rb_a, NA_SCOMPLEX);
  a = NA_PTR_TYPE(rb_a, complex*);
  if (!NA_IsNArray(rb_b))
    rb_raise(rb_eArgError, "b (4th argument) must be NArray");
  if (NA_RANK(rb_b) != 2)
    rb_raise(rb_eArgError, "rank of b (4th argument) must be %d", 2);
  if (NA_SHAPE1(rb_b) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of b must be the same as shape 1 of a");
  ldb = NA_SHAPE0(rb_b);
  if (NA_TYPE(rb_b) != NA_SCOMPLEX)
    rb_b = na_change_type(rb_b, NA_SCOMPLEX);
  b = NA_PTR_TYPE(rb_b, complex*);
  lwork = NUM2INT(rb_lwork);
  jobvsr = StringValueCStr(rb_jobvsr)[0];
  ldvsr = lsame_(&jobvsr,"V") ? n : 1;
  ldvsl = lsame_(&jobvsl,"V") ? n : 1;
  {
    int shape[1];
    shape[0] = n;
    rb_alpha = na_make_object(NA_SCOMPLEX, 1, shape, cNArray);
  }
  alpha = NA_PTR_TYPE(rb_alpha, complex*);
  {
    int shape[1];
    shape[0] = n;
    rb_beta = na_make_object(NA_SCOMPLEX, 1, shape, cNArray);
  }
  beta = NA_PTR_TYPE(rb_beta, complex*);
  {
    int shape[2];
    shape[0] = ldvsl;
    shape[1] = n;
    rb_vsl = na_make_object(NA_SCOMPLEX, 2, shape, cNArray);
  }
  vsl = NA_PTR_TYPE(rb_vsl, complex*);
  {
    int shape[2];
    shape[0] = ldvsr;
    shape[1] = n;
    rb_vsr = na_make_object(NA_SCOMPLEX, 2, shape, cNArray);
  }
  vsr = NA_PTR_TYPE(rb_vsr, complex*);
  {
    int shape[1];
    shape[0] = MAX(1,lwork);
    rb_work = na_make_object(NA_SCOMPLEX, 1, shape, cNArray);
  }
  work = NA_PTR_TYPE(rb_work, complex*);
  {
    int shape[2];
    shape[0] = lda;
    shape[1] = n;
    rb_a_out__ = na_make_object(NA_SCOMPLEX, 2, shape, cNArray);
  }
  a_out__ = NA_PTR_TYPE(rb_a_out__, complex*);
  MEMCPY(a_out__, a, complex, NA_TOTAL(rb_a));
  rb_a = rb_a_out__;
  a = a_out__;
  {
    int shape[2];
    shape[0] = ldb;
    shape[1] = n;
    rb_b_out__ = na_make_object(NA_SCOMPLEX, 2, shape, cNArray);
  }
  b_out__ = NA_PTR_TYPE(rb_b_out__, complex*);
  MEMCPY(b_out__, b, complex, NA_TOTAL(rb_b));
  rb_b = rb_b_out__;
  b = b_out__;
  rwork = ALLOC_N(real, (3*n));

  cgegs_(&jobvsl, &jobvsr, &n, a, &lda, b, &ldb, alpha, beta, vsl, &ldvsl, vsr, &ldvsr, work, &lwork, rwork, &info);

  free(rwork);
  rb_info = INT2NUM(info);
  return rb_ary_new3(8, rb_alpha, rb_beta, rb_vsl, rb_vsr, rb_work, rb_info, rb_a, rb_b);
}

void
init_lapack_cgegs(VALUE mLapack){
  rb_define_module_function(mLapack, "cgegs", rb_cgegs, -1);
}

#include "rb_lapack.h"

static VALUE
rb_dgegs(int argc, VALUE *argv, VALUE self){
  VALUE rb_jobvsl;
  char jobvsl; 
  VALUE rb_jobvsr;
  char jobvsr; 
  VALUE rb_a;
  doublereal *a; 
  VALUE rb_b;
  doublereal *b; 
  VALUE rb_lwork;
  integer lwork; 
  VALUE rb_alphar;
  doublereal *alphar; 
  VALUE rb_alphai;
  doublereal *alphai; 
  VALUE rb_beta;
  doublereal *beta; 
  VALUE rb_vsl;
  doublereal *vsl; 
  VALUE rb_vsr;
  doublereal *vsr; 
  VALUE rb_work;
  doublereal *work; 
  VALUE rb_info;
  integer info; 
  VALUE rb_a_out__;
  doublereal *a_out__;
  VALUE rb_b_out__;
  doublereal *b_out__;

  integer lda;
  integer n;
  integer ldb;
  integer ldvsl;
  integer ldvsr;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  alphar, alphai, beta, vsl, vsr, work, info, a, b = NumRu::Lapack.dgegs( jobvsl, jobvsr, a, b, lwork)\n    or\n  NumRu::Lapack.dgegs  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE DGEGS( JOBVSL, JOBVSR, N, A, LDA, B, LDB, ALPHAR, ALPHAI, BETA, VSL, LDVSL, VSR, LDVSR, WORK, LWORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  This routine is deprecated and has been replaced by routine DGGES.\n*\n*  DGEGS computes the eigenvalues, real Schur form, and, optionally,\n*  left and or/right Schur vectors of a real matrix pair (A,B).\n*  Given two square matrices A and B, the generalized real Schur\n*  factorization has the form\n*\n*    A = Q*S*Z**T,  B = Q*T*Z**T\n*\n*  where Q and Z are orthogonal matrices, T is upper triangular, and S\n*  is an upper quasi-triangular matrix with 1-by-1 and 2-by-2 diagonal\n*  blocks, the 2-by-2 blocks corresponding to complex conjugate pairs\n*  of eigenvalues of (A,B).  The columns of Q are the left Schur vectors\n*  and the columns of Z are the right Schur vectors.\n*\n*  If only the eigenvalues of (A,B) are needed, the driver routine\n*  DGEGV should be used instead.  See DGEGV for a description of the\n*  eigenvalues of the generalized nonsymmetric eigenvalue problem\n*  (GNEP).\n*\n\n*  Arguments\n*  =========\n*\n*  JOBVSL  (input) CHARACTER*1\n*          = 'N':  do not compute the left Schur vectors;\n*          = 'V':  compute the left Schur vectors (returned in VSL).\n*\n*  JOBVSR  (input) CHARACTER*1\n*          = 'N':  do not compute the right Schur vectors;\n*          = 'V':  compute the right Schur vectors (returned in VSR).\n*\n*  N       (input) INTEGER\n*          The order of the matrices A, B, VSL, and VSR.  N >= 0.\n*\n*  A       (input/output) DOUBLE PRECISION array, dimension (LDA, N)\n*          On entry, the matrix A.\n*          On exit, the upper quasi-triangular matrix S from the\n*          generalized real Schur factorization.\n*\n*  LDA     (input) INTEGER\n*          The leading dimension of A.  LDA >= max(1,N).\n*\n*  B       (input/output) DOUBLE PRECISION array, dimension (LDB, N)\n*          On entry, the matrix B.\n*          On exit, the upper triangular matrix T from the generalized\n*          real Schur factorization.\n*\n*  LDB     (input) INTEGER\n*          The leading dimension of B.  LDB >= max(1,N).\n*\n*  ALPHAR  (output) DOUBLE PRECISION array, dimension (N)\n*          The real parts of each scalar alpha defining an eigenvalue\n*          of GNEP.\n*\n*  ALPHAI  (output) DOUBLE PRECISION array, dimension (N)\n*          The imaginary parts of each scalar alpha defining an\n*          eigenvalue of GNEP.  If ALPHAI(j) is zero, then the j-th\n*          eigenvalue is real; if positive, then the j-th and (j+1)-st\n*          eigenvalues are a complex conjugate pair, with\n*          ALPHAI(j+1) = -ALPHAI(j).\n*\n*  BETA    (output) DOUBLE PRECISION array, dimension (N)\n*          The scalars beta that define the eigenvalues of GNEP.\n*          Together, the quantities alpha = (ALPHAR(j),ALPHAI(j)) and\n*          beta = BETA(j) represent the j-th eigenvalue of the matrix\n*          pair (A,B), in one of the forms lambda = alpha/beta or\n*          mu = beta/alpha.  Since either lambda or mu may overflow,\n*          they should not, in general, be computed.\n*\n*  VSL     (output) DOUBLE PRECISION array, dimension (LDVSL,N)\n*          If JOBVSL = 'V', the matrix of left Schur vectors Q.\n*          Not referenced if JOBVSL = 'N'.\n*\n*  LDVSL   (input) INTEGER\n*          The leading dimension of the matrix VSL. LDVSL >=1, and\n*          if JOBVSL = 'V', LDVSL >= N.\n*\n*  VSR     (output) DOUBLE PRECISION array, dimension (LDVSR,N)\n*          If JOBVSR = 'V', the matrix of right Schur vectors Z.\n*          Not referenced if JOBVSR = 'N'.\n*\n*  LDVSR   (input) INTEGER\n*          The leading dimension of the matrix VSR. LDVSR >= 1, and\n*          if JOBVSR = 'V', LDVSR >= N.\n*\n*  WORK    (workspace/output) DOUBLE PRECISION array, dimension (MAX(1,LWORK))\n*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.\n*\n*  LWORK   (input) INTEGER\n*          The dimension of the array WORK.  LWORK >= max(1,4*N).\n*          For good performance, LWORK must generally be larger.\n*          To compute the optimal value of LWORK, call ILAENV to get\n*          blocksizes (for DGEQRF, DORMQR, and DORGQR.)  Then compute:\n*          NB  -- MAX of the blocksizes for DGEQRF, DORMQR, and DORGQR\n*          The optimal LWORK is  2*N + N*(NB+1).\n*\n*          If LWORK = -1, then a workspace query is assumed; the routine\n*          only calculates the optimal size of the WORK array, returns\n*          this value as the first entry of the WORK array, and no error\n*          message related to LWORK is issued by XERBLA.\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit\n*          < 0:  if INFO = -i, the i-th argument had an illegal value.\n*          = 1,...,N:\n*                The QZ iteration failed.  (A,B) are not in Schur\n*                form, but ALPHAR(j), ALPHAI(j), and BETA(j) should\n*                be correct for j=INFO+1,...,N.\n*          > N:  errors that usually indicate LAPACK problems:\n*                =N+1: error return from DGGBAL\n*                =N+2: error return from DGEQRF\n*                =N+3: error return from DORMQR\n*                =N+4: error return from DORGQR\n*                =N+5: error return from DGGHRD\n*                =N+6: error return from DHGEQZ (other than failed\n*                                                iteration)\n*                =N+7: error return from DGGBAK (computing VSL)\n*                =N+8: error return from DGGBAK (computing VSR)\n*                =N+9: error return from DLASCL (various places)\n*\n\n*  =====================================================================\n*\n\n");
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
  jobvsr = StringValueCStr(rb_jobvsr)[0];
  lwork = NUM2INT(rb_lwork);
  if (!NA_IsNArray(rb_a))
    rb_raise(rb_eArgError, "a (3th argument) must be NArray");
  if (NA_RANK(rb_a) != 2)
    rb_raise(rb_eArgError, "rank of a (3th argument) must be %d", 2);
  lda = NA_SHAPE0(rb_a);
  n = NA_SHAPE1(rb_a);
  if (NA_TYPE(rb_a) != NA_DFLOAT)
    rb_a = na_change_type(rb_a, NA_DFLOAT);
  a = NA_PTR_TYPE(rb_a, doublereal*);
  if (!NA_IsNArray(rb_b))
    rb_raise(rb_eArgError, "b (4th argument) must be NArray");
  if (NA_RANK(rb_b) != 2)
    rb_raise(rb_eArgError, "rank of b (4th argument) must be %d", 2);
  ldb = NA_SHAPE0(rb_b);
  if (NA_SHAPE1(rb_b) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of b must be the same as shape 1 of a");
  if (NA_TYPE(rb_b) != NA_DFLOAT)
    rb_b = na_change_type(rb_b, NA_DFLOAT);
  b = NA_PTR_TYPE(rb_b, doublereal*);
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
  ldvsl = lsame_(&jobvsl,"V") ? n : 1;
  {
    int shape[2];
    shape[0] = ldvsl;
    shape[1] = n;
    rb_vsl = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  vsl = NA_PTR_TYPE(rb_vsl, doublereal*);
  ldvsr = lsame_(&jobvsr,"V") ? n : 1;
  {
    int shape[2];
    shape[0] = ldvsr;
    shape[1] = n;
    rb_vsr = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  vsr = NA_PTR_TYPE(rb_vsr, doublereal*);
  {
    int shape[1];
    shape[0] = MAX(1,lwork);
    rb_work = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  work = NA_PTR_TYPE(rb_work, doublereal*);
  {
    int shape[2];
    shape[0] = lda;
    shape[1] = n;
    rb_a_out__ = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  a_out__ = NA_PTR_TYPE(rb_a_out__, doublereal*);
  MEMCPY(a_out__, a, doublereal, NA_TOTAL(rb_a));
  rb_a = rb_a_out__;
  a = a_out__;
  {
    int shape[2];
    shape[0] = ldb;
    shape[1] = n;
    rb_b_out__ = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  b_out__ = NA_PTR_TYPE(rb_b_out__, doublereal*);
  MEMCPY(b_out__, b, doublereal, NA_TOTAL(rb_b));
  rb_b = rb_b_out__;
  b = b_out__;

  dgegs_(&jobvsl, &jobvsr, &n, a, &lda, b, &ldb, alphar, alphai, beta, vsl, &ldvsl, vsr, &ldvsr, work, &lwork, &info);

  rb_info = INT2NUM(info);
  return rb_ary_new3(9, rb_alphar, rb_alphai, rb_beta, rb_vsl, rb_vsr, rb_work, rb_info, rb_a, rb_b);
}

void
init_lapack_dgegs(VALUE mLapack){
  rb_define_module_function(mLapack, "dgegs", rb_dgegs, -1);
}

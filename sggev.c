#include "rb_lapack.h"

static VALUE
rb_sggev(int argc, VALUE *argv, VALUE self){
  VALUE rb_jobvl;
  char jobvl; 
  VALUE rb_jobvr;
  char jobvr; 
  VALUE rb_a;
  real *a; 
  VALUE rb_b;
  real *b; 
  VALUE rb_lwork;
  integer lwork; 
  VALUE rb_alphar;
  real *alphar; 
  VALUE rb_alphai;
  real *alphai; 
  VALUE rb_beta;
  real *beta; 
  VALUE rb_vl;
  real *vl; 
  VALUE rb_vr;
  real *vr; 
  VALUE rb_work;
  real *work; 
  VALUE rb_info;
  integer info; 
  VALUE rb_a_out__;
  real *a_out__;
  VALUE rb_b_out__;
  real *b_out__;

  integer lda;
  integer n;
  integer ldb;
  integer ldvl;
  integer ldvr;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  alphar, alphai, beta, vl, vr, work, info, a, b = NumRu::Lapack.sggev( jobvl, jobvr, a, b, lwork)\n    or\n  NumRu::Lapack.sggev  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE SGGEV( JOBVL, JOBVR, N, A, LDA, B, LDB, ALPHAR, ALPHAI, BETA, VL, LDVL, VR, LDVR, WORK, LWORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  SGGEV computes for a pair of N-by-N real nonsymmetric matrices (A,B)\n*  the generalized eigenvalues, and optionally, the left and/or right\n*  generalized eigenvectors.\n*\n*  A generalized eigenvalue for a pair of matrices (A,B) is a scalar\n*  lambda or a ratio alpha/beta = lambda, such that A - lambda*B is\n*  singular. It is usually represented as the pair (alpha,beta), as\n*  there is a reasonable interpretation for beta=0, and even for both\n*  being zero.\n*\n*  The right eigenvector v(j) corresponding to the eigenvalue lambda(j)\n*  of (A,B) satisfies\n*\n*                   A * v(j) = lambda(j) * B * v(j).\n*\n*  The left eigenvector u(j) corresponding to the eigenvalue lambda(j)\n*  of (A,B) satisfies\n*\n*                   u(j)**H * A  = lambda(j) * u(j)**H * B .\n*\n*  where u(j)**H is the conjugate-transpose of u(j).\n*\n*\n\n*  Arguments\n*  =========\n*\n*  JOBVL   (input) CHARACTER*1\n*          = 'N':  do not compute the left generalized eigenvectors;\n*          = 'V':  compute the left generalized eigenvectors.\n*\n*  JOBVR   (input) CHARACTER*1\n*          = 'N':  do not compute the right generalized eigenvectors;\n*          = 'V':  compute the right generalized eigenvectors.\n*\n*  N       (input) INTEGER\n*          The order of the matrices A, B, VL, and VR.  N >= 0.\n*\n*  A       (input/output) REAL array, dimension (LDA, N)\n*          On entry, the matrix A in the pair (A,B).\n*          On exit, A has been overwritten.\n*\n*  LDA     (input) INTEGER\n*          The leading dimension of A.  LDA >= max(1,N).\n*\n*  B       (input/output) REAL array, dimension (LDB, N)\n*          On entry, the matrix B in the pair (A,B).\n*          On exit, B has been overwritten.\n*\n*  LDB     (input) INTEGER\n*          The leading dimension of B.  LDB >= max(1,N).\n*\n*  ALPHAR  (output) REAL array, dimension (N)\n*  ALPHAI  (output) REAL array, dimension (N)\n*  BETA    (output) REAL array, dimension (N)\n*          On exit, (ALPHAR(j) + ALPHAI(j)*i)/BETA(j), j=1,...,N, will\n*          be the generalized eigenvalues.  If ALPHAI(j) is zero, then\n*          the j-th eigenvalue is real; if positive, then the j-th and\n*          (j+1)-st eigenvalues are a complex conjugate pair, with\n*          ALPHAI(j+1) negative.\n*\n*          Note: the quotients ALPHAR(j)/BETA(j) and ALPHAI(j)/BETA(j)\n*          may easily over- or underflow, and BETA(j) may even be zero.\n*          Thus, the user should avoid naively computing the ratio\n*          alpha/beta.  However, ALPHAR and ALPHAI will be always less\n*          than and usually comparable with norm(A) in magnitude, and\n*          BETA always less than and usually comparable with norm(B).\n*\n*  VL      (output) REAL array, dimension (LDVL,N)\n*          If JOBVL = 'V', the left eigenvectors u(j) are stored one\n*          after another in the columns of VL, in the same order as\n*          their eigenvalues. If the j-th eigenvalue is real, then\n*          u(j) = VL(:,j), the j-th column of VL. If the j-th and\n*          (j+1)-th eigenvalues form a complex conjugate pair, then\n*          u(j) = VL(:,j)+i*VL(:,j+1) and u(j+1) = VL(:,j)-i*VL(:,j+1).\n*          Each eigenvector is scaled so the largest component has\n*          abs(real part)+abs(imag. part)=1.\n*          Not referenced if JOBVL = 'N'.\n*\n*  LDVL    (input) INTEGER\n*          The leading dimension of the matrix VL. LDVL >= 1, and\n*          if JOBVL = 'V', LDVL >= N.\n*\n*  VR      (output) REAL array, dimension (LDVR,N)\n*          If JOBVR = 'V', the right eigenvectors v(j) are stored one\n*          after another in the columns of VR, in the same order as\n*          their eigenvalues. If the j-th eigenvalue is real, then\n*          v(j) = VR(:,j), the j-th column of VR. If the j-th and\n*          (j+1)-th eigenvalues form a complex conjugate pair, then\n*          v(j) = VR(:,j)+i*VR(:,j+1) and v(j+1) = VR(:,j)-i*VR(:,j+1).\n*          Each eigenvector is scaled so the largest component has\n*          abs(real part)+abs(imag. part)=1.\n*          Not referenced if JOBVR = 'N'.\n*\n*  LDVR    (input) INTEGER\n*          The leading dimension of the matrix VR. LDVR >= 1, and\n*          if JOBVR = 'V', LDVR >= N.\n*\n*  WORK    (workspace/output) REAL array, dimension (MAX(1,LWORK))\n*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.\n*\n*  LWORK   (input) INTEGER\n*          The dimension of the array WORK.  LWORK >= max(1,8*N).\n*          For good performance, LWORK must generally be larger.\n*\n*          If LWORK = -1, then a workspace query is assumed; the routine\n*          only calculates the optimal size of the WORK array, returns\n*          this value as the first entry of the WORK array, and no error\n*          message related to LWORK is issued by XERBLA.\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit\n*          < 0:  if INFO = -i, the i-th argument had an illegal value.\n*          = 1,...,N:\n*                The QZ iteration failed.  No eigenvectors have been\n*                calculated, but ALPHAR(j), ALPHAI(j), and BETA(j)\n*                should be correct for j=INFO+1,...,N.\n*          > N:  =N+1: other than QZ iteration failed in SHGEQZ.\n*                =N+2: error return from STGEVC.\n*\n\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 5)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 5)", argc);
  rb_jobvl = argv[0];
  rb_jobvr = argv[1];
  rb_a = argv[2];
  rb_b = argv[3];
  rb_lwork = argv[4];

  jobvl = StringValueCStr(rb_jobvl)[0];
  jobvr = StringValueCStr(rb_jobvr)[0];
  lwork = NUM2INT(rb_lwork);
  if (!NA_IsNArray(rb_a))
    rb_raise(rb_eArgError, "a (3th argument) must be NArray");
  if (NA_RANK(rb_a) != 2)
    rb_raise(rb_eArgError, "rank of a (3th argument) must be %d", 2);
  lda = NA_SHAPE0(rb_a);
  n = NA_SHAPE1(rb_a);
  if (NA_TYPE(rb_a) != NA_SFLOAT)
    rb_a = na_change_type(rb_a, NA_SFLOAT);
  a = NA_PTR_TYPE(rb_a, real*);
  if (!NA_IsNArray(rb_b))
    rb_raise(rb_eArgError, "b (4th argument) must be NArray");
  if (NA_RANK(rb_b) != 2)
    rb_raise(rb_eArgError, "rank of b (4th argument) must be %d", 2);
  ldb = NA_SHAPE0(rb_b);
  if (NA_SHAPE1(rb_b) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of b must be the same as shape 1 of a");
  if (NA_TYPE(rb_b) != NA_SFLOAT)
    rb_b = na_change_type(rb_b, NA_SFLOAT);
  b = NA_PTR_TYPE(rb_b, real*);
  {
    int shape[1];
    shape[0] = n;
    rb_alphar = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  alphar = NA_PTR_TYPE(rb_alphar, real*);
  {
    int shape[1];
    shape[0] = n;
    rb_alphai = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  alphai = NA_PTR_TYPE(rb_alphai, real*);
  {
    int shape[1];
    shape[0] = n;
    rb_beta = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  beta = NA_PTR_TYPE(rb_beta, real*);
  ldvl = lsame_(&jobvl,"V") ? n : 1;
  {
    int shape[2];
    shape[0] = ldvl;
    shape[1] = n;
    rb_vl = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  vl = NA_PTR_TYPE(rb_vl, real*);
  ldvr = lsame_(&jobvr,"V") ? n : 1;
  {
    int shape[2];
    shape[0] = ldvr;
    shape[1] = n;
    rb_vr = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  vr = NA_PTR_TYPE(rb_vr, real*);
  {
    int shape[1];
    shape[0] = MAX(1,lwork);
    rb_work = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  work = NA_PTR_TYPE(rb_work, real*);
  {
    int shape[2];
    shape[0] = lda;
    shape[1] = n;
    rb_a_out__ = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  a_out__ = NA_PTR_TYPE(rb_a_out__, real*);
  MEMCPY(a_out__, a, real, NA_TOTAL(rb_a));
  rb_a = rb_a_out__;
  a = a_out__;
  {
    int shape[2];
    shape[0] = ldb;
    shape[1] = n;
    rb_b_out__ = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  b_out__ = NA_PTR_TYPE(rb_b_out__, real*);
  MEMCPY(b_out__, b, real, NA_TOTAL(rb_b));
  rb_b = rb_b_out__;
  b = b_out__;

  sggev_(&jobvl, &jobvr, &n, a, &lda, b, &ldb, alphar, alphai, beta, vl, &ldvl, vr, &ldvr, work, &lwork, &info);

  rb_info = INT2NUM(info);
  return rb_ary_new3(9, rb_alphar, rb_alphai, rb_beta, rb_vl, rb_vr, rb_work, rb_info, rb_a, rb_b);
}

void
init_lapack_sggev(VALUE mLapack){
  rb_define_module_function(mLapack, "sggev", rb_sggev, -1);
}

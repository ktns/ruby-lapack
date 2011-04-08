#include "rb_lapack.h"

extern VOID zggev_(char *jobvl, char *jobvr, integer *n, doublecomplex *a, integer *lda, doublecomplex *b, integer *ldb, doublecomplex *alpha, doublecomplex *beta, doublecomplex *vl, integer *ldvl, doublecomplex *vr, integer *ldvr, doublecomplex *work, integer *lwork, doublereal *rwork, integer *info);

static VALUE sHelp, sUsage;

static VALUE
rblapack_zggev(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_jobvl;
  char jobvl; 
  VALUE rblapack_jobvr;
  char jobvr; 
  VALUE rblapack_a;
  doublecomplex *a; 
  VALUE rblapack_b;
  doublecomplex *b; 
  VALUE rblapack_lwork;
  integer lwork; 
  VALUE rblapack_alpha;
  doublecomplex *alpha; 
  VALUE rblapack_beta;
  doublecomplex *beta; 
  VALUE rblapack_vl;
  doublecomplex *vl; 
  VALUE rblapack_vr;
  doublecomplex *vr; 
  VALUE rblapack_work;
  doublecomplex *work; 
  VALUE rblapack_rwork;
  doublereal *rwork; 
  VALUE rblapack_info;
  integer info; 
  VALUE rblapack_a_out__;
  doublecomplex *a_out__;
  VALUE rblapack_b_out__;
  doublecomplex *b_out__;

  integer lda;
  integer n;
  integer ldb;
  integer ldvl;
  integer ldvr;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  alpha, beta, vl, vr, work, rwork, info, a, b = NumRu::Lapack.zggev( jobvl, jobvr, a, b, lwork, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE ZGGEV( JOBVL, JOBVR, N, A, LDA, B, LDB, ALPHA, BETA, VL, LDVL, VR, LDVR, WORK, LWORK, RWORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  ZGGEV computes for a pair of N-by-N complex nonsymmetric matrices\n*  (A,B), the generalized eigenvalues, and optionally, the left and/or\n*  right generalized eigenvectors.\n*\n*  A generalized eigenvalue for a pair of matrices (A,B) is a scalar\n*  lambda or a ratio alpha/beta = lambda, such that A - lambda*B is\n*  singular. It is usually represented as the pair (alpha,beta), as\n*  there is a reasonable interpretation for beta=0, and even for both\n*  being zero.\n*\n*  The right generalized eigenvector v(j) corresponding to the\n*  generalized eigenvalue lambda(j) of (A,B) satisfies\n*\n*               A * v(j) = lambda(j) * B * v(j).\n*\n*  The left generalized eigenvector u(j) corresponding to the\n*  generalized eigenvalues lambda(j) of (A,B) satisfies\n*\n*               u(j)**H * A = lambda(j) * u(j)**H * B\n*\n*  where u(j)**H is the conjugate-transpose of u(j).\n*\n\n*  Arguments\n*  =========\n*\n*  JOBVL   (input) CHARACTER*1\n*          = 'N':  do not compute the left generalized eigenvectors;\n*          = 'V':  compute the left generalized eigenvectors.\n*\n*  JOBVR   (input) CHARACTER*1\n*          = 'N':  do not compute the right generalized eigenvectors;\n*          = 'V':  compute the right generalized eigenvectors.\n*\n*  N       (input) INTEGER\n*          The order of the matrices A, B, VL, and VR.  N >= 0.\n*\n*  A       (input/output) COMPLEX*16 array, dimension (LDA, N)\n*          On entry, the matrix A in the pair (A,B).\n*          On exit, A has been overwritten.\n*\n*  LDA     (input) INTEGER\n*          The leading dimension of A.  LDA >= max(1,N).\n*\n*  B       (input/output) COMPLEX*16 array, dimension (LDB, N)\n*          On entry, the matrix B in the pair (A,B).\n*          On exit, B has been overwritten.\n*\n*  LDB     (input) INTEGER\n*          The leading dimension of B.  LDB >= max(1,N).\n*\n*  ALPHA   (output) COMPLEX*16 array, dimension (N)\n*  BETA    (output) COMPLEX*16 array, dimension (N)\n*          On exit, ALPHA(j)/BETA(j), j=1,...,N, will be the\n*          generalized eigenvalues.\n*\n*          Note: the quotients ALPHA(j)/BETA(j) may easily over- or\n*          underflow, and BETA(j) may even be zero.  Thus, the user\n*          should avoid naively computing the ratio alpha/beta.\n*          However, ALPHA will be always less than and usually\n*          comparable with norm(A) in magnitude, and BETA always less\n*          than and usually comparable with norm(B).\n*\n*  VL      (output) COMPLEX*16 array, dimension (LDVL,N)\n*          If JOBVL = 'V', the left generalized eigenvectors u(j) are\n*          stored one after another in the columns of VL, in the same\n*          order as their eigenvalues.\n*          Each eigenvector is scaled so the largest component has\n*          abs(real part) + abs(imag. part) = 1.\n*          Not referenced if JOBVL = 'N'.\n*\n*  LDVL    (input) INTEGER\n*          The leading dimension of the matrix VL. LDVL >= 1, and\n*          if JOBVL = 'V', LDVL >= N.\n*\n*  VR      (output) COMPLEX*16 array, dimension (LDVR,N)\n*          If JOBVR = 'V', the right generalized eigenvectors v(j) are\n*          stored one after another in the columns of VR, in the same\n*          order as their eigenvalues.\n*          Each eigenvector is scaled so the largest component has\n*          abs(real part) + abs(imag. part) = 1.\n*          Not referenced if JOBVR = 'N'.\n*\n*  LDVR    (input) INTEGER\n*          The leading dimension of the matrix VR. LDVR >= 1, and\n*          if JOBVR = 'V', LDVR >= N.\n*\n*  WORK    (workspace/output) COMPLEX*16 array, dimension (MAX(1,LWORK))\n*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.\n*\n*  LWORK   (input) INTEGER\n*          The dimension of the array WORK.  LWORK >= max(1,2*N).\n*          For good performance, LWORK must generally be larger.\n*\n*          If LWORK = -1, then a workspace query is assumed; the routine\n*          only calculates the optimal size of the WORK array, returns\n*          this value as the first entry of the WORK array, and no error\n*          message related to LWORK is issued by XERBLA.\n*\n*  RWORK   (workspace/output) DOUBLE PRECISION array, dimension (8*N)\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit\n*          < 0:  if INFO = -i, the i-th argument had an illegal value.\n*          =1,...,N:\n*                The QZ iteration failed.  No eigenvectors have been\n*                calculated, but ALPHA(j) and BETA(j) should be\n*                correct for j=INFO+1,...,N.\n*          > N:  =N+1: other then QZ iteration failed in DHGEQZ,\n*                =N+2: error return from DTGEVC.\n*\n\n*  =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  alpha, beta, vl, vr, work, rwork, info, a, b = NumRu::Lapack.zggev( jobvl, jobvr, a, b, lwork, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 5)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 5)", argc);
  rblapack_jobvl = argv[0];
  rblapack_jobvr = argv[1];
  rblapack_a = argv[2];
  rblapack_b = argv[3];
  rblapack_lwork = argv[4];
  if (rb_options != Qnil) {
  }

  if (!NA_IsNArray(rblapack_a))
    rb_raise(rb_eArgError, "a (3th argument) must be NArray");
  if (NA_RANK(rblapack_a) != 2)
    rb_raise(rb_eArgError, "rank of a (3th argument) must be %d", 2);
  n = NA_SHAPE1(rblapack_a);
  lda = NA_SHAPE0(rblapack_a);
  if (NA_TYPE(rblapack_a) != NA_DCOMPLEX)
    rblapack_a = na_change_type(rblapack_a, NA_DCOMPLEX);
  a = NA_PTR_TYPE(rblapack_a, doublecomplex*);
  jobvr = StringValueCStr(rblapack_jobvr)[0];
  jobvl = StringValueCStr(rblapack_jobvl)[0];
  lwork = NUM2INT(rblapack_lwork);
  if (!NA_IsNArray(rblapack_b))
    rb_raise(rb_eArgError, "b (4th argument) must be NArray");
  if (NA_RANK(rblapack_b) != 2)
    rb_raise(rb_eArgError, "rank of b (4th argument) must be %d", 2);
  if (NA_SHAPE1(rblapack_b) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of b must be the same as shape 1 of a");
  ldb = NA_SHAPE0(rblapack_b);
  if (NA_TYPE(rblapack_b) != NA_DCOMPLEX)
    rblapack_b = na_change_type(rblapack_b, NA_DCOMPLEX);
  b = NA_PTR_TYPE(rblapack_b, doublecomplex*);
  ldvr = lsame_(&jobvr,"V") ? n : 1;
  ldvl = lsame_(&jobvl,"V") ? n : 1;
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
    int shape[2];
    shape[0] = ldvl;
    shape[1] = n;
    rblapack_vl = na_make_object(NA_DCOMPLEX, 2, shape, cNArray);
  }
  vl = NA_PTR_TYPE(rblapack_vl, doublecomplex*);
  {
    int shape[2];
    shape[0] = ldvr;
    shape[1] = n;
    rblapack_vr = na_make_object(NA_DCOMPLEX, 2, shape, cNArray);
  }
  vr = NA_PTR_TYPE(rblapack_vr, doublecomplex*);
  {
    int shape[1];
    shape[0] = MAX(1,lwork);
    rblapack_work = na_make_object(NA_DCOMPLEX, 1, shape, cNArray);
  }
  work = NA_PTR_TYPE(rblapack_work, doublecomplex*);
  {
    int shape[1];
    shape[0] = 8*n;
    rblapack_rwork = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  rwork = NA_PTR_TYPE(rblapack_rwork, doublereal*);
  {
    int shape[2];
    shape[0] = lda;
    shape[1] = n;
    rblapack_a_out__ = na_make_object(NA_DCOMPLEX, 2, shape, cNArray);
  }
  a_out__ = NA_PTR_TYPE(rblapack_a_out__, doublecomplex*);
  MEMCPY(a_out__, a, doublecomplex, NA_TOTAL(rblapack_a));
  rblapack_a = rblapack_a_out__;
  a = a_out__;
  {
    int shape[2];
    shape[0] = ldb;
    shape[1] = n;
    rblapack_b_out__ = na_make_object(NA_DCOMPLEX, 2, shape, cNArray);
  }
  b_out__ = NA_PTR_TYPE(rblapack_b_out__, doublecomplex*);
  MEMCPY(b_out__, b, doublecomplex, NA_TOTAL(rblapack_b));
  rblapack_b = rblapack_b_out__;
  b = b_out__;

  zggev_(&jobvl, &jobvr, &n, a, &lda, b, &ldb, alpha, beta, vl, &ldvl, vr, &ldvr, work, &lwork, rwork, &info);

  rblapack_info = INT2NUM(info);
  return rb_ary_new3(9, rblapack_alpha, rblapack_beta, rblapack_vl, rblapack_vr, rblapack_work, rblapack_rwork, rblapack_info, rblapack_a, rblapack_b);
}

void
init_lapack_zggev(VALUE mLapack){
  rb_define_module_function(mLapack, "zggev", rblapack_zggev, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}

#include "rb_lapack.h"

extern VOID cgegv_(char *jobvl, char *jobvr, integer *n, complex *a, integer *lda, complex *b, integer *ldb, complex *alpha, complex *beta, complex *vl, integer *ldvl, complex *vr, integer *ldvr, complex *work, integer *lwork, real *rwork, integer *info);

static VALUE sHelp, sUsage;

static VALUE
rblapack_cgegv(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_jobvl;
  char jobvl; 
  VALUE rblapack_jobvr;
  char jobvr; 
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
  VALUE rblapack_vl;
  complex *vl; 
  VALUE rblapack_vr;
  complex *vr; 
  VALUE rblapack_work;
  complex *work; 
  VALUE rblapack_rwork;
  real *rwork; 
  VALUE rblapack_info;
  integer info; 
  VALUE rblapack_a_out__;
  complex *a_out__;
  VALUE rblapack_b_out__;
  complex *b_out__;

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
      printf("%s\n", "USAGE:\n  alpha, beta, vl, vr, work, rwork, info, a, b = NumRu::Lapack.cgegv( jobvl, jobvr, a, b, lwork, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE CGEGV( JOBVL, JOBVR, N, A, LDA, B, LDB, ALPHA, BETA, VL, LDVL, VR, LDVR, WORK, LWORK, RWORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  This routine is deprecated and has been replaced by routine CGGEV.\n*\n*  CGEGV computes the eigenvalues and, optionally, the left and/or right\n*  eigenvectors of a complex matrix pair (A,B).\n*  Given two square matrices A and B,\n*  the generalized nonsymmetric eigenvalue problem (GNEP) is to find the\n*  eigenvalues lambda and corresponding (non-zero) eigenvectors x such\n*  that\n*     A*x = lambda*B*x.\n*\n*  An alternate form is to find the eigenvalues mu and corresponding\n*  eigenvectors y such that\n*     mu*A*y = B*y.\n*\n*  These two forms are equivalent with mu = 1/lambda and x = y if\n*  neither lambda nor mu is zero.  In order to deal with the case that\n*  lambda or mu is zero or small, two values alpha and beta are returned\n*  for each eigenvalue, such that lambda = alpha/beta and\n*  mu = beta/alpha.\n*  \n*  The vectors x and y in the above equations are right eigenvectors of\n*  the matrix pair (A,B).  Vectors u and v satisfying\n*     u**H*A = lambda*u**H*B  or  mu*v**H*A = v**H*B\n*  are left eigenvectors of (A,B).\n*\n*  Note: this routine performs \"full balancing\" on A and B -- see\n*  \"Further Details\", below.\n*\n\n*  Arguments\n*  =========\n*\n*  JOBVL   (input) CHARACTER*1\n*          = 'N':  do not compute the left generalized eigenvectors;\n*          = 'V':  compute the left generalized eigenvectors (returned\n*                  in VL).\n*\n*  JOBVR   (input) CHARACTER*1\n*          = 'N':  do not compute the right generalized eigenvectors;\n*          = 'V':  compute the right generalized eigenvectors (returned\n*                  in VR).\n*\n*  N       (input) INTEGER\n*          The order of the matrices A, B, VL, and VR.  N >= 0.\n*\n*  A       (input/output) COMPLEX array, dimension (LDA, N)\n*          On entry, the matrix A.\n*          If JOBVL = 'V' or JOBVR = 'V', then on exit A\n*          contains the Schur form of A from the generalized Schur\n*          factorization of the pair (A,B) after balancing.  If no\n*          eigenvectors were computed, then only the diagonal elements\n*          of the Schur form will be correct.  See CGGHRD and CHGEQZ\n*          for details.\n*\n*  LDA     (input) INTEGER\n*          The leading dimension of A.  LDA >= max(1,N).\n*\n*  B       (input/output) COMPLEX array, dimension (LDB, N)\n*          On entry, the matrix B.\n*          If JOBVL = 'V' or JOBVR = 'V', then on exit B contains the\n*          upper triangular matrix obtained from B in the generalized\n*          Schur factorization of the pair (A,B) after balancing.\n*          If no eigenvectors were computed, then only the diagonal\n*          elements of B will be correct.  See CGGHRD and CHGEQZ for\n*          details.\n*\n*  LDB     (input) INTEGER\n*          The leading dimension of B.  LDB >= max(1,N).\n*\n*  ALPHA   (output) COMPLEX array, dimension (N)\n*          The complex scalars alpha that define the eigenvalues of\n*          GNEP.\n*\n*  BETA    (output) COMPLEX array, dimension (N)\n*          The complex scalars beta that define the eigenvalues of GNEP.\n*          \n*          Together, the quantities alpha = ALPHA(j) and beta = BETA(j)\n*          represent the j-th eigenvalue of the matrix pair (A,B), in\n*          one of the forms lambda = alpha/beta or mu = beta/alpha.\n*          Since either lambda or mu may overflow, they should not,\n*          in general, be computed.\n\n*\n*  VL      (output) COMPLEX array, dimension (LDVL,N)\n*          If JOBVL = 'V', the left eigenvectors u(j) are stored\n*          in the columns of VL, in the same order as their eigenvalues.\n*          Each eigenvector is scaled so that its largest component has\n*          abs(real part) + abs(imag. part) = 1, except for eigenvectors\n*          corresponding to an eigenvalue with alpha = beta = 0, which\n*          are set to zero.\n*          Not referenced if JOBVL = 'N'.\n*\n*  LDVL    (input) INTEGER\n*          The leading dimension of the matrix VL. LDVL >= 1, and\n*          if JOBVL = 'V', LDVL >= N.\n*\n*  VR      (output) COMPLEX array, dimension (LDVR,N)\n*          If JOBVR = 'V', the right eigenvectors x(j) are stored\n*          in the columns of VR, in the same order as their eigenvalues.\n*          Each eigenvector is scaled so that its largest component has\n*          abs(real part) + abs(imag. part) = 1, except for eigenvectors\n*          corresponding to an eigenvalue with alpha = beta = 0, which\n*          are set to zero.\n*          Not referenced if JOBVR = 'N'.\n*\n*  LDVR    (input) INTEGER\n*          The leading dimension of the matrix VR. LDVR >= 1, and\n*          if JOBVR = 'V', LDVR >= N.\n*\n*  WORK    (workspace/output) COMPLEX array, dimension (MAX(1,LWORK))\n*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.\n*\n*  LWORK   (input) INTEGER\n*          The dimension of the array WORK.  LWORK >= max(1,2*N).\n*          For good performance, LWORK must generally be larger.\n*          To compute the optimal value of LWORK, call ILAENV to get\n*          blocksizes (for CGEQRF, CUNMQR, and CUNGQR.)  Then compute:\n*          NB  -- MAX of the blocksizes for CGEQRF, CUNMQR, and CUNGQR;\n*          The optimal LWORK is  MAX( 2*N, N*(NB+1) ).\n*\n*          If LWORK = -1, then a workspace query is assumed; the routine\n*          only calculates the optimal size of the WORK array, returns\n*          this value as the first entry of the WORK array, and no error\n*          message related to LWORK is issued by XERBLA.\n*\n*  RWORK   (workspace/output) REAL array, dimension (8*N)\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit\n*          < 0:  if INFO = -i, the i-th argument had an illegal value.\n*          =1,...,N:\n*                The QZ iteration failed.  No eigenvectors have been\n*                calculated, but ALPHA(j) and BETA(j) should be\n*                correct for j=INFO+1,...,N.\n*          > N:  errors that usually indicate LAPACK problems:\n*                =N+1: error return from CGGBAL\n*                =N+2: error return from CGEQRF\n*                =N+3: error return from CUNMQR\n*                =N+4: error return from CUNGQR\n*                =N+5: error return from CGGHRD\n*                =N+6: error return from CHGEQZ (other than failed\n*                                               iteration)\n*                =N+7: error return from CTGEVC\n*                =N+8: error return from CGGBAK (computing VL)\n*                =N+9: error return from CGGBAK (computing VR)\n*                =N+10: error return from CLASCL (various calls)\n*\n\n*  Further Details\n*  ===============\n*\n*  Balancing\n*  ---------\n*\n*  This driver calls CGGBAL to both permute and scale rows and columns\n*  of A and B.  The permutations PL and PR are chosen so that PL*A*PR\n*  and PL*B*R will be upper triangular except for the diagonal blocks\n*  A(i:j,i:j) and B(i:j,i:j), with i and j as close together as\n*  possible.  The diagonal scaling matrices DL and DR are chosen so\n*  that the pair  DL*PL*A*PR*DR, DL*PL*B*PR*DR have elements close to\n*  one (except for the elements that start out zero.)\n*\n*  After the eigenvalues and eigenvectors of the balanced matrices\n*  have been computed, CGGBAK transforms the eigenvectors back to what\n*  they would have been (in perfect arithmetic) if they had not been\n*  balanced.\n*\n*  Contents of A and B on Exit\n*  -------- -- - --- - -- ----\n*\n*  If any eigenvectors are computed (either JOBVL='V' or JOBVR='V' or\n*  both), then on exit the arrays A and B will contain the complex Schur\n*  form[*] of the \"balanced\" versions of A and B.  If no eigenvectors\n*  are computed, then only the diagonal blocks will be correct.\n*\n*  [*] In other words, upper triangular form.\n*\n*  =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  alpha, beta, vl, vr, work, rwork, info, a, b = NumRu::Lapack.cgegv( jobvl, jobvr, a, b, lwork, [:usage => usage, :help => help])\n");
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
  if (NA_TYPE(rblapack_a) != NA_SCOMPLEX)
    rblapack_a = na_change_type(rblapack_a, NA_SCOMPLEX);
  a = NA_PTR_TYPE(rblapack_a, complex*);
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
  if (NA_TYPE(rblapack_b) != NA_SCOMPLEX)
    rblapack_b = na_change_type(rblapack_b, NA_SCOMPLEX);
  b = NA_PTR_TYPE(rblapack_b, complex*);
  ldvr = lsame_(&jobvr,"V") ? n : 1;
  ldvl = lsame_(&jobvl,"V") ? n : 1;
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
    shape[0] = ldvl;
    shape[1] = n;
    rblapack_vl = na_make_object(NA_SCOMPLEX, 2, shape, cNArray);
  }
  vl = NA_PTR_TYPE(rblapack_vl, complex*);
  {
    int shape[2];
    shape[0] = ldvr;
    shape[1] = n;
    rblapack_vr = na_make_object(NA_SCOMPLEX, 2, shape, cNArray);
  }
  vr = NA_PTR_TYPE(rblapack_vr, complex*);
  {
    int shape[1];
    shape[0] = MAX(1,lwork);
    rblapack_work = na_make_object(NA_SCOMPLEX, 1, shape, cNArray);
  }
  work = NA_PTR_TYPE(rblapack_work, complex*);
  {
    int shape[1];
    shape[0] = 8*n;
    rblapack_rwork = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  rwork = NA_PTR_TYPE(rblapack_rwork, real*);
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

  cgegv_(&jobvl, &jobvr, &n, a, &lda, b, &ldb, alpha, beta, vl, &ldvl, vr, &ldvr, work, &lwork, rwork, &info);

  rblapack_info = INT2NUM(info);
  return rb_ary_new3(9, rblapack_alpha, rblapack_beta, rblapack_vl, rblapack_vr, rblapack_work, rblapack_rwork, rblapack_info, rblapack_a, rblapack_b);
}

void
init_lapack_cgegv(VALUE mLapack){
  rb_define_module_function(mLapack, "cgegv", rblapack_cgegv, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}

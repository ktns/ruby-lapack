#include "rb_lapack.h"

extern VOID dgegv_(char *jobvl, char *jobvr, integer *n, doublereal *a, integer *lda, doublereal *b, integer *ldb, doublereal *alphar, doublereal *alphai, doublereal *beta, doublereal *vl, integer *ldvl, doublereal *vr, integer *ldvr, doublereal *work, integer *lwork, integer *info);

static VALUE sHelp, sUsage;

static VALUE
rblapack_dgegv(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_jobvl;
  char jobvl; 
  VALUE rblapack_jobvr;
  char jobvr; 
  VALUE rblapack_a;
  doublereal *a; 
  VALUE rblapack_b;
  doublereal *b; 
  VALUE rblapack_lwork;
  integer lwork; 
  VALUE rblapack_alphar;
  doublereal *alphar; 
  VALUE rblapack_alphai;
  doublereal *alphai; 
  VALUE rblapack_beta;
  doublereal *beta; 
  VALUE rblapack_vl;
  doublereal *vl; 
  VALUE rblapack_vr;
  doublereal *vr; 
  VALUE rblapack_work;
  doublereal *work; 
  VALUE rblapack_info;
  integer info; 
  VALUE rblapack_a_out__;
  doublereal *a_out__;
  VALUE rblapack_b_out__;
  doublereal *b_out__;

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
      printf("%s\n", "USAGE:\n  alphar, alphai, beta, vl, vr, work, info, a, b = NumRu::Lapack.dgegv( jobvl, jobvr, a, b, lwork, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE DGEGV( JOBVL, JOBVR, N, A, LDA, B, LDB, ALPHAR, ALPHAI, BETA, VL, LDVL, VR, LDVR, WORK, LWORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  This routine is deprecated and has been replaced by routine DGGEV.\n*\n*  DGEGV computes the eigenvalues and, optionally, the left and/or right\n*  eigenvectors of a real matrix pair (A,B).\n*  Given two square matrices A and B,\n*  the generalized nonsymmetric eigenvalue problem (GNEP) is to find the\n*  eigenvalues lambda and corresponding (non-zero) eigenvectors x such\n*  that\n*\n*     A*x = lambda*B*x.\n*\n*  An alternate form is to find the eigenvalues mu and corresponding\n*  eigenvectors y such that\n*\n*     mu*A*y = B*y.\n*\n*  These two forms are equivalent with mu = 1/lambda and x = y if\n*  neither lambda nor mu is zero.  In order to deal with the case that\n*  lambda or mu is zero or small, two values alpha and beta are returned\n*  for each eigenvalue, such that lambda = alpha/beta and\n*  mu = beta/alpha.\n*\n*  The vectors x and y in the above equations are right eigenvectors of\n*  the matrix pair (A,B).  Vectors u and v satisfying\n*\n*     u**H*A = lambda*u**H*B  or  mu*v**H*A = v**H*B\n*\n*  are left eigenvectors of (A,B).\n*\n*  Note: this routine performs \"full balancing\" on A and B -- see\n*  \"Further Details\", below.\n*\n\n*  Arguments\n*  =========\n*\n*  JOBVL   (input) CHARACTER*1\n*          = 'N':  do not compute the left generalized eigenvectors;\n*          = 'V':  compute the left generalized eigenvectors (returned\n*                  in VL).\n*\n*  JOBVR   (input) CHARACTER*1\n*          = 'N':  do not compute the right generalized eigenvectors;\n*          = 'V':  compute the right generalized eigenvectors (returned\n*                  in VR).\n*\n*  N       (input) INTEGER\n*          The order of the matrices A, B, VL, and VR.  N >= 0.\n*\n*  A       (input/output) DOUBLE PRECISION array, dimension (LDA, N)\n*          On entry, the matrix A.\n*          If JOBVL = 'V' or JOBVR = 'V', then on exit A\n*          contains the real Schur form of A from the generalized Schur\n*          factorization of the pair (A,B) after balancing.\n*          If no eigenvectors were computed, then only the diagonal\n*          blocks from the Schur form will be correct.  See DGGHRD and\n*          DHGEQZ for details.\n*\n*  LDA     (input) INTEGER\n*          The leading dimension of A.  LDA >= max(1,N).\n*\n*  B       (input/output) DOUBLE PRECISION array, dimension (LDB, N)\n*          On entry, the matrix B.\n*          If JOBVL = 'V' or JOBVR = 'V', then on exit B contains the\n*          upper triangular matrix obtained from B in the generalized\n*          Schur factorization of the pair (A,B) after balancing.\n*          If no eigenvectors were computed, then only those elements of\n*          B corresponding to the diagonal blocks from the Schur form of\n*          A will be correct.  See DGGHRD and DHGEQZ for details.\n*\n*  LDB     (input) INTEGER\n*          The leading dimension of B.  LDB >= max(1,N).\n*\n*  ALPHAR  (output) DOUBLE PRECISION array, dimension (N)\n*          The real parts of each scalar alpha defining an eigenvalue of\n*          GNEP.\n*\n*  ALPHAI  (output) DOUBLE PRECISION array, dimension (N)\n*          The imaginary parts of each scalar alpha defining an\n*          eigenvalue of GNEP.  If ALPHAI(j) is zero, then the j-th\n*          eigenvalue is real; if positive, then the j-th and\n*          (j+1)-st eigenvalues are a complex conjugate pair, with\n*          ALPHAI(j+1) = -ALPHAI(j).\n*\n*  BETA    (output) DOUBLE PRECISION array, dimension (N)\n*          The scalars beta that define the eigenvalues of GNEP.\n*          \n*          Together, the quantities alpha = (ALPHAR(j),ALPHAI(j)) and\n*          beta = BETA(j) represent the j-th eigenvalue of the matrix\n*          pair (A,B), in one of the forms lambda = alpha/beta or\n*          mu = beta/alpha.  Since either lambda or mu may overflow,\n*          they should not, in general, be computed.\n*\n*  VL      (output) DOUBLE PRECISION array, dimension (LDVL,N)\n*          If JOBVL = 'V', the left eigenvectors u(j) are stored\n*          in the columns of VL, in the same order as their eigenvalues.\n*          If the j-th eigenvalue is real, then u(j) = VL(:,j).\n*          If the j-th and (j+1)-st eigenvalues form a complex conjugate\n*          pair, then\n*             u(j) = VL(:,j) + i*VL(:,j+1)\n*          and\n*            u(j+1) = VL(:,j) - i*VL(:,j+1).\n*\n*          Each eigenvector is scaled so that its largest component has\n*          abs(real part) + abs(imag. part) = 1, except for eigenvectors\n*          corresponding to an eigenvalue with alpha = beta = 0, which\n*          are set to zero.\n*          Not referenced if JOBVL = 'N'.\n*\n*  LDVL    (input) INTEGER\n*          The leading dimension of the matrix VL. LDVL >= 1, and\n*          if JOBVL = 'V', LDVL >= N.\n*\n*  VR      (output) DOUBLE PRECISION array, dimension (LDVR,N)\n*          If JOBVR = 'V', the right eigenvectors x(j) are stored\n*          in the columns of VR, in the same order as their eigenvalues.\n*          If the j-th eigenvalue is real, then x(j) = VR(:,j).\n*          If the j-th and (j+1)-st eigenvalues form a complex conjugate\n*          pair, then\n*            x(j) = VR(:,j) + i*VR(:,j+1)\n*          and\n*            x(j+1) = VR(:,j) - i*VR(:,j+1).\n*\n*          Each eigenvector is scaled so that its largest component has\n*          abs(real part) + abs(imag. part) = 1, except for eigenvalues\n*          corresponding to an eigenvalue with alpha = beta = 0, which\n*          are set to zero.\n*          Not referenced if JOBVR = 'N'.\n*\n*  LDVR    (input) INTEGER\n*          The leading dimension of the matrix VR. LDVR >= 1, and\n*          if JOBVR = 'V', LDVR >= N.\n*\n*  WORK    (workspace/output) DOUBLE PRECISION array, dimension (MAX(1,LWORK))\n*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.\n*\n*  LWORK   (input) INTEGER\n*          The dimension of the array WORK.  LWORK >= max(1,8*N).\n*          For good performance, LWORK must generally be larger.\n*          To compute the optimal value of LWORK, call ILAENV to get\n*          blocksizes (for DGEQRF, DORMQR, and DORGQR.)  Then compute:\n*          NB  -- MAX of the blocksizes for DGEQRF, DORMQR, and DORGQR;\n*          The optimal LWORK is:\n*              2*N + MAX( 6*N, N*(NB+1) ).\n*\n*          If LWORK = -1, then a workspace query is assumed; the routine\n*          only calculates the optimal size of the WORK array, returns\n*          this value as the first entry of the WORK array, and no error\n*          message related to LWORK is issued by XERBLA.\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit\n*          < 0:  if INFO = -i, the i-th argument had an illegal value.\n*          = 1,...,N:\n*                The QZ iteration failed.  No eigenvectors have been\n*                calculated, but ALPHAR(j), ALPHAI(j), and BETA(j)\n*                should be correct for j=INFO+1,...,N.\n*          > N:  errors that usually indicate LAPACK problems:\n*                =N+1: error return from DGGBAL\n*                =N+2: error return from DGEQRF\n*                =N+3: error return from DORMQR\n*                =N+4: error return from DORGQR\n*                =N+5: error return from DGGHRD\n*                =N+6: error return from DHGEQZ (other than failed\n*                                                iteration)\n*                =N+7: error return from DTGEVC\n*                =N+8: error return from DGGBAK (computing VL)\n*                =N+9: error return from DGGBAK (computing VR)\n*                =N+10: error return from DLASCL (various calls)\n*\n\n*  Further Details\n*  ===============\n*\n*  Balancing\n*  ---------\n*\n*  This driver calls DGGBAL to both permute and scale rows and columns\n*  of A and B.  The permutations PL and PR are chosen so that PL*A*PR\n*  and PL*B*R will be upper triangular except for the diagonal blocks\n*  A(i:j,i:j) and B(i:j,i:j), with i and j as close together as\n*  possible.  The diagonal scaling matrices DL and DR are chosen so\n*  that the pair  DL*PL*A*PR*DR, DL*PL*B*PR*DR have elements close to\n*  one (except for the elements that start out zero.)\n*\n*  After the eigenvalues and eigenvectors of the balanced matrices\n*  have been computed, DGGBAK transforms the eigenvectors back to what\n*  they would have been (in perfect arithmetic) if they had not been\n*  balanced.\n*\n*  Contents of A and B on Exit\n*  -------- -- - --- - -- ----\n*\n*  If any eigenvectors are computed (either JOBVL='V' or JOBVR='V' or\n*  both), then on exit the arrays A and B will contain the real Schur\n*  form[*] of the \"balanced\" versions of A and B.  If no eigenvectors\n*  are computed, then only the diagonal blocks will be correct.\n*\n*  [*] See DHGEQZ, DGEGS, or read the book \"Matrix Computations\",\n*      by Golub & van Loan, pub. by Johns Hopkins U. Press.\n*\n*  =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  alphar, alphai, beta, vl, vr, work, info, a, b = NumRu::Lapack.dgegv( jobvl, jobvr, a, b, lwork, [:usage => usage, :help => help])\n");
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
  if (NA_TYPE(rblapack_a) != NA_DFLOAT)
    rblapack_a = na_change_type(rblapack_a, NA_DFLOAT);
  a = NA_PTR_TYPE(rblapack_a, doublereal*);
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
  if (NA_TYPE(rblapack_b) != NA_DFLOAT)
    rblapack_b = na_change_type(rblapack_b, NA_DFLOAT);
  b = NA_PTR_TYPE(rblapack_b, doublereal*);
  ldvr = lsame_(&jobvr,"V") ? n : 1;
  ldvl = lsame_(&jobvl,"V") ? n : 1;
  {
    int shape[1];
    shape[0] = n;
    rblapack_alphar = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  alphar = NA_PTR_TYPE(rblapack_alphar, doublereal*);
  {
    int shape[1];
    shape[0] = n;
    rblapack_alphai = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  alphai = NA_PTR_TYPE(rblapack_alphai, doublereal*);
  {
    int shape[1];
    shape[0] = n;
    rblapack_beta = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  beta = NA_PTR_TYPE(rblapack_beta, doublereal*);
  {
    int shape[2];
    shape[0] = ldvl;
    shape[1] = n;
    rblapack_vl = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  vl = NA_PTR_TYPE(rblapack_vl, doublereal*);
  {
    int shape[2];
    shape[0] = ldvr;
    shape[1] = n;
    rblapack_vr = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  vr = NA_PTR_TYPE(rblapack_vr, doublereal*);
  {
    int shape[1];
    shape[0] = MAX(1,lwork);
    rblapack_work = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  work = NA_PTR_TYPE(rblapack_work, doublereal*);
  {
    int shape[2];
    shape[0] = lda;
    shape[1] = n;
    rblapack_a_out__ = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  a_out__ = NA_PTR_TYPE(rblapack_a_out__, doublereal*);
  MEMCPY(a_out__, a, doublereal, NA_TOTAL(rblapack_a));
  rblapack_a = rblapack_a_out__;
  a = a_out__;
  {
    int shape[2];
    shape[0] = ldb;
    shape[1] = n;
    rblapack_b_out__ = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  b_out__ = NA_PTR_TYPE(rblapack_b_out__, doublereal*);
  MEMCPY(b_out__, b, doublereal, NA_TOTAL(rblapack_b));
  rblapack_b = rblapack_b_out__;
  b = b_out__;

  dgegv_(&jobvl, &jobvr, &n, a, &lda, b, &ldb, alphar, alphai, beta, vl, &ldvl, vr, &ldvr, work, &lwork, &info);

  rblapack_info = INT2NUM(info);
  return rb_ary_new3(9, rblapack_alphar, rblapack_alphai, rblapack_beta, rblapack_vl, rblapack_vr, rblapack_work, rblapack_info, rblapack_a, rblapack_b);
}

void
init_lapack_dgegv(VALUE mLapack){
  rb_define_module_function(mLapack, "dgegv", rblapack_dgegv, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}

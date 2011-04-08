#include "rb_lapack.h"

extern VOID cggevx_(char *balanc, char *jobvl, char *jobvr, char *sense, integer *n, complex *a, integer *lda, complex *b, integer *ldb, complex *alpha, complex *beta, complex *vl, integer *ldvl, complex *vr, integer *ldvr, integer *ilo, integer *ihi, real *lscale, real *rscale, real *abnrm, real *bbnrm, real *rconde, real *rcondv, complex *work, integer *lwork, real *rwork, integer *iwork, logical *bwork, integer *info);

static VALUE sHelp, sUsage;

static VALUE
rblapack_cggevx(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_balanc;
  char balanc; 
  VALUE rblapack_jobvl;
  char jobvl; 
  VALUE rblapack_jobvr;
  char jobvr; 
  VALUE rblapack_sense;
  char sense; 
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
  VALUE rblapack_ilo;
  integer ilo; 
  VALUE rblapack_ihi;
  integer ihi; 
  VALUE rblapack_lscale;
  real *lscale; 
  VALUE rblapack_rscale;
  real *rscale; 
  VALUE rblapack_abnrm;
  real abnrm; 
  VALUE rblapack_bbnrm;
  real bbnrm; 
  VALUE rblapack_rconde;
  real *rconde; 
  VALUE rblapack_rcondv;
  real *rcondv; 
  VALUE rblapack_work;
  complex *work; 
  VALUE rblapack_info;
  integer info; 
  VALUE rblapack_a_out__;
  complex *a_out__;
  VALUE rblapack_b_out__;
  complex *b_out__;
  real *rwork;
  integer *iwork;
  logical *bwork;

  integer lda;
  integer n;
  integer ldb;
  integer ldvl;
  integer ldvr;
  integer lrwork;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  alpha, beta, vl, vr, ilo, ihi, lscale, rscale, abnrm, bbnrm, rconde, rcondv, work, info, a, b = NumRu::Lapack.cggevx( balanc, jobvl, jobvr, sense, a, b, lwork, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE CGGEVX( BALANC, JOBVL, JOBVR, SENSE, N, A, LDA, B, LDB, ALPHA, BETA, VL, LDVL, VR, LDVR, ILO, IHI, LSCALE, RSCALE, ABNRM, BBNRM, RCONDE, RCONDV, WORK, LWORK, RWORK, IWORK, BWORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  CGGEVX computes for a pair of N-by-N complex nonsymmetric matrices\n*  (A,B) the generalized eigenvalues, and optionally, the left and/or\n*  right generalized eigenvectors.\n*\n*  Optionally, it also computes a balancing transformation to improve\n*  the conditioning of the eigenvalues and eigenvectors (ILO, IHI,\n*  LSCALE, RSCALE, ABNRM, and BBNRM), reciprocal condition numbers for\n*  the eigenvalues (RCONDE), and reciprocal condition numbers for the\n*  right eigenvectors (RCONDV).\n*\n*  A generalized eigenvalue for a pair of matrices (A,B) is a scalar\n*  lambda or a ratio alpha/beta = lambda, such that A - lambda*B is\n*  singular. It is usually represented as the pair (alpha,beta), as\n*  there is a reasonable interpretation for beta=0, and even for both\n*  being zero.\n*\n*  The right eigenvector v(j) corresponding to the eigenvalue lambda(j)\n*  of (A,B) satisfies\n*                   A * v(j) = lambda(j) * B * v(j) .\n*  The left eigenvector u(j) corresponding to the eigenvalue lambda(j)\n*  of (A,B) satisfies\n*                   u(j)**H * A  = lambda(j) * u(j)**H * B.\n*  where u(j)**H is the conjugate-transpose of u(j).\n*\n*\n\n*  Arguments\n*  =========\n*\n*  BALANC  (input) CHARACTER*1\n*          Specifies the balance option to be performed:\n*          = 'N':  do not diagonally scale or permute;\n*          = 'P':  permute only;\n*          = 'S':  scale only;\n*          = 'B':  both permute and scale.\n*          Computed reciprocal condition numbers will be for the\n*          matrices after permuting and/or balancing. Permuting does\n*          not change condition numbers (in exact arithmetic), but\n*          balancing does.\n*\n*  JOBVL   (input) CHARACTER*1\n*          = 'N':  do not compute the left generalized eigenvectors;\n*          = 'V':  compute the left generalized eigenvectors.\n*\n*  JOBVR   (input) CHARACTER*1\n*          = 'N':  do not compute the right generalized eigenvectors;\n*          = 'V':  compute the right generalized eigenvectors.\n*\n*  SENSE   (input) CHARACTER*1\n*          Determines which reciprocal condition numbers are computed.\n*          = 'N': none are computed;\n*          = 'E': computed for eigenvalues only;\n*          = 'V': computed for eigenvectors only;\n*          = 'B': computed for eigenvalues and eigenvectors.\n*\n*  N       (input) INTEGER\n*          The order of the matrices A, B, VL, and VR.  N >= 0.\n*\n*  A       (input/output) COMPLEX array, dimension (LDA, N)\n*          On entry, the matrix A in the pair (A,B).\n*          On exit, A has been overwritten. If JOBVL='V' or JOBVR='V'\n*          or both, then A contains the first part of the complex Schur\n*          form of the \"balanced\" versions of the input A and B.\n*\n*  LDA     (input) INTEGER\n*          The leading dimension of A.  LDA >= max(1,N).\n*\n*  B       (input/output) COMPLEX array, dimension (LDB, N)\n*          On entry, the matrix B in the pair (A,B).\n*          On exit, B has been overwritten. If JOBVL='V' or JOBVR='V'\n*          or both, then B contains the second part of the complex\n*          Schur form of the \"balanced\" versions of the input A and B.\n*\n*  LDB     (input) INTEGER\n*          The leading dimension of B.  LDB >= max(1,N).\n*\n*  ALPHA   (output) COMPLEX array, dimension (N)\n*  BETA    (output) COMPLEX array, dimension (N)\n*          On exit, ALPHA(j)/BETA(j), j=1,...,N, will be the generalized\n*          eigenvalues.\n*\n*          Note: the quotient ALPHA(j)/BETA(j) ) may easily over- or\n*          underflow, and BETA(j) may even be zero.  Thus, the user\n*          should avoid naively computing the ratio ALPHA/BETA.\n*          However, ALPHA will be always less than and usually\n*          comparable with norm(A) in magnitude, and BETA always less\n*          than and usually comparable with norm(B).\n*\n*  VL      (output) COMPLEX array, dimension (LDVL,N)\n*          If JOBVL = 'V', the left generalized eigenvectors u(j) are\n*          stored one after another in the columns of VL, in the same\n*          order as their eigenvalues.\n*          Each eigenvector will be scaled so the largest component\n*          will have abs(real part) + abs(imag. part) = 1.\n*          Not referenced if JOBVL = 'N'.\n*\n*  LDVL    (input) INTEGER\n*          The leading dimension of the matrix VL. LDVL >= 1, and\n*          if JOBVL = 'V', LDVL >= N.\n*\n*  VR      (output) COMPLEX array, dimension (LDVR,N)\n*          If JOBVR = 'V', the right generalized eigenvectors v(j) are\n*          stored one after another in the columns of VR, in the same\n*          order as their eigenvalues.\n*          Each eigenvector will be scaled so the largest component\n*          will have abs(real part) + abs(imag. part) = 1.\n*          Not referenced if JOBVR = 'N'.\n*\n*  LDVR    (input) INTEGER\n*          The leading dimension of the matrix VR. LDVR >= 1, and\n*          if JOBVR = 'V', LDVR >= N.\n*\n*  ILO     (output) INTEGER\n*  IHI     (output) INTEGER\n*          ILO and IHI are integer values such that on exit\n*          A(i,j) = 0 and B(i,j) = 0 if i > j and\n*          j = 1,...,ILO-1 or i = IHI+1,...,N.\n*          If BALANC = 'N' or 'S', ILO = 1 and IHI = N.\n*\n*  LSCALE  (output) REAL array, dimension (N)\n*          Details of the permutations and scaling factors applied\n*          to the left side of A and B.  If PL(j) is the index of the\n*          row interchanged with row j, and DL(j) is the scaling\n*          factor applied to row j, then\n*            LSCALE(j) = PL(j)  for j = 1,...,ILO-1\n*                      = DL(j)  for j = ILO,...,IHI\n*                      = PL(j)  for j = IHI+1,...,N.\n*          The order in which the interchanges are made is N to IHI+1,\n*          then 1 to ILO-1.\n*\n*  RSCALE  (output) REAL array, dimension (N)\n*          Details of the permutations and scaling factors applied\n*          to the right side of A and B.  If PR(j) is the index of the\n*          column interchanged with column j, and DR(j) is the scaling\n*          factor applied to column j, then\n*            RSCALE(j) = PR(j)  for j = 1,...,ILO-1\n*                      = DR(j)  for j = ILO,...,IHI\n*                      = PR(j)  for j = IHI+1,...,N\n*          The order in which the interchanges are made is N to IHI+1,\n*          then 1 to ILO-1.\n*\n*  ABNRM   (output) REAL\n*          The one-norm of the balanced matrix A.\n*\n*  BBNRM   (output) REAL\n*          The one-norm of the balanced matrix B.\n*\n*  RCONDE  (output) REAL array, dimension (N)\n*          If SENSE = 'E' or 'B', the reciprocal condition numbers of\n*          the eigenvalues, stored in consecutive elements of the array.\n*          If SENSE = 'N' or 'V', RCONDE is not referenced.\n*\n*  RCONDV  (output) REAL array, dimension (N)\n*          If SENSE = 'V' or 'B', the estimated reciprocal condition\n*          numbers of the eigenvectors, stored in consecutive elements\n*          of the array. If the eigenvalues cannot be reordered to\n*          compute RCONDV(j), RCONDV(j) is set to 0; this can only occur\n*          when the true value would be very small anyway. \n*          If SENSE = 'N' or 'E', RCONDV is not referenced.\n*\n*  WORK    (workspace/output) COMPLEX array, dimension (MAX(1,LWORK))\n*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.\n*\n*  LWORK   (input) INTEGER\n*          The dimension of the array WORK. LWORK >= max(1,2*N).\n*          If SENSE = 'E', LWORK >= max(1,4*N).\n*          If SENSE = 'V' or 'B', LWORK >= max(1,2*N*N+2*N).\n*\n*          If LWORK = -1, then a workspace query is assumed; the routine\n*          only calculates the optimal size of the WORK array, returns\n*          this value as the first entry of the WORK array, and no error\n*          message related to LWORK is issued by XERBLA.\n*\n*  RWORK   (workspace) REAL array, dimension (lrwork)\n*          lrwork must be at least max(1,6*N) if BALANC = 'S' or 'B',\n*          and at least max(1,2*N) otherwise.\n*          Real workspace.\n*\n*  IWORK   (workspace) INTEGER array, dimension (N+2)\n*          If SENSE = 'E', IWORK is not referenced.\n*\n*  BWORK   (workspace) LOGICAL array, dimension (N)\n*          If SENSE = 'N', BWORK is not referenced.\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit\n*          < 0:  if INFO = -i, the i-th argument had an illegal value.\n*          = 1,...,N:\n*                The QZ iteration failed.  No eigenvectors have been\n*                calculated, but ALPHA(j) and BETA(j) should be correct\n*                for j=INFO+1,...,N.\n*          > N:  =N+1: other than QZ iteration failed in CHGEQZ.\n*                =N+2: error return from CTGEVC.\n*\n\n*  Further Details\n*  ===============\n*\n*  Balancing a matrix pair (A,B) includes, first, permuting rows and\n*  columns to isolate eigenvalues, second, applying diagonal similarity\n*  transformation to the rows and columns to make the rows and columns\n*  as close in norm as possible. The computed reciprocal condition\n*  numbers correspond to the balanced matrix. Permuting rows and columns\n*  will not change the condition numbers (in exact arithmetic) but\n*  diagonal scaling will.  For further explanation of balancing, see\n*  section 4.11.1.2 of LAPACK Users' Guide.\n*\n*  An approximate error bound on the chordal distance between the i-th\n*  computed generalized eigenvalue w and the corresponding exact\n*  eigenvalue lambda is\n*\n*       chord(w, lambda) <= EPS * norm(ABNRM, BBNRM) / RCONDE(I)\n*\n*  An approximate error bound for the angle between the i-th computed\n*  eigenvector VL(i) or VR(i) is given by\n*\n*       EPS * norm(ABNRM, BBNRM) / DIF(i).\n*\n*  For further explanation of the reciprocal condition numbers RCONDE\n*  and RCONDV, see section 4.11 of LAPACK User's Guide.\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  alpha, beta, vl, vr, ilo, ihi, lscale, rscale, abnrm, bbnrm, rconde, rcondv, work, info, a, b = NumRu::Lapack.cggevx( balanc, jobvl, jobvr, sense, a, b, lwork, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 7)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 7)", argc);
  rblapack_balanc = argv[0];
  rblapack_jobvl = argv[1];
  rblapack_jobvr = argv[2];
  rblapack_sense = argv[3];
  rblapack_a = argv[4];
  rblapack_b = argv[5];
  rblapack_lwork = argv[6];
  if (rb_options != Qnil) {
  }

  if (!NA_IsNArray(rblapack_a))
    rb_raise(rb_eArgError, "a (5th argument) must be NArray");
  if (NA_RANK(rblapack_a) != 2)
    rb_raise(rb_eArgError, "rank of a (5th argument) must be %d", 2);
  n = NA_SHAPE1(rblapack_a);
  lda = NA_SHAPE0(rblapack_a);
  if (NA_TYPE(rblapack_a) != NA_SCOMPLEX)
    rblapack_a = na_change_type(rblapack_a, NA_SCOMPLEX);
  a = NA_PTR_TYPE(rblapack_a, complex*);
  jobvr = StringValueCStr(rblapack_jobvr)[0];
  balanc = StringValueCStr(rblapack_balanc)[0];
  if (!NA_IsNArray(rblapack_b))
    rb_raise(rb_eArgError, "b (6th argument) must be NArray");
  if (NA_RANK(rblapack_b) != 2)
    rb_raise(rb_eArgError, "rank of b (6th argument) must be %d", 2);
  if (NA_SHAPE1(rblapack_b) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of b must be the same as shape 1 of a");
  ldb = NA_SHAPE0(rblapack_b);
  if (NA_TYPE(rblapack_b) != NA_SCOMPLEX)
    rblapack_b = na_change_type(rblapack_b, NA_SCOMPLEX);
  b = NA_PTR_TYPE(rblapack_b, complex*);
  jobvl = StringValueCStr(rblapack_jobvl)[0];
  lwork = NUM2INT(rblapack_lwork);
  sense = StringValueCStr(rblapack_sense)[0];
  lrwork = ((lsame_(&balanc,"S")) || (lsame_(&balanc,"B"))) ? MAX(1,6*n) : MAX(1,2*n);
  ldvl = lsame_(&jobvl,"V") ? n : 1;
  ldvr = lsame_(&jobvr,"V") ? n : 1;
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
    shape[0] = n;
    rblapack_lscale = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  lscale = NA_PTR_TYPE(rblapack_lscale, real*);
  {
    int shape[1];
    shape[0] = n;
    rblapack_rscale = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  rscale = NA_PTR_TYPE(rblapack_rscale, real*);
  {
    int shape[1];
    shape[0] = n;
    rblapack_rconde = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  rconde = NA_PTR_TYPE(rblapack_rconde, real*);
  {
    int shape[1];
    shape[0] = n;
    rblapack_rcondv = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  rcondv = NA_PTR_TYPE(rblapack_rcondv, real*);
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
  rwork = ALLOC_N(real, (lrwork));
  iwork = ALLOC_N(integer, (lsame_(&sense,"E") ? 0 : n+2));
  bwork = ALLOC_N(logical, (lsame_(&sense,"N") ? 0 : n));

  cggevx_(&balanc, &jobvl, &jobvr, &sense, &n, a, &lda, b, &ldb, alpha, beta, vl, &ldvl, vr, &ldvr, &ilo, &ihi, lscale, rscale, &abnrm, &bbnrm, rconde, rcondv, work, &lwork, rwork, iwork, bwork, &info);

  free(rwork);
  free(iwork);
  free(bwork);
  rblapack_ilo = INT2NUM(ilo);
  rblapack_ihi = INT2NUM(ihi);
  rblapack_abnrm = rb_float_new((double)abnrm);
  rblapack_bbnrm = rb_float_new((double)bbnrm);
  rblapack_info = INT2NUM(info);
  return rb_ary_new3(16, rblapack_alpha, rblapack_beta, rblapack_vl, rblapack_vr, rblapack_ilo, rblapack_ihi, rblapack_lscale, rblapack_rscale, rblapack_abnrm, rblapack_bbnrm, rblapack_rconde, rblapack_rcondv, rblapack_work, rblapack_info, rblapack_a, rblapack_b);
}

void
init_lapack_cggevx(VALUE mLapack){
  rb_define_module_function(mLapack, "cggevx", rblapack_cggevx, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}

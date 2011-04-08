#include "rb_lapack.h"

extern VOID sgesvx_(char *fact, char *trans, integer *n, integer *nrhs, real *a, integer *lda, real *af, integer *ldaf, integer *ipiv, char *equed, real *r, real *c, real *b, integer *ldb, real *x, integer *ldx, real *rcond, real *ferr, real *berr, real *work, integer *iwork, integer *info);

static VALUE sHelp, sUsage;

static VALUE
rblapack_sgesvx(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_fact;
  char fact; 
  VALUE rblapack_trans;
  char trans; 
  VALUE rblapack_a;
  real *a; 
  VALUE rblapack_af;
  real *af; 
  VALUE rblapack_ipiv;
  integer *ipiv; 
  VALUE rblapack_equed;
  char equed; 
  VALUE rblapack_r;
  real *r; 
  VALUE rblapack_c;
  real *c; 
  VALUE rblapack_b;
  real *b; 
  VALUE rblapack_x;
  real *x; 
  VALUE rblapack_rcond;
  real rcond; 
  VALUE rblapack_ferr;
  real *ferr; 
  VALUE rblapack_berr;
  real *berr; 
  VALUE rblapack_work;
  real *work; 
  VALUE rblapack_info;
  integer info; 
  VALUE rblapack_a_out__;
  real *a_out__;
  VALUE rblapack_af_out__;
  real *af_out__;
  VALUE rblapack_ipiv_out__;
  integer *ipiv_out__;
  VALUE rblapack_r_out__;
  real *r_out__;
  VALUE rblapack_c_out__;
  real *c_out__;
  VALUE rblapack_b_out__;
  real *b_out__;
  integer *iwork;

  integer lda;
  integer n;
  integer ldaf;
  integer ldb;
  integer nrhs;
  integer ldx;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  x, rcond, ferr, berr, work, info, a, af, ipiv, equed, r, c, b = NumRu::Lapack.sgesvx( fact, trans, a, af, ipiv, equed, r, c, b, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE SGESVX( FACT, TRANS, N, NRHS, A, LDA, AF, LDAF, IPIV, EQUED, R, C, B, LDB, X, LDX, RCOND, FERR, BERR, WORK, IWORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  SGESVX uses the LU factorization to compute the solution to a real\n*  system of linear equations\n*     A * X = B,\n*  where A is an N-by-N matrix and X and B are N-by-NRHS matrices.\n*\n*  Error bounds on the solution and a condition estimate are also\n*  provided.\n*\n*  Description\n*  ===========\n*\n*  The following steps are performed:\n*\n*  1. If FACT = 'E', real scaling factors are computed to equilibrate\n*     the system:\n*        TRANS = 'N':  diag(R)*A*diag(C)     *inv(diag(C))*X = diag(R)*B\n*        TRANS = 'T': (diag(R)*A*diag(C))**T *inv(diag(R))*X = diag(C)*B\n*        TRANS = 'C': (diag(R)*A*diag(C))**H *inv(diag(R))*X = diag(C)*B\n*     Whether or not the system will be equilibrated depends on the\n*     scaling of the matrix A, but if equilibration is used, A is\n*     overwritten by diag(R)*A*diag(C) and B by diag(R)*B (if TRANS='N')\n*     or diag(C)*B (if TRANS = 'T' or 'C').\n*\n*  2. If FACT = 'N' or 'E', the LU decomposition is used to factor the\n*     matrix A (after equilibration if FACT = 'E') as\n*        A = P * L * U,\n*     where P is a permutation matrix, L is a unit lower triangular\n*     matrix, and U is upper triangular.\n*\n*  3. If some U(i,i)=0, so that U is exactly singular, then the routine\n*     returns with INFO = i. Otherwise, the factored form of A is used\n*     to estimate the condition number of the matrix A.  If the\n*     reciprocal of the condition number is less than machine precision,\n*     INFO = N+1 is returned as a warning, but the routine still goes on\n*     to solve for X and compute error bounds as described below.\n*\n*  4. The system of equations is solved for X using the factored form\n*     of A.\n*\n*  5. Iterative refinement is applied to improve the computed solution\n*     matrix and calculate error bounds and backward error estimates\n*     for it.\n*\n*  6. If equilibration was used, the matrix X is premultiplied by\n*     diag(C) (if TRANS = 'N') or diag(R) (if TRANS = 'T' or 'C') so\n*     that it solves the original system before equilibration.\n*\n\n*  Arguments\n*  =========\n*\n*  FACT    (input) CHARACTER*1\n*          Specifies whether or not the factored form of the matrix A is\n*          supplied on entry, and if not, whether the matrix A should be\n*          equilibrated before it is factored.\n*          = 'F':  On entry, AF and IPIV contain the factored form of A.\n*                  If EQUED is not 'N', the matrix A has been\n*                  equilibrated with scaling factors given by R and C.\n*                  A, AF, and IPIV are not modified.\n*          = 'N':  The matrix A will be copied to AF and factored.\n*          = 'E':  The matrix A will be equilibrated if necessary, then\n*                  copied to AF and factored.\n*\n*  TRANS   (input) CHARACTER*1\n*          Specifies the form of the system of equations:\n*          = 'N':  A * X = B     (No transpose)\n*          = 'T':  A**T * X = B  (Transpose)\n*          = 'C':  A**H * X = B  (Transpose)\n*\n*  N       (input) INTEGER\n*          The number of linear equations, i.e., the order of the\n*          matrix A.  N >= 0.\n*\n*  NRHS    (input) INTEGER\n*          The number of right hand sides, i.e., the number of columns\n*          of the matrices B and X.  NRHS >= 0.\n*\n*  A       (input/output) REAL array, dimension (LDA,N)\n*          On entry, the N-by-N matrix A.  If FACT = 'F' and EQUED is\n*          not 'N', then A must have been equilibrated by the scaling\n*          factors in R and/or C.  A is not modified if FACT = 'F' or\n*          'N', or if FACT = 'E' and EQUED = 'N' on exit.\n*\n*          On exit, if EQUED .ne. 'N', A is scaled as follows:\n*          EQUED = 'R':  A := diag(R) * A\n*          EQUED = 'C':  A := A * diag(C)\n*          EQUED = 'B':  A := diag(R) * A * diag(C).\n*\n*  LDA     (input) INTEGER\n*          The leading dimension of the array A.  LDA >= max(1,N).\n*\n*  AF      (input or output) REAL array, dimension (LDAF,N)\n*          If FACT = 'F', then AF is an input argument and on entry\n*          contains the factors L and U from the factorization\n*          A = P*L*U as computed by SGETRF.  If EQUED .ne. 'N', then\n*          AF is the factored form of the equilibrated matrix A.\n*\n*          If FACT = 'N', then AF is an output argument and on exit\n*          returns the factors L and U from the factorization A = P*L*U\n*          of the original matrix A.\n*\n*          If FACT = 'E', then AF is an output argument and on exit\n*          returns the factors L and U from the factorization A = P*L*U\n*          of the equilibrated matrix A (see the description of A for\n*          the form of the equilibrated matrix).\n*\n*  LDAF    (input) INTEGER\n*          The leading dimension of the array AF.  LDAF >= max(1,N).\n*\n*  IPIV    (input or output) INTEGER array, dimension (N)\n*          If FACT = 'F', then IPIV is an input argument and on entry\n*          contains the pivot indices from the factorization A = P*L*U\n*          as computed by SGETRF; row i of the matrix was interchanged\n*          with row IPIV(i).\n*\n*          If FACT = 'N', then IPIV is an output argument and on exit\n*          contains the pivot indices from the factorization A = P*L*U\n*          of the original matrix A.\n*\n*          If FACT = 'E', then IPIV is an output argument and on exit\n*          contains the pivot indices from the factorization A = P*L*U\n*          of the equilibrated matrix A.\n*\n*  EQUED   (input or output) CHARACTER*1\n*          Specifies the form of equilibration that was done.\n*          = 'N':  No equilibration (always true if FACT = 'N').\n*          = 'R':  Row equilibration, i.e., A has been premultiplied by\n*                  diag(R).\n*          = 'C':  Column equilibration, i.e., A has been postmultiplied\n*                  by diag(C).\n*          = 'B':  Both row and column equilibration, i.e., A has been\n*                  replaced by diag(R) * A * diag(C).\n*          EQUED is an input argument if FACT = 'F'; otherwise, it is an\n*          output argument.\n*\n*  R       (input or output) REAL array, dimension (N)\n*          The row scale factors for A.  If EQUED = 'R' or 'B', A is\n*          multiplied on the left by diag(R); if EQUED = 'N' or 'C', R\n*          is not accessed.  R is an input argument if FACT = 'F';\n*          otherwise, R is an output argument.  If FACT = 'F' and\n*          EQUED = 'R' or 'B', each element of R must be positive.\n*\n*  C       (input or output) REAL array, dimension (N)\n*          The column scale factors for A.  If EQUED = 'C' or 'B', A is\n*          multiplied on the right by diag(C); if EQUED = 'N' or 'R', C\n*          is not accessed.  C is an input argument if FACT = 'F';\n*          otherwise, C is an output argument.  If FACT = 'F' and\n*          EQUED = 'C' or 'B', each element of C must be positive.\n*\n*  B       (input/output) REAL array, dimension (LDB,NRHS)\n*          On entry, the N-by-NRHS right hand side matrix B.\n*          On exit,\n*          if EQUED = 'N', B is not modified;\n*          if TRANS = 'N' and EQUED = 'R' or 'B', B is overwritten by\n*          diag(R)*B;\n*          if TRANS = 'T' or 'C' and EQUED = 'C' or 'B', B is\n*          overwritten by diag(C)*B.\n*\n*  LDB     (input) INTEGER\n*          The leading dimension of the array B.  LDB >= max(1,N).\n*\n*  X       (output) REAL array, dimension (LDX,NRHS)\n*          If INFO = 0 or INFO = N+1, the N-by-NRHS solution matrix X\n*          to the original system of equations.  Note that A and B are\n*          modified on exit if EQUED .ne. 'N', and the solution to the\n*          equilibrated system is inv(diag(C))*X if TRANS = 'N' and\n*          EQUED = 'C' or 'B', or inv(diag(R))*X if TRANS = 'T' or 'C'\n*          and EQUED = 'R' or 'B'.\n*\n*  LDX     (input) INTEGER\n*          The leading dimension of the array X.  LDX >= max(1,N).\n*\n*  RCOND   (output) REAL\n*          The estimate of the reciprocal condition number of the matrix\n*          A after equilibration (if done).  If RCOND is less than the\n*          machine precision (in particular, if RCOND = 0), the matrix\n*          is singular to working precision.  This condition is\n*          indicated by a return code of INFO > 0.\n*\n*  FERR    (output) REAL array, dimension (NRHS)\n*          The estimated forward error bound for each solution vector\n*          X(j) (the j-th column of the solution matrix X).\n*          If XTRUE is the true solution corresponding to X(j), FERR(j)\n*          is an estimated upper bound for the magnitude of the largest\n*          element in (X(j) - XTRUE) divided by the magnitude of the\n*          largest element in X(j).  The estimate is as reliable as\n*          the estimate for RCOND, and is almost always a slight\n*          overestimate of the true error.\n*\n*  BERR    (output) REAL array, dimension (NRHS)\n*          The componentwise relative backward error of each solution\n*          vector X(j) (i.e., the smallest relative change in\n*          any element of A or B that makes X(j) an exact solution).\n*\n*  WORK    (workspace/output) REAL array, dimension (4*N)\n*          On exit, WORK(1) contains the reciprocal pivot growth\n*          factor norm(A)/norm(U). The \"max absolute element\" norm is\n*          used. If WORK(1) is much less than 1, then the stability\n*          of the LU factorization of the (equilibrated) matrix A\n*          could be poor. This also means that the solution X, condition\n*          estimator RCOND, and forward error bound FERR could be\n*          unreliable. If factorization fails with 0<INFO<=N, then\n*          WORK(1) contains the reciprocal pivot growth factor for the\n*          leading INFO columns of A.\n*\n*  IWORK   (workspace) INTEGER array, dimension (N)\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit\n*          < 0:  if INFO = -i, the i-th argument had an illegal value\n*          > 0:  if INFO = i, and i is\n*                <= N:  U(i,i) is exactly zero.  The factorization has\n*                       been completed, but the factor U is exactly\n*                       singular, so the solution and error bounds\n*                       could not be computed. RCOND = 0 is returned.\n*                = N+1: U is nonsingular, but RCOND is less than machine\n*                       precision, meaning that the matrix is singular\n*                       to working precision.  Nevertheless, the\n*                       solution and error bounds are computed because\n*                       there are a number of situations where the\n*                       computed solution can be more accurate than the\n*                       value of RCOND would suggest.\n*\n\n*  =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  x, rcond, ferr, berr, work, info, a, af, ipiv, equed, r, c, b = NumRu::Lapack.sgesvx( fact, trans, a, af, ipiv, equed, r, c, b, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 9)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 9)", argc);
  rblapack_fact = argv[0];
  rblapack_trans = argv[1];
  rblapack_a = argv[2];
  rblapack_af = argv[3];
  rblapack_ipiv = argv[4];
  rblapack_equed = argv[5];
  rblapack_r = argv[6];
  rblapack_c = argv[7];
  rblapack_b = argv[8];
  if (rb_options != Qnil) {
  }

  trans = StringValueCStr(rblapack_trans)[0];
  if (!NA_IsNArray(rblapack_ipiv))
    rb_raise(rb_eArgError, "ipiv (5th argument) must be NArray");
  if (NA_RANK(rblapack_ipiv) != 1)
    rb_raise(rb_eArgError, "rank of ipiv (5th argument) must be %d", 1);
  n = NA_SHAPE0(rblapack_ipiv);
  if (NA_TYPE(rblapack_ipiv) != NA_LINT)
    rblapack_ipiv = na_change_type(rblapack_ipiv, NA_LINT);
  ipiv = NA_PTR_TYPE(rblapack_ipiv, integer*);
  if (!NA_IsNArray(rblapack_a))
    rb_raise(rb_eArgError, "a (3th argument) must be NArray");
  if (NA_RANK(rblapack_a) != 2)
    rb_raise(rb_eArgError, "rank of a (3th argument) must be %d", 2);
  if (NA_SHAPE1(rblapack_a) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of a must be the same as shape 0 of ipiv");
  lda = NA_SHAPE0(rblapack_a);
  if (NA_TYPE(rblapack_a) != NA_SFLOAT)
    rblapack_a = na_change_type(rblapack_a, NA_SFLOAT);
  a = NA_PTR_TYPE(rblapack_a, real*);
  if (!NA_IsNArray(rblapack_b))
    rb_raise(rb_eArgError, "b (9th argument) must be NArray");
  if (NA_RANK(rblapack_b) != 2)
    rb_raise(rb_eArgError, "rank of b (9th argument) must be %d", 2);
  nrhs = NA_SHAPE1(rblapack_b);
  ldb = NA_SHAPE0(rblapack_b);
  if (NA_TYPE(rblapack_b) != NA_SFLOAT)
    rblapack_b = na_change_type(rblapack_b, NA_SFLOAT);
  b = NA_PTR_TYPE(rblapack_b, real*);
  if (!NA_IsNArray(rblapack_c))
    rb_raise(rb_eArgError, "c (8th argument) must be NArray");
  if (NA_RANK(rblapack_c) != 1)
    rb_raise(rb_eArgError, "rank of c (8th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_c) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of c must be the same as shape 0 of ipiv");
  if (NA_TYPE(rblapack_c) != NA_SFLOAT)
    rblapack_c = na_change_type(rblapack_c, NA_SFLOAT);
  c = NA_PTR_TYPE(rblapack_c, real*);
  equed = StringValueCStr(rblapack_equed)[0];
  if (!NA_IsNArray(rblapack_af))
    rb_raise(rb_eArgError, "af (4th argument) must be NArray");
  if (NA_RANK(rblapack_af) != 2)
    rb_raise(rb_eArgError, "rank of af (4th argument) must be %d", 2);
  if (NA_SHAPE1(rblapack_af) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of af must be the same as shape 0 of ipiv");
  ldaf = NA_SHAPE0(rblapack_af);
  if (NA_TYPE(rblapack_af) != NA_SFLOAT)
    rblapack_af = na_change_type(rblapack_af, NA_SFLOAT);
  af = NA_PTR_TYPE(rblapack_af, real*);
  if (!NA_IsNArray(rblapack_r))
    rb_raise(rb_eArgError, "r (7th argument) must be NArray");
  if (NA_RANK(rblapack_r) != 1)
    rb_raise(rb_eArgError, "rank of r (7th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_r) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of r must be the same as shape 0 of ipiv");
  if (NA_TYPE(rblapack_r) != NA_SFLOAT)
    rblapack_r = na_change_type(rblapack_r, NA_SFLOAT);
  r = NA_PTR_TYPE(rblapack_r, real*);
  fact = StringValueCStr(rblapack_fact)[0];
  ldx = MAX(1,n);
  {
    int shape[2];
    shape[0] = ldx;
    shape[1] = nrhs;
    rblapack_x = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  x = NA_PTR_TYPE(rblapack_x, real*);
  {
    int shape[1];
    shape[0] = nrhs;
    rblapack_ferr = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  ferr = NA_PTR_TYPE(rblapack_ferr, real*);
  {
    int shape[1];
    shape[0] = nrhs;
    rblapack_berr = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  berr = NA_PTR_TYPE(rblapack_berr, real*);
  {
    int shape[1];
    shape[0] = 4*n;
    rblapack_work = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  work = NA_PTR_TYPE(rblapack_work, real*);
  {
    int shape[2];
    shape[0] = lda;
    shape[1] = n;
    rblapack_a_out__ = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  a_out__ = NA_PTR_TYPE(rblapack_a_out__, real*);
  MEMCPY(a_out__, a, real, NA_TOTAL(rblapack_a));
  rblapack_a = rblapack_a_out__;
  a = a_out__;
  {
    int shape[2];
    shape[0] = ldaf;
    shape[1] = n;
    rblapack_af_out__ = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  af_out__ = NA_PTR_TYPE(rblapack_af_out__, real*);
  MEMCPY(af_out__, af, real, NA_TOTAL(rblapack_af));
  rblapack_af = rblapack_af_out__;
  af = af_out__;
  {
    int shape[1];
    shape[0] = n;
    rblapack_ipiv_out__ = na_make_object(NA_LINT, 1, shape, cNArray);
  }
  ipiv_out__ = NA_PTR_TYPE(rblapack_ipiv_out__, integer*);
  MEMCPY(ipiv_out__, ipiv, integer, NA_TOTAL(rblapack_ipiv));
  rblapack_ipiv = rblapack_ipiv_out__;
  ipiv = ipiv_out__;
  {
    int shape[1];
    shape[0] = n;
    rblapack_r_out__ = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  r_out__ = NA_PTR_TYPE(rblapack_r_out__, real*);
  MEMCPY(r_out__, r, real, NA_TOTAL(rblapack_r));
  rblapack_r = rblapack_r_out__;
  r = r_out__;
  {
    int shape[1];
    shape[0] = n;
    rblapack_c_out__ = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  c_out__ = NA_PTR_TYPE(rblapack_c_out__, real*);
  MEMCPY(c_out__, c, real, NA_TOTAL(rblapack_c));
  rblapack_c = rblapack_c_out__;
  c = c_out__;
  {
    int shape[2];
    shape[0] = ldb;
    shape[1] = nrhs;
    rblapack_b_out__ = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  b_out__ = NA_PTR_TYPE(rblapack_b_out__, real*);
  MEMCPY(b_out__, b, real, NA_TOTAL(rblapack_b));
  rblapack_b = rblapack_b_out__;
  b = b_out__;
  iwork = ALLOC_N(integer, (n));

  sgesvx_(&fact, &trans, &n, &nrhs, a, &lda, af, &ldaf, ipiv, &equed, r, c, b, &ldb, x, &ldx, &rcond, ferr, berr, work, iwork, &info);

  free(iwork);
  rblapack_rcond = rb_float_new((double)rcond);
  rblapack_info = INT2NUM(info);
  rblapack_equed = rb_str_new(&equed,1);
  return rb_ary_new3(13, rblapack_x, rblapack_rcond, rblapack_ferr, rblapack_berr, rblapack_work, rblapack_info, rblapack_a, rblapack_af, rblapack_ipiv, rblapack_equed, rblapack_r, rblapack_c, rblapack_b);
}

void
init_lapack_sgesvx(VALUE mLapack){
  rb_define_module_function(mLapack, "sgesvx", rblapack_sgesvx, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}

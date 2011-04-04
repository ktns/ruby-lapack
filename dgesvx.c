#include "rb_lapack.h"

extern VOID dgesvx_(char *fact, char *trans, integer *n, integer *nrhs, doublereal *a, integer *lda, doublereal *af, integer *ldaf, integer *ipiv, char *equed, doublereal *r, doublereal *c, doublereal *b, integer *ldb, doublereal *x, integer *ldx, doublereal *rcond, doublereal *ferr, doublereal *berr, doublereal *work, integer *iwork, integer *info);

static VALUE
rb_dgesvx(int argc, VALUE *argv, VALUE self){
  VALUE rb_fact;
  char fact; 
  VALUE rb_trans;
  char trans; 
  VALUE rb_a;
  doublereal *a; 
  VALUE rb_af;
  doublereal *af; 
  VALUE rb_ipiv;
  integer *ipiv; 
  VALUE rb_equed;
  char equed; 
  VALUE rb_r;
  doublereal *r; 
  VALUE rb_c;
  doublereal *c; 
  VALUE rb_b;
  doublereal *b; 
  VALUE rb_x;
  doublereal *x; 
  VALUE rb_rcond;
  doublereal rcond; 
  VALUE rb_ferr;
  doublereal *ferr; 
  VALUE rb_berr;
  doublereal *berr; 
  VALUE rb_work;
  doublereal *work; 
  VALUE rb_info;
  integer info; 
  VALUE rb_a_out__;
  doublereal *a_out__;
  VALUE rb_af_out__;
  doublereal *af_out__;
  VALUE rb_ipiv_out__;
  integer *ipiv_out__;
  VALUE rb_r_out__;
  doublereal *r_out__;
  VALUE rb_c_out__;
  doublereal *c_out__;
  VALUE rb_b_out__;
  doublereal *b_out__;
  integer *iwork;

  integer lda;
  integer n;
  integer ldaf;
  integer ldb;
  integer nrhs;
  integer ldx;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  x, rcond, ferr, berr, work, info, a, af, ipiv, equed, r, c, b = NumRu::Lapack.dgesvx( fact, trans, a, af, ipiv, equed, r, c, b)\n    or\n  NumRu::Lapack.dgesvx  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE DGESVX( FACT, TRANS, N, NRHS, A, LDA, AF, LDAF, IPIV, EQUED, R, C, B, LDB, X, LDX, RCOND, FERR, BERR, WORK, IWORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  DGESVX uses the LU factorization to compute the solution to a real\n*  system of linear equations\n*     A * X = B,\n*  where A is an N-by-N matrix and X and B are N-by-NRHS matrices.\n*\n*  Error bounds on the solution and a condition estimate are also\n*  provided.\n*\n*  Description\n*  ===========\n*\n*  The following steps are performed:\n*\n*  1. If FACT = 'E', real scaling factors are computed to equilibrate\n*     the system:\n*        TRANS = 'N':  diag(R)*A*diag(C)     *inv(diag(C))*X = diag(R)*B\n*        TRANS = 'T': (diag(R)*A*diag(C))**T *inv(diag(R))*X = diag(C)*B\n*        TRANS = 'C': (diag(R)*A*diag(C))**H *inv(diag(R))*X = diag(C)*B\n*     Whether or not the system will be equilibrated depends on the\n*     scaling of the matrix A, but if equilibration is used, A is\n*     overwritten by diag(R)*A*diag(C) and B by diag(R)*B (if TRANS='N')\n*     or diag(C)*B (if TRANS = 'T' or 'C').\n*\n*  2. If FACT = 'N' or 'E', the LU decomposition is used to factor the\n*     matrix A (after equilibration if FACT = 'E') as\n*        A = P * L * U,\n*     where P is a permutation matrix, L is a unit lower triangular\n*     matrix, and U is upper triangular.\n*\n*  3. If some U(i,i)=0, so that U is exactly singular, then the routine\n*     returns with INFO = i. Otherwise, the factored form of A is used\n*     to estimate the condition number of the matrix A.  If the\n*     reciprocal of the condition number is less than machine precision,\n*     INFO = N+1 is returned as a warning, but the routine still goes on\n*     to solve for X and compute error bounds as described below.\n*\n*  4. The system of equations is solved for X using the factored form\n*     of A.\n*\n*  5. Iterative refinement is applied to improve the computed solution\n*     matrix and calculate error bounds and backward error estimates\n*     for it.\n*\n*  6. If equilibration was used, the matrix X is premultiplied by\n*     diag(C) (if TRANS = 'N') or diag(R) (if TRANS = 'T' or 'C') so\n*     that it solves the original system before equilibration.\n*\n\n*  Arguments\n*  =========\n*\n*  FACT    (input) CHARACTER*1\n*          Specifies whether or not the factored form of the matrix A is\n*          supplied on entry, and if not, whether the matrix A should be\n*          equilibrated before it is factored.\n*          = 'F':  On entry, AF and IPIV contain the factored form of A.\n*                  If EQUED is not 'N', the matrix A has been\n*                  equilibrated with scaling factors given by R and C.\n*                  A, AF, and IPIV are not modified.\n*          = 'N':  The matrix A will be copied to AF and factored.\n*          = 'E':  The matrix A will be equilibrated if necessary, then\n*                  copied to AF and factored.\n*\n*  TRANS   (input) CHARACTER*1\n*          Specifies the form of the system of equations:\n*          = 'N':  A * X = B     (No transpose)\n*          = 'T':  A**T * X = B  (Transpose)\n*          = 'C':  A**H * X = B  (Transpose)\n*\n*  N       (input) INTEGER\n*          The number of linear equations, i.e., the order of the\n*          matrix A.  N >= 0.\n*\n*  NRHS    (input) INTEGER\n*          The number of right hand sides, i.e., the number of columns\n*          of the matrices B and X.  NRHS >= 0.\n*\n*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)\n*          On entry, the N-by-N matrix A.  If FACT = 'F' and EQUED is\n*          not 'N', then A must have been equilibrated by the scaling\n*          factors in R and/or C.  A is not modified if FACT = 'F' or\n*          'N', or if FACT = 'E' and EQUED = 'N' on exit.\n*\n*          On exit, if EQUED .ne. 'N', A is scaled as follows:\n*          EQUED = 'R':  A := diag(R) * A\n*          EQUED = 'C':  A := A * diag(C)\n*          EQUED = 'B':  A := diag(R) * A * diag(C).\n*\n*  LDA     (input) INTEGER\n*          The leading dimension of the array A.  LDA >= max(1,N).\n*\n*  AF      (input or output) DOUBLE PRECISION array, dimension (LDAF,N)\n*          If FACT = 'F', then AF is an input argument and on entry\n*          contains the factors L and U from the factorization\n*          A = P*L*U as computed by DGETRF.  If EQUED .ne. 'N', then\n*          AF is the factored form of the equilibrated matrix A.\n*\n*          If FACT = 'N', then AF is an output argument and on exit\n*          returns the factors L and U from the factorization A = P*L*U\n*          of the original matrix A.\n*\n*          If FACT = 'E', then AF is an output argument and on exit\n*          returns the factors L and U from the factorization A = P*L*U\n*          of the equilibrated matrix A (see the description of A for\n*          the form of the equilibrated matrix).\n*\n*  LDAF    (input) INTEGER\n*          The leading dimension of the array AF.  LDAF >= max(1,N).\n*\n*  IPIV    (input or output) INTEGER array, dimension (N)\n*          If FACT = 'F', then IPIV is an input argument and on entry\n*          contains the pivot indices from the factorization A = P*L*U\n*          as computed by DGETRF; row i of the matrix was interchanged\n*          with row IPIV(i).\n*\n*          If FACT = 'N', then IPIV is an output argument and on exit\n*          contains the pivot indices from the factorization A = P*L*U\n*          of the original matrix A.\n*\n*          If FACT = 'E', then IPIV is an output argument and on exit\n*          contains the pivot indices from the factorization A = P*L*U\n*          of the equilibrated matrix A.\n*\n*  EQUED   (input or output) CHARACTER*1\n*          Specifies the form of equilibration that was done.\n*          = 'N':  No equilibration (always true if FACT = 'N').\n*          = 'R':  Row equilibration, i.e., A has been premultiplied by\n*                  diag(R).\n*          = 'C':  Column equilibration, i.e., A has been postmultiplied\n*                  by diag(C).\n*          = 'B':  Both row and column equilibration, i.e., A has been\n*                  replaced by diag(R) * A * diag(C).\n*          EQUED is an input argument if FACT = 'F'; otherwise, it is an\n*          output argument.\n*\n*  R       (input or output) DOUBLE PRECISION array, dimension (N)\n*          The row scale factors for A.  If EQUED = 'R' or 'B', A is\n*          multiplied on the left by diag(R); if EQUED = 'N' or 'C', R\n*          is not accessed.  R is an input argument if FACT = 'F';\n*          otherwise, R is an output argument.  If FACT = 'F' and\n*          EQUED = 'R' or 'B', each element of R must be positive.\n*\n*  C       (input or output) DOUBLE PRECISION array, dimension (N)\n*          The column scale factors for A.  If EQUED = 'C' or 'B', A is\n*          multiplied on the right by diag(C); if EQUED = 'N' or 'R', C\n*          is not accessed.  C is an input argument if FACT = 'F';\n*          otherwise, C is an output argument.  If FACT = 'F' and\n*          EQUED = 'C' or 'B', each element of C must be positive.\n*\n*  B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)\n*          On entry, the N-by-NRHS right hand side matrix B.\n*          On exit,\n*          if EQUED = 'N', B is not modified;\n*          if TRANS = 'N' and EQUED = 'R' or 'B', B is overwritten by\n*          diag(R)*B;\n*          if TRANS = 'T' or 'C' and EQUED = 'C' or 'B', B is\n*          overwritten by diag(C)*B.\n*\n*  LDB     (input) INTEGER\n*          The leading dimension of the array B.  LDB >= max(1,N).\n*\n*  X       (output) DOUBLE PRECISION array, dimension (LDX,NRHS)\n*          If INFO = 0 or INFO = N+1, the N-by-NRHS solution matrix X\n*          to the original system of equations.  Note that A and B are\n*          modified on exit if EQUED .ne. 'N', and the solution to the\n*          equilibrated system is inv(diag(C))*X if TRANS = 'N' and\n*          EQUED = 'C' or 'B', or inv(diag(R))*X if TRANS = 'T' or 'C'\n*          and EQUED = 'R' or 'B'.\n*\n*  LDX     (input) INTEGER\n*          The leading dimension of the array X.  LDX >= max(1,N).\n*\n*  RCOND   (output) DOUBLE PRECISION\n*          The estimate of the reciprocal condition number of the matrix\n*          A after equilibration (if done).  If RCOND is less than the\n*          machine precision (in particular, if RCOND = 0), the matrix\n*          is singular to working precision.  This condition is\n*          indicated by a return code of INFO > 0.\n*\n*  FERR    (output) DOUBLE PRECISION array, dimension (NRHS)\n*          The estimated forward error bound for each solution vector\n*          X(j) (the j-th column of the solution matrix X).\n*          If XTRUE is the true solution corresponding to X(j), FERR(j)\n*          is an estimated upper bound for the magnitude of the largest\n*          element in (X(j) - XTRUE) divided by the magnitude of the\n*          largest element in X(j).  The estimate is as reliable as\n*          the estimate for RCOND, and is almost always a slight\n*          overestimate of the true error.\n*\n*  BERR    (output) DOUBLE PRECISION array, dimension (NRHS)\n*          The componentwise relative backward error of each solution\n*          vector X(j) (i.e., the smallest relative change in\n*          any element of A or B that makes X(j) an exact solution).\n*\n*  WORK    (workspace/output) DOUBLE PRECISION array, dimension (4*N)\n*          On exit, WORK(1) contains the reciprocal pivot growth\n*          factor norm(A)/norm(U). The \"max absolute element\" norm is\n*          used. If WORK(1) is much less than 1, then the stability\n*          of the LU factorization of the (equilibrated) matrix A\n*          could be poor. This also means that the solution X, condition\n*          estimator RCOND, and forward error bound FERR could be\n*          unreliable. If factorization fails with 0<INFO<=N, then\n*          WORK(1) contains the reciprocal pivot growth factor for the\n*          leading INFO columns of A.\n*\n*  IWORK   (workspace) INTEGER array, dimension (N)\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit\n*          < 0:  if INFO = -i, the i-th argument had an illegal value\n*          > 0:  if INFO = i, and i is\n*                <= N:  U(i,i) is exactly zero.  The factorization has\n*                       been completed, but the factor U is exactly\n*                       singular, so the solution and error bounds\n*                       could not be computed. RCOND = 0 is returned.\n*                = N+1: U is nonsingular, but RCOND is less than machine\n*                       precision, meaning that the matrix is singular\n*                       to working precision.  Nevertheless, the\n*                       solution and error bounds are computed because\n*                       there are a number of situations where the\n*                       computed solution can be more accurate than the\n*                       value of RCOND would suggest.\n*\n\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 9)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 9)", argc);
  rb_fact = argv[0];
  rb_trans = argv[1];
  rb_a = argv[2];
  rb_af = argv[3];
  rb_ipiv = argv[4];
  rb_equed = argv[5];
  rb_r = argv[6];
  rb_c = argv[7];
  rb_b = argv[8];

  trans = StringValueCStr(rb_trans)[0];
  if (!NA_IsNArray(rb_ipiv))
    rb_raise(rb_eArgError, "ipiv (5th argument) must be NArray");
  if (NA_RANK(rb_ipiv) != 1)
    rb_raise(rb_eArgError, "rank of ipiv (5th argument) must be %d", 1);
  n = NA_SHAPE0(rb_ipiv);
  if (NA_TYPE(rb_ipiv) != NA_LINT)
    rb_ipiv = na_change_type(rb_ipiv, NA_LINT);
  ipiv = NA_PTR_TYPE(rb_ipiv, integer*);
  if (!NA_IsNArray(rb_a))
    rb_raise(rb_eArgError, "a (3th argument) must be NArray");
  if (NA_RANK(rb_a) != 2)
    rb_raise(rb_eArgError, "rank of a (3th argument) must be %d", 2);
  if (NA_SHAPE1(rb_a) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of a must be the same as shape 0 of ipiv");
  lda = NA_SHAPE0(rb_a);
  if (NA_TYPE(rb_a) != NA_DFLOAT)
    rb_a = na_change_type(rb_a, NA_DFLOAT);
  a = NA_PTR_TYPE(rb_a, doublereal*);
  if (!NA_IsNArray(rb_b))
    rb_raise(rb_eArgError, "b (9th argument) must be NArray");
  if (NA_RANK(rb_b) != 2)
    rb_raise(rb_eArgError, "rank of b (9th argument) must be %d", 2);
  nrhs = NA_SHAPE1(rb_b);
  ldb = NA_SHAPE0(rb_b);
  if (NA_TYPE(rb_b) != NA_DFLOAT)
    rb_b = na_change_type(rb_b, NA_DFLOAT);
  b = NA_PTR_TYPE(rb_b, doublereal*);
  if (!NA_IsNArray(rb_c))
    rb_raise(rb_eArgError, "c (8th argument) must be NArray");
  if (NA_RANK(rb_c) != 1)
    rb_raise(rb_eArgError, "rank of c (8th argument) must be %d", 1);
  if (NA_SHAPE0(rb_c) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of c must be the same as shape 0 of ipiv");
  if (NA_TYPE(rb_c) != NA_DFLOAT)
    rb_c = na_change_type(rb_c, NA_DFLOAT);
  c = NA_PTR_TYPE(rb_c, doublereal*);
  equed = StringValueCStr(rb_equed)[0];
  if (!NA_IsNArray(rb_af))
    rb_raise(rb_eArgError, "af (4th argument) must be NArray");
  if (NA_RANK(rb_af) != 2)
    rb_raise(rb_eArgError, "rank of af (4th argument) must be %d", 2);
  if (NA_SHAPE1(rb_af) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of af must be the same as shape 0 of ipiv");
  ldaf = NA_SHAPE0(rb_af);
  if (NA_TYPE(rb_af) != NA_DFLOAT)
    rb_af = na_change_type(rb_af, NA_DFLOAT);
  af = NA_PTR_TYPE(rb_af, doublereal*);
  if (!NA_IsNArray(rb_r))
    rb_raise(rb_eArgError, "r (7th argument) must be NArray");
  if (NA_RANK(rb_r) != 1)
    rb_raise(rb_eArgError, "rank of r (7th argument) must be %d", 1);
  if (NA_SHAPE0(rb_r) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of r must be the same as shape 0 of ipiv");
  if (NA_TYPE(rb_r) != NA_DFLOAT)
    rb_r = na_change_type(rb_r, NA_DFLOAT);
  r = NA_PTR_TYPE(rb_r, doublereal*);
  fact = StringValueCStr(rb_fact)[0];
  ldx = MAX(1,n);
  {
    int shape[2];
    shape[0] = ldx;
    shape[1] = nrhs;
    rb_x = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  x = NA_PTR_TYPE(rb_x, doublereal*);
  {
    int shape[1];
    shape[0] = nrhs;
    rb_ferr = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  ferr = NA_PTR_TYPE(rb_ferr, doublereal*);
  {
    int shape[1];
    shape[0] = nrhs;
    rb_berr = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  berr = NA_PTR_TYPE(rb_berr, doublereal*);
  {
    int shape[1];
    shape[0] = 4*n;
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
    shape[0] = ldaf;
    shape[1] = n;
    rb_af_out__ = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  af_out__ = NA_PTR_TYPE(rb_af_out__, doublereal*);
  MEMCPY(af_out__, af, doublereal, NA_TOTAL(rb_af));
  rb_af = rb_af_out__;
  af = af_out__;
  {
    int shape[1];
    shape[0] = n;
    rb_ipiv_out__ = na_make_object(NA_LINT, 1, shape, cNArray);
  }
  ipiv_out__ = NA_PTR_TYPE(rb_ipiv_out__, integer*);
  MEMCPY(ipiv_out__, ipiv, integer, NA_TOTAL(rb_ipiv));
  rb_ipiv = rb_ipiv_out__;
  ipiv = ipiv_out__;
  {
    int shape[1];
    shape[0] = n;
    rb_r_out__ = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  r_out__ = NA_PTR_TYPE(rb_r_out__, doublereal*);
  MEMCPY(r_out__, r, doublereal, NA_TOTAL(rb_r));
  rb_r = rb_r_out__;
  r = r_out__;
  {
    int shape[1];
    shape[0] = n;
    rb_c_out__ = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  c_out__ = NA_PTR_TYPE(rb_c_out__, doublereal*);
  MEMCPY(c_out__, c, doublereal, NA_TOTAL(rb_c));
  rb_c = rb_c_out__;
  c = c_out__;
  {
    int shape[2];
    shape[0] = ldb;
    shape[1] = nrhs;
    rb_b_out__ = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  b_out__ = NA_PTR_TYPE(rb_b_out__, doublereal*);
  MEMCPY(b_out__, b, doublereal, NA_TOTAL(rb_b));
  rb_b = rb_b_out__;
  b = b_out__;
  iwork = ALLOC_N(integer, (n));

  dgesvx_(&fact, &trans, &n, &nrhs, a, &lda, af, &ldaf, ipiv, &equed, r, c, b, &ldb, x, &ldx, &rcond, ferr, berr, work, iwork, &info);

  free(iwork);
  rb_rcond = rb_float_new((double)rcond);
  rb_info = INT2NUM(info);
  rb_equed = rb_str_new(&equed,1);
  return rb_ary_new3(13, rb_x, rb_rcond, rb_ferr, rb_berr, rb_work, rb_info, rb_a, rb_af, rb_ipiv, rb_equed, rb_r, rb_c, rb_b);
}

void
init_lapack_dgesvx(VALUE mLapack){
  rb_define_module_function(mLapack, "dgesvx", rb_dgesvx, -1);
}

#include "rb_lapack.h"

static VALUE
rb_sgerfsx(int argc, VALUE *argv, VALUE self){
  VALUE rb_trans;
  char trans; 
  VALUE rb_equed;
  char equed; 
  VALUE rb_a;
  real *a; 
  VALUE rb_af;
  real *af; 
  VALUE rb_ipiv;
  integer *ipiv; 
  VALUE rb_r;
  real *r; 
  VALUE rb_c;
  real *c; 
  VALUE rb_b;
  real *b; 
  VALUE rb_x;
  real *x; 
  VALUE rb_nparams;
  integer nparams; 
  VALUE rb_params;
  real *params; 
  VALUE rb_rcond;
  real rcond; 
  VALUE rb_berr;
  real *berr; 
  VALUE rb_err_bnds_norm;
  real *err_bnds_norm; 
  VALUE rb_err_bnds_comp;
  real *err_bnds_comp; 
  VALUE rb_info;
  integer info; 
  VALUE rb_r_out__;
  real *r_out__;
  VALUE rb_c_out__;
  real *c_out__;
  VALUE rb_x_out__;
  real *x_out__;
  VALUE rb_params_out__;
  real *params_out__;
  real *work;
  integer *iwork;

  integer lda;
  integer n;
  integer ldaf;
  integer ldb;
  integer nrhs;
  integer ldx;
  integer n_err_bnds;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  rcond, berr, err_bnds_norm, err_bnds_comp, info, r, c, x, params = NumRu::Lapack.sgerfsx( trans, equed, a, af, ipiv, r, c, b, x, nparams, params)\n    or\n  NumRu::Lapack.sgerfsx  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE SGERFSX( TRANS, EQUED, N, NRHS, A, LDA, AF, LDAF, IPIV, R, C, B, LDB, X, LDX, RCOND, BERR, N_ERR_BNDS, ERR_BNDS_NORM, ERR_BNDS_COMP, NPARAMS, PARAMS, WORK, IWORK, INFO )\n\n*     Purpose\n*     =======\n*\n*     SGERFSX improves the computed solution to a system of linear\n*     equations and provides error bounds and backward error estimates\n*     for the solution.  In addition to normwise error bound, the code\n*     provides maximum componentwise error bound if possible.  See\n*     comments for ERR_BNDS_NORM and ERR_BNDS_COMP for details of the\n*     error bounds.\n*\n*     The original system of linear equations may have been equilibrated\n*     before calling this routine, as described by arguments EQUED, R\n*     and C below. In this case, the solution and error bounds returned\n*     are for the original unequilibrated system.\n*\n\n*     Arguments\n*     =========\n*\n*     Some optional parameters are bundled in the PARAMS array.  These\n*     settings determine how refinement is performed, but often the\n*     defaults are acceptable.  If the defaults are acceptable, users\n*     can pass NPARAMS = 0 which prevents the source code from accessing\n*     the PARAMS argument.\n*\n*     TRANS   (input) CHARACTER*1\n*     Specifies the form of the system of equations:\n*       = 'N':  A * X = B     (No transpose)\n*       = 'T':  A**T * X = B  (Transpose)\n*       = 'C':  A**H * X = B  (Conjugate transpose = Transpose)\n*\n*     EQUED   (input) CHARACTER*1\n*     Specifies the form of equilibration that was done to A\n*     before calling this routine. This is needed to compute\n*     the solution and error bounds correctly.\n*       = 'N':  No equilibration\n*       = 'R':  Row equilibration, i.e., A has been premultiplied by\n*               diag(R).\n*       = 'C':  Column equilibration, i.e., A has been postmultiplied\n*               by diag(C).\n*       = 'B':  Both row and column equilibration, i.e., A has been\n*               replaced by diag(R) * A * diag(C).\n*               The right hand side B has been changed accordingly.\n*\n*     N       (input) INTEGER\n*     The order of the matrix A.  N >= 0.\n*\n*     NRHS    (input) INTEGER\n*     The number of right hand sides, i.e., the number of columns\n*     of the matrices B and X.  NRHS >= 0.\n*\n*     A       (input) REAL array, dimension (LDA,N)\n*     The original N-by-N matrix A.\n*\n*     LDA     (input) INTEGER\n*     The leading dimension of the array A.  LDA >= max(1,N).\n*\n*     AF      (input) REAL array, dimension (LDAF,N)\n*     The factors L and U from the factorization A = P*L*U\n*     as computed by SGETRF.\n*\n*     LDAF    (input) INTEGER\n*     The leading dimension of the array AF.  LDAF >= max(1,N).\n*\n*     IPIV    (input) INTEGER array, dimension (N)\n*     The pivot indices from SGETRF; for 1<=i<=N, row i of the\n*     matrix was interchanged with row IPIV(i).\n*\n*     R       (input or output) REAL array, dimension (N)\n*     The row scale factors for A.  If EQUED = 'R' or 'B', A is\n*     multiplied on the left by diag(R); if EQUED = 'N' or 'C', R\n*     is not accessed.  R is an input argument if FACT = 'F';\n*     otherwise, R is an output argument.  If FACT = 'F' and\n*     EQUED = 'R' or 'B', each element of R must be positive.\n*     If R is output, each element of R is a power of the radix.\n*     If R is input, each element of R should be a power of the radix\n*     to ensure a reliable solution and error estimates. Scaling by\n*     powers of the radix does not cause rounding errors unless the\n*     result underflows or overflows. Rounding errors during scaling\n*     lead to refining with a matrix that is not equivalent to the\n*     input matrix, producing error estimates that may not be\n*     reliable.\n*\n*     C       (input or output) REAL array, dimension (N)\n*     The column scale factors for A.  If EQUED = 'C' or 'B', A is\n*     multiplied on the right by diag(C); if EQUED = 'N' or 'R', C\n*     is not accessed.  C is an input argument if FACT = 'F';\n*     otherwise, C is an output argument.  If FACT = 'F' and\n*     EQUED = 'C' or 'B', each element of C must be positive.\n*     If C is output, each element of C is a power of the radix.\n*     If C is input, each element of C should be a power of the radix\n*     to ensure a reliable solution and error estimates. Scaling by\n*     powers of the radix does not cause rounding errors unless the\n*     result underflows or overflows. Rounding errors during scaling\n*     lead to refining with a matrix that is not equivalent to the\n*     input matrix, producing error estimates that may not be\n*     reliable.\n*\n*     B       (input) REAL array, dimension (LDB,NRHS)\n*     The right hand side matrix B.\n*\n*     LDB     (input) INTEGER\n*     The leading dimension of the array B.  LDB >= max(1,N).\n*\n*     X       (input/output) REAL array, dimension (LDX,NRHS)\n*     On entry, the solution matrix X, as computed by SGETRS.\n*     On exit, the improved solution matrix X.\n*\n*     LDX     (input) INTEGER\n*     The leading dimension of the array X.  LDX >= max(1,N).\n*\n*     RCOND   (output) REAL\n*     Reciprocal scaled condition number.  This is an estimate of the\n*     reciprocal Skeel condition number of the matrix A after\n*     equilibration (if done).  If this is less than the machine\n*     precision (in particular, if it is zero), the matrix is singular\n*     to working precision.  Note that the error may still be small even\n*     if this number is very small and the matrix appears ill-\n*     conditioned.\n*\n*     BERR    (output) REAL array, dimension (NRHS)\n*     Componentwise relative backward error.  This is the\n*     componentwise relative backward error of each solution vector X(j)\n*     (i.e., the smallest relative change in any element of A or B that\n*     makes X(j) an exact solution).\n*\n*     N_ERR_BNDS (input) INTEGER\n*     Number of error bounds to return for each right hand side\n*     and each type (normwise or componentwise).  See ERR_BNDS_NORM and\n*     ERR_BNDS_COMP below.\n*\n*     ERR_BNDS_NORM  (output) REAL array, dimension (NRHS, N_ERR_BNDS)\n*     For each right-hand side, this array contains information about\n*     various error bounds and condition numbers corresponding to the\n*     normwise relative error, which is defined as follows:\n*\n*     Normwise relative error in the ith solution vector:\n*             max_j (abs(XTRUE(j,i) - X(j,i)))\n*            ------------------------------\n*                  max_j abs(X(j,i))\n*\n*     The array is indexed by the type of error information as described\n*     below. There currently are up to three pieces of information\n*     returned.\n*\n*     The first index in ERR_BNDS_NORM(i,:) corresponds to the ith\n*     right-hand side.\n*\n*     The second index in ERR_BNDS_NORM(:,err) contains the following\n*     three fields:\n*     err = 1 \"Trust/don't trust\" boolean. Trust the answer if the\n*              reciprocal condition number is less than the threshold\n*              sqrt(n) * slamch('Epsilon').\n*\n*     err = 2 \"Guaranteed\" error bound: The estimated forward error,\n*              almost certainly within a factor of 10 of the true error\n*              so long as the next entry is greater than the threshold\n*              sqrt(n) * slamch('Epsilon'). This error bound should only\n*              be trusted if the previous boolean is true.\n*\n*     err = 3  Reciprocal condition number: Estimated normwise\n*              reciprocal condition number.  Compared with the threshold\n*              sqrt(n) * slamch('Epsilon') to determine if the error\n*              estimate is \"guaranteed\". These reciprocal condition\n*              numbers are 1 / (norm(Z^{-1},inf) * norm(Z,inf)) for some\n*              appropriately scaled matrix Z.\n*              Let Z = S*A, where S scales each row by a power of the\n*              radix so all absolute row sums of Z are approximately 1.\n*\n*     See Lapack Working Note 165 for further details and extra\n*     cautions.\n*\n*     ERR_BNDS_COMP  (output) REAL array, dimension (NRHS, N_ERR_BNDS)\n*     For each right-hand side, this array contains information about\n*     various error bounds and condition numbers corresponding to the\n*     componentwise relative error, which is defined as follows:\n*\n*     Componentwise relative error in the ith solution vector:\n*                    abs(XTRUE(j,i) - X(j,i))\n*             max_j ----------------------\n*                         abs(X(j,i))\n*\n*     The array is indexed by the right-hand side i (on which the\n*     componentwise relative error depends), and the type of error\n*     information as described below. There currently are up to three\n*     pieces of information returned for each right-hand side. If\n*     componentwise accuracy is not requested (PARAMS(3) = 0.0), then\n*     ERR_BNDS_COMP is not accessed.  If N_ERR_BNDS .LT. 3, then at most\n*     the first (:,N_ERR_BNDS) entries are returned.\n*\n*     The first index in ERR_BNDS_COMP(i,:) corresponds to the ith\n*     right-hand side.\n*\n*     The second index in ERR_BNDS_COMP(:,err) contains the following\n*     three fields:\n*     err = 1 \"Trust/don't trust\" boolean. Trust the answer if the\n*              reciprocal condition number is less than the threshold\n*              sqrt(n) * slamch('Epsilon').\n*\n*     err = 2 \"Guaranteed\" error bound: The estimated forward error,\n*              almost certainly within a factor of 10 of the true error\n*              so long as the next entry is greater than the threshold\n*              sqrt(n) * slamch('Epsilon'). This error bound should only\n*              be trusted if the previous boolean is true.\n*\n*     err = 3  Reciprocal condition number: Estimated componentwise\n*              reciprocal condition number.  Compared with the threshold\n*              sqrt(n) * slamch('Epsilon') to determine if the error\n*              estimate is \"guaranteed\". These reciprocal condition\n*              numbers are 1 / (norm(Z^{-1},inf) * norm(Z,inf)) for some\n*              appropriately scaled matrix Z.\n*              Let Z = S*(A*diag(x)), where x is the solution for the\n*              current right-hand side and S scales each row of\n*              A*diag(x) by a power of the radix so all absolute row\n*              sums of Z are approximately 1.\n*\n*     See Lapack Working Note 165 for further details and extra\n*     cautions.\n*\n*     NPARAMS (input) INTEGER\n*     Specifies the number of parameters set in PARAMS.  If .LE. 0, the\n*     PARAMS array is never referenced and default values are used.\n*\n*     PARAMS  (input / output) REAL array, dimension NPARAMS\n*     Specifies algorithm parameters.  If an entry is .LT. 0.0, then\n*     that entry will be filled with default value used for that\n*     parameter.  Only positions up to NPARAMS are accessed; defaults\n*     are used for higher-numbered parameters.\n*\n*       PARAMS(LA_LINRX_ITREF_I = 1) : Whether to perform iterative\n*            refinement or not.\n*         Default: 1.0\n*            = 0.0 : No refinement is performed, and no error bounds are\n*                    computed.\n*            = 1.0 : Use the double-precision refinement algorithm,\n*                    possibly with doubled-single computations if the\n*                    compilation environment does not support DOUBLE\n*                    PRECISION.\n*              (other values are reserved for future use)\n*\n*       PARAMS(LA_LINRX_ITHRESH_I = 2) : Maximum number of residual\n*            computations allowed for refinement.\n*         Default: 10\n*         Aggressive: Set to 100 to permit convergence using approximate\n*                     factorizations or factorizations other than LU. If\n*                     the factorization uses a technique other than\n*                     Gaussian elimination, the guarantees in\n*                     err_bnds_norm and err_bnds_comp may no longer be\n*                     trustworthy.\n*\n*       PARAMS(LA_LINRX_CWISE_I = 3) : Flag determining if the code\n*            will attempt to find a solution with small componentwise\n*            relative error in the double-precision algorithm.  Positive\n*            is true, 0.0 is false.\n*         Default: 1.0 (attempt componentwise convergence)\n*\n*     WORK    (workspace) REAL array, dimension (4*N)\n*\n*     IWORK   (workspace) INTEGER array, dimension (N)\n*\n*     INFO    (output) INTEGER\n*       = 0:  Successful exit. The solution to every right-hand side is\n*         guaranteed.\n*       < 0:  If INFO = -i, the i-th argument had an illegal value\n*       > 0 and <= N:  U(INFO,INFO) is exactly zero.  The factorization\n*         has been completed, but the factor U is exactly singular, so\n*         the solution and error bounds could not be computed. RCOND = 0\n*         is returned.\n*       = N+J: The solution corresponding to the Jth right-hand side is\n*         not guaranteed. The solutions corresponding to other right-\n*         hand sides K with K > J may not be guaranteed as well, but\n*         only the first such right-hand side is reported. If a small\n*         componentwise error is not requested (PARAMS(3) = 0.0) then\n*         the Jth right-hand side is the first with a normwise error\n*         bound that is not guaranteed (the smallest J such\n*         that ERR_BNDS_NORM(J,1) = 0.0). By default (PARAMS(3) = 1.0)\n*         the Jth right-hand side is the first with either a normwise or\n*         componentwise error bound that is not guaranteed (the smallest\n*         J such that either ERR_BNDS_NORM(J,1) = 0.0 or\n*         ERR_BNDS_COMP(J,1) = 0.0). See the definition of\n*         ERR_BNDS_NORM(:,1) and ERR_BNDS_COMP(:,1). To get information\n*         about all of the right-hand sides check ERR_BNDS_NORM or\n*         ERR_BNDS_COMP.\n*\n\n*     ==================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 11)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 11)", argc);
  rb_trans = argv[0];
  rb_equed = argv[1];
  rb_a = argv[2];
  rb_af = argv[3];
  rb_ipiv = argv[4];
  rb_r = argv[5];
  rb_c = argv[6];
  rb_b = argv[7];
  rb_x = argv[8];
  rb_nparams = argv[9];
  rb_params = argv[10];

  trans = StringValueCStr(rb_trans)[0];
  equed = StringValueCStr(rb_equed)[0];
  nparams = NUM2INT(rb_nparams);
  if (!NA_IsNArray(rb_a))
    rb_raise(rb_eArgError, "a (3th argument) must be NArray");
  if (NA_RANK(rb_a) != 2)
    rb_raise(rb_eArgError, "rank of a (3th argument) must be %d", 2);
  lda = NA_SHAPE0(rb_a);
  n = NA_SHAPE1(rb_a);
  if (NA_TYPE(rb_a) != NA_SFLOAT)
    rb_a = na_change_type(rb_a, NA_SFLOAT);
  a = NA_PTR_TYPE(rb_a, real*);
  if (!NA_IsNArray(rb_af))
    rb_raise(rb_eArgError, "af (4th argument) must be NArray");
  if (NA_RANK(rb_af) != 2)
    rb_raise(rb_eArgError, "rank of af (4th argument) must be %d", 2);
  ldaf = NA_SHAPE0(rb_af);
  if (NA_SHAPE1(rb_af) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of af must be the same as shape 1 of a");
  if (NA_TYPE(rb_af) != NA_SFLOAT)
    rb_af = na_change_type(rb_af, NA_SFLOAT);
  af = NA_PTR_TYPE(rb_af, real*);
  if (!NA_IsNArray(rb_ipiv))
    rb_raise(rb_eArgError, "ipiv (5th argument) must be NArray");
  if (NA_RANK(rb_ipiv) != 1)
    rb_raise(rb_eArgError, "rank of ipiv (5th argument) must be %d", 1);
  if (NA_SHAPE0(rb_ipiv) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of ipiv must be the same as shape 1 of a");
  if (NA_TYPE(rb_ipiv) != NA_LINT)
    rb_ipiv = na_change_type(rb_ipiv, NA_LINT);
  ipiv = NA_PTR_TYPE(rb_ipiv, integer*);
  if (!NA_IsNArray(rb_r))
    rb_raise(rb_eArgError, "r (6th argument) must be NArray");
  if (NA_RANK(rb_r) != 1)
    rb_raise(rb_eArgError, "rank of r (6th argument) must be %d", 1);
  if (NA_SHAPE0(rb_r) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of r must be the same as shape 1 of a");
  if (NA_TYPE(rb_r) != NA_SFLOAT)
    rb_r = na_change_type(rb_r, NA_SFLOAT);
  r = NA_PTR_TYPE(rb_r, real*);
  if (!NA_IsNArray(rb_c))
    rb_raise(rb_eArgError, "c (7th argument) must be NArray");
  if (NA_RANK(rb_c) != 1)
    rb_raise(rb_eArgError, "rank of c (7th argument) must be %d", 1);
  if (NA_SHAPE0(rb_c) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of c must be the same as shape 1 of a");
  if (NA_TYPE(rb_c) != NA_SFLOAT)
    rb_c = na_change_type(rb_c, NA_SFLOAT);
  c = NA_PTR_TYPE(rb_c, real*);
  if (!NA_IsNArray(rb_b))
    rb_raise(rb_eArgError, "b (8th argument) must be NArray");
  if (NA_RANK(rb_b) != 2)
    rb_raise(rb_eArgError, "rank of b (8th argument) must be %d", 2);
  ldb = NA_SHAPE0(rb_b);
  nrhs = NA_SHAPE1(rb_b);
  if (NA_TYPE(rb_b) != NA_SFLOAT)
    rb_b = na_change_type(rb_b, NA_SFLOAT);
  b = NA_PTR_TYPE(rb_b, real*);
  if (!NA_IsNArray(rb_x))
    rb_raise(rb_eArgError, "x (9th argument) must be NArray");
  if (NA_RANK(rb_x) != 2)
    rb_raise(rb_eArgError, "rank of x (9th argument) must be %d", 2);
  ldx = NA_SHAPE0(rb_x);
  if (NA_SHAPE1(rb_x) != nrhs)
    rb_raise(rb_eRuntimeError, "shape 1 of x must be the same as shape 1 of b");
  if (NA_TYPE(rb_x) != NA_SFLOAT)
    rb_x = na_change_type(rb_x, NA_SFLOAT);
  x = NA_PTR_TYPE(rb_x, real*);
  if (!NA_IsNArray(rb_params))
    rb_raise(rb_eArgError, "params (11th argument) must be NArray");
  if (NA_RANK(rb_params) != 0)
    rb_raise(rb_eArgError, "rank of params (11th argument) must be %d", 0);
  if (NA_TYPE(rb_params) != NA_SFLOAT)
    rb_params = na_change_type(rb_params, NA_SFLOAT);
  params = NA_PTR_TYPE(rb_params, real*);
  {
    int shape[1];
    shape[0] = nrhs;
    rb_berr = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  berr = NA_PTR_TYPE(rb_berr, real*);
  {
    int shape[2];
    shape[0] = nrhs;
    shape[1] = n_err_bnds;
    rb_err_bnds_norm = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  err_bnds_norm = NA_PTR_TYPE(rb_err_bnds_norm, real*);
  {
    int shape[2];
    shape[0] = nrhs;
    shape[1] = n_err_bnds;
    rb_err_bnds_comp = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  err_bnds_comp = NA_PTR_TYPE(rb_err_bnds_comp, real*);
  {
    int shape[1];
    shape[0] = n;
    rb_r_out__ = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  r_out__ = NA_PTR_TYPE(rb_r_out__, real*);
  MEMCPY(r_out__, r, real, NA_TOTAL(rb_r));
  rb_r = rb_r_out__;
  r = r_out__;
  {
    int shape[1];
    shape[0] = n;
    rb_c_out__ = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  c_out__ = NA_PTR_TYPE(rb_c_out__, real*);
  MEMCPY(c_out__, c, real, NA_TOTAL(rb_c));
  rb_c = rb_c_out__;
  c = c_out__;
  {
    int shape[2];
    shape[0] = ldx;
    shape[1] = nrhs;
    rb_x_out__ = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  x_out__ = NA_PTR_TYPE(rb_x_out__, real*);
  MEMCPY(x_out__, x, real, NA_TOTAL(rb_x));
  rb_x = rb_x_out__;
  x = x_out__;
  {
    int shape[0];
    rb_params_out__ = na_make_object(NA_SFLOAT, 0, shape, cNArray);
  }
  params_out__ = NA_PTR_TYPE(rb_params_out__, real*);
  MEMCPY(params_out__, params, real, NA_TOTAL(rb_params));
  rb_params = rb_params_out__;
  params = params_out__;
  work = ALLOC_N(real, (4*n));
  iwork = ALLOC_N(integer, (n));

  sgerfsx_(&trans, &equed, &n, &nrhs, a, &lda, af, &ldaf, ipiv, r, c, b, &ldb, x, &ldx, &rcond, berr, &n_err_bnds, err_bnds_norm, err_bnds_comp, &nparams, params, work, iwork, &info);

  free(work);
  free(iwork);
  rb_rcond = rb_float_new((double)rcond);
  rb_info = INT2NUM(info);
  return rb_ary_new3(9, rb_rcond, rb_berr, rb_err_bnds_norm, rb_err_bnds_comp, rb_info, rb_r, rb_c, rb_x, rb_params);
}

void
init_lapack_sgerfsx(VALUE mLapack){
  rb_define_module_function(mLapack, "sgerfsx", rb_sgerfsx, -1);
}

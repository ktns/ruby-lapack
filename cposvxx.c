#include "rb_lapack.h"

extern VOID cposvxx_(char *fact, char *uplo, integer *n, integer *nrhs, complex *a, integer *lda, complex *af, integer *ldaf, char *equed, real *s, complex *b, integer *ldb, complex *x, integer *ldx, real *rcond, real *rpvgrw, real *berr, integer *n_err_bnds, real *err_bnds_norm, real *err_bnds_comp, integer *nparams, real *params, complex *work, real *rwork, integer *info);

static VALUE
rb_cposvxx(int argc, VALUE *argv, VALUE self){
  VALUE rb_fact;
  char fact; 
  VALUE rb_uplo;
  char uplo; 
  VALUE rb_a;
  complex *a; 
  VALUE rb_af;
  complex *af; 
  VALUE rb_equed;
  char equed; 
  VALUE rb_s;
  real *s; 
  VALUE rb_b;
  complex *b; 
  VALUE rb_params;
  real *params; 
  VALUE rb_x;
  complex *x; 
  VALUE rb_rcond;
  real rcond; 
  VALUE rb_rpvgrw;
  real rpvgrw; 
  VALUE rb_berr;
  real *berr; 
  VALUE rb_err_bnds_norm;
  real *err_bnds_norm; 
  VALUE rb_err_bnds_comp;
  real *err_bnds_comp; 
  VALUE rb_info;
  integer info; 
  VALUE rb_a_out__;
  complex *a_out__;
  VALUE rb_af_out__;
  complex *af_out__;
  VALUE rb_s_out__;
  real *s_out__;
  VALUE rb_b_out__;
  complex *b_out__;
  VALUE rb_params_out__;
  real *params_out__;
  complex *work;
  real *rwork;

  integer lda;
  integer n;
  integer ldaf;
  integer ldb;
  integer nrhs;
  integer nparams;
  integer ldx;
  integer n_err_bnds;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  x, rcond, rpvgrw, berr, err_bnds_norm, err_bnds_comp, info, a, af, equed, s, b, params = NumRu::Lapack.cposvxx( fact, uplo, a, af, equed, s, b, params)\n    or\n  NumRu::Lapack.cposvxx  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE CPOSVXX( FACT, UPLO, N, NRHS, A, LDA, AF, LDAF, EQUED, S, B, LDB, X, LDX, RCOND, RPVGRW, BERR, N_ERR_BNDS, ERR_BNDS_NORM, ERR_BNDS_COMP, NPARAMS, PARAMS, WORK, RWORK, INFO )\n\n*     Purpose\n*     =======\n*\n*     CPOSVXX uses the Cholesky factorization A = U**T*U or A = L*L**T\n*     to compute the solution to a complex system of linear equations\n*     A * X = B, where A is an N-by-N symmetric positive definite matrix\n*     and X and B are N-by-NRHS matrices.\n*\n*     If requested, both normwise and maximum componentwise error bounds\n*     are returned. CPOSVXX will return a solution with a tiny\n*     guaranteed error (O(eps) where eps is the working machine\n*     precision) unless the matrix is very ill-conditioned, in which\n*     case a warning is returned. Relevant condition numbers also are\n*     calculated and returned.\n*\n*     CPOSVXX accepts user-provided factorizations and equilibration\n*     factors; see the definitions of the FACT and EQUED options.\n*     Solving with refinement and using a factorization from a previous\n*     CPOSVXX call will also produce a solution with either O(eps)\n*     errors or warnings, but we cannot make that claim for general\n*     user-provided factorizations and equilibration factors if they\n*     differ from what CPOSVXX would itself produce.\n*\n*     Description\n*     ===========\n*\n*     The following steps are performed:\n*\n*     1. If FACT = 'E', real scaling factors are computed to equilibrate\n*     the system:\n*\n*       diag(S)*A*diag(S)     *inv(diag(S))*X = diag(S)*B\n*\n*     Whether or not the system will be equilibrated depends on the\n*     scaling of the matrix A, but if equilibration is used, A is\n*     overwritten by diag(S)*A*diag(S) and B by diag(S)*B.\n*\n*     2. If FACT = 'N' or 'E', the Cholesky decomposition is used to\n*     factor the matrix A (after equilibration if FACT = 'E') as\n*        A = U**T* U,  if UPLO = 'U', or\n*        A = L * L**T,  if UPLO = 'L',\n*     where U is an upper triangular matrix and L is a lower triangular\n*     matrix.\n*\n*     3. If the leading i-by-i principal minor is not positive definite,\n*     then the routine returns with INFO = i. Otherwise, the factored\n*     form of A is used to estimate the condition number of the matrix\n*     A (see argument RCOND).  If the reciprocal of the condition number\n*     is less than machine precision, the routine still goes on to solve\n*     for X and compute error bounds as described below.\n*\n*     4. The system of equations is solved for X using the factored form\n*     of A.\n*\n*     5. By default (unless PARAMS(LA_LINRX_ITREF_I) is set to zero),\n*     the routine will use iterative refinement to try to get a small\n*     error and error bounds.  Refinement calculates the residual to at\n*     least twice the working precision.\n*\n*     6. If equilibration was used, the matrix X is premultiplied by\n*     diag(S) so that it solves the original system before\n*     equilibration.\n*\n\n*     Arguments\n*     =========\n*\n*     Some optional parameters are bundled in the PARAMS array.  These\n*     settings determine how refinement is performed, but often the\n*     defaults are acceptable.  If the defaults are acceptable, users\n*     can pass NPARAMS = 0 which prevents the source code from accessing\n*     the PARAMS argument.\n*\n*     FACT    (input) CHARACTER*1\n*     Specifies whether or not the factored form of the matrix A is\n*     supplied on entry, and if not, whether the matrix A should be\n*     equilibrated before it is factored.\n*       = 'F':  On entry, AF contains the factored form of A.\n*               If EQUED is not 'N', the matrix A has been\n*               equilibrated with scaling factors given by S.\n*               A and AF are not modified.\n*       = 'N':  The matrix A will be copied to AF and factored.\n*       = 'E':  The matrix A will be equilibrated if necessary, then\n*               copied to AF and factored.\n*\n*     UPLO    (input) CHARACTER*1\n*       = 'U':  Upper triangle of A is stored;\n*       = 'L':  Lower triangle of A is stored.\n*\n*     N       (input) INTEGER\n*     The number of linear equations, i.e., the order of the\n*     matrix A.  N >= 0.\n*\n*     NRHS    (input) INTEGER\n*     The number of right hand sides, i.e., the number of columns\n*     of the matrices B and X.  NRHS >= 0.\n*\n*     A       (input/output) COMPLEX array, dimension (LDA,N)\n*     On entry, the symmetric matrix A, except if FACT = 'F' and EQUED =\n*     'Y', then A must contain the equilibrated matrix\n*     diag(S)*A*diag(S).  If UPLO = 'U', the leading N-by-N upper\n*     triangular part of A contains the upper triangular part of the\n*     matrix A, and the strictly lower triangular part of A is not\n*     referenced.  If UPLO = 'L', the leading N-by-N lower triangular\n*     part of A contains the lower triangular part of the matrix A, and\n*     the strictly upper triangular part of A is not referenced.  A is\n*     not modified if FACT = 'F' or 'N', or if FACT = 'E' and EQUED =\n*     'N' on exit.\n*\n*     On exit, if FACT = 'E' and EQUED = 'Y', A is overwritten by\n*     diag(S)*A*diag(S).\n*\n*     LDA     (input) INTEGER\n*     The leading dimension of the array A.  LDA >= max(1,N).\n*\n*     AF      (input or output) COMPLEX array, dimension (LDAF,N)\n*     If FACT = 'F', then AF is an input argument and on entry\n*     contains the triangular factor U or L from the Cholesky\n*     factorization A = U**T*U or A = L*L**T, in the same storage\n*     format as A.  If EQUED .ne. 'N', then AF is the factored\n*     form of the equilibrated matrix diag(S)*A*diag(S).\n*\n*     If FACT = 'N', then AF is an output argument and on exit\n*     returns the triangular factor U or L from the Cholesky\n*     factorization A = U**T*U or A = L*L**T of the original\n*     matrix A.\n*\n*     If FACT = 'E', then AF is an output argument and on exit\n*     returns the triangular factor U or L from the Cholesky\n*     factorization A = U**T*U or A = L*L**T of the equilibrated\n*     matrix A (see the description of A for the form of the\n*     equilibrated matrix).\n*\n*     LDAF    (input) INTEGER\n*     The leading dimension of the array AF.  LDAF >= max(1,N).\n*\n*     EQUED   (input or output) CHARACTER*1\n*     Specifies the form of equilibration that was done.\n*       = 'N':  No equilibration (always true if FACT = 'N').\n*       = 'Y':  Both row and column equilibration, i.e., A has been\n*               replaced by diag(S) * A * diag(S).\n*     EQUED is an input argument if FACT = 'F'; otherwise, it is an\n*     output argument.\n*\n*     S       (input or output) REAL array, dimension (N)\n*     The row scale factors for A.  If EQUED = 'Y', A is multiplied on\n*     the left and right by diag(S).  S is an input argument if FACT =\n*     'F'; otherwise, S is an output argument.  If FACT = 'F' and EQUED\n*     = 'Y', each element of S must be positive.  If S is output, each\n*     element of S is a power of the radix. If S is input, each element\n*     of S should be a power of the radix to ensure a reliable solution\n*     and error estimates. Scaling by powers of the radix does not cause\n*     rounding errors unless the result underflows or overflows.\n*     Rounding errors during scaling lead to refining with a matrix that\n*     is not equivalent to the input matrix, producing error estimates\n*     that may not be reliable.\n*\n*     B       (input/output) COMPLEX array, dimension (LDB,NRHS)\n*     On entry, the N-by-NRHS right hand side matrix B.\n*     On exit,\n*     if EQUED = 'N', B is not modified;\n*     if EQUED = 'Y', B is overwritten by diag(S)*B;\n*\n*     LDB     (input) INTEGER\n*     The leading dimension of the array B.  LDB >= max(1,N).\n*\n*     X       (output) COMPLEX array, dimension (LDX,NRHS)\n*     If INFO = 0, the N-by-NRHS solution matrix X to the original\n*     system of equations.  Note that A and B are modified on exit if\n*     EQUED .ne. 'N', and the solution to the equilibrated system is\n*     inv(diag(S))*X.\n*\n*     LDX     (input) INTEGER\n*     The leading dimension of the array X.  LDX >= max(1,N).\n*\n*     RCOND   (output) REAL\n*     Reciprocal scaled condition number.  This is an estimate of the\n*     reciprocal Skeel condition number of the matrix A after\n*     equilibration (if done).  If this is less than the machine\n*     precision (in particular, if it is zero), the matrix is singular\n*     to working precision.  Note that the error may still be small even\n*     if this number is very small and the matrix appears ill-\n*     conditioned.\n*\n*     RPVGRW  (output) REAL\n*     Reciprocal pivot growth.  On exit, this contains the reciprocal\n*     pivot growth factor norm(A)/norm(U). The \"max absolute element\"\n*     norm is used.  If this is much less than 1, then the stability of\n*     the LU factorization of the (equilibrated) matrix A could be poor.\n*     This also means that the solution X, estimated condition numbers,\n*     and error bounds could be unreliable. If factorization fails with\n*     0<INFO<=N, then this contains the reciprocal pivot growth factor\n*     for the leading INFO columns of A.\n*\n*     BERR    (output) REAL array, dimension (NRHS)\n*     Componentwise relative backward error.  This is the\n*     componentwise relative backward error of each solution vector X(j)\n*     (i.e., the smallest relative change in any element of A or B that\n*     makes X(j) an exact solution).\n*\n*     N_ERR_BNDS (input) INTEGER\n*     Number of error bounds to return for each right hand side\n*     and each type (normwise or componentwise).  See ERR_BNDS_NORM and\n*     ERR_BNDS_COMP below.\n*\n*     ERR_BNDS_NORM  (output) REAL array, dimension (NRHS, N_ERR_BNDS)\n*     For each right-hand side, this array contains information about\n*     various error bounds and condition numbers corresponding to the\n*     normwise relative error, which is defined as follows:\n*\n*     Normwise relative error in the ith solution vector:\n*             max_j (abs(XTRUE(j,i) - X(j,i)))\n*            ------------------------------\n*                  max_j abs(X(j,i))\n*\n*     The array is indexed by the type of error information as described\n*     below. There currently are up to three pieces of information\n*     returned.\n*\n*     The first index in ERR_BNDS_NORM(i,:) corresponds to the ith\n*     right-hand side.\n*\n*     The second index in ERR_BNDS_NORM(:,err) contains the following\n*     three fields:\n*     err = 1 \"Trust/don't trust\" boolean. Trust the answer if the\n*              reciprocal condition number is less than the threshold\n*              sqrt(n) * slamch('Epsilon').\n*\n*     err = 2 \"Guaranteed\" error bound: The estimated forward error,\n*              almost certainly within a factor of 10 of the true error\n*              so long as the next entry is greater than the threshold\n*              sqrt(n) * slamch('Epsilon'). This error bound should only\n*              be trusted if the previous boolean is true.\n*\n*     err = 3  Reciprocal condition number: Estimated normwise\n*              reciprocal condition number.  Compared with the threshold\n*              sqrt(n) * slamch('Epsilon') to determine if the error\n*              estimate is \"guaranteed\". These reciprocal condition\n*              numbers are 1 / (norm(Z^{-1},inf) * norm(Z,inf)) for some\n*              appropriately scaled matrix Z.\n*              Let Z = S*A, where S scales each row by a power of the\n*              radix so all absolute row sums of Z are approximately 1.\n*\n*     See Lapack Working Note 165 for further details and extra\n*     cautions.\n*\n*     ERR_BNDS_COMP  (output) REAL array, dimension (NRHS, N_ERR_BNDS)\n*     For each right-hand side, this array contains information about\n*     various error bounds and condition numbers corresponding to the\n*     componentwise relative error, which is defined as follows:\n*\n*     Componentwise relative error in the ith solution vector:\n*                    abs(XTRUE(j,i) - X(j,i))\n*             max_j ----------------------\n*                         abs(X(j,i))\n*\n*     The array is indexed by the right-hand side i (on which the\n*     componentwise relative error depends), and the type of error\n*     information as described below. There currently are up to three\n*     pieces of information returned for each right-hand side. If\n*     componentwise accuracy is not requested (PARAMS(3) = 0.0), then\n*     ERR_BNDS_COMP is not accessed.  If N_ERR_BNDS .LT. 3, then at most\n*     the first (:,N_ERR_BNDS) entries are returned.\n*\n*     The first index in ERR_BNDS_COMP(i,:) corresponds to the ith\n*     right-hand side.\n*\n*     The second index in ERR_BNDS_COMP(:,err) contains the following\n*     three fields:\n*     err = 1 \"Trust/don't trust\" boolean. Trust the answer if the\n*              reciprocal condition number is less than the threshold\n*              sqrt(n) * slamch('Epsilon').\n*\n*     err = 2 \"Guaranteed\" error bound: The estimated forward error,\n*              almost certainly within a factor of 10 of the true error\n*              so long as the next entry is greater than the threshold\n*              sqrt(n) * slamch('Epsilon'). This error bound should only\n*              be trusted if the previous boolean is true.\n*\n*     err = 3  Reciprocal condition number: Estimated componentwise\n*              reciprocal condition number.  Compared with the threshold\n*              sqrt(n) * slamch('Epsilon') to determine if the error\n*              estimate is \"guaranteed\". These reciprocal condition\n*              numbers are 1 / (norm(Z^{-1},inf) * norm(Z,inf)) for some\n*              appropriately scaled matrix Z.\n*              Let Z = S*(A*diag(x)), where x is the solution for the\n*              current right-hand side and S scales each row of\n*              A*diag(x) by a power of the radix so all absolute row\n*              sums of Z are approximately 1.\n*\n*     See Lapack Working Note 165 for further details and extra\n*     cautions.\n*\n*     NPARAMS (input) INTEGER\n*     Specifies the number of parameters set in PARAMS.  If .LE. 0, the\n*     PARAMS array is never referenced and default values are used.\n*\n*     PARAMS  (input / output) REAL array, dimension NPARAMS\n*     Specifies algorithm parameters.  If an entry is .LT. 0.0, then\n*     that entry will be filled with default value used for that\n*     parameter.  Only positions up to NPARAMS are accessed; defaults\n*     are used for higher-numbered parameters.\n*\n*       PARAMS(LA_LINRX_ITREF_I = 1) : Whether to perform iterative\n*            refinement or not.\n*         Default: 1.0\n*            = 0.0 : No refinement is performed, and no error bounds are\n*                    computed.\n*            = 1.0 : Use the double-precision refinement algorithm,\n*                    possibly with doubled-single computations if the\n*                    compilation environment does not support DOUBLE\n*                    PRECISION.\n*              (other values are reserved for future use)\n*\n*       PARAMS(LA_LINRX_ITHRESH_I = 2) : Maximum number of residual\n*            computations allowed for refinement.\n*         Default: 10\n*         Aggressive: Set to 100 to permit convergence using approximate\n*                     factorizations or factorizations other than LU. If\n*                     the factorization uses a technique other than\n*                     Gaussian elimination, the guarantees in\n*                     err_bnds_norm and err_bnds_comp may no longer be\n*                     trustworthy.\n*\n*       PARAMS(LA_LINRX_CWISE_I = 3) : Flag determining if the code\n*            will attempt to find a solution with small componentwise\n*            relative error in the double-precision algorithm.  Positive\n*            is true, 0.0 is false.\n*         Default: 1.0 (attempt componentwise convergence)\n*\n*     WORK    (workspace) COMPLEX array, dimension (2*N)\n*\n*     RWORK   (workspace) REAL array, dimension (2*N)\n*\n*     INFO    (output) INTEGER\n*       = 0:  Successful exit. The solution to every right-hand side is\n*         guaranteed.\n*       < 0:  If INFO = -i, the i-th argument had an illegal value\n*       > 0 and <= N:  U(INFO,INFO) is exactly zero.  The factorization\n*         has been completed, but the factor U is exactly singular, so\n*         the solution and error bounds could not be computed. RCOND = 0\n*         is returned.\n*       = N+J: The solution corresponding to the Jth right-hand side is\n*         not guaranteed. The solutions corresponding to other right-\n*         hand sides K with K > J may not be guaranteed as well, but\n*         only the first such right-hand side is reported. If a small\n*         componentwise error is not requested (PARAMS(3) = 0.0) then\n*         the Jth right-hand side is the first with a normwise error\n*         bound that is not guaranteed (the smallest J such\n*         that ERR_BNDS_NORM(J,1) = 0.0). By default (PARAMS(3) = 1.0)\n*         the Jth right-hand side is the first with either a normwise or\n*         componentwise error bound that is not guaranteed (the smallest\n*         J such that either ERR_BNDS_NORM(J,1) = 0.0 or\n*         ERR_BNDS_COMP(J,1) = 0.0). See the definition of\n*         ERR_BNDS_NORM(:,1) and ERR_BNDS_COMP(:,1). To get information\n*         about all of the right-hand sides check ERR_BNDS_NORM or\n*         ERR_BNDS_COMP.\n*\n\n*     ==================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 8)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 8)", argc);
  rb_fact = argv[0];
  rb_uplo = argv[1];
  rb_a = argv[2];
  rb_af = argv[3];
  rb_equed = argv[4];
  rb_s = argv[5];
  rb_b = argv[6];
  rb_params = argv[7];

  uplo = StringValueCStr(rb_uplo)[0];
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
    rb_raise(rb_eArgError, "b (7th argument) must be NArray");
  if (NA_RANK(rb_b) != 2)
    rb_raise(rb_eArgError, "rank of b (7th argument) must be %d", 2);
  nrhs = NA_SHAPE1(rb_b);
  ldb = NA_SHAPE0(rb_b);
  if (NA_TYPE(rb_b) != NA_SCOMPLEX)
    rb_b = na_change_type(rb_b, NA_SCOMPLEX);
  b = NA_PTR_TYPE(rb_b, complex*);
  fact = StringValueCStr(rb_fact)[0];
  if (!NA_IsNArray(rb_af))
    rb_raise(rb_eArgError, "af (4th argument) must be NArray");
  if (NA_RANK(rb_af) != 2)
    rb_raise(rb_eArgError, "rank of af (4th argument) must be %d", 2);
  if (NA_SHAPE1(rb_af) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of af must be the same as shape 1 of a");
  ldaf = NA_SHAPE0(rb_af);
  if (NA_TYPE(rb_af) != NA_SCOMPLEX)
    rb_af = na_change_type(rb_af, NA_SCOMPLEX);
  af = NA_PTR_TYPE(rb_af, complex*);
  equed = StringValueCStr(rb_equed)[0];
  if (!NA_IsNArray(rb_params))
    rb_raise(rb_eArgError, "params (8th argument) must be NArray");
  if (NA_RANK(rb_params) != 1)
    rb_raise(rb_eArgError, "rank of params (8th argument) must be %d", 1);
  nparams = NA_SHAPE0(rb_params);
  if (NA_TYPE(rb_params) != NA_SFLOAT)
    rb_params = na_change_type(rb_params, NA_SFLOAT);
  params = NA_PTR_TYPE(rb_params, real*);
  if (!NA_IsNArray(rb_s))
    rb_raise(rb_eArgError, "s (6th argument) must be NArray");
  if (NA_RANK(rb_s) != 1)
    rb_raise(rb_eArgError, "rank of s (6th argument) must be %d", 1);
  if (NA_SHAPE0(rb_s) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of s must be the same as shape 1 of a");
  if (NA_TYPE(rb_s) != NA_SFLOAT)
    rb_s = na_change_type(rb_s, NA_SFLOAT);
  s = NA_PTR_TYPE(rb_s, real*);
  n_err_bnds = 3;
  ldx = MAX(1,n);
  {
    int shape[2];
    shape[0] = ldx;
    shape[1] = nrhs;
    rb_x = na_make_object(NA_SCOMPLEX, 2, shape, cNArray);
  }
  x = NA_PTR_TYPE(rb_x, complex*);
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
    shape[0] = ldaf;
    shape[1] = n;
    rb_af_out__ = na_make_object(NA_SCOMPLEX, 2, shape, cNArray);
  }
  af_out__ = NA_PTR_TYPE(rb_af_out__, complex*);
  MEMCPY(af_out__, af, complex, NA_TOTAL(rb_af));
  rb_af = rb_af_out__;
  af = af_out__;
  {
    int shape[1];
    shape[0] = n;
    rb_s_out__ = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  s_out__ = NA_PTR_TYPE(rb_s_out__, real*);
  MEMCPY(s_out__, s, real, NA_TOTAL(rb_s));
  rb_s = rb_s_out__;
  s = s_out__;
  {
    int shape[2];
    shape[0] = ldb;
    shape[1] = nrhs;
    rb_b_out__ = na_make_object(NA_SCOMPLEX, 2, shape, cNArray);
  }
  b_out__ = NA_PTR_TYPE(rb_b_out__, complex*);
  MEMCPY(b_out__, b, complex, NA_TOTAL(rb_b));
  rb_b = rb_b_out__;
  b = b_out__;
  {
    int shape[1];
    shape[0] = nparams;
    rb_params_out__ = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  params_out__ = NA_PTR_TYPE(rb_params_out__, real*);
  MEMCPY(params_out__, params, real, NA_TOTAL(rb_params));
  rb_params = rb_params_out__;
  params = params_out__;
  work = ALLOC_N(complex, (2*n));
  rwork = ALLOC_N(real, (2*n));

  cposvxx_(&fact, &uplo, &n, &nrhs, a, &lda, af, &ldaf, &equed, s, b, &ldb, x, &ldx, &rcond, &rpvgrw, berr, &n_err_bnds, err_bnds_norm, err_bnds_comp, &nparams, params, work, rwork, &info);

  free(work);
  free(rwork);
  rb_rcond = rb_float_new((double)rcond);
  rb_rpvgrw = rb_float_new((double)rpvgrw);
  rb_info = INT2NUM(info);
  rb_equed = rb_str_new(&equed,1);
  return rb_ary_new3(13, rb_x, rb_rcond, rb_rpvgrw, rb_berr, rb_err_bnds_norm, rb_err_bnds_comp, rb_info, rb_a, rb_af, rb_equed, rb_s, rb_b, rb_params);
}

void
init_lapack_cposvxx(VALUE mLapack){
  rb_define_module_function(mLapack, "cposvxx", rb_cposvxx, -1);
}

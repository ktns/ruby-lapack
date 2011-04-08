#include "rb_lapack.h"

extern VOID dsysvxx_(char *fact, char *uplo, integer *n, integer *nrhs, doublereal *a, integer *lda, doublereal *af, integer *ldaf, integer *ipiv, char *equed, doublereal *s, doublereal *b, integer *ldb, doublereal *x, integer *ldx, doublereal *rcond, doublereal *rpvgrw, doublereal *berr, integer *n_err_bnds, doublereal *err_bnds_norm, doublereal *err_bnds_comp, integer *nparams, doublereal *params, doublereal *work, integer *iwork, integer *info);

static VALUE sHelp, sUsage;

static VALUE
rblapack_dsysvxx(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_fact;
  char fact; 
  VALUE rblapack_uplo;
  char uplo; 
  VALUE rblapack_a;
  doublereal *a; 
  VALUE rblapack_af;
  doublereal *af; 
  VALUE rblapack_ipiv;
  integer *ipiv; 
  VALUE rblapack_equed;
  char equed; 
  VALUE rblapack_s;
  doublereal *s; 
  VALUE rblapack_b;
  doublereal *b; 
  VALUE rblapack_params;
  doublereal *params; 
  VALUE rblapack_x;
  doublereal *x; 
  VALUE rblapack_rcond;
  doublereal rcond; 
  VALUE rblapack_rpvgrw;
  doublereal rpvgrw; 
  VALUE rblapack_berr;
  doublereal *berr; 
  VALUE rblapack_err_bnds_norm;
  doublereal *err_bnds_norm; 
  VALUE rblapack_err_bnds_comp;
  doublereal *err_bnds_comp; 
  VALUE rblapack_info;
  integer info; 
  VALUE rblapack_a_out__;
  doublereal *a_out__;
  VALUE rblapack_af_out__;
  doublereal *af_out__;
  VALUE rblapack_ipiv_out__;
  integer *ipiv_out__;
  VALUE rblapack_s_out__;
  doublereal *s_out__;
  VALUE rblapack_b_out__;
  doublereal *b_out__;
  VALUE rblapack_params_out__;
  doublereal *params_out__;
  doublereal *work;
  integer *iwork;

  integer lda;
  integer n;
  integer ldaf;
  integer ldb;
  integer nrhs;
  integer nparams;
  integer ldx;
  integer n_err_bnds;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  x, rcond, rpvgrw, berr, err_bnds_norm, err_bnds_comp, info, a, af, ipiv, equed, s, b, params = NumRu::Lapack.dsysvxx( fact, uplo, a, af, ipiv, equed, s, b, params, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE DSYSVXX( FACT, UPLO, N, NRHS, A, LDA, AF, LDAF, IPIV, EQUED, S, B, LDB, X, LDX, RCOND, RPVGRW, BERR, N_ERR_BNDS, ERR_BNDS_NORM, ERR_BNDS_COMP, NPARAMS, PARAMS, WORK, IWORK, INFO )\n\n*     Purpose\n*     =======\n*\n*     DSYSVXX uses the diagonal pivoting factorization to compute the\n*     solution to a double precision system of linear equations A * X = B, where A\n*     is an N-by-N symmetric matrix and X and B are N-by-NRHS matrices.\n*\n*     If requested, both normwise and maximum componentwise error bounds\n*     are returned. DSYSVXX will return a solution with a tiny\n*     guaranteed error (O(eps) where eps is the working machine\n*     precision) unless the matrix is very ill-conditioned, in which\n*     case a warning is returned. Relevant condition numbers also are\n*     calculated and returned.\n*\n*     DSYSVXX accepts user-provided factorizations and equilibration\n*     factors; see the definitions of the FACT and EQUED options.\n*     Solving with refinement and using a factorization from a previous\n*     DSYSVXX call will also produce a solution with either O(eps)\n*     errors or warnings, but we cannot make that claim for general\n*     user-provided factorizations and equilibration factors if they\n*     differ from what DSYSVXX would itself produce.\n*\n*     Description\n*     ===========\n*\n*     The following steps are performed:\n*\n*     1. If FACT = 'E', double precision scaling factors are computed to equilibrate\n*     the system:\n*\n*       diag(S)*A*diag(S)     *inv(diag(S))*X = diag(S)*B\n*\n*     Whether or not the system will be equilibrated depends on the\n*     scaling of the matrix A, but if equilibration is used, A is\n*     overwritten by diag(S)*A*diag(S) and B by diag(S)*B.\n*\n*     2. If FACT = 'N' or 'E', the LU decomposition is used to factor\n*     the matrix A (after equilibration if FACT = 'E') as\n*\n*        A = U * D * U**T,  if UPLO = 'U', or\n*        A = L * D * L**T,  if UPLO = 'L',\n*\n*     where U (or L) is a product of permutation and unit upper (lower)\n*     triangular matrices, and D is symmetric and block diagonal with\n*     1-by-1 and 2-by-2 diagonal blocks.\n*\n*     3. If some D(i,i)=0, so that D is exactly singular, then the\n*     routine returns with INFO = i. Otherwise, the factored form of A\n*     is used to estimate the condition number of the matrix A (see\n*     argument RCOND).  If the reciprocal of the condition number is\n*     less than machine precision, the routine still goes on to solve\n*     for X and compute error bounds as described below.\n*\n*     4. The system of equations is solved for X using the factored form\n*     of A.\n*\n*     5. By default (unless PARAMS(LA_LINRX_ITREF_I) is set to zero),\n*     the routine will use iterative refinement to try to get a small\n*     error and error bounds.  Refinement calculates the residual to at\n*     least twice the working precision.\n*\n*     6. If equilibration was used, the matrix X is premultiplied by\n*     diag(R) so that it solves the original system before\n*     equilibration.\n*\n\n*     Arguments\n*     =========\n*\n*     Some optional parameters are bundled in the PARAMS array.  These\n*     settings determine how refinement is performed, but often the\n*     defaults are acceptable.  If the defaults are acceptable, users\n*     can pass NPARAMS = 0 which prevents the source code from accessing\n*     the PARAMS argument.\n*\n*     FACT    (input) CHARACTER*1\n*     Specifies whether or not the factored form of the matrix A is\n*     supplied on entry, and if not, whether the matrix A should be\n*     equilibrated before it is factored.\n*       = 'F':  On entry, AF and IPIV contain the factored form of A.\n*               If EQUED is not 'N', the matrix A has been\n*               equilibrated with scaling factors given by S.\n*               A, AF, and IPIV are not modified.\n*       = 'N':  The matrix A will be copied to AF and factored.\n*       = 'E':  The matrix A will be equilibrated if necessary, then\n*               copied to AF and factored.\n*\n*     UPLO    (input) CHARACTER*1\n*       = 'U':  Upper triangle of A is stored;\n*       = 'L':  Lower triangle of A is stored.\n*\n*     N       (input) INTEGER\n*     The number of linear equations, i.e., the order of the\n*     matrix A.  N >= 0.\n*\n*     NRHS    (input) INTEGER\n*     The number of right hand sides, i.e., the number of columns\n*     of the matrices B and X.  NRHS >= 0.\n*\n*     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)\n*     The symmetric matrix A.  If UPLO = 'U', the leading N-by-N\n*     upper triangular part of A contains the upper triangular\n*     part of the matrix A, and the strictly lower triangular\n*     part of A is not referenced.  If UPLO = 'L', the leading\n*     N-by-N lower triangular part of A contains the lower\n*     triangular part of the matrix A, and the strictly upper\n*     triangular part of A is not referenced.\n*\n*     On exit, if FACT = 'E' and EQUED = 'Y', A is overwritten by\n*     diag(S)*A*diag(S).\n*\n*     LDA     (input) INTEGER\n*     The leading dimension of the array A.  LDA >= max(1,N).\n*\n*     AF      (input or output) DOUBLE PRECISION array, dimension (LDAF,N)\n*     If FACT = 'F', then AF is an input argument and on entry\n*     contains the block diagonal matrix D and the multipliers\n*     used to obtain the factor U or L from the factorization A =\n*     U*D*U**T or A = L*D*L**T as computed by DSYTRF.\n*\n*     If FACT = 'N', then AF is an output argument and on exit\n*     returns the block diagonal matrix D and the multipliers\n*     used to obtain the factor U or L from the factorization A =\n*     U*D*U**T or A = L*D*L**T.\n*\n*     LDAF    (input) INTEGER\n*     The leading dimension of the array AF.  LDAF >= max(1,N).\n*\n*     IPIV    (input or output) INTEGER array, dimension (N)\n*     If FACT = 'F', then IPIV is an input argument and on entry\n*     contains details of the interchanges and the block\n*     structure of D, as determined by DSYTRF.  If IPIV(k) > 0,\n*     then rows and columns k and IPIV(k) were interchanged and\n*     D(k,k) is a 1-by-1 diagonal block.  If UPLO = 'U' and\n*     IPIV(k) = IPIV(k-1) < 0, then rows and columns k-1 and\n*     -IPIV(k) were interchanged and D(k-1:k,k-1:k) is a 2-by-2\n*     diagonal block.  If UPLO = 'L' and IPIV(k) = IPIV(k+1) < 0,\n*     then rows and columns k+1 and -IPIV(k) were interchanged\n*     and D(k:k+1,k:k+1) is a 2-by-2 diagonal block.\n*\n*     If FACT = 'N', then IPIV is an output argument and on exit\n*     contains details of the interchanges and the block\n*     structure of D, as determined by DSYTRF.\n*\n*     EQUED   (input or output) CHARACTER*1\n*     Specifies the form of equilibration that was done.\n*       = 'N':  No equilibration (always true if FACT = 'N').\n*       = 'Y':  Both row and column equilibration, i.e., A has been\n*               replaced by diag(S) * A * diag(S).\n*     EQUED is an input argument if FACT = 'F'; otherwise, it is an\n*     output argument.\n*\n*     S       (input or output) DOUBLE PRECISION array, dimension (N)\n*     The scale factors for A.  If EQUED = 'Y', A is multiplied on\n*     the left and right by diag(S).  S is an input argument if FACT =\n*     'F'; otherwise, S is an output argument.  If FACT = 'F' and EQUED\n*     = 'Y', each element of S must be positive.  If S is output, each\n*     element of S is a power of the radix. If S is input, each element\n*     of S should be a power of the radix to ensure a reliable solution\n*     and error estimates. Scaling by powers of the radix does not cause\n*     rounding errors unless the result underflows or overflows.\n*     Rounding errors during scaling lead to refining with a matrix that\n*     is not equivalent to the input matrix, producing error estimates\n*     that may not be reliable.\n*\n*     B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)\n*     On entry, the N-by-NRHS right hand side matrix B.\n*     On exit,\n*     if EQUED = 'N', B is not modified;\n*     if EQUED = 'Y', B is overwritten by diag(S)*B;\n*\n*     LDB     (input) INTEGER\n*     The leading dimension of the array B.  LDB >= max(1,N).\n*\n*     X       (output) DOUBLE PRECISION array, dimension (LDX,NRHS)\n*     If INFO = 0, the N-by-NRHS solution matrix X to the original\n*     system of equations.  Note that A and B are modified on exit if\n*     EQUED .ne. 'N', and the solution to the equilibrated system is\n*     inv(diag(S))*X.\n*\n*     LDX     (input) INTEGER\n*     The leading dimension of the array X.  LDX >= max(1,N).\n*\n*     RCOND   (output) DOUBLE PRECISION\n*     Reciprocal scaled condition number.  This is an estimate of the\n*     reciprocal Skeel condition number of the matrix A after\n*     equilibration (if done).  If this is less than the machine\n*     precision (in particular, if it is zero), the matrix is singular\n*     to working precision.  Note that the error may still be small even\n*     if this number is very small and the matrix appears ill-\n*     conditioned.\n*\n*     RPVGRW  (output) DOUBLE PRECISION\n*     Reciprocal pivot growth.  On exit, this contains the reciprocal\n*     pivot growth factor norm(A)/norm(U). The \"max absolute element\"\n*     norm is used.  If this is much less than 1, then the stability of\n*     the LU factorization of the (equilibrated) matrix A could be poor.\n*     This also means that the solution X, estimated condition numbers,\n*     and error bounds could be unreliable. If factorization fails with\n*     0<INFO<=N, then this contains the reciprocal pivot growth factor\n*     for the leading INFO columns of A.\n*\n*     BERR    (output) DOUBLE PRECISION array, dimension (NRHS)\n*     Componentwise relative backward error.  This is the\n*     componentwise relative backward error of each solution vector X(j)\n*     (i.e., the smallest relative change in any element of A or B that\n*     makes X(j) an exact solution).\n*\n*     N_ERR_BNDS (input) INTEGER\n*     Number of error bounds to return for each right hand side\n*     and each type (normwise or componentwise).  See ERR_BNDS_NORM and\n*     ERR_BNDS_COMP below.\n*\n*     ERR_BNDS_NORM  (output) DOUBLE PRECISION array, dimension (NRHS, N_ERR_BNDS)\n*     For each right-hand side, this array contains information about\n*     various error bounds and condition numbers corresponding to the\n*     normwise relative error, which is defined as follows:\n*\n*     Normwise relative error in the ith solution vector:\n*             max_j (abs(XTRUE(j,i) - X(j,i)))\n*            ------------------------------\n*                  max_j abs(X(j,i))\n*\n*     The array is indexed by the type of error information as described\n*     below. There currently are up to three pieces of information\n*     returned.\n*\n*     The first index in ERR_BNDS_NORM(i,:) corresponds to the ith\n*     right-hand side.\n*\n*     The second index in ERR_BNDS_NORM(:,err) contains the following\n*     three fields:\n*     err = 1 \"Trust/don't trust\" boolean. Trust the answer if the\n*              reciprocal condition number is less than the threshold\n*              sqrt(n) * dlamch('Epsilon').\n*\n*     err = 2 \"Guaranteed\" error bound: The estimated forward error,\n*              almost certainly within a factor of 10 of the true error\n*              so long as the next entry is greater than the threshold\n*              sqrt(n) * dlamch('Epsilon'). This error bound should only\n*              be trusted if the previous boolean is true.\n*\n*     err = 3  Reciprocal condition number: Estimated normwise\n*              reciprocal condition number.  Compared with the threshold\n*              sqrt(n) * dlamch('Epsilon') to determine if the error\n*              estimate is \"guaranteed\". These reciprocal condition\n*              numbers are 1 / (norm(Z^{-1},inf) * norm(Z,inf)) for some\n*              appropriately scaled matrix Z.\n*              Let Z = S*A, where S scales each row by a power of the\n*              radix so all absolute row sums of Z are approximately 1.\n*\n*     See Lapack Working Note 165 for further details and extra\n*     cautions.\n*\n*     ERR_BNDS_COMP  (output) DOUBLE PRECISION array, dimension (NRHS, N_ERR_BNDS)\n*     For each right-hand side, this array contains information about\n*     various error bounds and condition numbers corresponding to the\n*     componentwise relative error, which is defined as follows:\n*\n*     Componentwise relative error in the ith solution vector:\n*                    abs(XTRUE(j,i) - X(j,i))\n*             max_j ----------------------\n*                         abs(X(j,i))\n*\n*     The array is indexed by the right-hand side i (on which the\n*     componentwise relative error depends), and the type of error\n*     information as described below. There currently are up to three\n*     pieces of information returned for each right-hand side. If\n*     componentwise accuracy is not requested (PARAMS(3) = 0.0), then\n*     ERR_BNDS_COMP is not accessed.  If N_ERR_BNDS .LT. 3, then at most\n*     the first (:,N_ERR_BNDS) entries are returned.\n*\n*     The first index in ERR_BNDS_COMP(i,:) corresponds to the ith\n*     right-hand side.\n*\n*     The second index in ERR_BNDS_COMP(:,err) contains the following\n*     three fields:\n*     err = 1 \"Trust/don't trust\" boolean. Trust the answer if the\n*              reciprocal condition number is less than the threshold\n*              sqrt(n) * dlamch('Epsilon').\n*\n*     err = 2 \"Guaranteed\" error bound: The estimated forward error,\n*              almost certainly within a factor of 10 of the true error\n*              so long as the next entry is greater than the threshold\n*              sqrt(n) * dlamch('Epsilon'). This error bound should only\n*              be trusted if the previous boolean is true.\n*\n*     err = 3  Reciprocal condition number: Estimated componentwise\n*              reciprocal condition number.  Compared with the threshold\n*              sqrt(n) * dlamch('Epsilon') to determine if the error\n*              estimate is \"guaranteed\". These reciprocal condition\n*              numbers are 1 / (norm(Z^{-1},inf) * norm(Z,inf)) for some\n*              appropriately scaled matrix Z.\n*              Let Z = S*(A*diag(x)), where x is the solution for the\n*              current right-hand side and S scales each row of\n*              A*diag(x) by a power of the radix so all absolute row\n*              sums of Z are approximately 1.\n*\n*     See Lapack Working Note 165 for further details and extra\n*     cautions.\n*\n*     NPARAMS (input) INTEGER\n*     Specifies the number of parameters set in PARAMS.  If .LE. 0, the\n*     PARAMS array is never referenced and default values are used.\n*\n*     PARAMS  (input / output) DOUBLE PRECISION array, dimension (NPARAMS)\n*     Specifies algorithm parameters.  If an entry is .LT. 0.0, then\n*     that entry will be filled with default value used for that\n*     parameter.  Only positions up to NPARAMS are accessed; defaults\n*     are used for higher-numbered parameters.\n*\n*       PARAMS(LA_LINRX_ITREF_I = 1) : Whether to perform iterative\n*            refinement or not.\n*         Default: 1.0D+0\n*            = 0.0 : No refinement is performed, and no error bounds are\n*                    computed.\n*            = 1.0 : Use the extra-precise refinement algorithm.\n*              (other values are reserved for future use)\n*\n*       PARAMS(LA_LINRX_ITHRESH_I = 2) : Maximum number of residual\n*            computations allowed for refinement.\n*         Default: 10\n*         Aggressive: Set to 100 to permit convergence using approximate\n*                     factorizations or factorizations other than LU. If\n*                     the factorization uses a technique other than\n*                     Gaussian elimination, the guarantees in\n*                     err_bnds_norm and err_bnds_comp may no longer be\n*                     trustworthy.\n*\n*       PARAMS(LA_LINRX_CWISE_I = 3) : Flag determining if the code\n*            will attempt to find a solution with small componentwise\n*            relative error in the double-precision algorithm.  Positive\n*            is true, 0.0 is false.\n*         Default: 1.0 (attempt componentwise convergence)\n*\n*     WORK    (workspace) DOUBLE PRECISION array, dimension (4*N)\n*\n*     IWORK   (workspace) INTEGER array, dimension (N)\n*\n*     INFO    (output) INTEGER\n*       = 0:  Successful exit. The solution to every right-hand side is\n*         guaranteed.\n*       < 0:  If INFO = -i, the i-th argument had an illegal value\n*       > 0 and <= N:  U(INFO,INFO) is exactly zero.  The factorization\n*         has been completed, but the factor U is exactly singular, so\n*         the solution and error bounds could not be computed. RCOND = 0\n*         is returned.\n*       = N+J: The solution corresponding to the Jth right-hand side is\n*         not guaranteed. The solutions corresponding to other right-\n*         hand sides K with K > J may not be guaranteed as well, but\n*         only the first such right-hand side is reported. If a small\n*         componentwise error is not requested (PARAMS(3) = 0.0) then\n*         the Jth right-hand side is the first with a normwise error\n*         bound that is not guaranteed (the smallest J such\n*         that ERR_BNDS_NORM(J,1) = 0.0). By default (PARAMS(3) = 1.0)\n*         the Jth right-hand side is the first with either a normwise or\n*         componentwise error bound that is not guaranteed (the smallest\n*         J such that either ERR_BNDS_NORM(J,1) = 0.0 or\n*         ERR_BNDS_COMP(J,1) = 0.0). See the definition of\n*         ERR_BNDS_NORM(:,1) and ERR_BNDS_COMP(:,1). To get information\n*         about all of the right-hand sides check ERR_BNDS_NORM or\n*         ERR_BNDS_COMP.\n*\n\n*     ==================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  x, rcond, rpvgrw, berr, err_bnds_norm, err_bnds_comp, info, a, af, ipiv, equed, s, b, params = NumRu::Lapack.dsysvxx( fact, uplo, a, af, ipiv, equed, s, b, params, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 9)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 9)", argc);
  rblapack_fact = argv[0];
  rblapack_uplo = argv[1];
  rblapack_a = argv[2];
  rblapack_af = argv[3];
  rblapack_ipiv = argv[4];
  rblapack_equed = argv[5];
  rblapack_s = argv[6];
  rblapack_b = argv[7];
  rblapack_params = argv[8];
  if (rb_options != Qnil) {
  }

  uplo = StringValueCStr(rblapack_uplo)[0];
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
  if (NA_TYPE(rblapack_a) != NA_DFLOAT)
    rblapack_a = na_change_type(rblapack_a, NA_DFLOAT);
  a = NA_PTR_TYPE(rblapack_a, doublereal*);
  if (!NA_IsNArray(rblapack_b))
    rb_raise(rb_eArgError, "b (8th argument) must be NArray");
  if (NA_RANK(rblapack_b) != 2)
    rb_raise(rb_eArgError, "rank of b (8th argument) must be %d", 2);
  nrhs = NA_SHAPE1(rblapack_b);
  ldb = NA_SHAPE0(rblapack_b);
  if (NA_TYPE(rblapack_b) != NA_DFLOAT)
    rblapack_b = na_change_type(rblapack_b, NA_DFLOAT);
  b = NA_PTR_TYPE(rblapack_b, doublereal*);
  fact = StringValueCStr(rblapack_fact)[0];
  equed = StringValueCStr(rblapack_equed)[0];
  n_err_bnds = 3;
  if (!NA_IsNArray(rblapack_params))
    rb_raise(rb_eArgError, "params (9th argument) must be NArray");
  if (NA_RANK(rblapack_params) != 1)
    rb_raise(rb_eArgError, "rank of params (9th argument) must be %d", 1);
  nparams = NA_SHAPE0(rblapack_params);
  if (NA_TYPE(rblapack_params) != NA_DFLOAT)
    rblapack_params = na_change_type(rblapack_params, NA_DFLOAT);
  params = NA_PTR_TYPE(rblapack_params, doublereal*);
  if (!NA_IsNArray(rblapack_s))
    rb_raise(rb_eArgError, "s (7th argument) must be NArray");
  if (NA_RANK(rblapack_s) != 1)
    rb_raise(rb_eArgError, "rank of s (7th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_s) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of s must be the same as shape 0 of ipiv");
  if (NA_TYPE(rblapack_s) != NA_DFLOAT)
    rblapack_s = na_change_type(rblapack_s, NA_DFLOAT);
  s = NA_PTR_TYPE(rblapack_s, doublereal*);
  if (!NA_IsNArray(rblapack_af))
    rb_raise(rb_eArgError, "af (4th argument) must be NArray");
  if (NA_RANK(rblapack_af) != 2)
    rb_raise(rb_eArgError, "rank of af (4th argument) must be %d", 2);
  if (NA_SHAPE1(rblapack_af) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of af must be the same as shape 0 of ipiv");
  ldaf = NA_SHAPE0(rblapack_af);
  if (NA_TYPE(rblapack_af) != NA_DFLOAT)
    rblapack_af = na_change_type(rblapack_af, NA_DFLOAT);
  af = NA_PTR_TYPE(rblapack_af, doublereal*);
  ldx = MAX(1,n);
  {
    int shape[2];
    shape[0] = ldx;
    shape[1] = nrhs;
    rblapack_x = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  x = NA_PTR_TYPE(rblapack_x, doublereal*);
  {
    int shape[1];
    shape[0] = nrhs;
    rblapack_berr = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  berr = NA_PTR_TYPE(rblapack_berr, doublereal*);
  {
    int shape[2];
    shape[0] = nrhs;
    shape[1] = n_err_bnds;
    rblapack_err_bnds_norm = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  err_bnds_norm = NA_PTR_TYPE(rblapack_err_bnds_norm, doublereal*);
  {
    int shape[2];
    shape[0] = nrhs;
    shape[1] = n_err_bnds;
    rblapack_err_bnds_comp = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  err_bnds_comp = NA_PTR_TYPE(rblapack_err_bnds_comp, doublereal*);
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
    shape[0] = ldaf;
    shape[1] = n;
    rblapack_af_out__ = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  af_out__ = NA_PTR_TYPE(rblapack_af_out__, doublereal*);
  MEMCPY(af_out__, af, doublereal, NA_TOTAL(rblapack_af));
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
    rblapack_s_out__ = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  s_out__ = NA_PTR_TYPE(rblapack_s_out__, doublereal*);
  MEMCPY(s_out__, s, doublereal, NA_TOTAL(rblapack_s));
  rblapack_s = rblapack_s_out__;
  s = s_out__;
  {
    int shape[2];
    shape[0] = ldb;
    shape[1] = nrhs;
    rblapack_b_out__ = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  b_out__ = NA_PTR_TYPE(rblapack_b_out__, doublereal*);
  MEMCPY(b_out__, b, doublereal, NA_TOTAL(rblapack_b));
  rblapack_b = rblapack_b_out__;
  b = b_out__;
  {
    int shape[1];
    shape[0] = nparams;
    rblapack_params_out__ = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  params_out__ = NA_PTR_TYPE(rblapack_params_out__, doublereal*);
  MEMCPY(params_out__, params, doublereal, NA_TOTAL(rblapack_params));
  rblapack_params = rblapack_params_out__;
  params = params_out__;
  work = ALLOC_N(doublereal, (4*n));
  iwork = ALLOC_N(integer, (n));

  dsysvxx_(&fact, &uplo, &n, &nrhs, a, &lda, af, &ldaf, ipiv, &equed, s, b, &ldb, x, &ldx, &rcond, &rpvgrw, berr, &n_err_bnds, err_bnds_norm, err_bnds_comp, &nparams, params, work, iwork, &info);

  free(work);
  free(iwork);
  rblapack_rcond = rb_float_new((double)rcond);
  rblapack_rpvgrw = rb_float_new((double)rpvgrw);
  rblapack_info = INT2NUM(info);
  rblapack_equed = rb_str_new(&equed,1);
  return rb_ary_new3(14, rblapack_x, rblapack_rcond, rblapack_rpvgrw, rblapack_berr, rblapack_err_bnds_norm, rblapack_err_bnds_comp, rblapack_info, rblapack_a, rblapack_af, rblapack_ipiv, rblapack_equed, rblapack_s, rblapack_b, rblapack_params);
}

void
init_lapack_dsysvxx(VALUE mLapack){
  rb_define_module_function(mLapack, "dsysvxx", rblapack_dsysvxx, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}

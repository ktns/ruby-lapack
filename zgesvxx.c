#include "rb_lapack.h"

extern VOID zgesvxx_(char *fact, char *trans, integer *n, integer *nrhs, doublecomplex *a, integer *lda, doublecomplex *af, integer *ldaf, integer *ipiv, char *equed, doublereal *r, doublereal *c, doublecomplex *b, integer *ldb, doublecomplex *x, integer *ldx, doublereal *rcond, doublereal *rpvgrw, doublereal *berr, integer *n_err_bnds, doublereal *err_bnds_norm, doublereal *err_bnds_comp, integer *nparams, doublereal *params, doublecomplex *work, doublereal *rwork, integer *info);

static VALUE
rb_zgesvxx(int argc, VALUE *argv, VALUE self){
  VALUE rb_fact;
  char fact; 
  VALUE rb_trans;
  char trans; 
  VALUE rb_a;
  doublecomplex *a; 
  VALUE rb_af;
  doublecomplex *af; 
  VALUE rb_ipiv;
  integer *ipiv; 
  VALUE rb_equed;
  char equed; 
  VALUE rb_r;
  doublereal *r; 
  VALUE rb_c;
  doublereal *c; 
  VALUE rb_b;
  doublecomplex *b; 
  VALUE rb_params;
  doublereal *params; 
  VALUE rb_x;
  doublecomplex *x; 
  VALUE rb_rcond;
  doublereal rcond; 
  VALUE rb_rpvgrw;
  doublereal rpvgrw; 
  VALUE rb_berr;
  doublereal *berr; 
  VALUE rb_err_bnds_norm;
  doublereal *err_bnds_norm; 
  VALUE rb_err_bnds_comp;
  doublereal *err_bnds_comp; 
  VALUE rb_info;
  integer info; 
  VALUE rb_a_out__;
  doublecomplex *a_out__;
  VALUE rb_af_out__;
  doublecomplex *af_out__;
  VALUE rb_ipiv_out__;
  integer *ipiv_out__;
  VALUE rb_r_out__;
  doublereal *r_out__;
  VALUE rb_c_out__;
  doublereal *c_out__;
  VALUE rb_b_out__;
  doublecomplex *b_out__;
  VALUE rb_params_out__;
  doublereal *params_out__;
  doublecomplex *work;
  doublereal *rwork;

  integer lda;
  integer n;
  integer ldaf;
  integer ldb;
  integer nrhs;
  integer nparams;
  integer ldx;
  integer n_err_bnds;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  x, rcond, rpvgrw, berr, err_bnds_norm, err_bnds_comp, info, a, af, ipiv, equed, r, c, b, params = NumRu::Lapack.zgesvxx( fact, trans, a, af, ipiv, equed, r, c, b, params)\n    or\n  NumRu::Lapack.zgesvxx  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE ZGESVXX( FACT, TRANS, N, NRHS, A, LDA, AF, LDAF, IPIV, EQUED, R, C, B, LDB, X, LDX, RCOND, RPVGRW, BERR, N_ERR_BNDS, ERR_BNDS_NORM, ERR_BNDS_COMP, NPARAMS, PARAMS, WORK, RWORK, INFO )\n\n*     Purpose\n*     =======\n*\n*     ZGESVXX uses the LU factorization to compute the solution to a\n*     complex*16 system of linear equations  A * X = B,  where A is an\n*     N-by-N matrix and X and B are N-by-NRHS matrices.\n*\n*     If requested, both normwise and maximum componentwise error bounds\n*     are returned. ZGESVXX will return a solution with a tiny\n*     guaranteed error (O(eps) where eps is the working machine\n*     precision) unless the matrix is very ill-conditioned, in which\n*     case a warning is returned. Relevant condition numbers also are\n*     calculated and returned.\n*\n*     ZGESVXX accepts user-provided factorizations and equilibration\n*     factors; see the definitions of the FACT and EQUED options.\n*     Solving with refinement and using a factorization from a previous\n*     ZGESVXX call will also produce a solution with either O(eps)\n*     errors or warnings, but we cannot make that claim for general\n*     user-provided factorizations and equilibration factors if they\n*     differ from what ZGESVXX would itself produce.\n*\n*     Description\n*     ===========\n*\n*     The following steps are performed:\n*\n*     1. If FACT = 'E', double precision scaling factors are computed to equilibrate\n*     the system:\n*\n*       TRANS = 'N':  diag(R)*A*diag(C)     *inv(diag(C))*X = diag(R)*B\n*       TRANS = 'T': (diag(R)*A*diag(C))**T *inv(diag(R))*X = diag(C)*B\n*       TRANS = 'C': (diag(R)*A*diag(C))**H *inv(diag(R))*X = diag(C)*B\n*\n*     Whether or not the system will be equilibrated depends on the\n*     scaling of the matrix A, but if equilibration is used, A is\n*     overwritten by diag(R)*A*diag(C) and B by diag(R)*B (if TRANS='N')\n*     or diag(C)*B (if TRANS = 'T' or 'C').\n*\n*     2. If FACT = 'N' or 'E', the LU decomposition is used to factor\n*     the matrix A (after equilibration if FACT = 'E') as\n*\n*       A = P * L * U,\n*\n*     where P is a permutation matrix, L is a unit lower triangular\n*     matrix, and U is upper triangular.\n*\n*     3. If some U(i,i)=0, so that U is exactly singular, then the\n*     routine returns with INFO = i. Otherwise, the factored form of A\n*     is used to estimate the condition number of the matrix A (see\n*     argument RCOND). If the reciprocal of the condition number is less\n*     than machine precision, the routine still goes on to solve for X\n*     and compute error bounds as described below.\n*\n*     4. The system of equations is solved for X using the factored form\n*     of A.\n*\n*     5. By default (unless PARAMS(LA_LINRX_ITREF_I) is set to zero),\n*     the routine will use iterative refinement to try to get a small\n*     error and error bounds.  Refinement calculates the residual to at\n*     least twice the working precision.\n*\n*     6. If equilibration was used, the matrix X is premultiplied by\n*     diag(C) (if TRANS = 'N') or diag(R) (if TRANS = 'T' or 'C') so\n*     that it solves the original system before equilibration.\n*\n\n*     Arguments\n*     =========\n*\n*     Some optional parameters are bundled in the PARAMS array.  These\n*     settings determine how refinement is performed, but often the\n*     defaults are acceptable.  If the defaults are acceptable, users\n*     can pass NPARAMS = 0 which prevents the source code from accessing\n*     the PARAMS argument.\n*\n*     FACT    (input) CHARACTER*1\n*     Specifies whether or not the factored form of the matrix A is\n*     supplied on entry, and if not, whether the matrix A should be\n*     equilibrated before it is factored.\n*       = 'F':  On entry, AF and IPIV contain the factored form of A.\n*               If EQUED is not 'N', the matrix A has been\n*               equilibrated with scaling factors given by R and C.\n*               A, AF, and IPIV are not modified.\n*       = 'N':  The matrix A will be copied to AF and factored.\n*       = 'E':  The matrix A will be equilibrated if necessary, then\n*               copied to AF and factored.\n*\n*     TRANS   (input) CHARACTER*1\n*     Specifies the form of the system of equations:\n*       = 'N':  A * X = B     (No transpose)\n*       = 'T':  A**T * X = B  (Transpose)\n*       = 'C':  A**H * X = B  (Conjugate Transpose)\n*\n*     N       (input) INTEGER\n*     The number of linear equations, i.e., the order of the\n*     matrix A.  N >= 0.\n*\n*     NRHS    (input) INTEGER\n*     The number of right hand sides, i.e., the number of columns\n*     of the matrices B and X.  NRHS >= 0.\n*\n*     A       (input/output) COMPLEX*16 array, dimension (LDA,N)\n*     On entry, the N-by-N matrix A.  If FACT = 'F' and EQUED is\n*     not 'N', then A must have been equilibrated by the scaling\n*     factors in R and/or C.  A is not modified if FACT = 'F' or\n*     'N', or if FACT = 'E' and EQUED = 'N' on exit.\n*\n*     On exit, if EQUED .ne. 'N', A is scaled as follows:\n*     EQUED = 'R':  A := diag(R) * A\n*     EQUED = 'C':  A := A * diag(C)\n*     EQUED = 'B':  A := diag(R) * A * diag(C).\n*\n*     LDA     (input) INTEGER\n*     The leading dimension of the array A.  LDA >= max(1,N).\n*\n*     AF      (input or output) COMPLEX*16 array, dimension (LDAF,N)\n*     If FACT = 'F', then AF is an input argument and on entry\n*     contains the factors L and U from the factorization\n*     A = P*L*U as computed by ZGETRF.  If EQUED .ne. 'N', then\n*     AF is the factored form of the equilibrated matrix A.\n*\n*     If FACT = 'N', then AF is an output argument and on exit\n*     returns the factors L and U from the factorization A = P*L*U\n*     of the original matrix A.\n*\n*     If FACT = 'E', then AF is an output argument and on exit\n*     returns the factors L and U from the factorization A = P*L*U\n*     of the equilibrated matrix A (see the description of A for\n*     the form of the equilibrated matrix).\n*\n*     LDAF    (input) INTEGER\n*     The leading dimension of the array AF.  LDAF >= max(1,N).\n*\n*     IPIV    (input or output) INTEGER array, dimension (N)\n*     If FACT = 'F', then IPIV is an input argument and on entry\n*     contains the pivot indices from the factorization A = P*L*U\n*     as computed by ZGETRF; row i of the matrix was interchanged\n*     with row IPIV(i).\n*\n*     If FACT = 'N', then IPIV is an output argument and on exit\n*     contains the pivot indices from the factorization A = P*L*U\n*     of the original matrix A.\n*\n*     If FACT = 'E', then IPIV is an output argument and on exit\n*     contains the pivot indices from the factorization A = P*L*U\n*     of the equilibrated matrix A.\n*\n*     EQUED   (input or output) CHARACTER*1\n*     Specifies the form of equilibration that was done.\n*       = 'N':  No equilibration (always true if FACT = 'N').\n*       = 'R':  Row equilibration, i.e., A has been premultiplied by\n*               diag(R).\n*       = 'C':  Column equilibration, i.e., A has been postmultiplied\n*               by diag(C).\n*       = 'B':  Both row and column equilibration, i.e., A has been\n*               replaced by diag(R) * A * diag(C).\n*     EQUED is an input argument if FACT = 'F'; otherwise, it is an\n*     output argument.\n*\n*     R       (input or output) DOUBLE PRECISION array, dimension (N)\n*     The row scale factors for A.  If EQUED = 'R' or 'B', A is\n*     multiplied on the left by diag(R); if EQUED = 'N' or 'C', R\n*     is not accessed.  R is an input argument if FACT = 'F';\n*     otherwise, R is an output argument.  If FACT = 'F' and\n*     EQUED = 'R' or 'B', each element of R must be positive.\n*     If R is output, each element of R is a power of the radix.\n*     If R is input, each element of R should be a power of the radix\n*     to ensure a reliable solution and error estimates. Scaling by\n*     powers of the radix does not cause rounding errors unless the\n*     result underflows or overflows. Rounding errors during scaling\n*     lead to refining with a matrix that is not equivalent to the\n*     input matrix, producing error estimates that may not be\n*     reliable.\n*\n*     C       (input or output) DOUBLE PRECISION array, dimension (N)\n*     The column scale factors for A.  If EQUED = 'C' or 'B', A is\n*     multiplied on the right by diag(C); if EQUED = 'N' or 'R', C\n*     is not accessed.  C is an input argument if FACT = 'F';\n*     otherwise, C is an output argument.  If FACT = 'F' and\n*     EQUED = 'C' or 'B', each element of C must be positive.\n*     If C is output, each element of C is a power of the radix.\n*     If C is input, each element of C should be a power of the radix\n*     to ensure a reliable solution and error estimates. Scaling by\n*     powers of the radix does not cause rounding errors unless the\n*     result underflows or overflows. Rounding errors during scaling\n*     lead to refining with a matrix that is not equivalent to the\n*     input matrix, producing error estimates that may not be\n*     reliable.\n*\n*     B       (input/output) COMPLEX*16 array, dimension (LDB,NRHS)\n*     On entry, the N-by-NRHS right hand side matrix B.\n*     On exit,\n*     if EQUED = 'N', B is not modified;\n*     if TRANS = 'N' and EQUED = 'R' or 'B', B is overwritten by\n*        diag(R)*B;\n*     if TRANS = 'T' or 'C' and EQUED = 'C' or 'B', B is\n*        overwritten by diag(C)*B.\n*\n*     LDB     (input) INTEGER\n*     The leading dimension of the array B.  LDB >= max(1,N).\n*\n*     X       (output) COMPLEX*16 array, dimension (LDX,NRHS)\n*     If INFO = 0, the N-by-NRHS solution matrix X to the original\n*     system of equations.  Note that A and B are modified on exit\n*     if EQUED .ne. 'N', and the solution to the equilibrated system is\n*     inv(diag(C))*X if TRANS = 'N' and EQUED = 'C' or 'B', or\n*     inv(diag(R))*X if TRANS = 'T' or 'C' and EQUED = 'R' or 'B'.\n*\n*     LDX     (input) INTEGER\n*     The leading dimension of the array X.  LDX >= max(1,N).\n*\n*     RCOND   (output) DOUBLE PRECISION\n*     Reciprocal scaled condition number.  This is an estimate of the\n*     reciprocal Skeel condition number of the matrix A after\n*     equilibration (if done).  If this is less than the machine\n*     precision (in particular, if it is zero), the matrix is singular\n*     to working precision.  Note that the error may still be small even\n*     if this number is very small and the matrix appears ill-\n*     conditioned.\n*\n*     RPVGRW  (output) DOUBLE PRECISION\n*     Reciprocal pivot growth.  On exit, this contains the reciprocal\n*     pivot growth factor norm(A)/norm(U). The \"max absolute element\"\n*     norm is used.  If this is much less than 1, then the stability of\n*     the LU factorization of the (equilibrated) matrix A could be poor.\n*     This also means that the solution X, estimated condition numbers,\n*     and error bounds could be unreliable. If factorization fails with\n*     0<INFO<=N, then this contains the reciprocal pivot growth factor\n*     for the leading INFO columns of A.  In ZGESVX, this quantity is\n*     returned in WORK(1).\n*\n*     BERR    (output) DOUBLE PRECISION array, dimension (NRHS)\n*     Componentwise relative backward error.  This is the\n*     componentwise relative backward error of each solution vector X(j)\n*     (i.e., the smallest relative change in any element of A or B that\n*     makes X(j) an exact solution).\n*\n*     N_ERR_BNDS (input) INTEGER\n*     Number of error bounds to return for each right hand side\n*     and each type (normwise or componentwise).  See ERR_BNDS_NORM and\n*     ERR_BNDS_COMP below.\n*\n*     ERR_BNDS_NORM  (output) DOUBLE PRECISION array, dimension (NRHS, N_ERR_BNDS)\n*     For each right-hand side, this array contains information about\n*     various error bounds and condition numbers corresponding to the\n*     normwise relative error, which is defined as follows:\n*\n*     Normwise relative error in the ith solution vector:\n*             max_j (abs(XTRUE(j,i) - X(j,i)))\n*            ------------------------------\n*                  max_j abs(X(j,i))\n*\n*     The array is indexed by the type of error information as described\n*     below. There currently are up to three pieces of information\n*     returned.\n*\n*     The first index in ERR_BNDS_NORM(i,:) corresponds to the ith\n*     right-hand side.\n*\n*     The second index in ERR_BNDS_NORM(:,err) contains the following\n*     three fields:\n*     err = 1 \"Trust/don't trust\" boolean. Trust the answer if the\n*              reciprocal condition number is less than the threshold\n*              sqrt(n) * dlamch('Epsilon').\n*\n*     err = 2 \"Guaranteed\" error bound: The estimated forward error,\n*              almost certainly within a factor of 10 of the true error\n*              so long as the next entry is greater than the threshold\n*              sqrt(n) * dlamch('Epsilon'). This error bound should only\n*              be trusted if the previous boolean is true.\n*\n*     err = 3  Reciprocal condition number: Estimated normwise\n*              reciprocal condition number.  Compared with the threshold\n*              sqrt(n) * dlamch('Epsilon') to determine if the error\n*              estimate is \"guaranteed\". These reciprocal condition\n*              numbers are 1 / (norm(Z^{-1},inf) * norm(Z,inf)) for some\n*              appropriately scaled matrix Z.\n*              Let Z = S*A, where S scales each row by a power of the\n*              radix so all absolute row sums of Z are approximately 1.\n*\n*     See Lapack Working Note 165 for further details and extra\n*     cautions.\n*\n*     ERR_BNDS_COMP  (output) DOUBLE PRECISION array, dimension (NRHS, N_ERR_BNDS)\n*     For each right-hand side, this array contains information about\n*     various error bounds and condition numbers corresponding to the\n*     componentwise relative error, which is defined as follows:\n*\n*     Componentwise relative error in the ith solution vector:\n*                    abs(XTRUE(j,i) - X(j,i))\n*             max_j ----------------------\n*                         abs(X(j,i))\n*\n*     The array is indexed by the right-hand side i (on which the\n*     componentwise relative error depends), and the type of error\n*     information as described below. There currently are up to three\n*     pieces of information returned for each right-hand side. If\n*     componentwise accuracy is not requested (PARAMS(3) = 0.0), then\n*     ERR_BNDS_COMP is not accessed.  If N_ERR_BNDS .LT. 3, then at most\n*     the first (:,N_ERR_BNDS) entries are returned.\n*\n*     The first index in ERR_BNDS_COMP(i,:) corresponds to the ith\n*     right-hand side.\n*\n*     The second index in ERR_BNDS_COMP(:,err) contains the following\n*     three fields:\n*     err = 1 \"Trust/don't trust\" boolean. Trust the answer if the\n*              reciprocal condition number is less than the threshold\n*              sqrt(n) * dlamch('Epsilon').\n*\n*     err = 2 \"Guaranteed\" error bound: The estimated forward error,\n*              almost certainly within a factor of 10 of the true error\n*              so long as the next entry is greater than the threshold\n*              sqrt(n) * dlamch('Epsilon'). This error bound should only\n*              be trusted if the previous boolean is true.\n*\n*     err = 3  Reciprocal condition number: Estimated componentwise\n*              reciprocal condition number.  Compared with the threshold\n*              sqrt(n) * dlamch('Epsilon') to determine if the error\n*              estimate is \"guaranteed\". These reciprocal condition\n*              numbers are 1 / (norm(Z^{-1},inf) * norm(Z,inf)) for some\n*              appropriately scaled matrix Z.\n*              Let Z = S*(A*diag(x)), where x is the solution for the\n*              current right-hand side and S scales each row of\n*              A*diag(x) by a power of the radix so all absolute row\n*              sums of Z are approximately 1.\n*\n*     See Lapack Working Note 165 for further details and extra\n*     cautions.\n*\n*     NPARAMS (input) INTEGER\n*     Specifies the number of parameters set in PARAMS.  If .LE. 0, the\n*     PARAMS array is never referenced and default values are used.\n*\n*     PARAMS  (input / output) DOUBLE PRECISION array, dimension NPARAMS\n*     Specifies algorithm parameters.  If an entry is .LT. 0.0, then\n*     that entry will be filled with default value used for that\n*     parameter.  Only positions up to NPARAMS are accessed; defaults\n*     are used for higher-numbered parameters.\n*\n*       PARAMS(LA_LINRX_ITREF_I = 1) : Whether to perform iterative\n*            refinement or not.\n*         Default: 1.0D+0\n*            = 0.0 : No refinement is performed, and no error bounds are\n*                    computed.\n*            = 1.0 : Use the extra-precise refinement algorithm.\n*              (other values are reserved for future use)\n*\n*       PARAMS(LA_LINRX_ITHRESH_I = 2) : Maximum number of residual\n*            computations allowed for refinement.\n*         Default: 10\n*         Aggressive: Set to 100 to permit convergence using approximate\n*                     factorizations or factorizations other than LU. If\n*                     the factorization uses a technique other than\n*                     Gaussian elimination, the guarantees in\n*                     err_bnds_norm and err_bnds_comp may no longer be\n*                     trustworthy.\n*\n*       PARAMS(LA_LINRX_CWISE_I = 3) : Flag determining if the code\n*            will attempt to find a solution with small componentwise\n*            relative error in the double-precision algorithm.  Positive\n*            is true, 0.0 is false.\n*         Default: 1.0 (attempt componentwise convergence)\n*\n*     WORK    (workspace) COMPLEX*16 array, dimension (2*N)\n*\n*     RWORK   (workspace) DOUBLE PRECISION array, dimension (2*N)\n*\n*     INFO    (output) INTEGER\n*       = 0:  Successful exit. The solution to every right-hand side is\n*         guaranteed.\n*       < 0:  If INFO = -i, the i-th argument had an illegal value\n*       > 0 and <= N:  U(INFO,INFO) is exactly zero.  The factorization\n*         has been completed, but the factor U is exactly singular, so\n*         the solution and error bounds could not be computed. RCOND = 0\n*         is returned.\n*       = N+J: The solution corresponding to the Jth right-hand side is\n*         not guaranteed. The solutions corresponding to other right-\n*         hand sides K with K > J may not be guaranteed as well, but\n*         only the first such right-hand side is reported. If a small\n*         componentwise error is not requested (PARAMS(3) = 0.0) then\n*         the Jth right-hand side is the first with a normwise error\n*         bound that is not guaranteed (the smallest J such\n*         that ERR_BNDS_NORM(J,1) = 0.0). By default (PARAMS(3) = 1.0)\n*         the Jth right-hand side is the first with either a normwise or\n*         componentwise error bound that is not guaranteed (the smallest\n*         J such that either ERR_BNDS_NORM(J,1) = 0.0 or\n*         ERR_BNDS_COMP(J,1) = 0.0). See the definition of\n*         ERR_BNDS_NORM(:,1) and ERR_BNDS_COMP(:,1). To get information\n*         about all of the right-hand sides check ERR_BNDS_NORM or\n*         ERR_BNDS_COMP.\n*\n\n*     ==================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 10)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 10)", argc);
  rb_fact = argv[0];
  rb_trans = argv[1];
  rb_a = argv[2];
  rb_af = argv[3];
  rb_ipiv = argv[4];
  rb_equed = argv[5];
  rb_r = argv[6];
  rb_c = argv[7];
  rb_b = argv[8];
  rb_params = argv[9];

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
  if (NA_TYPE(rb_a) != NA_DCOMPLEX)
    rb_a = na_change_type(rb_a, NA_DCOMPLEX);
  a = NA_PTR_TYPE(rb_a, doublecomplex*);
  if (!NA_IsNArray(rb_b))
    rb_raise(rb_eArgError, "b (9th argument) must be NArray");
  if (NA_RANK(rb_b) != 2)
    rb_raise(rb_eArgError, "rank of b (9th argument) must be %d", 2);
  nrhs = NA_SHAPE1(rb_b);
  ldb = NA_SHAPE0(rb_b);
  if (NA_TYPE(rb_b) != NA_DCOMPLEX)
    rb_b = na_change_type(rb_b, NA_DCOMPLEX);
  b = NA_PTR_TYPE(rb_b, doublecomplex*);
  if (!NA_IsNArray(rb_c))
    rb_raise(rb_eArgError, "c (8th argument) must be NArray");
  if (NA_RANK(rb_c) != 1)
    rb_raise(rb_eArgError, "rank of c (8th argument) must be %d", 1);
  if (NA_SHAPE0(rb_c) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of c must be the same as shape 0 of ipiv");
  if (NA_TYPE(rb_c) != NA_DFLOAT)
    rb_c = na_change_type(rb_c, NA_DFLOAT);
  c = NA_PTR_TYPE(rb_c, doublereal*);
  fact = StringValueCStr(rb_fact)[0];
  equed = StringValueCStr(rb_equed)[0];
  if (!NA_IsNArray(rb_r))
    rb_raise(rb_eArgError, "r (7th argument) must be NArray");
  if (NA_RANK(rb_r) != 1)
    rb_raise(rb_eArgError, "rank of r (7th argument) must be %d", 1);
  if (NA_SHAPE0(rb_r) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of r must be the same as shape 0 of ipiv");
  if (NA_TYPE(rb_r) != NA_DFLOAT)
    rb_r = na_change_type(rb_r, NA_DFLOAT);
  r = NA_PTR_TYPE(rb_r, doublereal*);
  n_err_bnds = 3;
  if (!NA_IsNArray(rb_params))
    rb_raise(rb_eArgError, "params (10th argument) must be NArray");
  if (NA_RANK(rb_params) != 1)
    rb_raise(rb_eArgError, "rank of params (10th argument) must be %d", 1);
  nparams = NA_SHAPE0(rb_params);
  if (NA_TYPE(rb_params) != NA_DFLOAT)
    rb_params = na_change_type(rb_params, NA_DFLOAT);
  params = NA_PTR_TYPE(rb_params, doublereal*);
  if (!NA_IsNArray(rb_af))
    rb_raise(rb_eArgError, "af (4th argument) must be NArray");
  if (NA_RANK(rb_af) != 2)
    rb_raise(rb_eArgError, "rank of af (4th argument) must be %d", 2);
  if (NA_SHAPE1(rb_af) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of af must be the same as shape 0 of ipiv");
  ldaf = NA_SHAPE0(rb_af);
  if (NA_TYPE(rb_af) != NA_DCOMPLEX)
    rb_af = na_change_type(rb_af, NA_DCOMPLEX);
  af = NA_PTR_TYPE(rb_af, doublecomplex*);
  ldx = MAX(1,n);
  {
    int shape[2];
    shape[0] = ldx;
    shape[1] = nrhs;
    rb_x = na_make_object(NA_DCOMPLEX, 2, shape, cNArray);
  }
  x = NA_PTR_TYPE(rb_x, doublecomplex*);
  {
    int shape[1];
    shape[0] = nrhs;
    rb_berr = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  berr = NA_PTR_TYPE(rb_berr, doublereal*);
  {
    int shape[2];
    shape[0] = nrhs;
    shape[1] = n_err_bnds;
    rb_err_bnds_norm = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  err_bnds_norm = NA_PTR_TYPE(rb_err_bnds_norm, doublereal*);
  {
    int shape[2];
    shape[0] = nrhs;
    shape[1] = n_err_bnds;
    rb_err_bnds_comp = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  err_bnds_comp = NA_PTR_TYPE(rb_err_bnds_comp, doublereal*);
  {
    int shape[2];
    shape[0] = lda;
    shape[1] = n;
    rb_a_out__ = na_make_object(NA_DCOMPLEX, 2, shape, cNArray);
  }
  a_out__ = NA_PTR_TYPE(rb_a_out__, doublecomplex*);
  MEMCPY(a_out__, a, doublecomplex, NA_TOTAL(rb_a));
  rb_a = rb_a_out__;
  a = a_out__;
  {
    int shape[2];
    shape[0] = ldaf;
    shape[1] = n;
    rb_af_out__ = na_make_object(NA_DCOMPLEX, 2, shape, cNArray);
  }
  af_out__ = NA_PTR_TYPE(rb_af_out__, doublecomplex*);
  MEMCPY(af_out__, af, doublecomplex, NA_TOTAL(rb_af));
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
    rb_b_out__ = na_make_object(NA_DCOMPLEX, 2, shape, cNArray);
  }
  b_out__ = NA_PTR_TYPE(rb_b_out__, doublecomplex*);
  MEMCPY(b_out__, b, doublecomplex, NA_TOTAL(rb_b));
  rb_b = rb_b_out__;
  b = b_out__;
  {
    int shape[1];
    shape[0] = nparams;
    rb_params_out__ = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  params_out__ = NA_PTR_TYPE(rb_params_out__, doublereal*);
  MEMCPY(params_out__, params, doublereal, NA_TOTAL(rb_params));
  rb_params = rb_params_out__;
  params = params_out__;
  work = ALLOC_N(doublecomplex, (2*n));
  rwork = ALLOC_N(doublereal, (2*n));

  zgesvxx_(&fact, &trans, &n, &nrhs, a, &lda, af, &ldaf, ipiv, &equed, r, c, b, &ldb, x, &ldx, &rcond, &rpvgrw, berr, &n_err_bnds, err_bnds_norm, err_bnds_comp, &nparams, params, work, rwork, &info);

  free(work);
  free(rwork);
  rb_rcond = rb_float_new((double)rcond);
  rb_rpvgrw = rb_float_new((double)rpvgrw);
  rb_info = INT2NUM(info);
  rb_equed = rb_str_new(&equed,1);
  return rb_ary_new3(15, rb_x, rb_rcond, rb_rpvgrw, rb_berr, rb_err_bnds_norm, rb_err_bnds_comp, rb_info, rb_a, rb_af, rb_ipiv, rb_equed, rb_r, rb_c, rb_b, rb_params);
}

void
init_lapack_zgesvxx(VALUE mLapack){
  rb_define_module_function(mLapack, "zgesvxx", rb_zgesvxx, -1);
}

#include "rb_lapack.h"

extern VOID dgbrfsx_(char *trans, char *equed, integer *n, integer *kl, integer *ku, integer *nrhs, doublereal *ab, integer *ldab, doublereal *afb, integer *ldafb, integer *ipiv, doublereal *r, doublereal *c, doublereal *b, integer *ldb, doublereal *x, integer *ldx, doublereal *rcond, doublereal *berr, integer *n_err_bnds, doublereal *err_bnds_norm, doublereal *err_bnds_comp, integer *nparams, doublereal *params, doublereal *work, integer *iwork, integer *info);

static VALUE sHelp, sUsage;

static VALUE
rblapack_dgbrfsx(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_trans;
  char trans; 
  VALUE rblapack_equed;
  char equed; 
  VALUE rblapack_kl;
  integer kl; 
  VALUE rblapack_ku;
  integer ku; 
  VALUE rblapack_ab;
  doublereal *ab; 
  VALUE rblapack_afb;
  doublereal *afb; 
  VALUE rblapack_ipiv;
  integer *ipiv; 
  VALUE rblapack_r;
  doublereal *r; 
  VALUE rblapack_c;
  doublereal *c; 
  VALUE rblapack_b;
  doublereal *b; 
  VALUE rblapack_x;
  doublereal *x; 
  VALUE rblapack_params;
  doublereal *params; 
  VALUE rblapack_rcond;
  doublereal rcond; 
  VALUE rblapack_berr;
  doublereal *berr; 
  VALUE rblapack_err_bnds_norm;
  doublereal *err_bnds_norm; 
  VALUE rblapack_err_bnds_comp;
  doublereal *err_bnds_comp; 
  VALUE rblapack_info;
  integer info; 
  VALUE rblapack_r_out__;
  doublereal *r_out__;
  VALUE rblapack_c_out__;
  doublereal *c_out__;
  VALUE rblapack_x_out__;
  doublereal *x_out__;
  VALUE rblapack_params_out__;
  doublereal *params_out__;
  doublereal *work;
  integer *iwork;

  integer ldab;
  integer n;
  integer ldafb;
  integer ldb;
  integer nrhs;
  integer ldx;
  integer nparams;
  integer n_err_bnds;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  rcond, berr, err_bnds_norm, err_bnds_comp, info, r, c, x, params = NumRu::Lapack.dgbrfsx( trans, equed, kl, ku, ab, afb, ipiv, r, c, b, x, params, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE DGBRFSX( TRANS, EQUED, N, KL, KU, NRHS, AB, LDAB, AFB, LDAFB, IPIV, R, C, B, LDB, X, LDX, RCOND, BERR, N_ERR_BNDS, ERR_BNDS_NORM, ERR_BNDS_COMP, NPARAMS, PARAMS, WORK, IWORK, INFO )\n\n*     Purpose\n*     =======\n*\n*     DGBRFSX improves the computed solution to a system of linear\n*     equations and provides error bounds and backward error estimates\n*     for the solution.  In addition to normwise error bound, the code\n*     provides maximum componentwise error bound if possible.  See\n*     comments for ERR_BNDS_NORM and ERR_BNDS_COMP for details of the\n*     error bounds.\n*\n*     The original system of linear equations may have been equilibrated\n*     before calling this routine, as described by arguments EQUED, R\n*     and C below. In this case, the solution and error bounds returned\n*     are for the original unequilibrated system.\n*\n\n*     Arguments\n*     =========\n*\n*     Some optional parameters are bundled in the PARAMS array.  These\n*     settings determine how refinement is performed, but often the\n*     defaults are acceptable.  If the defaults are acceptable, users\n*     can pass NPARAMS = 0 which prevents the source code from accessing\n*     the PARAMS argument.\n*\n*     TRANS   (input) CHARACTER*1\n*     Specifies the form of the system of equations:\n*       = 'N':  A * X = B     (No transpose)\n*       = 'T':  A**T * X = B  (Transpose)\n*       = 'C':  A**H * X = B  (Conjugate transpose = Transpose)\n*\n*     EQUED   (input) CHARACTER*1\n*     Specifies the form of equilibration that was done to A\n*     before calling this routine. This is needed to compute\n*     the solution and error bounds correctly.\n*       = 'N':  No equilibration\n*       = 'R':  Row equilibration, i.e., A has been premultiplied by\n*               diag(R).\n*       = 'C':  Column equilibration, i.e., A has been postmultiplied\n*               by diag(C).\n*       = 'B':  Both row and column equilibration, i.e., A has been\n*               replaced by diag(R) * A * diag(C).\n*               The right hand side B has been changed accordingly.\n*\n*     N       (input) INTEGER\n*     The order of the matrix A.  N >= 0.\n*\n*     KL      (input) INTEGER\n*     The number of subdiagonals within the band of A.  KL >= 0.\n*\n*     KU      (input) INTEGER\n*     The number of superdiagonals within the band of A.  KU >= 0.\n*\n*     NRHS    (input) INTEGER\n*     The number of right hand sides, i.e., the number of columns\n*     of the matrices B and X.  NRHS >= 0.\n*\n*     AB      (input) DOUBLE PRECISION array, dimension (LDAB,N)\n*     The original band matrix A, stored in rows 1 to KL+KU+1.\n*     The j-th column of A is stored in the j-th column of the\n*     array AB as follows:\n*     AB(ku+1+i-j,j) = A(i,j) for max(1,j-ku)<=i<=min(n,j+kl).\n*\n*     LDAB    (input) INTEGER\n*     The leading dimension of the array AB.  LDAB >= KL+KU+1.\n*\n*     AFB     (input) DOUBLE PRECISION array, dimension (LDAFB,N)\n*     Details of the LU factorization of the band matrix A, as\n*     computed by DGBTRF.  U is stored as an upper triangular band\n*     matrix with KL+KU superdiagonals in rows 1 to KL+KU+1, and\n*     the multipliers used during the factorization are stored in\n*     rows KL+KU+2 to 2*KL+KU+1.\n*\n*     LDAFB   (input) INTEGER\n*     The leading dimension of the array AFB.  LDAFB >= 2*KL*KU+1.\n*\n*     IPIV    (input) INTEGER array, dimension (N)\n*     The pivot indices from DGETRF; for 1<=i<=N, row i of the\n*     matrix was interchanged with row IPIV(i).\n*\n*     R       (input or output) DOUBLE PRECISION array, dimension (N)\n*     The row scale factors for A.  If EQUED = 'R' or 'B', A is\n*     multiplied on the left by diag(R); if EQUED = 'N' or 'C', R\n*     is not accessed.  R is an input argument if FACT = 'F';\n*     otherwise, R is an output argument.  If FACT = 'F' and\n*     EQUED = 'R' or 'B', each element of R must be positive.\n*     If R is output, each element of R is a power of the radix.\n*     If R is input, each element of R should be a power of the radix\n*     to ensure a reliable solution and error estimates. Scaling by\n*     powers of the radix does not cause rounding errors unless the\n*     result underflows or overflows. Rounding errors during scaling\n*     lead to refining with a matrix that is not equivalent to the\n*     input matrix, producing error estimates that may not be\n*     reliable.\n*\n*     C       (input or output) DOUBLE PRECISION array, dimension (N)\n*     The column scale factors for A.  If EQUED = 'C' or 'B', A is\n*     multiplied on the right by diag(C); if EQUED = 'N' or 'R', C\n*     is not accessed.  C is an input argument if FACT = 'F';\n*     otherwise, C is an output argument.  If FACT = 'F' and\n*     EQUED = 'C' or 'B', each element of C must be positive.\n*     If C is output, each element of C is a power of the radix.\n*     If C is input, each element of C should be a power of the radix\n*     to ensure a reliable solution and error estimates. Scaling by\n*     powers of the radix does not cause rounding errors unless the\n*     result underflows or overflows. Rounding errors during scaling\n*     lead to refining with a matrix that is not equivalent to the\n*     input matrix, producing error estimates that may not be\n*     reliable.\n*\n*     B       (input) DOUBLE PRECISION array, dimension (LDB,NRHS)\n*     The right hand side matrix B.\n*\n*     LDB     (input) INTEGER\n*     The leading dimension of the array B.  LDB >= max(1,N).\n*\n*     X       (input/output) DOUBLE PRECISION array, dimension (LDX,NRHS)\n*     On entry, the solution matrix X, as computed by DGETRS.\n*     On exit, the improved solution matrix X.\n*\n*     LDX     (input) INTEGER\n*     The leading dimension of the array X.  LDX >= max(1,N).\n*\n*     RCOND   (output) DOUBLE PRECISION\n*     Reciprocal scaled condition number.  This is an estimate of the\n*     reciprocal Skeel condition number of the matrix A after\n*     equilibration (if done).  If this is less than the machine\n*     precision (in particular, if it is zero), the matrix is singular\n*     to working precision.  Note that the error may still be small even\n*     if this number is very small and the matrix appears ill-\n*     conditioned.\n*\n*     BERR    (output) DOUBLE PRECISION array, dimension (NRHS)\n*     Componentwise relative backward error.  This is the\n*     componentwise relative backward error of each solution vector X(j)\n*     (i.e., the smallest relative change in any element of A or B that\n*     makes X(j) an exact solution).\n*\n*     N_ERR_BNDS (input) INTEGER\n*     Number of error bounds to return for each right hand side\n*     and each type (normwise or componentwise).  See ERR_BNDS_NORM and\n*     ERR_BNDS_COMP below.\n*\n*     ERR_BNDS_NORM  (output) DOUBLE PRECISION array, dimension (NRHS, N_ERR_BNDS)\n*     For each right-hand side, this array contains information about\n*     various error bounds and condition numbers corresponding to the\n*     normwise relative error, which is defined as follows:\n*\n*     Normwise relative error in the ith solution vector:\n*             max_j (abs(XTRUE(j,i) - X(j,i)))\n*            ------------------------------\n*                  max_j abs(X(j,i))\n*\n*     The array is indexed by the type of error information as described\n*     below. There currently are up to three pieces of information\n*     returned.\n*\n*     The first index in ERR_BNDS_NORM(i,:) corresponds to the ith\n*     right-hand side.\n*\n*     The second index in ERR_BNDS_NORM(:,err) contains the following\n*     three fields:\n*     err = 1 \"Trust/don't trust\" boolean. Trust the answer if the\n*              reciprocal condition number is less than the threshold\n*              sqrt(n) * dlamch('Epsilon').\n*\n*     err = 2 \"Guaranteed\" error bound: The estimated forward error,\n*              almost certainly within a factor of 10 of the true error\n*              so long as the next entry is greater than the threshold\n*              sqrt(n) * dlamch('Epsilon'). This error bound should only\n*              be trusted if the previous boolean is true.\n*\n*     err = 3  Reciprocal condition number: Estimated normwise\n*              reciprocal condition number.  Compared with the threshold\n*              sqrt(n) * dlamch('Epsilon') to determine if the error\n*              estimate is \"guaranteed\". These reciprocal condition\n*              numbers are 1 / (norm(Z^{-1},inf) * norm(Z,inf)) for some\n*              appropriately scaled matrix Z.\n*              Let Z = S*A, where S scales each row by a power of the\n*              radix so all absolute row sums of Z are approximately 1.\n*\n*     See Lapack Working Note 165 for further details and extra\n*     cautions.\n*\n*     ERR_BNDS_COMP  (output) DOUBLE PRECISION array, dimension (NRHS, N_ERR_BNDS)\n*     For each right-hand side, this array contains information about\n*     various error bounds and condition numbers corresponding to the\n*     componentwise relative error, which is defined as follows:\n*\n*     Componentwise relative error in the ith solution vector:\n*                    abs(XTRUE(j,i) - X(j,i))\n*             max_j ----------------------\n*                         abs(X(j,i))\n*\n*     The array is indexed by the right-hand side i (on which the\n*     componentwise relative error depends), and the type of error\n*     information as described below. There currently are up to three\n*     pieces of information returned for each right-hand side. If\n*     componentwise accuracy is not requested (PARAMS(3) = 0.0), then\n*     ERR_BNDS_COMP is not accessed.  If N_ERR_BNDS .LT. 3, then at most\n*     the first (:,N_ERR_BNDS) entries are returned.\n*\n*     The first index in ERR_BNDS_COMP(i,:) corresponds to the ith\n*     right-hand side.\n*\n*     The second index in ERR_BNDS_COMP(:,err) contains the following\n*     three fields:\n*     err = 1 \"Trust/don't trust\" boolean. Trust the answer if the\n*              reciprocal condition number is less than the threshold\n*              sqrt(n) * dlamch('Epsilon').\n*\n*     err = 2 \"Guaranteed\" error bound: The estimated forward error,\n*              almost certainly within a factor of 10 of the true error\n*              so long as the next entry is greater than the threshold\n*              sqrt(n) * dlamch('Epsilon'). This error bound should only\n*              be trusted if the previous boolean is true.\n*\n*     err = 3  Reciprocal condition number: Estimated componentwise\n*              reciprocal condition number.  Compared with the threshold\n*              sqrt(n) * dlamch('Epsilon') to determine if the error\n*              estimate is \"guaranteed\". These reciprocal condition\n*              numbers are 1 / (norm(Z^{-1},inf) * norm(Z,inf)) for some\n*              appropriately scaled matrix Z.\n*              Let Z = S*(A*diag(x)), where x is the solution for the\n*              current right-hand side and S scales each row of\n*              A*diag(x) by a power of the radix so all absolute row\n*              sums of Z are approximately 1.\n*\n*     See Lapack Working Note 165 for further details and extra\n*     cautions.\n*\n*     NPARAMS (input) INTEGER\n*     Specifies the number of parameters set in PARAMS.  If .LE. 0, the\n*     PARAMS array is never referenced and default values are used.\n*\n*     PARAMS  (input / output) DOUBLE PRECISION array, dimension (NPARAMS)\n*     Specifies algorithm parameters.  If an entry is .LT. 0.0, then\n*     that entry will be filled with default value used for that\n*     parameter.  Only positions up to NPARAMS are accessed; defaults\n*     are used for higher-numbered parameters.\n*\n*       PARAMS(LA_LINRX_ITREF_I = 1) : Whether to perform iterative\n*            refinement or not.\n*         Default: 1.0D+0\n*            = 0.0 : No refinement is performed, and no error bounds are\n*                    computed.\n*            = 1.0 : Use the double-precision refinement algorithm,\n*                    possibly with doubled-single computations if the\n*                    compilation environment does not support DOUBLE\n*                    PRECISION.\n*              (other values are reserved for future use)\n*\n*       PARAMS(LA_LINRX_ITHRESH_I = 2) : Maximum number of residual\n*            computations allowed for refinement.\n*         Default: 10\n*         Aggressive: Set to 100 to permit convergence using approximate\n*                     factorizations or factorizations other than LU. If\n*                     the factorization uses a technique other than\n*                     Gaussian elimination, the guarantees in\n*                     err_bnds_norm and err_bnds_comp may no longer be\n*                     trustworthy.\n*\n*       PARAMS(LA_LINRX_CWISE_I = 3) : Flag determining if the code\n*            will attempt to find a solution with small componentwise\n*            relative error in the double-precision algorithm.  Positive\n*            is true, 0.0 is false.\n*         Default: 1.0 (attempt componentwise convergence)\n*\n*     WORK    (workspace) DOUBLE PRECISION array, dimension (4*N)\n*\n*     IWORK   (workspace) INTEGER array, dimension (N)\n*\n*     INFO    (output) INTEGER\n*       = 0:  Successful exit. The solution to every right-hand side is\n*         guaranteed.\n*       < 0:  If INFO = -i, the i-th argument had an illegal value\n*       > 0 and <= N:  U(INFO,INFO) is exactly zero.  The factorization\n*         has been completed, but the factor U is exactly singular, so\n*         the solution and error bounds could not be computed. RCOND = 0\n*         is returned.\n*       = N+J: The solution corresponding to the Jth right-hand side is\n*         not guaranteed. The solutions corresponding to other right-\n*         hand sides K with K > J may not be guaranteed as well, but\n*         only the first such right-hand side is reported. If a small\n*         componentwise error is not requested (PARAMS(3) = 0.0) then\n*         the Jth right-hand side is the first with a normwise error\n*         bound that is not guaranteed (the smallest J such\n*         that ERR_BNDS_NORM(J,1) = 0.0). By default (PARAMS(3) = 1.0)\n*         the Jth right-hand side is the first with either a normwise or\n*         componentwise error bound that is not guaranteed (the smallest\n*         J such that either ERR_BNDS_NORM(J,1) = 0.0 or\n*         ERR_BNDS_COMP(J,1) = 0.0). See the definition of\n*         ERR_BNDS_NORM(:,1) and ERR_BNDS_COMP(:,1). To get information\n*         about all of the right-hand sides check ERR_BNDS_NORM or\n*         ERR_BNDS_COMP.\n*\n\n*     ==================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  rcond, berr, err_bnds_norm, err_bnds_comp, info, r, c, x, params = NumRu::Lapack.dgbrfsx( trans, equed, kl, ku, ab, afb, ipiv, r, c, b, x, params, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 12)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 12)", argc);
  rblapack_trans = argv[0];
  rblapack_equed = argv[1];
  rblapack_kl = argv[2];
  rblapack_ku = argv[3];
  rblapack_ab = argv[4];
  rblapack_afb = argv[5];
  rblapack_ipiv = argv[6];
  rblapack_r = argv[7];
  rblapack_c = argv[8];
  rblapack_b = argv[9];
  rblapack_x = argv[10];
  rblapack_params = argv[11];
  if (rb_options != Qnil) {
  }

  if (!NA_IsNArray(rblapack_ipiv))
    rb_raise(rb_eArgError, "ipiv (7th argument) must be NArray");
  if (NA_RANK(rblapack_ipiv) != 1)
    rb_raise(rb_eArgError, "rank of ipiv (7th argument) must be %d", 1);
  n = NA_SHAPE0(rblapack_ipiv);
  if (NA_TYPE(rblapack_ipiv) != NA_LINT)
    rblapack_ipiv = na_change_type(rblapack_ipiv, NA_LINT);
  ipiv = NA_PTR_TYPE(rblapack_ipiv, integer*);
  if (!NA_IsNArray(rblapack_ab))
    rb_raise(rb_eArgError, "ab (5th argument) must be NArray");
  if (NA_RANK(rblapack_ab) != 2)
    rb_raise(rb_eArgError, "rank of ab (5th argument) must be %d", 2);
  if (NA_SHAPE1(rblapack_ab) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of ab must be the same as shape 0 of ipiv");
  ldab = NA_SHAPE0(rblapack_ab);
  if (NA_TYPE(rblapack_ab) != NA_DFLOAT)
    rblapack_ab = na_change_type(rblapack_ab, NA_DFLOAT);
  ab = NA_PTR_TYPE(rblapack_ab, doublereal*);
  kl = NUM2INT(rblapack_kl);
  if (!NA_IsNArray(rblapack_x))
    rb_raise(rb_eArgError, "x (11th argument) must be NArray");
  if (NA_RANK(rblapack_x) != 2)
    rb_raise(rb_eArgError, "rank of x (11th argument) must be %d", 2);
  nrhs = NA_SHAPE1(rblapack_x);
  ldx = NA_SHAPE0(rblapack_x);
  if (NA_TYPE(rblapack_x) != NA_DFLOAT)
    rblapack_x = na_change_type(rblapack_x, NA_DFLOAT);
  x = NA_PTR_TYPE(rblapack_x, doublereal*);
  if (!NA_IsNArray(rblapack_b))
    rb_raise(rb_eArgError, "b (10th argument) must be NArray");
  if (NA_RANK(rblapack_b) != 2)
    rb_raise(rb_eArgError, "rank of b (10th argument) must be %d", 2);
  if (NA_SHAPE1(rblapack_b) != nrhs)
    rb_raise(rb_eRuntimeError, "shape 1 of b must be the same as shape 1 of x");
  ldb = NA_SHAPE0(rblapack_b);
  if (NA_TYPE(rblapack_b) != NA_DFLOAT)
    rblapack_b = na_change_type(rblapack_b, NA_DFLOAT);
  b = NA_PTR_TYPE(rblapack_b, doublereal*);
  if (!NA_IsNArray(rblapack_c))
    rb_raise(rb_eArgError, "c (9th argument) must be NArray");
  if (NA_RANK(rblapack_c) != 1)
    rb_raise(rb_eArgError, "rank of c (9th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_c) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of c must be the same as shape 0 of ipiv");
  if (NA_TYPE(rblapack_c) != NA_DFLOAT)
    rblapack_c = na_change_type(rblapack_c, NA_DFLOAT);
  c = NA_PTR_TYPE(rblapack_c, doublereal*);
  if (!NA_IsNArray(rblapack_r))
    rb_raise(rb_eArgError, "r (8th argument) must be NArray");
  if (NA_RANK(rblapack_r) != 1)
    rb_raise(rb_eArgError, "rank of r (8th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_r) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of r must be the same as shape 0 of ipiv");
  if (NA_TYPE(rblapack_r) != NA_DFLOAT)
    rblapack_r = na_change_type(rblapack_r, NA_DFLOAT);
  r = NA_PTR_TYPE(rblapack_r, doublereal*);
  if (!NA_IsNArray(rblapack_afb))
    rb_raise(rb_eArgError, "afb (6th argument) must be NArray");
  if (NA_RANK(rblapack_afb) != 2)
    rb_raise(rb_eArgError, "rank of afb (6th argument) must be %d", 2);
  if (NA_SHAPE1(rblapack_afb) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of afb must be the same as shape 0 of ipiv");
  ldafb = NA_SHAPE0(rblapack_afb);
  if (NA_TYPE(rblapack_afb) != NA_DFLOAT)
    rblapack_afb = na_change_type(rblapack_afb, NA_DFLOAT);
  afb = NA_PTR_TYPE(rblapack_afb, doublereal*);
  n_err_bnds = 3;
  if (!NA_IsNArray(rblapack_params))
    rb_raise(rb_eArgError, "params (12th argument) must be NArray");
  if (NA_RANK(rblapack_params) != 1)
    rb_raise(rb_eArgError, "rank of params (12th argument) must be %d", 1);
  nparams = NA_SHAPE0(rblapack_params);
  if (NA_TYPE(rblapack_params) != NA_DFLOAT)
    rblapack_params = na_change_type(rblapack_params, NA_DFLOAT);
  params = NA_PTR_TYPE(rblapack_params, doublereal*);
  ku = NUM2INT(rblapack_ku);
  equed = StringValueCStr(rblapack_equed)[0];
  trans = StringValueCStr(rblapack_trans)[0];
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
    int shape[1];
    shape[0] = n;
    rblapack_r_out__ = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  r_out__ = NA_PTR_TYPE(rblapack_r_out__, doublereal*);
  MEMCPY(r_out__, r, doublereal, NA_TOTAL(rblapack_r));
  rblapack_r = rblapack_r_out__;
  r = r_out__;
  {
    int shape[1];
    shape[0] = n;
    rblapack_c_out__ = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  c_out__ = NA_PTR_TYPE(rblapack_c_out__, doublereal*);
  MEMCPY(c_out__, c, doublereal, NA_TOTAL(rblapack_c));
  rblapack_c = rblapack_c_out__;
  c = c_out__;
  {
    int shape[2];
    shape[0] = ldx;
    shape[1] = nrhs;
    rblapack_x_out__ = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  x_out__ = NA_PTR_TYPE(rblapack_x_out__, doublereal*);
  MEMCPY(x_out__, x, doublereal, NA_TOTAL(rblapack_x));
  rblapack_x = rblapack_x_out__;
  x = x_out__;
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

  dgbrfsx_(&trans, &equed, &n, &kl, &ku, &nrhs, ab, &ldab, afb, &ldafb, ipiv, r, c, b, &ldb, x, &ldx, &rcond, berr, &n_err_bnds, err_bnds_norm, err_bnds_comp, &nparams, params, work, iwork, &info);

  free(work);
  free(iwork);
  rblapack_rcond = rb_float_new((double)rcond);
  rblapack_info = INT2NUM(info);
  return rb_ary_new3(9, rblapack_rcond, rblapack_berr, rblapack_err_bnds_norm, rblapack_err_bnds_comp, rblapack_info, rblapack_r, rblapack_c, rblapack_x, rblapack_params);
}

void
init_lapack_dgbrfsx(VALUE mLapack){
  rb_define_module_function(mLapack, "dgbrfsx", rblapack_dgbrfsx, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}

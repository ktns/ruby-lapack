#include "rb_lapack.h"

extern VOID zherfsx_(char *uplo, char *equed, integer *n, integer *nrhs, doublecomplex *a, integer *lda, doublecomplex *af, integer *ldaf, integer *ipiv, doublereal *s, doublecomplex *b, integer *ldb, doublecomplex *x, integer *ldx, doublereal *rcond, doublereal *berr, integer *n_err_bnds, doublereal *err_bnds_norm, doublereal *err_bnds_comp, integer *nparams, doublereal *params, doublecomplex *work, doublereal *rwork, integer *info);

static VALUE sHelp, sUsage;

static VALUE
rblapack_zherfsx(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_uplo;
  char uplo; 
  VALUE rblapack_equed;
  char equed; 
  VALUE rblapack_a;
  doublecomplex *a; 
  VALUE rblapack_af;
  doublecomplex *af; 
  VALUE rblapack_ipiv;
  integer *ipiv; 
  VALUE rblapack_s;
  doublereal *s; 
  VALUE rblapack_b;
  doublecomplex *b; 
  VALUE rblapack_x;
  doublecomplex *x; 
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
  VALUE rblapack_s_out__;
  doublereal *s_out__;
  VALUE rblapack_x_out__;
  doublecomplex *x_out__;
  VALUE rblapack_params_out__;
  doublereal *params_out__;
  doublecomplex *work;
  doublereal *rwork;

  integer lda;
  integer n;
  integer ldaf;
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
      printf("%s\n", "USAGE:\n  rcond, berr, err_bnds_norm, err_bnds_comp, info, s, x, params = NumRu::Lapack.zherfsx( uplo, equed, a, af, ipiv, s, b, x, params, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE ZHERFSX( UPLO, EQUED, N, NRHS, A, LDA, AF, LDAF, IPIV, S, B, LDB, X, LDX, RCOND, BERR, N_ERR_BNDS, ERR_BNDS_NORM, ERR_BNDS_COMP, NPARAMS, PARAMS, WORK, RWORK, INFO )\n\n*     Purpose\n*     =======\n*\n*     ZHERFSX improves the computed solution to a system of linear\n*     equations when the coefficient matrix is Hermitian indefinite, and\n*     provides error bounds and backward error estimates for the\n*     solution.  In addition to normwise error bound, the code provides\n*     maximum componentwise error bound if possible.  See comments for\n*     ERR_BNDS_NORM and ERR_BNDS_COMP for details of the error bounds.\n*\n*     The original system of linear equations may have been equilibrated\n*     before calling this routine, as described by arguments EQUED and S\n*     below. In this case, the solution and error bounds returned are\n*     for the original unequilibrated system.\n*\n\n*     Arguments\n*     =========\n*\n*     Some optional parameters are bundled in the PARAMS array.  These\n*     settings determine how refinement is performed, but often the\n*     defaults are acceptable.  If the defaults are acceptable, users\n*     can pass NPARAMS = 0 which prevents the source code from accessing\n*     the PARAMS argument.\n*\n*     UPLO    (input) CHARACTER*1\n*       = 'U':  Upper triangle of A is stored;\n*       = 'L':  Lower triangle of A is stored.\n*\n*     EQUED   (input) CHARACTER*1\n*     Specifies the form of equilibration that was done to A\n*     before calling this routine. This is needed to compute\n*     the solution and error bounds correctly.\n*       = 'N':  No equilibration\n*       = 'Y':  Both row and column equilibration, i.e., A has been\n*               replaced by diag(S) * A * diag(S).\n*               The right hand side B has been changed accordingly.\n*\n*     N       (input) INTEGER\n*     The order of the matrix A.  N >= 0.\n*\n*     NRHS    (input) INTEGER\n*     The number of right hand sides, i.e., the number of columns\n*     of the matrices B and X.  NRHS >= 0.\n*\n*     A       (input) COMPLEX*16 array, dimension (LDA,N)\n*     The symmetric matrix A.  If UPLO = 'U', the leading N-by-N\n*     upper triangular part of A contains the upper triangular\n*     part of the matrix A, and the strictly lower triangular\n*     part of A is not referenced.  If UPLO = 'L', the leading\n*     N-by-N lower triangular part of A contains the lower\n*     triangular part of the matrix A, and the strictly upper\n*     triangular part of A is not referenced.\n*\n*     LDA     (input) INTEGER\n*     The leading dimension of the array A.  LDA >= max(1,N).\n*\n*     AF      (input) COMPLEX*16 array, dimension (LDAF,N)\n*     The factored form of the matrix A.  AF contains the block\n*     diagonal matrix D and the multipliers used to obtain the\n*     factor U or L from the factorization A = U*D*U**T or A =\n*     L*D*L**T as computed by DSYTRF.\n*\n*     LDAF    (input) INTEGER\n*     The leading dimension of the array AF.  LDAF >= max(1,N).\n*\n*     IPIV    (input) INTEGER array, dimension (N)\n*     Details of the interchanges and the block structure of D\n*     as determined by DSYTRF.\n*\n*     S       (input or output) DOUBLE PRECISION array, dimension (N)\n*     The scale factors for A.  If EQUED = 'Y', A is multiplied on\n*     the left and right by diag(S).  S is an input argument if FACT =\n*     'F'; otherwise, S is an output argument.  If FACT = 'F' and EQUED\n*     = 'Y', each element of S must be positive.  If S is output, each\n*     element of S is a power of the radix. If S is input, each element\n*     of S should be a power of the radix to ensure a reliable solution\n*     and error estimates. Scaling by powers of the radix does not cause\n*     rounding errors unless the result underflows or overflows.\n*     Rounding errors during scaling lead to refining with a matrix that\n*     is not equivalent to the input matrix, producing error estimates\n*     that may not be reliable.\n*\n*     B       (input) COMPLEX*16 array, dimension (LDB,NRHS)\n*     The right hand side matrix B.\n*\n*     LDB     (input) INTEGER\n*     The leading dimension of the array B.  LDB >= max(1,N).\n*\n*     X       (input/output) COMPLEX*16 array, dimension (LDX,NRHS)\n*     On entry, the solution matrix X, as computed by DGETRS.\n*     On exit, the improved solution matrix X.\n*\n*     LDX     (input) INTEGER\n*     The leading dimension of the array X.  LDX >= max(1,N).\n*\n*     RCOND   (output) DOUBLE PRECISION\n*     Reciprocal scaled condition number.  This is an estimate of the\n*     reciprocal Skeel condition number of the matrix A after\n*     equilibration (if done).  If this is less than the machine\n*     precision (in particular, if it is zero), the matrix is singular\n*     to working precision.  Note that the error may still be small even\n*     if this number is very small and the matrix appears ill-\n*     conditioned.\n*\n*     BERR    (output) DOUBLE PRECISION array, dimension (NRHS)\n*     Componentwise relative backward error.  This is the\n*     componentwise relative backward error of each solution vector X(j)\n*     (i.e., the smallest relative change in any element of A or B that\n*     makes X(j) an exact solution).\n*\n*     N_ERR_BNDS (input) INTEGER\n*     Number of error bounds to return for each right hand side\n*     and each type (normwise or componentwise).  See ERR_BNDS_NORM and\n*     ERR_BNDS_COMP below.\n*\n*     ERR_BNDS_NORM  (output) DOUBLE PRECISION array, dimension (NRHS, N_ERR_BNDS)\n*     For each right-hand side, this array contains information about\n*     various error bounds and condition numbers corresponding to the\n*     normwise relative error, which is defined as follows:\n*\n*     Normwise relative error in the ith solution vector:\n*             max_j (abs(XTRUE(j,i) - X(j,i)))\n*            ------------------------------\n*                  max_j abs(X(j,i))\n*\n*     The array is indexed by the type of error information as described\n*     below. There currently are up to three pieces of information\n*     returned.\n*\n*     The first index in ERR_BNDS_NORM(i,:) corresponds to the ith\n*     right-hand side.\n*\n*     The second index in ERR_BNDS_NORM(:,err) contains the following\n*     three fields:\n*     err = 1 \"Trust/don't trust\" boolean. Trust the answer if the\n*              reciprocal condition number is less than the threshold\n*              sqrt(n) * dlamch('Epsilon').\n*\n*     err = 2 \"Guaranteed\" error bound: The estimated forward error,\n*              almost certainly within a factor of 10 of the true error\n*              so long as the next entry is greater than the threshold\n*              sqrt(n) * dlamch('Epsilon'). This error bound should only\n*              be trusted if the previous boolean is true.\n*\n*     err = 3  Reciprocal condition number: Estimated normwise\n*              reciprocal condition number.  Compared with the threshold\n*              sqrt(n) * dlamch('Epsilon') to determine if the error\n*              estimate is \"guaranteed\". These reciprocal condition\n*              numbers are 1 / (norm(Z^{-1},inf) * norm(Z,inf)) for some\n*              appropriately scaled matrix Z.\n*              Let Z = S*A, where S scales each row by a power of the\n*              radix so all absolute row sums of Z are approximately 1.\n*\n*     See Lapack Working Note 165 for further details and extra\n*     cautions.\n*\n*     ERR_BNDS_COMP  (output) DOUBLE PRECISION array, dimension (NRHS, N_ERR_BNDS)\n*     For each right-hand side, this array contains information about\n*     various error bounds and condition numbers corresponding to the\n*     componentwise relative error, which is defined as follows:\n*\n*     Componentwise relative error in the ith solution vector:\n*                    abs(XTRUE(j,i) - X(j,i))\n*             max_j ----------------------\n*                         abs(X(j,i))\n*\n*     The array is indexed by the right-hand side i (on which the\n*     componentwise relative error depends), and the type of error\n*     information as described below. There currently are up to three\n*     pieces of information returned for each right-hand side. If\n*     componentwise accuracy is not requested (PARAMS(3) = 0.0), then\n*     ERR_BNDS_COMP is not accessed.  If N_ERR_BNDS .LT. 3, then at most\n*     the first (:,N_ERR_BNDS) entries are returned.\n*\n*     The first index in ERR_BNDS_COMP(i,:) corresponds to the ith\n*     right-hand side.\n*\n*     The second index in ERR_BNDS_COMP(:,err) contains the following\n*     three fields:\n*     err = 1 \"Trust/don't trust\" boolean. Trust the answer if the\n*              reciprocal condition number is less than the threshold\n*              sqrt(n) * dlamch('Epsilon').\n*\n*     err = 2 \"Guaranteed\" error bound: The estimated forward error,\n*              almost certainly within a factor of 10 of the true error\n*              so long as the next entry is greater than the threshold\n*              sqrt(n) * dlamch('Epsilon'). This error bound should only\n*              be trusted if the previous boolean is true.\n*\n*     err = 3  Reciprocal condition number: Estimated componentwise\n*              reciprocal condition number.  Compared with the threshold\n*              sqrt(n) * dlamch('Epsilon') to determine if the error\n*              estimate is \"guaranteed\". These reciprocal condition\n*              numbers are 1 / (norm(Z^{-1},inf) * norm(Z,inf)) for some\n*              appropriately scaled matrix Z.\n*              Let Z = S*(A*diag(x)), where x is the solution for the\n*              current right-hand side and S scales each row of\n*              A*diag(x) by a power of the radix so all absolute row\n*              sums of Z are approximately 1.\n*\n*     See Lapack Working Note 165 for further details and extra\n*     cautions.\n*\n*     NPARAMS (input) INTEGER\n*     Specifies the number of parameters set in PARAMS.  If .LE. 0, the\n*     PARAMS array is never referenced and default values are used.\n*\n*     PARAMS  (input / output) DOUBLE PRECISION array, dimension NPARAMS\n*     Specifies algorithm parameters.  If an entry is .LT. 0.0, then\n*     that entry will be filled with default value used for that\n*     parameter.  Only positions up to NPARAMS are accessed; defaults\n*     are used for higher-numbered parameters.\n*\n*       PARAMS(LA_LINRX_ITREF_I = 1) : Whether to perform iterative\n*            refinement or not.\n*         Default: 1.0D+0\n*            = 0.0 : No refinement is performed, and no error bounds are\n*                    computed.\n*            = 1.0 : Use the double-precision refinement algorithm,\n*                    possibly with doubled-single computations if the\n*                    compilation environment does not support DOUBLE\n*                    PRECISION.\n*              (other values are reserved for future use)\n*\n*       PARAMS(LA_LINRX_ITHRESH_I = 2) : Maximum number of residual\n*            computations allowed for refinement.\n*         Default: 10\n*         Aggressive: Set to 100 to permit convergence using approximate\n*                     factorizations or factorizations other than LU. If\n*                     the factorization uses a technique other than\n*                     Gaussian elimination, the guarantees in\n*                     err_bnds_norm and err_bnds_comp may no longer be\n*                     trustworthy.\n*\n*       PARAMS(LA_LINRX_CWISE_I = 3) : Flag determining if the code\n*            will attempt to find a solution with small componentwise\n*            relative error in the double-precision algorithm.  Positive\n*            is true, 0.0 is false.\n*         Default: 1.0 (attempt componentwise convergence)\n*\n*     WORK    (workspace) COMPLEX*16 array, dimension (2*N)\n*\n*     RWORK   (workspace) DOUBLE PRECISION array, dimension (2*N)\n*\n*     INFO    (output) INTEGER\n*       = 0:  Successful exit. The solution to every right-hand side is\n*         guaranteed.\n*       < 0:  If INFO = -i, the i-th argument had an illegal value\n*       > 0 and <= N:  U(INFO,INFO) is exactly zero.  The factorization\n*         has been completed, but the factor U is exactly singular, so\n*         the solution and error bounds could not be computed. RCOND = 0\n*         is returned.\n*       = N+J: The solution corresponding to the Jth right-hand side is\n*         not guaranteed. The solutions corresponding to other right-\n*         hand sides K with K > J may not be guaranteed as well, but\n*         only the first such right-hand side is reported. If a small\n*         componentwise error is not requested (PARAMS(3) = 0.0) then\n*         the Jth right-hand side is the first with a normwise error\n*         bound that is not guaranteed (the smallest J such\n*         that ERR_BNDS_NORM(J,1) = 0.0). By default (PARAMS(3) = 1.0)\n*         the Jth right-hand side is the first with either a normwise or\n*         componentwise error bound that is not guaranteed (the smallest\n*         J such that either ERR_BNDS_NORM(J,1) = 0.0 or\n*         ERR_BNDS_COMP(J,1) = 0.0). See the definition of\n*         ERR_BNDS_NORM(:,1) and ERR_BNDS_COMP(:,1). To get information\n*         about all of the right-hand sides check ERR_BNDS_NORM or\n*         ERR_BNDS_COMP.\n*\n\n*     ==================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  rcond, berr, err_bnds_norm, err_bnds_comp, info, s, x, params = NumRu::Lapack.zherfsx( uplo, equed, a, af, ipiv, s, b, x, params, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 9)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 9)", argc);
  rblapack_uplo = argv[0];
  rblapack_equed = argv[1];
  rblapack_a = argv[2];
  rblapack_af = argv[3];
  rblapack_ipiv = argv[4];
  rblapack_s = argv[5];
  rblapack_b = argv[6];
  rblapack_x = argv[7];
  rblapack_params = argv[8];
  if (rb_options != Qnil) {
  }

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
  if (NA_TYPE(rblapack_a) != NA_DCOMPLEX)
    rblapack_a = na_change_type(rblapack_a, NA_DCOMPLEX);
  a = NA_PTR_TYPE(rblapack_a, doublecomplex*);
  if (!NA_IsNArray(rblapack_x))
    rb_raise(rb_eArgError, "x (8th argument) must be NArray");
  if (NA_RANK(rblapack_x) != 2)
    rb_raise(rb_eArgError, "rank of x (8th argument) must be %d", 2);
  nrhs = NA_SHAPE1(rblapack_x);
  ldx = NA_SHAPE0(rblapack_x);
  if (NA_TYPE(rblapack_x) != NA_DCOMPLEX)
    rblapack_x = na_change_type(rblapack_x, NA_DCOMPLEX);
  x = NA_PTR_TYPE(rblapack_x, doublecomplex*);
  if (!NA_IsNArray(rblapack_b))
    rb_raise(rb_eArgError, "b (7th argument) must be NArray");
  if (NA_RANK(rblapack_b) != 2)
    rb_raise(rb_eArgError, "rank of b (7th argument) must be %d", 2);
  if (NA_SHAPE1(rblapack_b) != nrhs)
    rb_raise(rb_eRuntimeError, "shape 1 of b must be the same as shape 1 of x");
  ldb = NA_SHAPE0(rblapack_b);
  if (NA_TYPE(rblapack_b) != NA_DCOMPLEX)
    rblapack_b = na_change_type(rblapack_b, NA_DCOMPLEX);
  b = NA_PTR_TYPE(rblapack_b, doublecomplex*);
  if (!NA_IsNArray(rblapack_af))
    rb_raise(rb_eArgError, "af (4th argument) must be NArray");
  if (NA_RANK(rblapack_af) != 2)
    rb_raise(rb_eArgError, "rank of af (4th argument) must be %d", 2);
  if (NA_SHAPE1(rblapack_af) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of af must be the same as shape 0 of ipiv");
  ldaf = NA_SHAPE0(rblapack_af);
  if (NA_TYPE(rblapack_af) != NA_DCOMPLEX)
    rblapack_af = na_change_type(rblapack_af, NA_DCOMPLEX);
  af = NA_PTR_TYPE(rblapack_af, doublecomplex*);
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
    rb_raise(rb_eArgError, "s (6th argument) must be NArray");
  if (NA_RANK(rblapack_s) != 1)
    rb_raise(rb_eArgError, "rank of s (6th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_s) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of s must be the same as shape 0 of ipiv");
  if (NA_TYPE(rblapack_s) != NA_DFLOAT)
    rblapack_s = na_change_type(rblapack_s, NA_DFLOAT);
  s = NA_PTR_TYPE(rblapack_s, doublereal*);
  equed = StringValueCStr(rblapack_equed)[0];
  uplo = StringValueCStr(rblapack_uplo)[0];
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
    rblapack_s_out__ = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  s_out__ = NA_PTR_TYPE(rblapack_s_out__, doublereal*);
  MEMCPY(s_out__, s, doublereal, NA_TOTAL(rblapack_s));
  rblapack_s = rblapack_s_out__;
  s = s_out__;
  {
    int shape[2];
    shape[0] = ldx;
    shape[1] = nrhs;
    rblapack_x_out__ = na_make_object(NA_DCOMPLEX, 2, shape, cNArray);
  }
  x_out__ = NA_PTR_TYPE(rblapack_x_out__, doublecomplex*);
  MEMCPY(x_out__, x, doublecomplex, NA_TOTAL(rblapack_x));
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
  work = ALLOC_N(doublecomplex, (2*n));
  rwork = ALLOC_N(doublereal, (2*n));

  zherfsx_(&uplo, &equed, &n, &nrhs, a, &lda, af, &ldaf, ipiv, s, b, &ldb, x, &ldx, &rcond, berr, &n_err_bnds, err_bnds_norm, err_bnds_comp, &nparams, params, work, rwork, &info);

  free(work);
  free(rwork);
  rblapack_rcond = rb_float_new((double)rcond);
  rblapack_info = INT2NUM(info);
  return rb_ary_new3(8, rblapack_rcond, rblapack_berr, rblapack_err_bnds_norm, rblapack_err_bnds_comp, rblapack_info, rblapack_s, rblapack_x, rblapack_params);
}

void
init_lapack_zherfsx(VALUE mLapack){
  rb_define_module_function(mLapack, "zherfsx", rblapack_zherfsx, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}

#include "rb_lapack.h"

extern VOID sla_porfsx_extended_(integer *prec_type, char *uplo, integer *n, integer *nrhs, real *a, integer *lda, real *af, integer *ldaf, logical *colequ, real *c, real *b, integer *ldb, real *y, integer *ldy, real *berr_out, integer *n_norms, real *err_bnds_norm, real *err_bnds_comp, real *res, real *ayb, real *dy, real *y_tail, real *rcond, integer *ithresh, real *rthresh, real *dz_ub, logical *ignore_cwise, integer *info);

static VALUE
rb_sla_porfsx_extended(int argc, VALUE *argv, VALUE self){
  VALUE rb_prec_type;
  integer prec_type; 
  VALUE rb_uplo;
  char uplo; 
  VALUE rb_a;
  real *a; 
  VALUE rb_af;
  real *af; 
  VALUE rb_colequ;
  logical colequ; 
  VALUE rb_c;
  real *c; 
  VALUE rb_b;
  real *b; 
  VALUE rb_y;
  real *y; 
  VALUE rb_n_norms;
  integer n_norms; 
  VALUE rb_err_bnds_norm;
  real *err_bnds_norm; 
  VALUE rb_err_bnds_comp;
  real *err_bnds_comp; 
  VALUE rb_res;
  real *res; 
  VALUE rb_ayb;
  real *ayb; 
  VALUE rb_dy;
  real *dy; 
  VALUE rb_y_tail;
  real *y_tail; 
  VALUE rb_rcond;
  real rcond; 
  VALUE rb_ithresh;
  integer ithresh; 
  VALUE rb_rthresh;
  real rthresh; 
  VALUE rb_dz_ub;
  real dz_ub; 
  VALUE rb_ignore_cwise;
  logical ignore_cwise; 
  VALUE rb_berr_out;
  real *berr_out; 
  VALUE rb_info;
  integer info; 
  VALUE rb_y_out__;
  real *y_out__;
  VALUE rb_err_bnds_norm_out__;
  real *err_bnds_norm_out__;
  VALUE rb_err_bnds_comp_out__;
  real *err_bnds_comp_out__;

  integer lda;
  integer n;
  integer ldaf;
  integer ldb;
  integer nrhs;
  integer ldy;
  integer n_err_bnds;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  berr_out, info, y, err_bnds_norm, err_bnds_comp = NumRu::Lapack.sla_porfsx_extended( prec_type, uplo, a, af, colequ, c, b, y, n_norms, err_bnds_norm, err_bnds_comp, res, ayb, dy, y_tail, rcond, ithresh, rthresh, dz_ub, ignore_cwise)\n    or\n  NumRu::Lapack.sla_porfsx_extended  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE SLA_PORFSX_EXTENDED( PREC_TYPE, UPLO, N, NRHS, A, LDA, AF, LDAF, COLEQU, C, B, LDB, Y, LDY, BERR_OUT, N_NORMS, ERR_BNDS_NORM, ERR_BNDS_COMP, RES, AYB, DY, Y_TAIL, RCOND, ITHRESH, RTHRESH, DZ_UB, IGNORE_CWISE, INFO )\n\n*  Purpose\n*  =======\n*\n*  SLA_PORFSX_EXTENDED improves the computed solution to a system of\n*  linear equations by performing extra-precise iterative refinement\n*  and provides error bounds and backward error estimates for the solution.\n*  This subroutine is called by SPORFSX to perform iterative refinement.\n*  In addition to normwise error bound, the code provides maximum\n*  componentwise error bound if possible. See comments for ERR_BNDS_NORM\n*  and ERR_BNDS_COMP for details of the error bounds. Note that this\n*  subroutine is only resonsible for setting the second fields of\n*  ERR_BNDS_NORM and ERR_BNDS_COMP.\n*\n\n*  Arguments\n*  =========\n*\n*     PREC_TYPE      (input) INTEGER\n*     Specifies the intermediate precision to be used in refinement.\n*     The value is defined by ILAPREC(P) where P is a CHARACTER and\n*     P    = 'S':  Single\n*          = 'D':  Double\n*          = 'I':  Indigenous\n*          = 'X', 'E':  Extra\n*\n*     UPLO    (input) CHARACTER*1\n*       = 'U':  Upper triangle of A is stored;\n*       = 'L':  Lower triangle of A is stored.\n*\n*     N              (input) INTEGER\n*     The number of linear equations, i.e., the order of the\n*     matrix A.  N >= 0.\n*\n*     NRHS           (input) INTEGER\n*     The number of right-hand-sides, i.e., the number of columns of the\n*     matrix B.\n*\n*     A              (input) REAL array, dimension (LDA,N)\n*     On entry, the N-by-N matrix A.\n*\n*     LDA            (input) INTEGER\n*     The leading dimension of the array A.  LDA >= max(1,N).\n*\n*     AF             (input) REAL array, dimension (LDAF,N)\n*     The triangular factor U or L from the Cholesky factorization\n*     A = U**T*U or A = L*L**T, as computed by SPOTRF.\n*\n*     LDAF           (input) INTEGER\n*     The leading dimension of the array AF.  LDAF >= max(1,N).\n*\n*     COLEQU         (input) LOGICAL\n*     If .TRUE. then column equilibration was done to A before calling\n*     this routine. This is needed to compute the solution and error\n*     bounds correctly.\n*\n*     C              (input) REAL array, dimension (N)\n*     The column scale factors for A. If COLEQU = .FALSE., C\n*     is not accessed. If C is input, each element of C should be a power\n*     of the radix to ensure a reliable solution and error estimates.\n*     Scaling by powers of the radix does not cause rounding errors unless\n*     the result underflows or overflows. Rounding errors during scaling\n*     lead to refining with a matrix that is not equivalent to the\n*     input matrix, producing error estimates that may not be\n*     reliable.\n*\n*     B              (input) REAL array, dimension (LDB,NRHS)\n*     The right-hand-side matrix B.\n*\n*     LDB            (input) INTEGER\n*     The leading dimension of the array B.  LDB >= max(1,N).\n*\n*     Y              (input/output) REAL array, dimension (LDY,NRHS)\n*     On entry, the solution matrix X, as computed by SPOTRS.\n*     On exit, the improved solution matrix Y.\n*\n*     LDY            (input) INTEGER\n*     The leading dimension of the array Y.  LDY >= max(1,N).\n*\n*     BERR_OUT       (output) REAL array, dimension (NRHS)\n*     On exit, BERR_OUT(j) contains the componentwise relative backward\n*     error for right-hand-side j from the formula\n*         max(i) ( abs(RES(i)) / ( abs(op(A_s))*abs(Y) + abs(B_s) )(i) )\n*     where abs(Z) is the componentwise absolute value of the matrix\n*     or vector Z. This is computed by SLA_LIN_BERR.\n*\n*     N_NORMS        (input) INTEGER\n*     Determines which error bounds to return (see ERR_BNDS_NORM\n*     and ERR_BNDS_COMP).\n*     If N_NORMS >= 1 return normwise error bounds.\n*     If N_NORMS >= 2 return componentwise error bounds.\n*\n*     ERR_BNDS_NORM  (input/output) REAL array, dimension (NRHS, N_ERR_BNDS)\n*     For each right-hand side, this array contains information about\n*     various error bounds and condition numbers corresponding to the\n*     normwise relative error, which is defined as follows:\n*\n*     Normwise relative error in the ith solution vector:\n*             max_j (abs(XTRUE(j,i) - X(j,i)))\n*            ------------------------------\n*                  max_j abs(X(j,i))\n*\n*     The array is indexed by the type of error information as described\n*     below. There currently are up to three pieces of information\n*     returned.\n*\n*     The first index in ERR_BNDS_NORM(i,:) corresponds to the ith\n*     right-hand side.\n*\n*     The second index in ERR_BNDS_NORM(:,err) contains the following\n*     three fields:\n*     err = 1 \"Trust/don't trust\" boolean. Trust the answer if the\n*              reciprocal condition number is less than the threshold\n*              sqrt(n) * slamch('Epsilon').\n*\n*     err = 2 \"Guaranteed\" error bound: The estimated forward error,\n*              almost certainly within a factor of 10 of the true error\n*              so long as the next entry is greater than the threshold\n*              sqrt(n) * slamch('Epsilon'). This error bound should only\n*              be trusted if the previous boolean is true.\n*\n*     err = 3  Reciprocal condition number: Estimated normwise\n*              reciprocal condition number.  Compared with the threshold\n*              sqrt(n) * slamch('Epsilon') to determine if the error\n*              estimate is \"guaranteed\". These reciprocal condition\n*              numbers are 1 / (norm(Z^{-1},inf) * norm(Z,inf)) for some\n*              appropriately scaled matrix Z.\n*              Let Z = S*A, where S scales each row by a power of the\n*              radix so all absolute row sums of Z are approximately 1.\n*\n*     This subroutine is only responsible for setting the second field\n*     above.\n*     See Lapack Working Note 165 for further details and extra\n*     cautions.\n*\n*     ERR_BNDS_COMP  (input/output) REAL array, dimension (NRHS, N_ERR_BNDS)\n*     For each right-hand side, this array contains information about\n*     various error bounds and condition numbers corresponding to the\n*     componentwise relative error, which is defined as follows:\n*\n*     Componentwise relative error in the ith solution vector:\n*                    abs(XTRUE(j,i) - X(j,i))\n*             max_j ----------------------\n*                         abs(X(j,i))\n*\n*     The array is indexed by the right-hand side i (on which the\n*     componentwise relative error depends), and the type of error\n*     information as described below. There currently are up to three\n*     pieces of information returned for each right-hand side. If\n*     componentwise accuracy is not requested (PARAMS(3) = 0.0), then\n*     ERR_BNDS_COMP is not accessed.  If N_ERR_BNDS .LT. 3, then at most\n*     the first (:,N_ERR_BNDS) entries are returned.\n*\n*     The first index in ERR_BNDS_COMP(i,:) corresponds to the ith\n*     right-hand side.\n*\n*     The second index in ERR_BNDS_COMP(:,err) contains the following\n*     three fields:\n*     err = 1 \"Trust/don't trust\" boolean. Trust the answer if the\n*              reciprocal condition number is less than the threshold\n*              sqrt(n) * slamch('Epsilon').\n*\n*     err = 2 \"Guaranteed\" error bound: The estimated forward error,\n*              almost certainly within a factor of 10 of the true error\n*              so long as the next entry is greater than the threshold\n*              sqrt(n) * slamch('Epsilon'). This error bound should only\n*              be trusted if the previous boolean is true.\n*\n*     err = 3  Reciprocal condition number: Estimated componentwise\n*              reciprocal condition number.  Compared with the threshold\n*              sqrt(n) * slamch('Epsilon') to determine if the error\n*              estimate is \"guaranteed\". These reciprocal condition\n*              numbers are 1 / (norm(Z^{-1},inf) * norm(Z,inf)) for some\n*              appropriately scaled matrix Z.\n*              Let Z = S*(A*diag(x)), where x is the solution for the\n*              current right-hand side and S scales each row of\n*              A*diag(x) by a power of the radix so all absolute row\n*              sums of Z are approximately 1.\n*\n*     This subroutine is only responsible for setting the second field\n*     above.\n*     See Lapack Working Note 165 for further details and extra\n*     cautions.\n*\n*     RES            (input) REAL array, dimension (N)\n*     Workspace to hold the intermediate residual.\n*\n*     AYB            (input) REAL array, dimension (N)\n*     Workspace. This can be the same workspace passed for Y_TAIL.\n*\n*     DY             (input) REAL array, dimension (N)\n*     Workspace to hold the intermediate solution.\n*\n*     Y_TAIL         (input) REAL array, dimension (N)\n*     Workspace to hold the trailing bits of the intermediate solution.\n*\n*     RCOND          (input) REAL\n*     Reciprocal scaled condition number.  This is an estimate of the\n*     reciprocal Skeel condition number of the matrix A after\n*     equilibration (if done).  If this is less than the machine\n*     precision (in particular, if it is zero), the matrix is singular\n*     to working precision.  Note that the error may still be small even\n*     if this number is very small and the matrix appears ill-\n*     conditioned.\n*\n*     ITHRESH        (input) INTEGER\n*     The maximum number of residual computations allowed for\n*     refinement. The default is 10. For 'aggressive' set to 100 to\n*     permit convergence using approximate factorizations or\n*     factorizations other than LU. If the factorization uses a\n*     technique other than Gaussian elimination, the guarantees in\n*     ERR_BNDS_NORM and ERR_BNDS_COMP may no longer be trustworthy.\n*\n*     RTHRESH        (input) REAL\n*     Determines when to stop refinement if the error estimate stops\n*     decreasing. Refinement will stop when the next solution no longer\n*     satisfies norm(dx_{i+1}) < RTHRESH * norm(dx_i) where norm(Z) is\n*     the infinity norm of Z. RTHRESH satisfies 0 < RTHRESH <= 1. The\n*     default value is 0.5. For 'aggressive' set to 0.9 to permit\n*     convergence on extremely ill-conditioned matrices. See LAWN 165\n*     for more details.\n*\n*     DZ_UB          (input) REAL\n*     Determines when to start considering componentwise convergence.\n*     Componentwise convergence is only considered after each component\n*     of the solution Y is stable, which we definte as the relative\n*     change in each component being less than DZ_UB. The default value\n*     is 0.25, requiring the first bit to be stable. See LAWN 165 for\n*     more details.\n*\n*     IGNORE_CWISE   (input) LOGICAL\n*     If .TRUE. then ignore componentwise convergence. Default value\n*     is .FALSE..\n*\n*     INFO           (output) INTEGER\n*       = 0:  Successful exit.\n*       < 0:  if INFO = -i, the ith argument to SPOTRS had an illegal\n*             value\n*\n\n*  =====================================================================\n*\n*     .. Local Scalars ..\n      INTEGER            UPLO2, CNT, I, J, X_STATE, Z_STATE\n      REAL               YK, DYK, YMIN, NORMY, NORMX, NORMDX, DXRAT,\n     $                   DZRAT, PREVNORMDX, PREV_DZ_Z, DXRATMAX,\n     $                   DZRATMAX, DX_X, DZ_Z, FINAL_DX_X, FINAL_DZ_Z,\n     $                   EPS, HUGEVAL, INCR_THRESH\n      LOGICAL            INCR_PREC\n*     ..\n\n");
    return Qnil;
  }
  if (argc != 20)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 20)", argc);
  rb_prec_type = argv[0];
  rb_uplo = argv[1];
  rb_a = argv[2];
  rb_af = argv[3];
  rb_colequ = argv[4];
  rb_c = argv[5];
  rb_b = argv[6];
  rb_y = argv[7];
  rb_n_norms = argv[8];
  rb_err_bnds_norm = argv[9];
  rb_err_bnds_comp = argv[10];
  rb_res = argv[11];
  rb_ayb = argv[12];
  rb_dy = argv[13];
  rb_y_tail = argv[14];
  rb_rcond = argv[15];
  rb_ithresh = argv[16];
  rb_rthresh = argv[17];
  rb_dz_ub = argv[18];
  rb_ignore_cwise = argv[19];

  rcond = (real)NUM2DBL(rb_rcond);
  if (!NA_IsNArray(rb_res))
    rb_raise(rb_eArgError, "res (12th argument) must be NArray");
  if (NA_RANK(rb_res) != 1)
    rb_raise(rb_eArgError, "rank of res (12th argument) must be %d", 1);
  n = NA_SHAPE0(rb_res);
  if (NA_TYPE(rb_res) != NA_SFLOAT)
    rb_res = na_change_type(rb_res, NA_SFLOAT);
  res = NA_PTR_TYPE(rb_res, real*);
  if (!NA_IsNArray(rb_a))
    rb_raise(rb_eArgError, "a (3th argument) must be NArray");
  if (NA_RANK(rb_a) != 2)
    rb_raise(rb_eArgError, "rank of a (3th argument) must be %d", 2);
  if (NA_SHAPE1(rb_a) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of a must be the same as shape 0 of res");
  lda = NA_SHAPE0(rb_a);
  if (NA_TYPE(rb_a) != NA_SFLOAT)
    rb_a = na_change_type(rb_a, NA_SFLOAT);
  a = NA_PTR_TYPE(rb_a, real*);
  if (!NA_IsNArray(rb_b))
    rb_raise(rb_eArgError, "b (7th argument) must be NArray");
  if (NA_RANK(rb_b) != 2)
    rb_raise(rb_eArgError, "rank of b (7th argument) must be %d", 2);
  nrhs = NA_SHAPE1(rb_b);
  ldb = NA_SHAPE0(rb_b);
  if (NA_TYPE(rb_b) != NA_SFLOAT)
    rb_b = na_change_type(rb_b, NA_SFLOAT);
  b = NA_PTR_TYPE(rb_b, real*);
  colequ = (rb_colequ == Qtrue);
  if (!NA_IsNArray(rb_y))
    rb_raise(rb_eArgError, "y (8th argument) must be NArray");
  if (NA_RANK(rb_y) != 2)
    rb_raise(rb_eArgError, "rank of y (8th argument) must be %d", 2);
  if (NA_SHAPE1(rb_y) != nrhs)
    rb_raise(rb_eRuntimeError, "shape 1 of y must be the same as shape 1 of b");
  ldy = NA_SHAPE0(rb_y);
  if (NA_TYPE(rb_y) != NA_SFLOAT)
    rb_y = na_change_type(rb_y, NA_SFLOAT);
  y = NA_PTR_TYPE(rb_y, real*);
  if (!NA_IsNArray(rb_c))
    rb_raise(rb_eArgError, "c (6th argument) must be NArray");
  if (NA_RANK(rb_c) != 1)
    rb_raise(rb_eArgError, "rank of c (6th argument) must be %d", 1);
  if (NA_SHAPE0(rb_c) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of c must be the same as shape 0 of res");
  if (NA_TYPE(rb_c) != NA_SFLOAT)
    rb_c = na_change_type(rb_c, NA_SFLOAT);
  c = NA_PTR_TYPE(rb_c, real*);
  dz_ub = (real)NUM2DBL(rb_dz_ub);
  if (!NA_IsNArray(rb_y_tail))
    rb_raise(rb_eArgError, "y_tail (15th argument) must be NArray");
  if (NA_RANK(rb_y_tail) != 1)
    rb_raise(rb_eArgError, "rank of y_tail (15th argument) must be %d", 1);
  if (NA_SHAPE0(rb_y_tail) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of y_tail must be the same as shape 0 of res");
  if (NA_TYPE(rb_y_tail) != NA_SFLOAT)
    rb_y_tail = na_change_type(rb_y_tail, NA_SFLOAT);
  y_tail = NA_PTR_TYPE(rb_y_tail, real*);
  if (!NA_IsNArray(rb_af))
    rb_raise(rb_eArgError, "af (4th argument) must be NArray");
  if (NA_RANK(rb_af) != 2)
    rb_raise(rb_eArgError, "rank of af (4th argument) must be %d", 2);
  if (NA_SHAPE1(rb_af) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of af must be the same as shape 0 of res");
  ldaf = NA_SHAPE0(rb_af);
  if (NA_TYPE(rb_af) != NA_SFLOAT)
    rb_af = na_change_type(rb_af, NA_SFLOAT);
  af = NA_PTR_TYPE(rb_af, real*);
  prec_type = NUM2INT(rb_prec_type);
  rthresh = (real)NUM2DBL(rb_rthresh);
  ithresh = NUM2INT(rb_ithresh);
  if (!NA_IsNArray(rb_err_bnds_comp))
    rb_raise(rb_eArgError, "err_bnds_comp (11th argument) must be NArray");
  if (NA_RANK(rb_err_bnds_comp) != 2)
    rb_raise(rb_eArgError, "rank of err_bnds_comp (11th argument) must be %d", 2);
  n_err_bnds = NA_SHAPE1(rb_err_bnds_comp);
  if (NA_SHAPE0(rb_err_bnds_comp) != nrhs)
    rb_raise(rb_eRuntimeError, "shape 0 of err_bnds_comp must be the same as shape 1 of b");
  if (NA_TYPE(rb_err_bnds_comp) != NA_SFLOAT)
    rb_err_bnds_comp = na_change_type(rb_err_bnds_comp, NA_SFLOAT);
  err_bnds_comp = NA_PTR_TYPE(rb_err_bnds_comp, real*);
  ignore_cwise = (rb_ignore_cwise == Qtrue);
  if (!NA_IsNArray(rb_dy))
    rb_raise(rb_eArgError, "dy (14th argument) must be NArray");
  if (NA_RANK(rb_dy) != 1)
    rb_raise(rb_eArgError, "rank of dy (14th argument) must be %d", 1);
  if (NA_SHAPE0(rb_dy) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of dy must be the same as shape 0 of res");
  if (NA_TYPE(rb_dy) != NA_SFLOAT)
    rb_dy = na_change_type(rb_dy, NA_SFLOAT);
  dy = NA_PTR_TYPE(rb_dy, real*);
  if (!NA_IsNArray(rb_ayb))
    rb_raise(rb_eArgError, "ayb (13th argument) must be NArray");
  if (NA_RANK(rb_ayb) != 1)
    rb_raise(rb_eArgError, "rank of ayb (13th argument) must be %d", 1);
  if (NA_SHAPE0(rb_ayb) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of ayb must be the same as shape 0 of res");
  if (NA_TYPE(rb_ayb) != NA_SFLOAT)
    rb_ayb = na_change_type(rb_ayb, NA_SFLOAT);
  ayb = NA_PTR_TYPE(rb_ayb, real*);
  if (!NA_IsNArray(rb_err_bnds_norm))
    rb_raise(rb_eArgError, "err_bnds_norm (10th argument) must be NArray");
  if (NA_RANK(rb_err_bnds_norm) != 2)
    rb_raise(rb_eArgError, "rank of err_bnds_norm (10th argument) must be %d", 2);
  if (NA_SHAPE1(rb_err_bnds_norm) != n_err_bnds)
    rb_raise(rb_eRuntimeError, "shape 1 of err_bnds_norm must be the same as shape 1 of err_bnds_comp");
  if (NA_SHAPE0(rb_err_bnds_norm) != nrhs)
    rb_raise(rb_eRuntimeError, "shape 0 of err_bnds_norm must be the same as shape 1 of b");
  if (NA_TYPE(rb_err_bnds_norm) != NA_SFLOAT)
    rb_err_bnds_norm = na_change_type(rb_err_bnds_norm, NA_SFLOAT);
  err_bnds_norm = NA_PTR_TYPE(rb_err_bnds_norm, real*);
  n_norms = NUM2INT(rb_n_norms);
  uplo = StringValueCStr(rb_uplo)[0];
  {
    int shape[1];
    shape[0] = nrhs;
    rb_berr_out = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  berr_out = NA_PTR_TYPE(rb_berr_out, real*);
  {
    int shape[2];
    shape[0] = ldy;
    shape[1] = nrhs;
    rb_y_out__ = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  y_out__ = NA_PTR_TYPE(rb_y_out__, real*);
  MEMCPY(y_out__, y, real, NA_TOTAL(rb_y));
  rb_y = rb_y_out__;
  y = y_out__;
  {
    int shape[2];
    shape[0] = nrhs;
    shape[1] = n_err_bnds;
    rb_err_bnds_norm_out__ = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  err_bnds_norm_out__ = NA_PTR_TYPE(rb_err_bnds_norm_out__, real*);
  MEMCPY(err_bnds_norm_out__, err_bnds_norm, real, NA_TOTAL(rb_err_bnds_norm));
  rb_err_bnds_norm = rb_err_bnds_norm_out__;
  err_bnds_norm = err_bnds_norm_out__;
  {
    int shape[2];
    shape[0] = nrhs;
    shape[1] = n_err_bnds;
    rb_err_bnds_comp_out__ = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  err_bnds_comp_out__ = NA_PTR_TYPE(rb_err_bnds_comp_out__, real*);
  MEMCPY(err_bnds_comp_out__, err_bnds_comp, real, NA_TOTAL(rb_err_bnds_comp));
  rb_err_bnds_comp = rb_err_bnds_comp_out__;
  err_bnds_comp = err_bnds_comp_out__;

  sla_porfsx_extended_(&prec_type, &uplo, &n, &nrhs, a, &lda, af, &ldaf, &colequ, c, b, &ldb, y, &ldy, berr_out, &n_norms, err_bnds_norm, err_bnds_comp, res, ayb, dy, y_tail, &rcond, &ithresh, &rthresh, &dz_ub, &ignore_cwise, &info);

  rb_info = INT2NUM(info);
  return rb_ary_new3(5, rb_berr_out, rb_info, rb_y, rb_err_bnds_norm, rb_err_bnds_comp);
}

void
init_lapack_sla_porfsx_extended(VALUE mLapack){
  rb_define_module_function(mLapack, "sla_porfsx_extended", rb_sla_porfsx_extended, -1);
}

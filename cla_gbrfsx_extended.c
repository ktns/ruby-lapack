#include "rb_lapack.h"

extern VOID cla_gbrfsx_extended_(integer *prec_type, integer *trans_type, integer *n, integer *kl, integer *ku, integer *nrhs, complex *ab, integer *ldab, complex *afb, integer *ldafb, integer *ipiv, logical *colequ, real *c, complex *b, integer *ldb, complex *y, integer *ldy, real *berr_out, integer *n_norms, real *err_bnds_norm, real *err_bnds_comp, complex *res, real *ayb, complex *dy, complex *y_tail, real *rcond, integer *ithresh, real *rthresh, real *dz_ub, logical *ignore_cwise, integer *info);

static VALUE sHelp, sUsage;

static VALUE
rblapack_cla_gbrfsx_extended(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_prec_type;
  integer prec_type; 
  VALUE rblapack_trans_type;
  integer trans_type; 
  VALUE rblapack_kl;
  integer kl; 
  VALUE rblapack_ku;
  integer ku; 
  VALUE rblapack_ab;
  complex *ab; 
  VALUE rblapack_ldab;
  integer ldab; 
  VALUE rblapack_afb;
  complex *afb; 
  VALUE rblapack_ldafb;
  integer ldafb; 
  VALUE rblapack_ipiv;
  integer *ipiv; 
  VALUE rblapack_colequ;
  logical colequ; 
  VALUE rblapack_c;
  real *c; 
  VALUE rblapack_b;
  complex *b; 
  VALUE rblapack_y;
  complex *y; 
  VALUE rblapack_n_norms;
  integer n_norms; 
  VALUE rblapack_err_bnds_norm;
  real *err_bnds_norm; 
  VALUE rblapack_err_bnds_comp;
  real *err_bnds_comp; 
  VALUE rblapack_res;
  complex *res; 
  VALUE rblapack_ayb;
  real *ayb; 
  VALUE rblapack_dy;
  complex *dy; 
  VALUE rblapack_y_tail;
  complex *y_tail; 
  VALUE rblapack_rcond;
  real rcond; 
  VALUE rblapack_ithresh;
  integer ithresh; 
  VALUE rblapack_rthresh;
  real rthresh; 
  VALUE rblapack_dz_ub;
  real dz_ub; 
  VALUE rblapack_ignore_cwise;
  logical ignore_cwise; 
  VALUE rblapack_berr_out;
  real *berr_out; 
  VALUE rblapack_info;
  integer info; 
  VALUE rblapack_y_out__;
  complex *y_out__;
  VALUE rblapack_err_bnds_norm_out__;
  real *err_bnds_norm_out__;
  VALUE rblapack_err_bnds_comp_out__;
  real *err_bnds_comp_out__;

  integer lda;
  integer n;
  integer ldaf;
  integer ldb;
  integer nrhs;
  integer ldy;
  integer n_err_bnds;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  berr_out, info, y, err_bnds_norm, err_bnds_comp = NumRu::Lapack.cla_gbrfsx_extended( prec_type, trans_type, kl, ku, ab, ldab, afb, ldafb, ipiv, colequ, c, b, y, n_norms, err_bnds_norm, err_bnds_comp, res, ayb, dy, y_tail, rcond, ithresh, rthresh, dz_ub, ignore_cwise, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE CLA_GBRFSX_EXTENDED ( PREC_TYPE, TRANS_TYPE, N, KL, KU, NRHS, AB, LDAB, AFB, LDAFB, IPIV, COLEQU, C, B, LDB, Y, LDY, BERR_OUT, N_NORMS, ERR_BNDS_NORM, ERR_BNDS_COMP, RES, AYB, DY, Y_TAIL, RCOND, ITHRESH, RTHRESH, DZ_UB, IGNORE_CWISE, INFO )\n\n*  Purpose\n*  =======\n*\n*  CLA_GBRFSX_EXTENDED improves the computed solution to a system of\n*  linear equations by performing extra-precise iterative refinement\n*  and provides error bounds and backward error estimates for the solution.\n*  This subroutine is called by CGBRFSX to perform iterative refinement.\n*  In addition to normwise error bound, the code provides maximum\n*  componentwise error bound if possible. See comments for ERR_BNDS_NORM\n*  and ERR_BNDS_COMP for details of the error bounds. Note that this\n*  subroutine is only resonsible for setting the second fields of\n*  ERR_BNDS_NORM and ERR_BNDS_COMP.\n*\n\n*  Arguments\n*  =========\n*\n*     PREC_TYPE      (input) INTEGER\n*     Specifies the intermediate precision to be used in refinement.\n*     The value is defined by ILAPREC(P) where P is a CHARACTER and\n*     P    = 'S':  Single\n*          = 'D':  Double\n*          = 'I':  Indigenous\n*          = 'X', 'E':  Extra\n*\n*     TRANS_TYPE     (input) INTEGER\n*     Specifies the transposition operation on A.\n*     The value is defined by ILATRANS(T) where T is a CHARACTER and\n*     T    = 'N':  No transpose\n*          = 'T':  Transpose\n*          = 'C':  Conjugate transpose\n*\n*     N              (input) INTEGER\n*     The number of linear equations, i.e., the order of the\n*     matrix A.  N >= 0.\n*\n*     KL             (input) INTEGER\n*     The number of subdiagonals within the band of A.  KL >= 0.\n*\n*     KU             (input) INTEGER\n*     The number of superdiagonals within the band of A.  KU >= 0\n*\n*     NRHS           (input) INTEGER\n*     The number of right-hand-sides, i.e., the number of columns of the\n*     matrix B.\n*\n*     AB             (input) COMPLEX array, dimension (LDA,N)\n*     On entry, the N-by-N matrix A.\n*\n*     LDAB           (input) INTEGER\n*     The leading dimension of the array A.  LDA >= max(1,N).\n*\n*     AFB            (input) COMPLEX array, dimension (LDAF,N)\n*     The factors L and U from the factorization\n*     A = P*L*U as computed by CGBTRF.\n*\n*     LDAFB          (input) INTEGER\n*     The leading dimension of the array AF.  LDAF >= max(1,N).\n*\n*     IPIV           (input) INTEGER array, dimension (N)\n*     The pivot indices from the factorization A = P*L*U\n*     as computed by CGBTRF; row i of the matrix was interchanged\n*     with row IPIV(i).\n*\n*     COLEQU         (input) LOGICAL\n*     If .TRUE. then column equilibration was done to A before calling\n*     this routine. This is needed to compute the solution and error\n*     bounds correctly.\n*\n*     C              (input) REAL array, dimension (N)\n*     The column scale factors for A. If COLEQU = .FALSE., C\n*     is not accessed. If C is input, each element of C should be a power\n*     of the radix to ensure a reliable solution and error estimates.\n*     Scaling by powers of the radix does not cause rounding errors unless\n*     the result underflows or overflows. Rounding errors during scaling\n*     lead to refining with a matrix that is not equivalent to the\n*     input matrix, producing error estimates that may not be\n*     reliable.\n*\n*     B              (input) COMPLEX array, dimension (LDB,NRHS)\n*     The right-hand-side matrix B.\n*\n*     LDB            (input) INTEGER\n*     The leading dimension of the array B.  LDB >= max(1,N).\n*\n*     Y              (input/output) COMPLEX array, dimension (LDY,NRHS)\n*     On entry, the solution matrix X, as computed by CGBTRS.\n*     On exit, the improved solution matrix Y.\n*\n*     LDY            (input) INTEGER\n*     The leading dimension of the array Y.  LDY >= max(1,N).\n*\n*     BERR_OUT       (output) REAL array, dimension (NRHS)\n*     On exit, BERR_OUT(j) contains the componentwise relative backward\n*     error for right-hand-side j from the formula\n*         max(i) ( abs(RES(i)) / ( abs(op(A_s))*abs(Y) + abs(B_s) )(i) )\n*     where abs(Z) is the componentwise absolute value of the matrix\n*     or vector Z. This is computed by CLA_LIN_BERR.\n*\n*     N_NORMS        (input) INTEGER\n*     Determines which error bounds to return (see ERR_BNDS_NORM\n*     and ERR_BNDS_COMP).\n*     If N_NORMS >= 1 return normwise error bounds.\n*     If N_NORMS >= 2 return componentwise error bounds.\n*\n*     ERR_BNDS_NORM  (input/output) REAL array, dimension\n*                    (NRHS, N_ERR_BNDS)\n*     For each right-hand side, this array contains information about\n*     various error bounds and condition numbers corresponding to the\n*     normwise relative error, which is defined as follows:\n*\n*     Normwise relative error in the ith solution vector:\n*             max_j (abs(XTRUE(j,i) - X(j,i)))\n*            ------------------------------\n*                  max_j abs(X(j,i))\n*\n*     The array is indexed by the type of error information as described\n*     below. There currently are up to three pieces of information\n*     returned.\n*\n*     The first index in ERR_BNDS_NORM(i,:) corresponds to the ith\n*     right-hand side.\n*\n*     The second index in ERR_BNDS_NORM(:,err) contains the following\n*     three fields:\n*     err = 1 \"Trust/don't trust\" boolean. Trust the answer if the\n*              reciprocal condition number is less than the threshold\n*              sqrt(n) * slamch('Epsilon').\n*\n*     err = 2 \"Guaranteed\" error bound: The estimated forward error,\n*              almost certainly within a factor of 10 of the true error\n*              so long as the next entry is greater than the threshold\n*              sqrt(n) * slamch('Epsilon'). This error bound should only\n*              be trusted if the previous boolean is true.\n*\n*     err = 3  Reciprocal condition number: Estimated normwise\n*              reciprocal condition number.  Compared with the threshold\n*              sqrt(n) * slamch('Epsilon') to determine if the error\n*              estimate is \"guaranteed\". These reciprocal condition\n*              numbers are 1 / (norm(Z^{-1},inf) * norm(Z,inf)) for some\n*              appropriately scaled matrix Z.\n*              Let Z = S*A, where S scales each row by a power of the\n*              radix so all absolute row sums of Z are approximately 1.\n*\n*     This subroutine is only responsible for setting the second field\n*     above.\n*     See Lapack Working Note 165 for further details and extra\n*     cautions.\n*\n*     ERR_BNDS_COMP  (input/output) REAL array, dimension\n*                    (NRHS, N_ERR_BNDS)\n*     For each right-hand side, this array contains information about\n*     various error bounds and condition numbers corresponding to the\n*     componentwise relative error, which is defined as follows:\n*\n*     Componentwise relative error in the ith solution vector:\n*                    abs(XTRUE(j,i) - X(j,i))\n*             max_j ----------------------\n*                         abs(X(j,i))\n*\n*     The array is indexed by the right-hand side i (on which the\n*     componentwise relative error depends), and the type of error\n*     information as described below. There currently are up to three\n*     pieces of information returned for each right-hand side. If\n*     componentwise accuracy is not requested (PARAMS(3) = 0.0), then\n*     ERR_BNDS_COMP is not accessed.  If N_ERR_BNDS .LT. 3, then at most\n*     the first (:,N_ERR_BNDS) entries are returned.\n*\n*     The first index in ERR_BNDS_COMP(i,:) corresponds to the ith\n*     right-hand side.\n*\n*     The second index in ERR_BNDS_COMP(:,err) contains the following\n*     three fields:\n*     err = 1 \"Trust/don't trust\" boolean. Trust the answer if the\n*              reciprocal condition number is less than the threshold\n*              sqrt(n) * slamch('Epsilon').\n*\n*     err = 2 \"Guaranteed\" error bound: The estimated forward error,\n*              almost certainly within a factor of 10 of the true error\n*              so long as the next entry is greater than the threshold\n*              sqrt(n) * slamch('Epsilon'). This error bound should only\n*              be trusted if the previous boolean is true.\n*\n*     err = 3  Reciprocal condition number: Estimated componentwise\n*              reciprocal condition number.  Compared with the threshold\n*              sqrt(n) * slamch('Epsilon') to determine if the error\n*              estimate is \"guaranteed\". These reciprocal condition\n*              numbers are 1 / (norm(Z^{-1},inf) * norm(Z,inf)) for some\n*              appropriately scaled matrix Z.\n*              Let Z = S*(A*diag(x)), where x is the solution for the\n*              current right-hand side and S scales each row of\n*              A*diag(x) by a power of the radix so all absolute row\n*              sums of Z are approximately 1.\n*\n*     This subroutine is only responsible for setting the second field\n*     above.\n*     See Lapack Working Note 165 for further details and extra\n*     cautions.\n*\n*     RES            (input) COMPLEX array, dimension (N)\n*     Workspace to hold the intermediate residual.\n*\n*     AYB            (input) REAL array, dimension (N)\n*     Workspace.\n*\n*     DY             (input) COMPLEX array, dimension (N)\n*     Workspace to hold the intermediate solution.\n*\n*     Y_TAIL         (input) COMPLEX array, dimension (N)\n*     Workspace to hold the trailing bits of the intermediate solution.\n*\n*     RCOND          (input) REAL\n*     Reciprocal scaled condition number.  This is an estimate of the\n*     reciprocal Skeel condition number of the matrix A after\n*     equilibration (if done).  If this is less than the machine\n*     precision (in particular, if it is zero), the matrix is singular\n*     to working precision.  Note that the error may still be small even\n*     if this number is very small and the matrix appears ill-\n*     conditioned.\n*\n*     ITHRESH        (input) INTEGER\n*     The maximum number of residual computations allowed for\n*     refinement. The default is 10. For 'aggressive' set to 100 to\n*     permit convergence using approximate factorizations or\n*     factorizations other than LU. If the factorization uses a\n*     technique other than Gaussian elimination, the guarantees in\n*     ERR_BNDS_NORM and ERR_BNDS_COMP may no longer be trustworthy.\n*\n*     RTHRESH        (input) REAL\n*     Determines when to stop refinement if the error estimate stops\n*     decreasing. Refinement will stop when the next solution no longer\n*     satisfies norm(dx_{i+1}) < RTHRESH * norm(dx_i) where norm(Z) is\n*     the infinity norm of Z. RTHRESH satisfies 0 < RTHRESH <= 1. The\n*     default value is 0.5. For 'aggressive' set to 0.9 to permit\n*     convergence on extremely ill-conditioned matrices. See LAWN 165\n*     for more details.\n*\n*     DZ_UB          (input) REAL\n*     Determines when to start considering componentwise convergence.\n*     Componentwise convergence is only considered after each component\n*     of the solution Y is stable, which we definte as the relative\n*     change in each component being less than DZ_UB. The default value\n*     is 0.25, requiring the first bit to be stable. See LAWN 165 for\n*     more details.\n*\n*     IGNORE_CWISE   (input) LOGICAL\n*     If .TRUE. then ignore componentwise convergence. Default value\n*     is .FALSE..\n*\n*     INFO           (output) INTEGER\n*       = 0:  Successful exit.\n*       < 0:  if INFO = -i, the ith argument to CGBTRS had an illegal\n*             value\n*\n\n*  =====================================================================\n*\n*     .. Local Scalars ..\n      CHARACTER          TRANS\n      INTEGER            CNT, I, J, M, X_STATE, Z_STATE, Y_PREC_STATE\n      REAL               YK, DYK, YMIN, NORMY, NORMX, NORMDX, DXRAT,\n     $                   DZRAT, PREVNORMDX, PREV_DZ_Z, DXRATMAX,\n     $                   DZRATMAX, DX_X, DZ_Z, FINAL_DX_X, FINAL_DZ_Z,\n     $                   EPS, HUGEVAL, INCR_THRESH\n      LOGICAL            INCR_PREC\n      COMPLEX            ZDUM\n*     ..\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  berr_out, info, y, err_bnds_norm, err_bnds_comp = NumRu::Lapack.cla_gbrfsx_extended( prec_type, trans_type, kl, ku, ab, ldab, afb, ldafb, ipiv, colequ, c, b, y, n_norms, err_bnds_norm, err_bnds_comp, res, ayb, dy, y_tail, rcond, ithresh, rthresh, dz_ub, ignore_cwise, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 25)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 25)", argc);
  rblapack_prec_type = argv[0];
  rblapack_trans_type = argv[1];
  rblapack_kl = argv[2];
  rblapack_ku = argv[3];
  rblapack_ab = argv[4];
  rblapack_ldab = argv[5];
  rblapack_afb = argv[6];
  rblapack_ldafb = argv[7];
  rblapack_ipiv = argv[8];
  rblapack_colequ = argv[9];
  rblapack_c = argv[10];
  rblapack_b = argv[11];
  rblapack_y = argv[12];
  rblapack_n_norms = argv[13];
  rblapack_err_bnds_norm = argv[14];
  rblapack_err_bnds_comp = argv[15];
  rblapack_res = argv[16];
  rblapack_ayb = argv[17];
  rblapack_dy = argv[18];
  rblapack_y_tail = argv[19];
  rblapack_rcond = argv[20];
  rblapack_ithresh = argv[21];
  rblapack_rthresh = argv[22];
  rblapack_dz_ub = argv[23];
  rblapack_ignore_cwise = argv[24];
  if (rb_options != Qnil) {
  }

  rcond = (real)NUM2DBL(rblapack_rcond);
  if (!NA_IsNArray(rblapack_res))
    rb_raise(rb_eArgError, "res (17th argument) must be NArray");
  if (NA_RANK(rblapack_res) != 1)
    rb_raise(rb_eArgError, "rank of res (17th argument) must be %d", 1);
  n = NA_SHAPE0(rblapack_res);
  if (NA_TYPE(rblapack_res) != NA_SCOMPLEX)
    rblapack_res = na_change_type(rblapack_res, NA_SCOMPLEX);
  res = NA_PTR_TYPE(rblapack_res, complex*);
  if (!NA_IsNArray(rblapack_ipiv))
    rb_raise(rb_eArgError, "ipiv (9th argument) must be NArray");
  if (NA_RANK(rblapack_ipiv) != 1)
    rb_raise(rb_eArgError, "rank of ipiv (9th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_ipiv) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of ipiv must be the same as shape 0 of res");
  if (NA_TYPE(rblapack_ipiv) != NA_LINT)
    rblapack_ipiv = na_change_type(rblapack_ipiv, NA_LINT);
  ipiv = NA_PTR_TYPE(rblapack_ipiv, integer*);
  if (!NA_IsNArray(rblapack_ab))
    rb_raise(rb_eArgError, "ab (5th argument) must be NArray");
  if (NA_RANK(rblapack_ab) != 2)
    rb_raise(rb_eArgError, "rank of ab (5th argument) must be %d", 2);
  if (NA_SHAPE1(rblapack_ab) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of ab must be the same as shape 0 of res");
  lda = NA_SHAPE0(rblapack_ab);
  if (NA_TYPE(rblapack_ab) != NA_SCOMPLEX)
    rblapack_ab = na_change_type(rblapack_ab, NA_SCOMPLEX);
  ab = NA_PTR_TYPE(rblapack_ab, complex*);
  kl = NUM2INT(rblapack_kl);
  if (!NA_IsNArray(rblapack_b))
    rb_raise(rb_eArgError, "b (12th argument) must be NArray");
  if (NA_RANK(rblapack_b) != 2)
    rb_raise(rb_eArgError, "rank of b (12th argument) must be %d", 2);
  nrhs = NA_SHAPE1(rblapack_b);
  ldb = NA_SHAPE0(rblapack_b);
  if (NA_TYPE(rblapack_b) != NA_SCOMPLEX)
    rblapack_b = na_change_type(rblapack_b, NA_SCOMPLEX);
  b = NA_PTR_TYPE(rblapack_b, complex*);
  colequ = (rblapack_colequ == Qtrue);
  if (!NA_IsNArray(rblapack_y))
    rb_raise(rb_eArgError, "y (13th argument) must be NArray");
  if (NA_RANK(rblapack_y) != 2)
    rb_raise(rb_eArgError, "rank of y (13th argument) must be %d", 2);
  if (NA_SHAPE1(rblapack_y) != nrhs)
    rb_raise(rb_eRuntimeError, "shape 1 of y must be the same as shape 1 of b");
  ldy = NA_SHAPE0(rblapack_y);
  if (NA_TYPE(rblapack_y) != NA_SCOMPLEX)
    rblapack_y = na_change_type(rblapack_y, NA_SCOMPLEX);
  y = NA_PTR_TYPE(rblapack_y, complex*);
  if (!NA_IsNArray(rblapack_c))
    rb_raise(rb_eArgError, "c (11th argument) must be NArray");
  if (NA_RANK(rblapack_c) != 1)
    rb_raise(rb_eArgError, "rank of c (11th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_c) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of c must be the same as shape 0 of res");
  if (NA_TYPE(rblapack_c) != NA_SFLOAT)
    rblapack_c = na_change_type(rblapack_c, NA_SFLOAT);
  c = NA_PTR_TYPE(rblapack_c, real*);
  dz_ub = (real)NUM2DBL(rblapack_dz_ub);
  if (!NA_IsNArray(rblapack_y_tail))
    rb_raise(rb_eArgError, "y_tail (20th argument) must be NArray");
  if (NA_RANK(rblapack_y_tail) != 1)
    rb_raise(rb_eArgError, "rank of y_tail (20th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_y_tail) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of y_tail must be the same as shape 0 of res");
  if (NA_TYPE(rblapack_y_tail) != NA_SCOMPLEX)
    rblapack_y_tail = na_change_type(rblapack_y_tail, NA_SCOMPLEX);
  y_tail = NA_PTR_TYPE(rblapack_y_tail, complex*);
  ku = NUM2INT(rblapack_ku);
  if (!NA_IsNArray(rblapack_err_bnds_norm))
    rb_raise(rb_eArgError, "err_bnds_norm (15th argument) must be NArray");
  if (NA_RANK(rblapack_err_bnds_norm) != 2)
    rb_raise(rb_eArgError, "rank of err_bnds_norm (15th argument) must be %d", 2);
  n_err_bnds = NA_SHAPE1(rblapack_err_bnds_norm);
  if (NA_SHAPE0(rblapack_err_bnds_norm) != nrhs)
    rb_raise(rb_eRuntimeError, "shape 0 of err_bnds_norm must be the same as shape 1 of b");
  if (NA_TYPE(rblapack_err_bnds_norm) != NA_SFLOAT)
    rblapack_err_bnds_norm = na_change_type(rblapack_err_bnds_norm, NA_SFLOAT);
  err_bnds_norm = NA_PTR_TYPE(rblapack_err_bnds_norm, real*);
  rthresh = (real)NUM2DBL(rblapack_rthresh);
  ithresh = NUM2INT(rblapack_ithresh);
  if (!NA_IsNArray(rblapack_err_bnds_comp))
    rb_raise(rb_eArgError, "err_bnds_comp (16th argument) must be NArray");
  if (NA_RANK(rblapack_err_bnds_comp) != 2)
    rb_raise(rb_eArgError, "rank of err_bnds_comp (16th argument) must be %d", 2);
  if (NA_SHAPE1(rblapack_err_bnds_comp) != n_err_bnds)
    rb_raise(rb_eRuntimeError, "shape 1 of err_bnds_comp must be the same as shape 1 of err_bnds_norm");
  if (NA_SHAPE0(rblapack_err_bnds_comp) != nrhs)
    rb_raise(rb_eRuntimeError, "shape 0 of err_bnds_comp must be the same as shape 1 of b");
  if (NA_TYPE(rblapack_err_bnds_comp) != NA_SFLOAT)
    rblapack_err_bnds_comp = na_change_type(rblapack_err_bnds_comp, NA_SFLOAT);
  err_bnds_comp = NA_PTR_TYPE(rblapack_err_bnds_comp, real*);
  if (!NA_IsNArray(rblapack_afb))
    rb_raise(rb_eArgError, "afb (7th argument) must be NArray");
  if (NA_RANK(rblapack_afb) != 2)
    rb_raise(rb_eArgError, "rank of afb (7th argument) must be %d", 2);
  if (NA_SHAPE1(rblapack_afb) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of afb must be the same as shape 0 of res");
  ldaf = NA_SHAPE0(rblapack_afb);
  if (NA_TYPE(rblapack_afb) != NA_SCOMPLEX)
    rblapack_afb = na_change_type(rblapack_afb, NA_SCOMPLEX);
  afb = NA_PTR_TYPE(rblapack_afb, complex*);
  n_norms = NUM2INT(rblapack_n_norms);
  ignore_cwise = (rblapack_ignore_cwise == Qtrue);
  trans_type = NUM2INT(rblapack_trans_type);
  if (!NA_IsNArray(rblapack_dy))
    rb_raise(rb_eArgError, "dy (19th argument) must be NArray");
  if (NA_RANK(rblapack_dy) != 1)
    rb_raise(rb_eArgError, "rank of dy (19th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_dy) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of dy must be the same as shape 0 of res");
  if (NA_TYPE(rblapack_dy) != NA_SCOMPLEX)
    rblapack_dy = na_change_type(rblapack_dy, NA_SCOMPLEX);
  dy = NA_PTR_TYPE(rblapack_dy, complex*);
  if (!NA_IsNArray(rblapack_ayb))
    rb_raise(rb_eArgError, "ayb (18th argument) must be NArray");
  if (NA_RANK(rblapack_ayb) != 1)
    rb_raise(rb_eArgError, "rank of ayb (18th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_ayb) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of ayb must be the same as shape 0 of res");
  if (NA_TYPE(rblapack_ayb) != NA_SFLOAT)
    rblapack_ayb = na_change_type(rblapack_ayb, NA_SFLOAT);
  ayb = NA_PTR_TYPE(rblapack_ayb, real*);
  prec_type = NUM2INT(rblapack_prec_type);
  ldab = lda = MAX(1,n);
  ldafb = ldaf = MAX(1,n);
  {
    int shape[1];
    shape[0] = nrhs;
    rblapack_berr_out = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  berr_out = NA_PTR_TYPE(rblapack_berr_out, real*);
  {
    int shape[2];
    shape[0] = ldy;
    shape[1] = nrhs;
    rblapack_y_out__ = na_make_object(NA_SCOMPLEX, 2, shape, cNArray);
  }
  y_out__ = NA_PTR_TYPE(rblapack_y_out__, complex*);
  MEMCPY(y_out__, y, complex, NA_TOTAL(rblapack_y));
  rblapack_y = rblapack_y_out__;
  y = y_out__;
  {
    int shape[2];
    shape[0] = nrhs;
    shape[1] = n_err_bnds;
    rblapack_err_bnds_norm_out__ = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  err_bnds_norm_out__ = NA_PTR_TYPE(rblapack_err_bnds_norm_out__, real*);
  MEMCPY(err_bnds_norm_out__, err_bnds_norm, real, NA_TOTAL(rblapack_err_bnds_norm));
  rblapack_err_bnds_norm = rblapack_err_bnds_norm_out__;
  err_bnds_norm = err_bnds_norm_out__;
  {
    int shape[2];
    shape[0] = nrhs;
    shape[1] = n_err_bnds;
    rblapack_err_bnds_comp_out__ = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  err_bnds_comp_out__ = NA_PTR_TYPE(rblapack_err_bnds_comp_out__, real*);
  MEMCPY(err_bnds_comp_out__, err_bnds_comp, real, NA_TOTAL(rblapack_err_bnds_comp));
  rblapack_err_bnds_comp = rblapack_err_bnds_comp_out__;
  err_bnds_comp = err_bnds_comp_out__;

  cla_gbrfsx_extended_(&prec_type, &trans_type, &n, &kl, &ku, &nrhs, ab, &ldab, afb, &ldafb, ipiv, &colequ, c, b, &ldb, y, &ldy, berr_out, &n_norms, err_bnds_norm, err_bnds_comp, res, ayb, dy, y_tail, &rcond, &ithresh, &rthresh, &dz_ub, &ignore_cwise, &info);

  rblapack_info = INT2NUM(info);
  return rb_ary_new3(5, rblapack_berr_out, rblapack_info, rblapack_y, rblapack_err_bnds_norm, rblapack_err_bnds_comp);
}

void
init_lapack_cla_gbrfsx_extended(VALUE mLapack){
  rb_define_module_function(mLapack, "cla_gbrfsx_extended", rblapack_cla_gbrfsx_extended, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}

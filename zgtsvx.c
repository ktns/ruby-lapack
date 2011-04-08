#include "rb_lapack.h"

extern VOID zgtsvx_(char *fact, char *trans, integer *n, integer *nrhs, doublecomplex *dl, doublecomplex *d, doublecomplex *du, doublecomplex *dlf, doublecomplex *df, doublecomplex *duf, doublecomplex *du2, integer *ipiv, doublecomplex *b, integer *ldb, doublecomplex *x, integer *ldx, doublereal *rcond, doublereal *ferr, doublereal *berr, doublecomplex *work, doublereal *rwork, integer *info);

static VALUE sHelp, sUsage;

static VALUE
rblapack_zgtsvx(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_fact;
  char fact; 
  VALUE rblapack_trans;
  char trans; 
  VALUE rblapack_dl;
  doublecomplex *dl; 
  VALUE rblapack_d;
  doublecomplex *d; 
  VALUE rblapack_du;
  doublecomplex *du; 
  VALUE rblapack_dlf;
  doublecomplex *dlf; 
  VALUE rblapack_df;
  doublecomplex *df; 
  VALUE rblapack_duf;
  doublecomplex *duf; 
  VALUE rblapack_du2;
  doublecomplex *du2; 
  VALUE rblapack_ipiv;
  integer *ipiv; 
  VALUE rblapack_b;
  doublecomplex *b; 
  VALUE rblapack_x;
  doublecomplex *x; 
  VALUE rblapack_rcond;
  doublereal rcond; 
  VALUE rblapack_ferr;
  doublereal *ferr; 
  VALUE rblapack_berr;
  doublereal *berr; 
  VALUE rblapack_info;
  integer info; 
  VALUE rblapack_dlf_out__;
  doublecomplex *dlf_out__;
  VALUE rblapack_df_out__;
  doublecomplex *df_out__;
  VALUE rblapack_duf_out__;
  doublecomplex *duf_out__;
  VALUE rblapack_du2_out__;
  doublecomplex *du2_out__;
  VALUE rblapack_ipiv_out__;
  integer *ipiv_out__;
  doublecomplex *work;
  doublereal *rwork;

  integer n;
  integer ldb;
  integer nrhs;
  integer ldx;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  x, rcond, ferr, berr, info, dlf, df, duf, du2, ipiv = NumRu::Lapack.zgtsvx( fact, trans, dl, d, du, dlf, df, duf, du2, ipiv, b, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE ZGTSVX( FACT, TRANS, N, NRHS, DL, D, DU, DLF, DF, DUF, DU2, IPIV, B, LDB, X, LDX, RCOND, FERR, BERR, WORK, RWORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  ZGTSVX uses the LU factorization to compute the solution to a complex\n*  system of linear equations A * X = B, A**T * X = B, or A**H * X = B,\n*  where A is a tridiagonal matrix of order N and X and B are N-by-NRHS\n*  matrices.\n*\n*  Error bounds on the solution and a condition estimate are also\n*  provided.\n*\n*  Description\n*  ===========\n*\n*  The following steps are performed:\n*\n*  1. If FACT = 'N', the LU decomposition is used to factor the matrix A\n*     as A = L * U, where L is a product of permutation and unit lower\n*     bidiagonal matrices and U is upper triangular with nonzeros in\n*     only the main diagonal and first two superdiagonals.\n*\n*  2. If some U(i,i)=0, so that U is exactly singular, then the routine\n*     returns with INFO = i. Otherwise, the factored form of A is used\n*     to estimate the condition number of the matrix A.  If the\n*     reciprocal of the condition number is less than machine precision,\n*     INFO = N+1 is returned as a warning, but the routine still goes on\n*     to solve for X and compute error bounds as described below.\n*\n*  3. The system of equations is solved for X using the factored form\n*     of A.\n*\n*  4. Iterative refinement is applied to improve the computed solution\n*     matrix and calculate error bounds and backward error estimates\n*     for it.\n*\n\n*  Arguments\n*  =========\n*\n*  FACT    (input) CHARACTER*1\n*          Specifies whether or not the factored form of A has been\n*          supplied on entry.\n*          = 'F':  DLF, DF, DUF, DU2, and IPIV contain the factored form\n*                  of A; DL, D, DU, DLF, DF, DUF, DU2 and IPIV will not\n*                  be modified.\n*          = 'N':  The matrix will be copied to DLF, DF, and DUF\n*                  and factored.\n*\n*  TRANS   (input) CHARACTER*1\n*          Specifies the form of the system of equations:\n*          = 'N':  A * X = B     (No transpose)\n*          = 'T':  A**T * X = B  (Transpose)\n*          = 'C':  A**H * X = B  (Conjugate transpose)\n*\n*  N       (input) INTEGER\n*          The order of the matrix A.  N >= 0.\n*\n*  NRHS    (input) INTEGER\n*          The number of right hand sides, i.e., the number of columns\n*          of the matrix B.  NRHS >= 0.\n*\n*  DL      (input) COMPLEX*16 array, dimension (N-1)\n*          The (n-1) subdiagonal elements of A.\n*\n*  D       (input) COMPLEX*16 array, dimension (N)\n*          The n diagonal elements of A.\n*\n*  DU      (input) COMPLEX*16 array, dimension (N-1)\n*          The (n-1) superdiagonal elements of A.\n*\n*  DLF     (input or output) COMPLEX*16 array, dimension (N-1)\n*          If FACT = 'F', then DLF is an input argument and on entry\n*          contains the (n-1) multipliers that define the matrix L from\n*          the LU factorization of A as computed by ZGTTRF.\n*\n*          If FACT = 'N', then DLF is an output argument and on exit\n*          contains the (n-1) multipliers that define the matrix L from\n*          the LU factorization of A.\n*\n*  DF      (input or output) COMPLEX*16 array, dimension (N)\n*          If FACT = 'F', then DF is an input argument and on entry\n*          contains the n diagonal elements of the upper triangular\n*          matrix U from the LU factorization of A.\n*\n*          If FACT = 'N', then DF is an output argument and on exit\n*          contains the n diagonal elements of the upper triangular\n*          matrix U from the LU factorization of A.\n*\n*  DUF     (input or output) COMPLEX*16 array, dimension (N-1)\n*          If FACT = 'F', then DUF is an input argument and on entry\n*          contains the (n-1) elements of the first superdiagonal of U.\n*\n*          If FACT = 'N', then DUF is an output argument and on exit\n*          contains the (n-1) elements of the first superdiagonal of U.\n*\n*  DU2     (input or output) COMPLEX*16 array, dimension (N-2)\n*          If FACT = 'F', then DU2 is an input argument and on entry\n*          contains the (n-2) elements of the second superdiagonal of\n*          U.\n*\n*          If FACT = 'N', then DU2 is an output argument and on exit\n*          contains the (n-2) elements of the second superdiagonal of\n*          U.\n*\n*  IPIV    (input or output) INTEGER array, dimension (N)\n*          If FACT = 'F', then IPIV is an input argument and on entry\n*          contains the pivot indices from the LU factorization of A as\n*          computed by ZGTTRF.\n*\n*          If FACT = 'N', then IPIV is an output argument and on exit\n*          contains the pivot indices from the LU factorization of A;\n*          row i of the matrix was interchanged with row IPIV(i).\n*          IPIV(i) will always be either i or i+1; IPIV(i) = i indicates\n*          a row interchange was not required.\n*\n*  B       (input) COMPLEX*16 array, dimension (LDB,NRHS)\n*          The N-by-NRHS right hand side matrix B.\n*\n*  LDB     (input) INTEGER\n*          The leading dimension of the array B.  LDB >= max(1,N).\n*\n*  X       (output) COMPLEX*16 array, dimension (LDX,NRHS)\n*          If INFO = 0 or INFO = N+1, the N-by-NRHS solution matrix X.\n*\n*  LDX     (input) INTEGER\n*          The leading dimension of the array X.  LDX >= max(1,N).\n*\n*  RCOND   (output) DOUBLE PRECISION\n*          The estimate of the reciprocal condition number of the matrix\n*          A.  If RCOND is less than the machine precision (in\n*          particular, if RCOND = 0), the matrix is singular to working\n*          precision.  This condition is indicated by a return code of\n*          INFO > 0.\n*\n*  FERR    (output) DOUBLE PRECISION array, dimension (NRHS)\n*          The estimated forward error bound for each solution vector\n*          X(j) (the j-th column of the solution matrix X).\n*          If XTRUE is the true solution corresponding to X(j), FERR(j)\n*          is an estimated upper bound for the magnitude of the largest\n*          element in (X(j) - XTRUE) divided by the magnitude of the\n*          largest element in X(j).  The estimate is as reliable as\n*          the estimate for RCOND, and is almost always a slight\n*          overestimate of the true error.\n*\n*  BERR    (output) DOUBLE PRECISION array, dimension (NRHS)\n*          The componentwise relative backward error of each solution\n*          vector X(j) (i.e., the smallest relative change in\n*          any element of A or B that makes X(j) an exact solution).\n*\n*  WORK    (workspace) COMPLEX*16 array, dimension (2*N)\n*\n*  RWORK   (workspace) DOUBLE PRECISION array, dimension (N)\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit\n*          < 0:  if INFO = -i, the i-th argument had an illegal value\n*          > 0:  if INFO = i, and i is\n*                <= N:  U(i,i) is exactly zero.  The factorization\n*                       has not been completed unless i = N, but the\n*                       factor U is exactly singular, so the solution\n*                       and error bounds could not be computed.\n*                       RCOND = 0 is returned.\n*                = N+1: U is nonsingular, but RCOND is less than machine\n*                       precision, meaning that the matrix is singular\n*                       to working precision.  Nevertheless, the\n*                       solution and error bounds are computed because\n*                       there are a number of situations where the\n*                       computed solution can be more accurate than the\n*                       value of RCOND would suggest.\n*\n\n*  =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  x, rcond, ferr, berr, info, dlf, df, duf, du2, ipiv = NumRu::Lapack.zgtsvx( fact, trans, dl, d, du, dlf, df, duf, du2, ipiv, b, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 11)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 11)", argc);
  rblapack_fact = argv[0];
  rblapack_trans = argv[1];
  rblapack_dl = argv[2];
  rblapack_d = argv[3];
  rblapack_du = argv[4];
  rblapack_dlf = argv[5];
  rblapack_df = argv[6];
  rblapack_duf = argv[7];
  rblapack_du2 = argv[8];
  rblapack_ipiv = argv[9];
  rblapack_b = argv[10];
  if (rb_options != Qnil) {
  }

  trans = StringValueCStr(rblapack_trans)[0];
  if (!NA_IsNArray(rblapack_ipiv))
    rb_raise(rb_eArgError, "ipiv (10th argument) must be NArray");
  if (NA_RANK(rblapack_ipiv) != 1)
    rb_raise(rb_eArgError, "rank of ipiv (10th argument) must be %d", 1);
  n = NA_SHAPE0(rblapack_ipiv);
  if (NA_TYPE(rblapack_ipiv) != NA_LINT)
    rblapack_ipiv = na_change_type(rblapack_ipiv, NA_LINT);
  ipiv = NA_PTR_TYPE(rblapack_ipiv, integer*);
  fact = StringValueCStr(rblapack_fact)[0];
  if (!NA_IsNArray(rblapack_b))
    rb_raise(rb_eArgError, "b (11th argument) must be NArray");
  if (NA_RANK(rblapack_b) != 2)
    rb_raise(rb_eArgError, "rank of b (11th argument) must be %d", 2);
  nrhs = NA_SHAPE1(rblapack_b);
  ldb = NA_SHAPE0(rblapack_b);
  if (NA_TYPE(rblapack_b) != NA_DCOMPLEX)
    rblapack_b = na_change_type(rblapack_b, NA_DCOMPLEX);
  b = NA_PTR_TYPE(rblapack_b, doublecomplex*);
  if (!NA_IsNArray(rblapack_df))
    rb_raise(rb_eArgError, "df (7th argument) must be NArray");
  if (NA_RANK(rblapack_df) != 1)
    rb_raise(rb_eArgError, "rank of df (7th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_df) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of df must be the same as shape 0 of ipiv");
  if (NA_TYPE(rblapack_df) != NA_DCOMPLEX)
    rblapack_df = na_change_type(rblapack_df, NA_DCOMPLEX);
  df = NA_PTR_TYPE(rblapack_df, doublecomplex*);
  if (!NA_IsNArray(rblapack_d))
    rb_raise(rb_eArgError, "d (4th argument) must be NArray");
  if (NA_RANK(rblapack_d) != 1)
    rb_raise(rb_eArgError, "rank of d (4th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_d) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of d must be the same as shape 0 of ipiv");
  if (NA_TYPE(rblapack_d) != NA_DCOMPLEX)
    rblapack_d = na_change_type(rblapack_d, NA_DCOMPLEX);
  d = NA_PTR_TYPE(rblapack_d, doublecomplex*);
  if (!NA_IsNArray(rblapack_du2))
    rb_raise(rb_eArgError, "du2 (9th argument) must be NArray");
  if (NA_RANK(rblapack_du2) != 1)
    rb_raise(rb_eArgError, "rank of du2 (9th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_du2) != (n-2))
    rb_raise(rb_eRuntimeError, "shape 0 of du2 must be %d", n-2);
  if (NA_TYPE(rblapack_du2) != NA_DCOMPLEX)
    rblapack_du2 = na_change_type(rblapack_du2, NA_DCOMPLEX);
  du2 = NA_PTR_TYPE(rblapack_du2, doublecomplex*);
  if (!NA_IsNArray(rblapack_du))
    rb_raise(rb_eArgError, "du (5th argument) must be NArray");
  if (NA_RANK(rblapack_du) != 1)
    rb_raise(rb_eArgError, "rank of du (5th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_du) != (n-1))
    rb_raise(rb_eRuntimeError, "shape 0 of du must be %d", n-1);
  if (NA_TYPE(rblapack_du) != NA_DCOMPLEX)
    rblapack_du = na_change_type(rblapack_du, NA_DCOMPLEX);
  du = NA_PTR_TYPE(rblapack_du, doublecomplex*);
  if (!NA_IsNArray(rblapack_dlf))
    rb_raise(rb_eArgError, "dlf (6th argument) must be NArray");
  if (NA_RANK(rblapack_dlf) != 1)
    rb_raise(rb_eArgError, "rank of dlf (6th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_dlf) != (n-1))
    rb_raise(rb_eRuntimeError, "shape 0 of dlf must be %d", n-1);
  if (NA_TYPE(rblapack_dlf) != NA_DCOMPLEX)
    rblapack_dlf = na_change_type(rblapack_dlf, NA_DCOMPLEX);
  dlf = NA_PTR_TYPE(rblapack_dlf, doublecomplex*);
  if (!NA_IsNArray(rblapack_dl))
    rb_raise(rb_eArgError, "dl (3th argument) must be NArray");
  if (NA_RANK(rblapack_dl) != 1)
    rb_raise(rb_eArgError, "rank of dl (3th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_dl) != (n-1))
    rb_raise(rb_eRuntimeError, "shape 0 of dl must be %d", n-1);
  if (NA_TYPE(rblapack_dl) != NA_DCOMPLEX)
    rblapack_dl = na_change_type(rblapack_dl, NA_DCOMPLEX);
  dl = NA_PTR_TYPE(rblapack_dl, doublecomplex*);
  if (!NA_IsNArray(rblapack_duf))
    rb_raise(rb_eArgError, "duf (8th argument) must be NArray");
  if (NA_RANK(rblapack_duf) != 1)
    rb_raise(rb_eArgError, "rank of duf (8th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_duf) != (n-1))
    rb_raise(rb_eRuntimeError, "shape 0 of duf must be %d", n-1);
  if (NA_TYPE(rblapack_duf) != NA_DCOMPLEX)
    rblapack_duf = na_change_type(rblapack_duf, NA_DCOMPLEX);
  duf = NA_PTR_TYPE(rblapack_duf, doublecomplex*);
  ldx = MAX(1,n);
  {
    int shape[2];
    shape[0] = ldx;
    shape[1] = nrhs;
    rblapack_x = na_make_object(NA_DCOMPLEX, 2, shape, cNArray);
  }
  x = NA_PTR_TYPE(rblapack_x, doublecomplex*);
  {
    int shape[1];
    shape[0] = nrhs;
    rblapack_ferr = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  ferr = NA_PTR_TYPE(rblapack_ferr, doublereal*);
  {
    int shape[1];
    shape[0] = nrhs;
    rblapack_berr = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  berr = NA_PTR_TYPE(rblapack_berr, doublereal*);
  {
    int shape[1];
    shape[0] = n-1;
    rblapack_dlf_out__ = na_make_object(NA_DCOMPLEX, 1, shape, cNArray);
  }
  dlf_out__ = NA_PTR_TYPE(rblapack_dlf_out__, doublecomplex*);
  MEMCPY(dlf_out__, dlf, doublecomplex, NA_TOTAL(rblapack_dlf));
  rblapack_dlf = rblapack_dlf_out__;
  dlf = dlf_out__;
  {
    int shape[1];
    shape[0] = n;
    rblapack_df_out__ = na_make_object(NA_DCOMPLEX, 1, shape, cNArray);
  }
  df_out__ = NA_PTR_TYPE(rblapack_df_out__, doublecomplex*);
  MEMCPY(df_out__, df, doublecomplex, NA_TOTAL(rblapack_df));
  rblapack_df = rblapack_df_out__;
  df = df_out__;
  {
    int shape[1];
    shape[0] = n-1;
    rblapack_duf_out__ = na_make_object(NA_DCOMPLEX, 1, shape, cNArray);
  }
  duf_out__ = NA_PTR_TYPE(rblapack_duf_out__, doublecomplex*);
  MEMCPY(duf_out__, duf, doublecomplex, NA_TOTAL(rblapack_duf));
  rblapack_duf = rblapack_duf_out__;
  duf = duf_out__;
  {
    int shape[1];
    shape[0] = n-2;
    rblapack_du2_out__ = na_make_object(NA_DCOMPLEX, 1, shape, cNArray);
  }
  du2_out__ = NA_PTR_TYPE(rblapack_du2_out__, doublecomplex*);
  MEMCPY(du2_out__, du2, doublecomplex, NA_TOTAL(rblapack_du2));
  rblapack_du2 = rblapack_du2_out__;
  du2 = du2_out__;
  {
    int shape[1];
    shape[0] = n;
    rblapack_ipiv_out__ = na_make_object(NA_LINT, 1, shape, cNArray);
  }
  ipiv_out__ = NA_PTR_TYPE(rblapack_ipiv_out__, integer*);
  MEMCPY(ipiv_out__, ipiv, integer, NA_TOTAL(rblapack_ipiv));
  rblapack_ipiv = rblapack_ipiv_out__;
  ipiv = ipiv_out__;
  work = ALLOC_N(doublecomplex, (2*n));
  rwork = ALLOC_N(doublereal, (n));

  zgtsvx_(&fact, &trans, &n, &nrhs, dl, d, du, dlf, df, duf, du2, ipiv, b, &ldb, x, &ldx, &rcond, ferr, berr, work, rwork, &info);

  free(work);
  free(rwork);
  rblapack_rcond = rb_float_new((double)rcond);
  rblapack_info = INT2NUM(info);
  return rb_ary_new3(10, rblapack_x, rblapack_rcond, rblapack_ferr, rblapack_berr, rblapack_info, rblapack_dlf, rblapack_df, rblapack_duf, rblapack_du2, rblapack_ipiv);
}

void
init_lapack_zgtsvx(VALUE mLapack){
  rb_define_module_function(mLapack, "zgtsvx", rblapack_zgtsvx, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}

#include "rb_lapack.h"

static VALUE
rb_cgtsvx(int argc, VALUE *argv, VALUE self){
  VALUE rb_fact;
  char fact; 
  VALUE rb_trans;
  char trans; 
  VALUE rb_dl;
  complex *dl; 
  VALUE rb_d;
  complex *d; 
  VALUE rb_du;
  complex *du; 
  VALUE rb_dlf;
  complex *dlf; 
  VALUE rb_df;
  complex *df; 
  VALUE rb_duf;
  complex *duf; 
  VALUE rb_du2;
  complex *du2; 
  VALUE rb_ipiv;
  integer *ipiv; 
  VALUE rb_b;
  complex *b; 
  VALUE rb_x;
  complex *x; 
  VALUE rb_rcond;
  real rcond; 
  VALUE rb_ferr;
  real *ferr; 
  VALUE rb_berr;
  real *berr; 
  VALUE rb_info;
  integer info; 
  VALUE rb_dlf_out__;
  complex *dlf_out__;
  VALUE rb_df_out__;
  complex *df_out__;
  VALUE rb_duf_out__;
  complex *duf_out__;
  VALUE rb_du2_out__;
  complex *du2_out__;
  VALUE rb_ipiv_out__;
  integer *ipiv_out__;
  complex *work;
  real *rwork;

  integer n;
  integer ldb;
  integer nrhs;
  integer ldx;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  x, rcond, ferr, berr, info, dlf, df, duf, du2, ipiv = NumRu::Lapack.cgtsvx( fact, trans, dl, d, du, dlf, df, duf, du2, ipiv, b)\n    or\n  NumRu::Lapack.cgtsvx  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE CGTSVX( FACT, TRANS, N, NRHS, DL, D, DU, DLF, DF, DUF, DU2, IPIV, B, LDB, X, LDX, RCOND, FERR, BERR, WORK, RWORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  CGTSVX uses the LU factorization to compute the solution to a complex\n*  system of linear equations A * X = B, A**T * X = B, or A**H * X = B,\n*  where A is a tridiagonal matrix of order N and X and B are N-by-NRHS\n*  matrices.\n*\n*  Error bounds on the solution and a condition estimate are also\n*  provided.\n*\n*  Description\n*  ===========\n*\n*  The following steps are performed:\n*\n*  1. If FACT = 'N', the LU decomposition is used to factor the matrix A\n*     as A = L * U, where L is a product of permutation and unit lower\n*     bidiagonal matrices and U is upper triangular with nonzeros in\n*     only the main diagonal and first two superdiagonals.\n*\n*  2. If some U(i,i)=0, so that U is exactly singular, then the routine\n*     returns with INFO = i. Otherwise, the factored form of A is used\n*     to estimate the condition number of the matrix A.  If the\n*     reciprocal of the condition number is less than machine precision,\n*     INFO = N+1 is returned as a warning, but the routine still goes on\n*     to solve for X and compute error bounds as described below.\n*\n*  3. The system of equations is solved for X using the factored form\n*     of A.\n*\n*  4. Iterative refinement is applied to improve the computed solution\n*     matrix and calculate error bounds and backward error estimates\n*     for it.\n*\n\n*  Arguments\n*  =========\n*\n*  FACT    (input) CHARACTER*1\n*          Specifies whether or not the factored form of A has been\n*          supplied on entry.\n*          = 'F':  DLF, DF, DUF, DU2, and IPIV contain the factored form\n*                  of A; DL, D, DU, DLF, DF, DUF, DU2 and IPIV will not\n*                  be modified.\n*          = 'N':  The matrix will be copied to DLF, DF, and DUF\n*                  and factored.\n*\n*  TRANS   (input) CHARACTER*1\n*          Specifies the form of the system of equations:\n*          = 'N':  A * X = B     (No transpose)\n*          = 'T':  A**T * X = B  (Transpose)\n*          = 'C':  A**H * X = B  (Conjugate transpose)\n*\n*  N       (input) INTEGER\n*          The order of the matrix A.  N >= 0.\n*\n*  NRHS    (input) INTEGER\n*          The number of right hand sides, i.e., the number of columns\n*          of the matrix B.  NRHS >= 0.\n*\n*  DL      (input) COMPLEX array, dimension (N-1)\n*          The (n-1) subdiagonal elements of A.\n*\n*  D       (input) COMPLEX array, dimension (N)\n*          The n diagonal elements of A.\n*\n*  DU      (input) COMPLEX array, dimension (N-1)\n*          The (n-1) superdiagonal elements of A.\n*\n*  DLF     (input or output) COMPLEX array, dimension (N-1)\n*          If FACT = 'F', then DLF is an input argument and on entry\n*          contains the (n-1) multipliers that define the matrix L from\n*          the LU factorization of A as computed by CGTTRF.\n*\n*          If FACT = 'N', then DLF is an output argument and on exit\n*          contains the (n-1) multipliers that define the matrix L from\n*          the LU factorization of A.\n*\n*  DF      (input or output) COMPLEX array, dimension (N)\n*          If FACT = 'F', then DF is an input argument and on entry\n*          contains the n diagonal elements of the upper triangular\n*          matrix U from the LU factorization of A.\n*\n*          If FACT = 'N', then DF is an output argument and on exit\n*          contains the n diagonal elements of the upper triangular\n*          matrix U from the LU factorization of A.\n*\n*  DUF     (input or output) COMPLEX array, dimension (N-1)\n*          If FACT = 'F', then DUF is an input argument and on entry\n*          contains the (n-1) elements of the first superdiagonal of U.\n*\n*          If FACT = 'N', then DUF is an output argument and on exit\n*          contains the (n-1) elements of the first superdiagonal of U.\n*\n*  DU2     (input or output) COMPLEX array, dimension (N-2)\n*          If FACT = 'F', then DU2 is an input argument and on entry\n*          contains the (n-2) elements of the second superdiagonal of\n*          U.\n*\n*          If FACT = 'N', then DU2 is an output argument and on exit\n*          contains the (n-2) elements of the second superdiagonal of\n*          U.\n*\n*  IPIV    (input or output) INTEGER array, dimension (N)\n*          If FACT = 'F', then IPIV is an input argument and on entry\n*          contains the pivot indices from the LU factorization of A as\n*          computed by CGTTRF.\n*\n*          If FACT = 'N', then IPIV is an output argument and on exit\n*          contains the pivot indices from the LU factorization of A;\n*          row i of the matrix was interchanged with row IPIV(i).\n*          IPIV(i) will always be either i or i+1; IPIV(i) = i indicates\n*          a row interchange was not required.\n*\n*  B       (input) COMPLEX array, dimension (LDB,NRHS)\n*          The N-by-NRHS right hand side matrix B.\n*\n*  LDB     (input) INTEGER\n*          The leading dimension of the array B.  LDB >= max(1,N).\n*\n*  X       (output) COMPLEX array, dimension (LDX,NRHS)\n*          If INFO = 0 or INFO = N+1, the N-by-NRHS solution matrix X.\n*\n*  LDX     (input) INTEGER\n*          The leading dimension of the array X.  LDX >= max(1,N).\n*\n*  RCOND   (output) REAL\n*          The estimate of the reciprocal condition number of the matrix\n*          A.  If RCOND is less than the machine precision (in\n*          particular, if RCOND = 0), the matrix is singular to working\n*          precision.  This condition is indicated by a return code of\n*          INFO > 0.\n*\n*  FERR    (output) REAL array, dimension (NRHS)\n*          The estimated forward error bound for each solution vector\n*          X(j) (the j-th column of the solution matrix X).\n*          If XTRUE is the true solution corresponding to X(j), FERR(j)\n*          is an estimated upper bound for the magnitude of the largest\n*          element in (X(j) - XTRUE) divided by the magnitude of the\n*          largest element in X(j).  The estimate is as reliable as\n*          the estimate for RCOND, and is almost always a slight\n*          overestimate of the true error.\n*\n*  BERR    (output) REAL array, dimension (NRHS)\n*          The componentwise relative backward error of each solution\n*          vector X(j) (i.e., the smallest relative change in\n*          any element of A or B that makes X(j) an exact solution).\n*\n*  WORK    (workspace) COMPLEX array, dimension (2*N)\n*\n*  RWORK   (workspace) REAL array, dimension (N)\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit\n*          < 0:  if INFO = -i, the i-th argument had an illegal value\n*          > 0:  if INFO = i, and i is\n*                <= N:  U(i,i) is exactly zero.  The factorization\n*                       has not been completed unless i = N, but the\n*                       factor U is exactly singular, so the solution\n*                       and error bounds could not be computed.\n*                       RCOND = 0 is returned.\n*                = N+1: U is nonsingular, but RCOND is less than machine\n*                       precision, meaning that the matrix is singular\n*                       to working precision.  Nevertheless, the\n*                       solution and error bounds are computed because\n*                       there are a number of situations where the\n*                       computed solution can be more accurate than the\n*                       value of RCOND would suggest.\n*\n\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 11)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 11)", argc);
  rb_fact = argv[0];
  rb_trans = argv[1];
  rb_dl = argv[2];
  rb_d = argv[3];
  rb_du = argv[4];
  rb_dlf = argv[5];
  rb_df = argv[6];
  rb_duf = argv[7];
  rb_du2 = argv[8];
  rb_ipiv = argv[9];
  rb_b = argv[10];

  fact = StringValueCStr(rb_fact)[0];
  trans = StringValueCStr(rb_trans)[0];
  if (!NA_IsNArray(rb_ipiv))
    rb_raise(rb_eArgError, "ipiv (3th argument) must be NArray");
  if (NA_RANK(rb_ipiv) != 1)
    rb_raise(rb_eArgError, "rank of ipiv (3th argument) must be %d", 1);
  n = NA_SHAPE0(rb_ipiv);
  if (NA_TYPE(rb_ipiv) != NA_LINT)
    rb_ipiv = na_change_type(rb_ipiv, NA_LINT);
  ipiv = NA_PTR_TYPE(rb_ipiv, integer*);
  if (!NA_IsNArray(rb_duf))
    rb_raise(rb_eArgError, "duf (4th argument) must be NArray");
  if (NA_RANK(rb_duf) != 1)
    rb_raise(rb_eArgError, "rank of duf (4th argument) must be %d", 1);
  if (NA_SHAPE0(rb_duf) != (n-1))
    rb_raise(rb_eRuntimeError, "shape 0 of duf must be %d", n-1);
  if (NA_TYPE(rb_duf) != NA_SCOMPLEX)
    rb_duf = na_change_type(rb_duf, NA_SCOMPLEX);
  duf = NA_PTR_TYPE(rb_duf, complex*);
  if (!NA_IsNArray(rb_d))
    rb_raise(rb_eArgError, "d (5th argument) must be NArray");
  if (NA_RANK(rb_d) != 1)
    rb_raise(rb_eArgError, "rank of d (5th argument) must be %d", 1);
  if (NA_SHAPE0(rb_d) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of d must be the same as shape 0 of ipiv");
  if (NA_TYPE(rb_d) != NA_SCOMPLEX)
    rb_d = na_change_type(rb_d, NA_SCOMPLEX);
  d = NA_PTR_TYPE(rb_d, complex*);
  if (!NA_IsNArray(rb_dl))
    rb_raise(rb_eArgError, "dl (6th argument) must be NArray");
  if (NA_RANK(rb_dl) != 1)
    rb_raise(rb_eArgError, "rank of dl (6th argument) must be %d", 1);
  if (NA_SHAPE0(rb_dl) != (n-1))
    rb_raise(rb_eRuntimeError, "shape 0 of dl must be %d", n-1);
  if (NA_TYPE(rb_dl) != NA_SCOMPLEX)
    rb_dl = na_change_type(rb_dl, NA_SCOMPLEX);
  dl = NA_PTR_TYPE(rb_dl, complex*);
  if (!NA_IsNArray(rb_df))
    rb_raise(rb_eArgError, "df (7th argument) must be NArray");
  if (NA_RANK(rb_df) != 1)
    rb_raise(rb_eArgError, "rank of df (7th argument) must be %d", 1);
  if (NA_SHAPE0(rb_df) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of df must be the same as shape 0 of ipiv");
  if (NA_TYPE(rb_df) != NA_SCOMPLEX)
    rb_df = na_change_type(rb_df, NA_SCOMPLEX);
  df = NA_PTR_TYPE(rb_df, complex*);
  if (!NA_IsNArray(rb_du))
    rb_raise(rb_eArgError, "du (8th argument) must be NArray");
  if (NA_RANK(rb_du) != 1)
    rb_raise(rb_eArgError, "rank of du (8th argument) must be %d", 1);
  if (NA_SHAPE0(rb_du) != (n-1))
    rb_raise(rb_eRuntimeError, "shape 0 of du must be %d", n-1);
  if (NA_TYPE(rb_du) != NA_SCOMPLEX)
    rb_du = na_change_type(rb_du, NA_SCOMPLEX);
  du = NA_PTR_TYPE(rb_du, complex*);
  if (!NA_IsNArray(rb_du2))
    rb_raise(rb_eArgError, "du2 (9th argument) must be NArray");
  if (NA_RANK(rb_du2) != 1)
    rb_raise(rb_eArgError, "rank of du2 (9th argument) must be %d", 1);
  if (NA_SHAPE0(rb_du2) != (n-2))
    rb_raise(rb_eRuntimeError, "shape 0 of du2 must be %d", n-2);
  if (NA_TYPE(rb_du2) != NA_SCOMPLEX)
    rb_du2 = na_change_type(rb_du2, NA_SCOMPLEX);
  du2 = NA_PTR_TYPE(rb_du2, complex*);
  if (!NA_IsNArray(rb_dlf))
    rb_raise(rb_eArgError, "dlf (10th argument) must be NArray");
  if (NA_RANK(rb_dlf) != 1)
    rb_raise(rb_eArgError, "rank of dlf (10th argument) must be %d", 1);
  if (NA_SHAPE0(rb_dlf) != (n-1))
    rb_raise(rb_eRuntimeError, "shape 0 of dlf must be %d", n-1);
  if (NA_TYPE(rb_dlf) != NA_SCOMPLEX)
    rb_dlf = na_change_type(rb_dlf, NA_SCOMPLEX);
  dlf = NA_PTR_TYPE(rb_dlf, complex*);
  if (!NA_IsNArray(rb_b))
    rb_raise(rb_eArgError, "b (11th argument) must be NArray");
  if (NA_RANK(rb_b) != 2)
    rb_raise(rb_eArgError, "rank of b (11th argument) must be %d", 2);
  ldb = NA_SHAPE0(rb_b);
  nrhs = NA_SHAPE1(rb_b);
  if (NA_TYPE(rb_b) != NA_SCOMPLEX)
    rb_b = na_change_type(rb_b, NA_SCOMPLEX);
  b = NA_PTR_TYPE(rb_b, complex*);
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
    rb_ferr = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  ferr = NA_PTR_TYPE(rb_ferr, real*);
  {
    int shape[1];
    shape[0] = nrhs;
    rb_berr = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  berr = NA_PTR_TYPE(rb_berr, real*);
  {
    int shape[1];
    shape[0] = n-1;
    rb_dlf_out__ = na_make_object(NA_SCOMPLEX, 1, shape, cNArray);
  }
  dlf_out__ = NA_PTR_TYPE(rb_dlf_out__, complex*);
  MEMCPY(dlf_out__, dlf, complex, NA_TOTAL(rb_dlf));
  rb_dlf = rb_dlf_out__;
  dlf = dlf_out__;
  {
    int shape[1];
    shape[0] = n;
    rb_df_out__ = na_make_object(NA_SCOMPLEX, 1, shape, cNArray);
  }
  df_out__ = NA_PTR_TYPE(rb_df_out__, complex*);
  MEMCPY(df_out__, df, complex, NA_TOTAL(rb_df));
  rb_df = rb_df_out__;
  df = df_out__;
  {
    int shape[1];
    shape[0] = n-1;
    rb_duf_out__ = na_make_object(NA_SCOMPLEX, 1, shape, cNArray);
  }
  duf_out__ = NA_PTR_TYPE(rb_duf_out__, complex*);
  MEMCPY(duf_out__, duf, complex, NA_TOTAL(rb_duf));
  rb_duf = rb_duf_out__;
  duf = duf_out__;
  {
    int shape[1];
    shape[0] = n-2;
    rb_du2_out__ = na_make_object(NA_SCOMPLEX, 1, shape, cNArray);
  }
  du2_out__ = NA_PTR_TYPE(rb_du2_out__, complex*);
  MEMCPY(du2_out__, du2, complex, NA_TOTAL(rb_du2));
  rb_du2 = rb_du2_out__;
  du2 = du2_out__;
  {
    int shape[1];
    shape[0] = DIM_LEN(n);
    rb_ipiv_out__ = na_make_object(NA_LINT, 1, shape, cNArray);
  }
  ipiv_out__ = NA_PTR_TYPE(rb_ipiv_out__, integer*);
  MEMCPY(ipiv_out__, ipiv, integer, NA_TOTAL(rb_ipiv));
  rb_ipiv = rb_ipiv_out__;
  ipiv = ipiv_out__;
  work = ALLOC_N(complex, (2*n));
  rwork = ALLOC_N(real, (n));

  cgtsvx_(&fact, &trans, &n, &nrhs, dl, d, du, dlf, df, duf, du2, ipiv, b, &ldb, x, &ldx, &rcond, ferr, berr, work, rwork, &info);

  free(work);
  free(rwork);
  rb_rcond = rb_float_new((double)rcond);
  rb_info = INT2NUM(info);
  return rb_ary_new3(10, rb_x, rb_rcond, rb_ferr, rb_berr, rb_info, rb_dlf, rb_df, rb_duf, rb_du2, rb_ipiv);
}

void
init_lapack_cgtsvx(VALUE mLapack){
  rb_define_module_function(mLapack, "cgtsvx", rb_cgtsvx, -1);
}

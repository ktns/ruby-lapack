#include "rb_lapack.h"

static VALUE
rb_dptsvx(int argc, VALUE *argv, VALUE self){
  VALUE rb_fact;
  char fact; 
  VALUE rb_d;
  doublereal *d; 
  VALUE rb_e;
  doublereal *e; 
  VALUE rb_df;
  doublereal *df; 
  VALUE rb_ef;
  doublereal *ef; 
  VALUE rb_b;
  doublereal *b; 
  VALUE rb_x;
  doublereal *x; 
  VALUE rb_rcond;
  doublereal rcond; 
  VALUE rb_ferr;
  doublereal *ferr; 
  VALUE rb_berr;
  doublereal *berr; 
  VALUE rb_info;
  integer info; 
  VALUE rb_df_out__;
  doublereal *df_out__;
  VALUE rb_ef_out__;
  doublereal *ef_out__;
  doublereal *work;

  integer n;
  integer ldb;
  integer nrhs;
  integer ldx;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  x, rcond, ferr, berr, info, df, ef = NumRu::Lapack.dptsvx( fact, d, e, df, ef, b)\n    or\n  NumRu::Lapack.dptsvx  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE DPTSVX( FACT, N, NRHS, D, E, DF, EF, B, LDB, X, LDX, RCOND, FERR, BERR, WORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  DPTSVX uses the factorization A = L*D*L**T to compute the solution\n*  to a real system of linear equations A*X = B, where A is an N-by-N\n*  symmetric positive definite tridiagonal matrix and X and B are\n*  N-by-NRHS matrices.\n*\n*  Error bounds on the solution and a condition estimate are also\n*  provided.\n*\n*  Description\n*  ===========\n*\n*  The following steps are performed:\n*\n*  1. If FACT = 'N', the matrix A is factored as A = L*D*L**T, where L\n*     is a unit lower bidiagonal matrix and D is diagonal.  The\n*     factorization can also be regarded as having the form\n*     A = U**T*D*U.\n*\n*  2. If the leading i-by-i principal minor is not positive definite,\n*     then the routine returns with INFO = i. Otherwise, the factored\n*     form of A is used to estimate the condition number of the matrix\n*     A.  If the reciprocal of the condition number is less than machine\n*     precision, INFO = N+1 is returned as a warning, but the routine\n*     still goes on to solve for X and compute error bounds as\n*     described below.\n*\n*  3. The system of equations is solved for X using the factored form\n*     of A.\n*\n*  4. Iterative refinement is applied to improve the computed solution\n*     matrix and calculate error bounds and backward error estimates\n*     for it.\n*\n\n*  Arguments\n*  =========\n*\n*  FACT    (input) CHARACTER*1\n*          Specifies whether or not the factored form of A has been\n*          supplied on entry.\n*          = 'F':  On entry, DF and EF contain the factored form of A.\n*                  D, E, DF, and EF will not be modified.\n*          = 'N':  The matrix A will be copied to DF and EF and\n*                  factored.\n*\n*  N       (input) INTEGER\n*          The order of the matrix A.  N >= 0.\n*\n*  NRHS    (input) INTEGER\n*          The number of right hand sides, i.e., the number of columns\n*          of the matrices B and X.  NRHS >= 0.\n*\n*  D       (input) DOUBLE PRECISION array, dimension (N)\n*          The n diagonal elements of the tridiagonal matrix A.\n*\n*  E       (input) DOUBLE PRECISION array, dimension (N-1)\n*          The (n-1) subdiagonal elements of the tridiagonal matrix A.\n*\n*  DF      (input or output) DOUBLE PRECISION array, dimension (N)\n*          If FACT = 'F', then DF is an input argument and on entry\n*          contains the n diagonal elements of the diagonal matrix D\n*          from the L*D*L**T factorization of A.\n*          If FACT = 'N', then DF is an output argument and on exit\n*          contains the n diagonal elements of the diagonal matrix D\n*          from the L*D*L**T factorization of A.\n*\n*  EF      (input or output) DOUBLE PRECISION array, dimension (N-1)\n*          If FACT = 'F', then EF is an input argument and on entry\n*          contains the (n-1) subdiagonal elements of the unit\n*          bidiagonal factor L from the L*D*L**T factorization of A.\n*          If FACT = 'N', then EF is an output argument and on exit\n*          contains the (n-1) subdiagonal elements of the unit\n*          bidiagonal factor L from the L*D*L**T factorization of A.\n*\n*  B       (input) DOUBLE PRECISION array, dimension (LDB,NRHS)\n*          The N-by-NRHS right hand side matrix B.\n*\n*  LDB     (input) INTEGER\n*          The leading dimension of the array B.  LDB >= max(1,N).\n*\n*  X       (output) DOUBLE PRECISION array, dimension (LDX,NRHS)\n*          If INFO = 0 of INFO = N+1, the N-by-NRHS solution matrix X.\n*\n*  LDX     (input) INTEGER\n*          The leading dimension of the array X.  LDX >= max(1,N).\n*\n*  RCOND   (output) DOUBLE PRECISION\n*          The reciprocal condition number of the matrix A.  If RCOND\n*          is less than the machine precision (in particular, if\n*          RCOND = 0), the matrix is singular to working precision.\n*          This condition is indicated by a return code of INFO > 0.\n*\n*  FERR    (output) DOUBLE PRECISION array, dimension (NRHS)\n*          The forward error bound for each solution vector\n*          X(j) (the j-th column of the solution matrix X).\n*          If XTRUE is the true solution corresponding to X(j), FERR(j)\n*          is an estimated upper bound for the magnitude of the largest\n*          element in (X(j) - XTRUE) divided by the magnitude of the\n*          largest element in X(j).\n*\n*  BERR    (output) DOUBLE PRECISION array, dimension (NRHS)\n*          The componentwise relative backward error of each solution\n*          vector X(j) (i.e., the smallest relative change in any\n*          element of A or B that makes X(j) an exact solution).\n*\n*  WORK    (workspace) DOUBLE PRECISION array, dimension (2*N)\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit\n*          < 0:  if INFO = -i, the i-th argument had an illegal value\n*          > 0:  if INFO = i, and i is\n*                <= N:  the leading minor of order i of A is\n*                       not positive definite, so the factorization\n*                       could not be completed, and the solution has not\n*                       been computed. RCOND = 0 is returned.\n*                = N+1: U is nonsingular, but RCOND is less than machine\n*                       precision, meaning that the matrix is singular\n*                       to working precision.  Nevertheless, the\n*                       solution and error bounds are computed because\n*                       there are a number of situations where the\n*                       computed solution can be more accurate than the\n*                       value of RCOND would suggest.\n*\n\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 6)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 6)", argc);
  rb_fact = argv[0];
  rb_d = argv[1];
  rb_e = argv[2];
  rb_df = argv[3];
  rb_ef = argv[4];
  rb_b = argv[5];

  fact = StringValueCStr(rb_fact)[0];
  if (!NA_IsNArray(rb_d))
    rb_raise(rb_eArgError, "d (2th argument) must be NArray");
  if (NA_RANK(rb_d) != 1)
    rb_raise(rb_eArgError, "rank of d (2th argument) must be %d", 1);
  n = NA_SHAPE0(rb_d);
  if (NA_TYPE(rb_d) != NA_DFLOAT)
    rb_d = na_change_type(rb_d, NA_DFLOAT);
  d = NA_PTR_TYPE(rb_d, doublereal*);
  if (!NA_IsNArray(rb_e))
    rb_raise(rb_eArgError, "e (3th argument) must be NArray");
  if (NA_RANK(rb_e) != 1)
    rb_raise(rb_eArgError, "rank of e (3th argument) must be %d", 1);
  if (NA_SHAPE0(rb_e) != (n-1))
    rb_raise(rb_eRuntimeError, "shape 0 of e must be %d", n-1);
  if (NA_TYPE(rb_e) != NA_DFLOAT)
    rb_e = na_change_type(rb_e, NA_DFLOAT);
  e = NA_PTR_TYPE(rb_e, doublereal*);
  if (!NA_IsNArray(rb_df))
    rb_raise(rb_eArgError, "df (4th argument) must be NArray");
  if (NA_RANK(rb_df) != 1)
    rb_raise(rb_eArgError, "rank of df (4th argument) must be %d", 1);
  if (NA_SHAPE0(rb_df) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of df must be the same as shape 0 of d");
  if (NA_TYPE(rb_df) != NA_DFLOAT)
    rb_df = na_change_type(rb_df, NA_DFLOAT);
  df = NA_PTR_TYPE(rb_df, doublereal*);
  if (!NA_IsNArray(rb_ef))
    rb_raise(rb_eArgError, "ef (5th argument) must be NArray");
  if (NA_RANK(rb_ef) != 1)
    rb_raise(rb_eArgError, "rank of ef (5th argument) must be %d", 1);
  if (NA_SHAPE0(rb_ef) != (n-1))
    rb_raise(rb_eRuntimeError, "shape 0 of ef must be %d", n-1);
  if (NA_TYPE(rb_ef) != NA_DFLOAT)
    rb_ef = na_change_type(rb_ef, NA_DFLOAT);
  ef = NA_PTR_TYPE(rb_ef, doublereal*);
  if (!NA_IsNArray(rb_b))
    rb_raise(rb_eArgError, "b (6th argument) must be NArray");
  if (NA_RANK(rb_b) != 2)
    rb_raise(rb_eArgError, "rank of b (6th argument) must be %d", 2);
  ldb = NA_SHAPE0(rb_b);
  nrhs = NA_SHAPE1(rb_b);
  if (NA_TYPE(rb_b) != NA_DFLOAT)
    rb_b = na_change_type(rb_b, NA_DFLOAT);
  b = NA_PTR_TYPE(rb_b, doublereal*);
  ldx = MAX(1,n);
  {
    int shape[2];
    shape[0] = ldx;
    shape[1] = nrhs;
    rb_x = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  x = NA_PTR_TYPE(rb_x, doublereal*);
  {
    int shape[1];
    shape[0] = nrhs;
    rb_ferr = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  ferr = NA_PTR_TYPE(rb_ferr, doublereal*);
  {
    int shape[1];
    shape[0] = nrhs;
    rb_berr = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  berr = NA_PTR_TYPE(rb_berr, doublereal*);
  {
    int shape[1];
    shape[0] = n;
    rb_df_out__ = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  df_out__ = NA_PTR_TYPE(rb_df_out__, doublereal*);
  MEMCPY(df_out__, df, doublereal, NA_TOTAL(rb_df));
  rb_df = rb_df_out__;
  df = df_out__;
  {
    int shape[1];
    shape[0] = n-1;
    rb_ef_out__ = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  ef_out__ = NA_PTR_TYPE(rb_ef_out__, doublereal*);
  MEMCPY(ef_out__, ef, doublereal, NA_TOTAL(rb_ef));
  rb_ef = rb_ef_out__;
  ef = ef_out__;
  work = ALLOC_N(doublereal, (2*n));

  dptsvx_(&fact, &n, &nrhs, d, e, df, ef, b, &ldb, x, &ldx, &rcond, ferr, berr, work, &info);

  free(work);
  rb_rcond = rb_float_new((double)rcond);
  rb_info = INT2NUM(info);
  return rb_ary_new3(7, rb_x, rb_rcond, rb_ferr, rb_berr, rb_info, rb_df, rb_ef);
}

void
init_lapack_dptsvx(VALUE mLapack){
  rb_define_module_function(mLapack, "dptsvx", rb_dptsvx, -1);
}

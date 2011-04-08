#include "rb_lapack.h"

extern VOID zptsvx_(char *fact, integer *n, integer *nrhs, doublereal *d, doublecomplex *e, doublereal *df, doublecomplex *ef, doublecomplex *b, integer *ldb, doublecomplex *x, integer *ldx, doublereal *rcond, doublereal *ferr, doublereal *berr, doublecomplex *work, doublereal *rwork, integer *info);

static VALUE sHelp, sUsage;

static VALUE
rblapack_zptsvx(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_fact;
  char fact; 
  VALUE rblapack_d;
  doublereal *d; 
  VALUE rblapack_e;
  doublecomplex *e; 
  VALUE rblapack_df;
  doublereal *df; 
  VALUE rblapack_ef;
  doublecomplex *ef; 
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
  VALUE rblapack_df_out__;
  doublereal *df_out__;
  VALUE rblapack_ef_out__;
  doublecomplex *ef_out__;
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
      printf("%s\n", "USAGE:\n  x, rcond, ferr, berr, info, df, ef = NumRu::Lapack.zptsvx( fact, d, e, df, ef, b, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE ZPTSVX( FACT, N, NRHS, D, E, DF, EF, B, LDB, X, LDX, RCOND, FERR, BERR, WORK, RWORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  ZPTSVX uses the factorization A = L*D*L**H to compute the solution\n*  to a complex system of linear equations A*X = B, where A is an\n*  N-by-N Hermitian positive definite tridiagonal matrix and X and B\n*  are N-by-NRHS matrices.\n*\n*  Error bounds on the solution and a condition estimate are also\n*  provided.\n*\n*  Description\n*  ===========\n*\n*  The following steps are performed:\n*\n*  1. If FACT = 'N', the matrix A is factored as A = L*D*L**H, where L\n*     is a unit lower bidiagonal matrix and D is diagonal.  The\n*     factorization can also be regarded as having the form\n*     A = U**H*D*U.\n*\n*  2. If the leading i-by-i principal minor is not positive definite,\n*     then the routine returns with INFO = i. Otherwise, the factored\n*     form of A is used to estimate the condition number of the matrix\n*     A.  If the reciprocal of the condition number is less than machine\n*     precision, INFO = N+1 is returned as a warning, but the routine\n*     still goes on to solve for X and compute error bounds as\n*     described below.\n*\n*  3. The system of equations is solved for X using the factored form\n*     of A.\n*\n*  4. Iterative refinement is applied to improve the computed solution\n*     matrix and calculate error bounds and backward error estimates\n*     for it.\n*\n\n*  Arguments\n*  =========\n*\n*  FACT    (input) CHARACTER*1\n*          Specifies whether or not the factored form of the matrix\n*          A is supplied on entry.\n*          = 'F':  On entry, DF and EF contain the factored form of A.\n*                  D, E, DF, and EF will not be modified.\n*          = 'N':  The matrix A will be copied to DF and EF and\n*                  factored.\n*\n*  N       (input) INTEGER\n*          The order of the matrix A.  N >= 0.\n*\n*  NRHS    (input) INTEGER\n*          The number of right hand sides, i.e., the number of columns\n*          of the matrices B and X.  NRHS >= 0.\n*\n*  D       (input) DOUBLE PRECISION array, dimension (N)\n*          The n diagonal elements of the tridiagonal matrix A.\n*\n*  E       (input) COMPLEX*16 array, dimension (N-1)\n*          The (n-1) subdiagonal elements of the tridiagonal matrix A.\n*\n*  DF      (input or output) DOUBLE PRECISION array, dimension (N)\n*          If FACT = 'F', then DF is an input argument and on entry\n*          contains the n diagonal elements of the diagonal matrix D\n*          from the L*D*L**H factorization of A.\n*          If FACT = 'N', then DF is an output argument and on exit\n*          contains the n diagonal elements of the diagonal matrix D\n*          from the L*D*L**H factorization of A.\n*\n*  EF      (input or output) COMPLEX*16 array, dimension (N-1)\n*          If FACT = 'F', then EF is an input argument and on entry\n*          contains the (n-1) subdiagonal elements of the unit\n*          bidiagonal factor L from the L*D*L**H factorization of A.\n*          If FACT = 'N', then EF is an output argument and on exit\n*          contains the (n-1) subdiagonal elements of the unit\n*          bidiagonal factor L from the L*D*L**H factorization of A.\n*\n*  B       (input) COMPLEX*16 array, dimension (LDB,NRHS)\n*          The N-by-NRHS right hand side matrix B.\n*\n*  LDB     (input) INTEGER\n*          The leading dimension of the array B.  LDB >= max(1,N).\n*\n*  X       (output) COMPLEX*16 array, dimension (LDX,NRHS)\n*          If INFO = 0 or INFO = N+1, the N-by-NRHS solution matrix X.\n*\n*  LDX     (input) INTEGER\n*          The leading dimension of the array X.  LDX >= max(1,N).\n*\n*  RCOND   (output) DOUBLE PRECISION\n*          The reciprocal condition number of the matrix A.  If RCOND\n*          is less than the machine precision (in particular, if\n*          RCOND = 0), the matrix is singular to working precision.\n*          This condition is indicated by a return code of INFO > 0.\n*\n*  FERR    (output) DOUBLE PRECISION array, dimension (NRHS)\n*          The forward error bound for each solution vector\n*          X(j) (the j-th column of the solution matrix X).\n*          If XTRUE is the true solution corresponding to X(j), FERR(j)\n*          is an estimated upper bound for the magnitude of the largest\n*          element in (X(j) - XTRUE) divided by the magnitude of the\n*          largest element in X(j).\n*\n*  BERR    (output) DOUBLE PRECISION array, dimension (NRHS)\n*          The componentwise relative backward error of each solution\n*          vector X(j) (i.e., the smallest relative change in any\n*          element of A or B that makes X(j) an exact solution).\n*\n*  WORK    (workspace) COMPLEX*16 array, dimension (N)\n*\n*  RWORK   (workspace) DOUBLE PRECISION array, dimension (N)\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit\n*          < 0:  if INFO = -i, the i-th argument had an illegal value\n*          > 0:  if INFO = i, and i is\n*                <= N:  the leading minor of order i of A is\n*                       not positive definite, so the factorization\n*                       could not be completed, and the solution has not\n*                       been computed. RCOND = 0 is returned.\n*                = N+1: U is nonsingular, but RCOND is less than machine\n*                       precision, meaning that the matrix is singular\n*                       to working precision.  Nevertheless, the\n*                       solution and error bounds are computed because\n*                       there are a number of situations where the\n*                       computed solution can be more accurate than the\n*                       value of RCOND would suggest.\n*\n\n*  =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  x, rcond, ferr, berr, info, df, ef = NumRu::Lapack.zptsvx( fact, d, e, df, ef, b, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 6)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 6)", argc);
  rblapack_fact = argv[0];
  rblapack_d = argv[1];
  rblapack_e = argv[2];
  rblapack_df = argv[3];
  rblapack_ef = argv[4];
  rblapack_b = argv[5];
  if (rb_options != Qnil) {
  }

  if (!NA_IsNArray(rblapack_df))
    rb_raise(rb_eArgError, "df (4th argument) must be NArray");
  if (NA_RANK(rblapack_df) != 1)
    rb_raise(rb_eArgError, "rank of df (4th argument) must be %d", 1);
  n = NA_SHAPE0(rblapack_df);
  if (NA_TYPE(rblapack_df) != NA_DFLOAT)
    rblapack_df = na_change_type(rblapack_df, NA_DFLOAT);
  df = NA_PTR_TYPE(rblapack_df, doublereal*);
  if (!NA_IsNArray(rblapack_d))
    rb_raise(rb_eArgError, "d (2th argument) must be NArray");
  if (NA_RANK(rblapack_d) != 1)
    rb_raise(rb_eArgError, "rank of d (2th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_d) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of d must be the same as shape 0 of df");
  if (NA_TYPE(rblapack_d) != NA_DFLOAT)
    rblapack_d = na_change_type(rblapack_d, NA_DFLOAT);
  d = NA_PTR_TYPE(rblapack_d, doublereal*);
  if (!NA_IsNArray(rblapack_b))
    rb_raise(rb_eArgError, "b (6th argument) must be NArray");
  if (NA_RANK(rblapack_b) != 2)
    rb_raise(rb_eArgError, "rank of b (6th argument) must be %d", 2);
  nrhs = NA_SHAPE1(rblapack_b);
  ldb = NA_SHAPE0(rblapack_b);
  if (NA_TYPE(rblapack_b) != NA_DCOMPLEX)
    rblapack_b = na_change_type(rblapack_b, NA_DCOMPLEX);
  b = NA_PTR_TYPE(rblapack_b, doublecomplex*);
  fact = StringValueCStr(rblapack_fact)[0];
  ldx = MAX(1,n);
  if (!NA_IsNArray(rblapack_ef))
    rb_raise(rb_eArgError, "ef (5th argument) must be NArray");
  if (NA_RANK(rblapack_ef) != 1)
    rb_raise(rb_eArgError, "rank of ef (5th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_ef) != (n-1))
    rb_raise(rb_eRuntimeError, "shape 0 of ef must be %d", n-1);
  if (NA_TYPE(rblapack_ef) != NA_DCOMPLEX)
    rblapack_ef = na_change_type(rblapack_ef, NA_DCOMPLEX);
  ef = NA_PTR_TYPE(rblapack_ef, doublecomplex*);
  if (!NA_IsNArray(rblapack_e))
    rb_raise(rb_eArgError, "e (3th argument) must be NArray");
  if (NA_RANK(rblapack_e) != 1)
    rb_raise(rb_eArgError, "rank of e (3th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_e) != (n-1))
    rb_raise(rb_eRuntimeError, "shape 0 of e must be %d", n-1);
  if (NA_TYPE(rblapack_e) != NA_DCOMPLEX)
    rblapack_e = na_change_type(rblapack_e, NA_DCOMPLEX);
  e = NA_PTR_TYPE(rblapack_e, doublecomplex*);
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
    shape[0] = n;
    rblapack_df_out__ = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  df_out__ = NA_PTR_TYPE(rblapack_df_out__, doublereal*);
  MEMCPY(df_out__, df, doublereal, NA_TOTAL(rblapack_df));
  rblapack_df = rblapack_df_out__;
  df = df_out__;
  {
    int shape[1];
    shape[0] = n-1;
    rblapack_ef_out__ = na_make_object(NA_DCOMPLEX, 1, shape, cNArray);
  }
  ef_out__ = NA_PTR_TYPE(rblapack_ef_out__, doublecomplex*);
  MEMCPY(ef_out__, ef, doublecomplex, NA_TOTAL(rblapack_ef));
  rblapack_ef = rblapack_ef_out__;
  ef = ef_out__;
  work = ALLOC_N(doublecomplex, (n));
  rwork = ALLOC_N(doublereal, (n));

  zptsvx_(&fact, &n, &nrhs, d, e, df, ef, b, &ldb, x, &ldx, &rcond, ferr, berr, work, rwork, &info);

  free(work);
  free(rwork);
  rblapack_rcond = rb_float_new((double)rcond);
  rblapack_info = INT2NUM(info);
  return rb_ary_new3(7, rblapack_x, rblapack_rcond, rblapack_ferr, rblapack_berr, rblapack_info, rblapack_df, rblapack_ef);
}

void
init_lapack_zptsvx(VALUE mLapack){
  rb_define_module_function(mLapack, "zptsvx", rblapack_zptsvx, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}

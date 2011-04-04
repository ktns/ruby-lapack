#include "rb_lapack.h"

extern VOID dtprfs_(char *uplo, char *trans, char *diag, integer *n, integer *nrhs, doublereal *ap, doublereal *b, integer *ldb, doublereal *x, integer *ldx, doublereal *ferr, doublereal *berr, doublereal *work, integer *iwork, integer *info);

static VALUE
rb_dtprfs(int argc, VALUE *argv, VALUE self){
  VALUE rb_uplo;
  char uplo; 
  VALUE rb_trans;
  char trans; 
  VALUE rb_diag;
  char diag; 
  VALUE rb_ap;
  doublereal *ap; 
  VALUE rb_b;
  doublereal *b; 
  VALUE rb_x;
  doublereal *x; 
  VALUE rb_ferr;
  doublereal *ferr; 
  VALUE rb_berr;
  doublereal *berr; 
  VALUE rb_info;
  integer info; 
  doublereal *work;
  integer *iwork;

  integer ldb;
  integer nrhs;
  integer ldx;
  integer n;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  ferr, berr, info = NumRu::Lapack.dtprfs( uplo, trans, diag, ap, b, x)\n    or\n  NumRu::Lapack.dtprfs  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE DTPRFS( UPLO, TRANS, DIAG, N, NRHS, AP, B, LDB, X, LDX, FERR, BERR, WORK, IWORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  DTPRFS provides error bounds and backward error estimates for the\n*  solution to a system of linear equations with a triangular packed\n*  coefficient matrix.\n*\n*  The solution matrix X must be computed by DTPTRS or some other\n*  means before entering this routine.  DTPRFS does not do iterative\n*  refinement because doing so cannot improve the backward error.\n*\n\n*  Arguments\n*  =========\n*\n*  UPLO    (input) CHARACTER*1\n*          = 'U':  A is upper triangular;\n*          = 'L':  A is lower triangular.\n*\n*  TRANS   (input) CHARACTER*1\n*          Specifies the form of the system of equations:\n*          = 'N':  A * X = B  (No transpose)\n*          = 'T':  A**T * X = B  (Transpose)\n*          = 'C':  A**H * X = B  (Conjugate transpose = Transpose)\n*\n*  DIAG    (input) CHARACTER*1\n*          = 'N':  A is non-unit triangular;\n*          = 'U':  A is unit triangular.\n*\n*  N       (input) INTEGER\n*          The order of the matrix A.  N >= 0.\n*\n*  NRHS    (input) INTEGER\n*          The number of right hand sides, i.e., the number of columns\n*          of the matrices B and X.  NRHS >= 0.\n*\n*  AP      (input) DOUBLE PRECISION array, dimension (N*(N+1)/2)\n*          The upper or lower triangular matrix A, packed columnwise in\n*          a linear array.  The j-th column of A is stored in the array\n*          AP as follows:\n*          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;\n*          if UPLO = 'L', AP(i + (j-1)*(2*n-j)/2) = A(i,j) for j<=i<=n.\n*          If DIAG = 'U', the diagonal elements of A are not referenced\n*          and are assumed to be 1.\n*\n*  B       (input) DOUBLE PRECISION array, dimension (LDB,NRHS)\n*          The right hand side matrix B.\n*\n*  LDB     (input) INTEGER\n*          The leading dimension of the array B.  LDB >= max(1,N).\n*\n*  X       (input) DOUBLE PRECISION array, dimension (LDX,NRHS)\n*          The solution matrix X.\n*\n*  LDX     (input) INTEGER\n*          The leading dimension of the array X.  LDX >= max(1,N).\n*\n*  FERR    (output) DOUBLE PRECISION array, dimension (NRHS)\n*          The estimated forward error bound for each solution vector\n*          X(j) (the j-th column of the solution matrix X).\n*          If XTRUE is the true solution corresponding to X(j), FERR(j)\n*          is an estimated upper bound for the magnitude of the largest\n*          element in (X(j) - XTRUE) divided by the magnitude of the\n*          largest element in X(j).  The estimate is as reliable as\n*          the estimate for RCOND, and is almost always a slight\n*          overestimate of the true error.\n*\n*  BERR    (output) DOUBLE PRECISION array, dimension (NRHS)\n*          The componentwise relative backward error of each solution\n*          vector X(j) (i.e., the smallest relative change in\n*          any element of A or B that makes X(j) an exact solution).\n*\n*  WORK    (workspace) DOUBLE PRECISION array, dimension (3*N)\n*\n*  IWORK   (workspace) INTEGER array, dimension (N)\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit\n*          < 0:  if INFO = -i, the i-th argument had an illegal value\n*\n\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 6)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 6)", argc);
  rb_uplo = argv[0];
  rb_trans = argv[1];
  rb_diag = argv[2];
  rb_ap = argv[3];
  rb_b = argv[4];
  rb_x = argv[5];

  if (!NA_IsNArray(rb_x))
    rb_raise(rb_eArgError, "x (6th argument) must be NArray");
  if (NA_RANK(rb_x) != 2)
    rb_raise(rb_eArgError, "rank of x (6th argument) must be %d", 2);
  nrhs = NA_SHAPE1(rb_x);
  ldx = NA_SHAPE0(rb_x);
  if (NA_TYPE(rb_x) != NA_DFLOAT)
    rb_x = na_change_type(rb_x, NA_DFLOAT);
  x = NA_PTR_TYPE(rb_x, doublereal*);
  if (!NA_IsNArray(rb_b))
    rb_raise(rb_eArgError, "b (5th argument) must be NArray");
  if (NA_RANK(rb_b) != 2)
    rb_raise(rb_eArgError, "rank of b (5th argument) must be %d", 2);
  if (NA_SHAPE1(rb_b) != nrhs)
    rb_raise(rb_eRuntimeError, "shape 1 of b must be the same as shape 1 of x");
  ldb = NA_SHAPE0(rb_b);
  if (NA_TYPE(rb_b) != NA_DFLOAT)
    rb_b = na_change_type(rb_b, NA_DFLOAT);
  b = NA_PTR_TYPE(rb_b, doublereal*);
  diag = StringValueCStr(rb_diag)[0];
  trans = StringValueCStr(rb_trans)[0];
  uplo = StringValueCStr(rb_uplo)[0];
  n = ldb;
  if (!NA_IsNArray(rb_ap))
    rb_raise(rb_eArgError, "ap (4th argument) must be NArray");
  if (NA_RANK(rb_ap) != 1)
    rb_raise(rb_eArgError, "rank of ap (4th argument) must be %d", 1);
  if (NA_SHAPE0(rb_ap) != (n*(n+1)/2))
    rb_raise(rb_eRuntimeError, "shape 0 of ap must be %d", n*(n+1)/2);
  if (NA_TYPE(rb_ap) != NA_DFLOAT)
    rb_ap = na_change_type(rb_ap, NA_DFLOAT);
  ap = NA_PTR_TYPE(rb_ap, doublereal*);
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
  work = ALLOC_N(doublereal, (3*n));
  iwork = ALLOC_N(integer, (n));

  dtprfs_(&uplo, &trans, &diag, &n, &nrhs, ap, b, &ldb, x, &ldx, ferr, berr, work, iwork, &info);

  free(work);
  free(iwork);
  rb_info = INT2NUM(info);
  return rb_ary_new3(3, rb_ferr, rb_berr, rb_info);
}

void
init_lapack_dtprfs(VALUE mLapack){
  rb_define_module_function(mLapack, "dtprfs", rb_dtprfs, -1);
}

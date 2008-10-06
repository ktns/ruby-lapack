#include "rb_lapack.h"

static VALUE
rb_ctrrfs(int argc, VALUE *argv, VALUE self){
  VALUE rb_uplo;
  char uplo; 
  VALUE rb_trans;
  char trans; 
  VALUE rb_diag;
  char diag; 
  VALUE rb_a;
  complex *a; 
  VALUE rb_b;
  complex *b; 
  VALUE rb_x;
  complex *x; 
  VALUE rb_ferr;
  real *ferr; 
  VALUE rb_berr;
  real *berr; 
  VALUE rb_info;
  integer info; 
  complex *work;
  real *rwork;

  integer lda;
  integer n;
  integer ldb;
  integer nrhs;
  integer ldx;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  ferr, berr, info = NumRu::Lapack.ctrrfs( uplo, trans, diag, a, b, x)\n    or\n  NumRu::Lapack.ctrrfs  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE CTRRFS( UPLO, TRANS, DIAG, N, NRHS, A, LDA, B, LDB, X, LDX, FERR, BERR, WORK, RWORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  CTRRFS provides error bounds and backward error estimates for the\n*  solution to a system of linear equations with a triangular\n*  coefficient matrix.\n*\n*  The solution matrix X must be computed by CTRTRS or some other\n*  means before entering this routine.  CTRRFS does not do iterative\n*  refinement because doing so cannot improve the backward error.\n*\n\n*  Arguments\n*  =========\n*\n*  UPLO    (input) CHARACTER*1\n*          = 'U':  A is upper triangular;\n*          = 'L':  A is lower triangular.\n*\n*  TRANS   (input) CHARACTER*1\n*          Specifies the form of the system of equations:\n*          = 'N':  A * X = B     (No transpose)\n*          = 'T':  A**T * X = B  (Transpose)\n*          = 'C':  A**H * X = B  (Conjugate transpose)\n*\n*  DIAG    (input) CHARACTER*1\n*          = 'N':  A is non-unit triangular;\n*          = 'U':  A is unit triangular.\n*\n*  N       (input) INTEGER\n*          The order of the matrix A.  N >= 0.\n*\n*  NRHS    (input) INTEGER\n*          The number of right hand sides, i.e., the number of columns\n*          of the matrices B and X.  NRHS >= 0.\n*\n*  A       (input) COMPLEX array, dimension (LDA,N)\n*          The triangular matrix A.  If UPLO = 'U', the leading N-by-N\n*          upper triangular part of the array A contains the upper\n*          triangular matrix, and the strictly lower triangular part of\n*          A is not referenced.  If UPLO = 'L', the leading N-by-N lower\n*          triangular part of the array A contains the lower triangular\n*          matrix, and the strictly upper triangular part of A is not\n*          referenced.  If DIAG = 'U', the diagonal elements of A are\n*          also not referenced and are assumed to be 1.\n*\n*  LDA     (input) INTEGER\n*          The leading dimension of the array A.  LDA >= max(1,N).\n*\n*  B       (input) COMPLEX array, dimension (LDB,NRHS)\n*          The right hand side matrix B.\n*\n*  LDB     (input) INTEGER\n*          The leading dimension of the array B.  LDB >= max(1,N).\n*\n*  X       (input) COMPLEX array, dimension (LDX,NRHS)\n*          The solution matrix X.\n*\n*  LDX     (input) INTEGER\n*          The leading dimension of the array X.  LDX >= max(1,N).\n*\n*  FERR    (output) REAL array, dimension (NRHS)\n*          The estimated forward error bound for each solution vector\n*          X(j) (the j-th column of the solution matrix X).\n*          If XTRUE is the true solution corresponding to X(j), FERR(j)\n*          is an estimated upper bound for the magnitude of the largest\n*          element in (X(j) - XTRUE) divided by the magnitude of the\n*          largest element in X(j).  The estimate is as reliable as\n*          the estimate for RCOND, and is almost always a slight\n*          overestimate of the true error.\n*\n*  BERR    (output) REAL array, dimension (NRHS)\n*          The componentwise relative backward error of each solution\n*          vector X(j) (i.e., the smallest relative change in\n*          any element of A or B that makes X(j) an exact solution).\n*\n*  WORK    (workspace) COMPLEX array, dimension (2*N)\n*\n*  RWORK   (workspace) REAL array, dimension (N)\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit\n*          < 0:  if INFO = -i, the i-th argument had an illegal value\n*\n\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 6)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 6)", argc);
  rb_uplo = argv[0];
  rb_trans = argv[1];
  rb_diag = argv[2];
  rb_a = argv[3];
  rb_b = argv[4];
  rb_x = argv[5];

  uplo = StringValueCStr(rb_uplo)[0];
  trans = StringValueCStr(rb_trans)[0];
  diag = StringValueCStr(rb_diag)[0];
  if (!NA_IsNArray(rb_a))
    rb_raise(rb_eArgError, "a (4th argument) must be NArray");
  if (NA_RANK(rb_a) != 2)
    rb_raise(rb_eArgError, "rank of a (4th argument) must be %d", 2);
  lda = NA_SHAPE0(rb_a);
  n = NA_SHAPE1(rb_a);
  if (NA_TYPE(rb_a) != NA_SCOMPLEX)
    rb_a = na_change_type(rb_a, NA_SCOMPLEX);
  a = NA_PTR_TYPE(rb_a, complex*);
  if (!NA_IsNArray(rb_b))
    rb_raise(rb_eArgError, "b (5th argument) must be NArray");
  if (NA_RANK(rb_b) != 2)
    rb_raise(rb_eArgError, "rank of b (5th argument) must be %d", 2);
  ldb = NA_SHAPE0(rb_b);
  nrhs = NA_SHAPE1(rb_b);
  if (NA_TYPE(rb_b) != NA_SCOMPLEX)
    rb_b = na_change_type(rb_b, NA_SCOMPLEX);
  b = NA_PTR_TYPE(rb_b, complex*);
  if (!NA_IsNArray(rb_x))
    rb_raise(rb_eArgError, "x (6th argument) must be NArray");
  if (NA_RANK(rb_x) != 2)
    rb_raise(rb_eArgError, "rank of x (6th argument) must be %d", 2);
  ldx = NA_SHAPE0(rb_x);
  if (NA_SHAPE1(rb_x) != nrhs)
    rb_raise(rb_eRuntimeError, "shape 1 of x must be the same as shape 1 of b");
  if (NA_TYPE(rb_x) != NA_SCOMPLEX)
    rb_x = na_change_type(rb_x, NA_SCOMPLEX);
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
  work = ALLOC_N(complex, (2*n));
  rwork = ALLOC_N(real, (n));

  ctrrfs_(&uplo, &trans, &diag, &n, &nrhs, a, &lda, b, &ldb, x, &ldx, ferr, berr, work, rwork, &info);

  free(work);
  free(rwork);
  rb_info = INT2NUM(info);
  return rb_ary_new3(3, rb_ferr, rb_berr, rb_info);
}

void
init_lapack_ctrrfs(VALUE mLapack){
  rb_define_module_function(mLapack, "ctrrfs", rb_ctrrfs, -1);
}

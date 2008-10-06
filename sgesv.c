#include "rb_lapack.h"

static VALUE
rb_sgesv(int argc, VALUE *argv, VALUE self){
  VALUE rb_a;
  real *a; 
  VALUE rb_b;
  real *b; 
  VALUE rb_ipiv;
  integer *ipiv; 
  VALUE rb_info;
  integer info; 
  VALUE rb_a_out__;
  real *a_out__;
  VALUE rb_b_out__;
  real *b_out__;

  integer lda;
  integer n;
  integer ldb;
  integer nrhs;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  ipiv, info, a, b = NumRu::Lapack.sgesv( a, b)\n    or\n  NumRu::Lapack.sgesv  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE SGESV( N, NRHS, A, LDA, IPIV, B, LDB, INFO )\n\n*  Purpose\n*  =======\n*\n*  SGESV computes the solution to a real system of linear equations\n*     A * X = B,\n*  where A is an N-by-N matrix and X and B are N-by-NRHS matrices.\n*\n*  The LU decomposition with partial pivoting and row interchanges is\n*  used to factor A as\n*     A = P * L * U,\n*  where P is a permutation matrix, L is unit lower triangular, and U is\n*  upper triangular.  The factored form of A is then used to solve the\n*  system of equations A * X = B.\n*\n\n*  Arguments\n*  =========\n*\n*  N       (input) INTEGER\n*          The number of linear equations, i.e., the order of the\n*          matrix A.  N >= 0.\n*\n*  NRHS    (input) INTEGER\n*          The number of right hand sides, i.e., the number of columns\n*          of the matrix B.  NRHS >= 0.\n*\n*  A       (input/output) REAL array, dimension (LDA,N)\n*          On entry, the N-by-N coefficient matrix A.\n*          On exit, the factors L and U from the factorization\n*          A = P*L*U; the unit diagonal elements of L are not stored.\n*\n*  LDA     (input) INTEGER\n*          The leading dimension of the array A.  LDA >= max(1,N).\n*\n*  IPIV    (output) INTEGER array, dimension (N)\n*          The pivot indices that define the permutation matrix P;\n*          row i of the matrix was interchanged with row IPIV(i).\n*\n*  B       (input/output) REAL array, dimension (LDB,NRHS)\n*          On entry, the N-by-NRHS matrix of right hand side matrix B.\n*          On exit, if INFO = 0, the N-by-NRHS solution matrix X.\n*\n*  LDB     (input) INTEGER\n*          The leading dimension of the array B.  LDB >= max(1,N).\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit\n*          < 0:  if INFO = -i, the i-th argument had an illegal value\n*          > 0:  if INFO = i, U(i,i) is exactly zero.  The factorization\n*                has been completed, but the factor U is exactly\n*                singular, so the solution could not be computed.\n*\n\n*  =====================================================================\n*\n*     .. External Subroutines ..\n      EXTERNAL           SGETRF, SGETRS, XERBLA\n*     ..\n*     .. Intrinsic Functions ..\n      INTRINSIC          MAX\n*     ..\n\n");
    return Qnil;
  }
  if (argc != 2)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 2)", argc);
  rb_a = argv[0];
  rb_b = argv[1];

  if (!NA_IsNArray(rb_a))
    rb_raise(rb_eArgError, "a (1th argument) must be NArray");
  if (NA_RANK(rb_a) != 2)
    rb_raise(rb_eArgError, "rank of a (1th argument) must be %d", 2);
  lda = NA_SHAPE0(rb_a);
  n = NA_SHAPE1(rb_a);
  if (NA_TYPE(rb_a) != NA_SFLOAT)
    rb_a = na_change_type(rb_a, NA_SFLOAT);
  a = NA_PTR_TYPE(rb_a, real*);
  if (!NA_IsNArray(rb_b))
    rb_raise(rb_eArgError, "b (2th argument) must be NArray");
  if (NA_RANK(rb_b) != 2)
    rb_raise(rb_eArgError, "rank of b (2th argument) must be %d", 2);
  ldb = NA_SHAPE0(rb_b);
  nrhs = NA_SHAPE1(rb_b);
  if (NA_TYPE(rb_b) != NA_SFLOAT)
    rb_b = na_change_type(rb_b, NA_SFLOAT);
  b = NA_PTR_TYPE(rb_b, real*);
  {
    int shape[1];
    shape[0] = n;
    rb_ipiv = na_make_object(NA_LINT, 1, shape, cNArray);
  }
  ipiv = NA_PTR_TYPE(rb_ipiv, integer*);
  {
    int shape[2];
    shape[0] = lda;
    shape[1] = n;
    rb_a_out__ = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  a_out__ = NA_PTR_TYPE(rb_a_out__, real*);
  MEMCPY(a_out__, a, real, NA_TOTAL(rb_a));
  rb_a = rb_a_out__;
  a = a_out__;
  {
    int shape[2];
    shape[0] = ldb;
    shape[1] = nrhs;
    rb_b_out__ = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  b_out__ = NA_PTR_TYPE(rb_b_out__, real*);
  MEMCPY(b_out__, b, real, NA_TOTAL(rb_b));
  rb_b = rb_b_out__;
  b = b_out__;

  sgesv_(&n, &nrhs, a, &lda, ipiv, b, &ldb, &info);

  rb_info = INT2NUM(info);
  return rb_ary_new3(4, rb_ipiv, rb_info, rb_a, rb_b);
}

void
init_lapack_sgesv(VALUE mLapack){
  rb_define_module_function(mLapack, "sgesv", rb_sgesv, -1);
}

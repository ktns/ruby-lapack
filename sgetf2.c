#include "rb_lapack.h"

static VALUE
rb_sgetf2(int argc, VALUE *argv, VALUE self){
  VALUE rb_m;
  integer m; 
  VALUE rb_a;
  real *a; 
  VALUE rb_ipiv;
  integer *ipiv; 
  VALUE rb_info;
  integer info; 
  VALUE rb_a_out__;
  real *a_out__;

  integer lda;
  integer n;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  ipiv, info, a = NumRu::Lapack.sgetf2( m, a)\n    or\n  NumRu::Lapack.sgetf2  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE SGETF2( M, N, A, LDA, IPIV, INFO )\n\n*  Purpose\n*  =======\n*\n*  SGETF2 computes an LU factorization of a general m-by-n matrix A\n*  using partial pivoting with row interchanges.\n*\n*  The factorization has the form\n*     A = P * L * U\n*  where P is a permutation matrix, L is lower triangular with unit\n*  diagonal elements (lower trapezoidal if m > n), and U is upper\n*  triangular (upper trapezoidal if m < n).\n*\n*  This is the right-looking Level 2 BLAS version of the algorithm.\n*\n\n*  Arguments\n*  =========\n*\n*  M       (input) INTEGER\n*          The number of rows of the matrix A.  M >= 0.\n*\n*  N       (input) INTEGER\n*          The number of columns of the matrix A.  N >= 0.\n*\n*  A       (input/output) REAL array, dimension (LDA,N)\n*          On entry, the m by n matrix to be factored.\n*          On exit, the factors L and U from the factorization\n*          A = P*L*U; the unit diagonal elements of L are not stored.\n*\n*  LDA     (input) INTEGER\n*          The leading dimension of the array A.  LDA >= max(1,M).\n*\n*  IPIV    (output) INTEGER array, dimension (min(M,N))\n*          The pivot indices; for 1 <= i <= min(M,N), row i of the\n*          matrix was interchanged with row IPIV(i).\n*\n*  INFO    (output) INTEGER\n*          = 0: successful exit\n*          < 0: if INFO = -k, the k-th argument had an illegal value\n*          > 0: if INFO = k, U(k,k) is exactly zero. The factorization\n*               has been completed, but the factor U is exactly\n*               singular, and division by zero will occur if it is used\n*               to solve a system of equations.\n*\n\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 2)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 2)", argc);
  rb_m = argv[0];
  rb_a = argv[1];

  m = NUM2INT(rb_m);
  if (!NA_IsNArray(rb_a))
    rb_raise(rb_eArgError, "a (2th argument) must be NArray");
  if (NA_RANK(rb_a) != 2)
    rb_raise(rb_eArgError, "rank of a (2th argument) must be %d", 2);
  lda = NA_SHAPE0(rb_a);
  n = NA_SHAPE1(rb_a);
  if (NA_TYPE(rb_a) != NA_SFLOAT)
    rb_a = na_change_type(rb_a, NA_SFLOAT);
  a = NA_PTR_TYPE(rb_a, real*);
  {
    int shape[1];
    shape[0] = DIM_LEN(MIN(m,n));
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

  sgetf2_(&m, &n, a, &lda, ipiv, &info);

  rb_info = INT2NUM(info);
  return rb_ary_new3(3, rb_ipiv, rb_info, rb_a);
}

void
init_lapack_sgetf2(VALUE mLapack){
  rb_define_module_function(mLapack, "sgetf2", rb_sgetf2, -1);
}

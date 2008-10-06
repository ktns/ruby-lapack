#include "rb_lapack.h"

static VALUE
rb_cgetc2(int argc, VALUE *argv, VALUE self){
  VALUE rb_a;
  complex *a; 
  VALUE rb_ipiv;
  integer *ipiv; 
  VALUE rb_jpiv;
  integer *jpiv; 
  VALUE rb_info;
  integer info; 
  VALUE rb_a_out__;
  complex *a_out__;

  integer lda;
  integer n;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  ipiv, jpiv, info, a = NumRu::Lapack.cgetc2( a)\n    or\n  NumRu::Lapack.cgetc2  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE CGETC2( N, A, LDA, IPIV, JPIV, INFO )\n\n*  Purpose\n*  =======\n*\n*  CGETC2 computes an LU factorization, using complete pivoting, of the\n*  n-by-n matrix A. The factorization has the form A = P * L * U * Q,\n*  where P and Q are permutation matrices, L is lower triangular with\n*  unit diagonal elements and U is upper triangular.\n*\n*  This is a level 1 BLAS version of the algorithm.\n*\n\n*  Arguments\n*  =========\n*\n*  N       (input) INTEGER\n*          The order of the matrix A. N >= 0.\n*\n*  A       (input/output) COMPLEX array, dimension (LDA, N)\n*          On entry, the n-by-n matrix to be factored.\n*          On exit, the factors L and U from the factorization\n*          A = P*L*U*Q; the unit diagonal elements of L are not stored.\n*          If U(k, k) appears to be less than SMIN, U(k, k) is given the\n*          value of SMIN, giving a nonsingular perturbed system.\n*\n*  LDA     (input) INTEGER\n*          The leading dimension of the array A.  LDA >= max(1, N).\n*\n*  IPIV    (output) INTEGER array, dimension (N).\n*          The pivot indices; for 1 <= i <= N, row i of the\n*          matrix has been interchanged with row IPIV(i).\n*\n*  JPIV    (output) INTEGER array, dimension (N).\n*          The pivot indices; for 1 <= j <= N, column j of the\n*          matrix has been interchanged with column JPIV(j).\n*\n*  INFO    (output) INTEGER\n*           = 0: successful exit\n*           > 0: if INFO = k, U(k, k) is likely to produce overflow if\n*                one tries to solve for x in Ax = b. So U is perturbed\n*                to avoid the overflow.\n*\n\n*  Further Details\n*  ===============\n*\n*  Based on contributions by\n*     Bo Kagstrom and Peter Poromaa, Department of Computing Science,\n*     Umea University, S-901 87 Umea, Sweden.\n*\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 1)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 1)", argc);
  rb_a = argv[0];

  if (!NA_IsNArray(rb_a))
    rb_raise(rb_eArgError, "a (1th argument) must be NArray");
  if (NA_RANK(rb_a) != 2)
    rb_raise(rb_eArgError, "rank of a (1th argument) must be %d", 2);
  lda = NA_SHAPE0(rb_a);
  n = NA_SHAPE1(rb_a);
  if (NA_TYPE(rb_a) != NA_SCOMPLEX)
    rb_a = na_change_type(rb_a, NA_SCOMPLEX);
  a = NA_PTR_TYPE(rb_a, complex*);
  {
    int shape[1];
    shape[0] = n;
    rb_ipiv = na_make_object(NA_LINT, 1, shape, cNArray);
  }
  ipiv = NA_PTR_TYPE(rb_ipiv, integer*);
  {
    int shape[1];
    shape[0] = n;
    rb_jpiv = na_make_object(NA_LINT, 1, shape, cNArray);
  }
  jpiv = NA_PTR_TYPE(rb_jpiv, integer*);
  {
    int shape[2];
    shape[0] = lda;
    shape[1] = n;
    rb_a_out__ = na_make_object(NA_SCOMPLEX, 2, shape, cNArray);
  }
  a_out__ = NA_PTR_TYPE(rb_a_out__, complex*);
  MEMCPY(a_out__, a, complex, NA_TOTAL(rb_a));
  rb_a = rb_a_out__;
  a = a_out__;

  cgetc2_(&n, a, &lda, ipiv, jpiv, &info);

  rb_info = INT2NUM(info);
  return rb_ary_new3(4, rb_ipiv, rb_jpiv, rb_info, rb_a);
}

void
init_lapack_cgetc2(VALUE mLapack){
  rb_define_module_function(mLapack, "cgetc2", rb_cgetc2, -1);
}

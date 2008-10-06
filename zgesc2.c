#include "rb_lapack.h"

static VALUE
rb_zgesc2(int argc, VALUE *argv, VALUE self){
  VALUE rb_a;
  doublecomplex *a; 
  VALUE rb_rhs;
  doublecomplex *rhs; 
  VALUE rb_ipiv;
  integer *ipiv; 
  VALUE rb_jpiv;
  integer *jpiv; 
  VALUE rb_scale;
  doublereal scale; 
  VALUE rb_rhs_out__;
  doublecomplex *rhs_out__;

  integer lda;
  integer n;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  scale, rhs = NumRu::Lapack.zgesc2( a, rhs, ipiv, jpiv)\n    or\n  NumRu::Lapack.zgesc2  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE ZGESC2( N, A, LDA, RHS, IPIV, JPIV, SCALE )\n\n*  Purpose\n*  =======\n*\n*  ZGESC2 solves a system of linear equations\n*\n*            A * X = scale* RHS\n*\n*  with a general N-by-N matrix A using the LU factorization with\n*  complete pivoting computed by ZGETC2.\n*\n*\n\n*  Arguments\n*  =========\n*\n*  N       (input) INTEGER\n*          The number of columns of the matrix A.\n*\n*  A       (input) COMPLEX*16 array, dimension (LDA, N)\n*          On entry, the  LU part of the factorization of the n-by-n\n*          matrix A computed by ZGETC2:  A = P * L * U * Q\n*\n*  LDA     (input) INTEGER\n*          The leading dimension of the array A.  LDA >= max(1, N).\n*\n*  RHS     (input/output) COMPLEX*16 array, dimension N.\n*          On entry, the right hand side vector b.\n*          On exit, the solution vector X.\n*\n*  IPIV    (input) INTEGER array, dimension (N).\n*          The pivot indices; for 1 <= i <= N, row i of the\n*          matrix has been interchanged with row IPIV(i).\n*\n*  JPIV    (input) INTEGER array, dimension (N).\n*          The pivot indices; for 1 <= j <= N, column j of the\n*          matrix has been interchanged with column JPIV(j).\n*\n*  SCALE    (output) DOUBLE PRECISION\n*           On exit, SCALE contains the scale factor. SCALE is chosen\n*           0 <= SCALE <= 1 to prevent owerflow in the solution.\n*\n\n*  Further Details\n*  ===============\n*\n*  Based on contributions by\n*     Bo Kagstrom and Peter Poromaa, Department of Computing Science,\n*     Umea University, S-901 87 Umea, Sweden.\n*\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 4)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 4)", argc);
  rb_a = argv[0];
  rb_rhs = argv[1];
  rb_ipiv = argv[2];
  rb_jpiv = argv[3];

  if (!NA_IsNArray(rb_a))
    rb_raise(rb_eArgError, "a (1th argument) must be NArray");
  if (NA_RANK(rb_a) != 2)
    rb_raise(rb_eArgError, "rank of a (1th argument) must be %d", 2);
  lda = NA_SHAPE0(rb_a);
  n = NA_SHAPE1(rb_a);
  if (NA_TYPE(rb_a) != NA_DCOMPLEX)
    rb_a = na_change_type(rb_a, NA_DCOMPLEX);
  a = NA_PTR_TYPE(rb_a, doublecomplex*);
  if (!NA_IsNArray(rb_rhs))
    rb_raise(rb_eArgError, "rhs (2th argument) must be NArray");
  if (NA_RANK(rb_rhs) != 1)
    rb_raise(rb_eArgError, "rank of rhs (2th argument) must be %d", 1);
  if (NA_SHAPE0(rb_rhs) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of rhs must be the same as shape 1 of a");
  if (NA_TYPE(rb_rhs) != NA_DCOMPLEX)
    rb_rhs = na_change_type(rb_rhs, NA_DCOMPLEX);
  rhs = NA_PTR_TYPE(rb_rhs, doublecomplex*);
  if (!NA_IsNArray(rb_ipiv))
    rb_raise(rb_eArgError, "ipiv (3th argument) must be NArray");
  if (NA_RANK(rb_ipiv) != 1)
    rb_raise(rb_eArgError, "rank of ipiv (3th argument) must be %d", 1);
  if (NA_SHAPE0(rb_ipiv) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of ipiv must be the same as shape 1 of a");
  if (NA_TYPE(rb_ipiv) != NA_LINT)
    rb_ipiv = na_change_type(rb_ipiv, NA_LINT);
  ipiv = NA_PTR_TYPE(rb_ipiv, integer*);
  if (!NA_IsNArray(rb_jpiv))
    rb_raise(rb_eArgError, "jpiv (4th argument) must be NArray");
  if (NA_RANK(rb_jpiv) != 1)
    rb_raise(rb_eArgError, "rank of jpiv (4th argument) must be %d", 1);
  if (NA_SHAPE0(rb_jpiv) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of jpiv must be the same as shape 1 of a");
  if (NA_TYPE(rb_jpiv) != NA_LINT)
    rb_jpiv = na_change_type(rb_jpiv, NA_LINT);
  jpiv = NA_PTR_TYPE(rb_jpiv, integer*);
  {
    int shape[1];
    shape[0] = n;
    rb_rhs_out__ = na_make_object(NA_DCOMPLEX, 1, shape, cNArray);
  }
  rhs_out__ = NA_PTR_TYPE(rb_rhs_out__, doublecomplex*);
  MEMCPY(rhs_out__, rhs, doublecomplex, NA_TOTAL(rb_rhs));
  rb_rhs = rb_rhs_out__;
  rhs = rhs_out__;

  zgesc2_(&n, a, &lda, rhs, ipiv, jpiv, &scale);

  rb_scale = rb_float_new((double)scale);
  return rb_ary_new3(2, rb_scale, rb_rhs);
}

void
init_lapack_zgesc2(VALUE mLapack){
  rb_define_module_function(mLapack, "zgesc2", rb_zgesc2, -1);
}

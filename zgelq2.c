#include "rb_lapack.h"

static VALUE
rb_zgelq2(int argc, VALUE *argv, VALUE self){
  VALUE rb_a;
  doublecomplex *a; 
  VALUE rb_tau;
  doublecomplex *tau; 
  VALUE rb_info;
  integer info; 
  VALUE rb_a_out__;
  doublecomplex *a_out__;
  doublecomplex *work;

  integer lda;
  integer n;
  integer m;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  tau, info, a = NumRu::Lapack.zgelq2( a)\n    or\n  NumRu::Lapack.zgelq2  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE ZGELQ2( M, N, A, LDA, TAU, WORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  ZGELQ2 computes an LQ factorization of a complex m by n matrix A:\n*  A = L * Q.\n*\n\n*  Arguments\n*  =========\n*\n*  M       (input) INTEGER\n*          The number of rows of the matrix A.  M >= 0.\n*\n*  N       (input) INTEGER\n*          The number of columns of the matrix A.  N >= 0.\n*\n*  A       (input/output) COMPLEX*16 array, dimension (LDA,N)\n*          On entry, the m by n matrix A.\n*          On exit, the elements on and below the diagonal of the array\n*          contain the m by min(m,n) lower trapezoidal matrix L (L is\n*          lower triangular if m <= n); the elements above the diagonal,\n*          with the array TAU, represent the unitary matrix Q as a\n*          product of elementary reflectors (see Further Details).\n*\n*  LDA     (input) INTEGER\n*          The leading dimension of the array A.  LDA >= max(1,M).\n*\n*  TAU     (output) COMPLEX*16 array, dimension (min(M,N))\n*          The scalar factors of the elementary reflectors (see Further\n*          Details).\n*\n*  WORK    (workspace) COMPLEX*16 array, dimension (M)\n*\n*  INFO    (output) INTEGER\n*          = 0: successful exit\n*          < 0: if INFO = -i, the i-th argument had an illegal value\n*\n\n*  Further Details\n*  ===============\n*\n*  The matrix Q is represented as a product of elementary reflectors\n*\n*     Q = H(k)' . . . H(2)' H(1)', where k = min(m,n).\n*\n*  Each H(i) has the form\n*\n*     H(i) = I - tau * v * v'\n*\n*  where tau is a complex scalar, and v is a complex vector with\n*  v(1:i-1) = 0 and v(i) = 1; conjg(v(i+1:n)) is stored on exit in\n*  A(i,i+1:n), and tau in TAU(i).\n*\n*  =====================================================================\n*\n\n");
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
  if (NA_TYPE(rb_a) != NA_DCOMPLEX)
    rb_a = na_change_type(rb_a, NA_DCOMPLEX);
  a = NA_PTR_TYPE(rb_a, doublecomplex*);
  m = lda;
  {
    int shape[1];
    shape[0] = MIN(m,n);
    rb_tau = na_make_object(NA_DCOMPLEX, 1, shape, cNArray);
  }
  tau = NA_PTR_TYPE(rb_tau, doublecomplex*);
  {
    int shape[2];
    shape[0] = lda;
    shape[1] = n;
    rb_a_out__ = na_make_object(NA_DCOMPLEX, 2, shape, cNArray);
  }
  a_out__ = NA_PTR_TYPE(rb_a_out__, doublecomplex*);
  MEMCPY(a_out__, a, doublecomplex, NA_TOTAL(rb_a));
  rb_a = rb_a_out__;
  a = a_out__;
  work = ALLOC_N(doublecomplex, (m));

  zgelq2_(&m, &n, a, &lda, tau, work, &info);

  free(work);
  rb_info = INT2NUM(info);
  return rb_ary_new3(3, rb_tau, rb_info, rb_a);
}

void
init_lapack_zgelq2(VALUE mLapack){
  rb_define_module_function(mLapack, "zgelq2", rb_zgelq2, -1);
}

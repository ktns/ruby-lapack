#include "rb_lapack.h"

extern VOID sgeqr2_(integer *m, integer *n, real *a, integer *lda, real *tau, real *work, integer *info);

static VALUE
rb_sgeqr2(int argc, VALUE *argv, VALUE self){
  VALUE rb_m;
  integer m; 
  VALUE rb_a;
  real *a; 
  VALUE rb_tau;
  real *tau; 
  VALUE rb_info;
  integer info; 
  VALUE rb_a_out__;
  real *a_out__;
  real *work;

  integer lda;
  integer n;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  tau, info, a = NumRu::Lapack.sgeqr2( m, a)\n    or\n  NumRu::Lapack.sgeqr2  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE SGEQR2( M, N, A, LDA, TAU, WORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  SGEQR2 computes a QR factorization of a real m by n matrix A:\n*  A = Q * R.\n*\n\n*  Arguments\n*  =========\n*\n*  M       (input) INTEGER\n*          The number of rows of the matrix A.  M >= 0.\n*\n*  N       (input) INTEGER\n*          The number of columns of the matrix A.  N >= 0.\n*\n*  A       (input/output) REAL array, dimension (LDA,N)\n*          On entry, the m by n matrix A.\n*          On exit, the elements on and above the diagonal of the array\n*          contain the min(m,n) by n upper trapezoidal matrix R (R is\n*          upper triangular if m >= n); the elements below the diagonal,\n*          with the array TAU, represent the orthogonal matrix Q as a\n*          product of elementary reflectors (see Further Details).\n*\n*  LDA     (input) INTEGER\n*          The leading dimension of the array A.  LDA >= max(1,M).\n*\n*  TAU     (output) REAL array, dimension (min(M,N))\n*          The scalar factors of the elementary reflectors (see Further\n*          Details).\n*\n*  WORK    (workspace) REAL array, dimension (N)\n*\n*  INFO    (output) INTEGER\n*          = 0: successful exit\n*          < 0: if INFO = -i, the i-th argument had an illegal value\n*\n\n*  Further Details\n*  ===============\n*\n*  The matrix Q is represented as a product of elementary reflectors\n*\n*     Q = H(1) H(2) . . . H(k), where k = min(m,n).\n*\n*  Each H(i) has the form\n*\n*     H(i) = I - tau * v * v'\n*\n*  where tau is a real scalar, and v is a real vector with\n*  v(1:i-1) = 0 and v(i) = 1; v(i+1:m) is stored on exit in A(i+1:m,i),\n*  and tau in TAU(i).\n*\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 2)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 2)", argc);
  rb_m = argv[0];
  rb_a = argv[1];

  if (!NA_IsNArray(rb_a))
    rb_raise(rb_eArgError, "a (2th argument) must be NArray");
  if (NA_RANK(rb_a) != 2)
    rb_raise(rb_eArgError, "rank of a (2th argument) must be %d", 2);
  n = NA_SHAPE1(rb_a);
  lda = NA_SHAPE0(rb_a);
  if (NA_TYPE(rb_a) != NA_SFLOAT)
    rb_a = na_change_type(rb_a, NA_SFLOAT);
  a = NA_PTR_TYPE(rb_a, real*);
  m = NUM2INT(rb_m);
  {
    int shape[1];
    shape[0] = MIN(m,n);
    rb_tau = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  tau = NA_PTR_TYPE(rb_tau, real*);
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
  work = ALLOC_N(real, (n));

  sgeqr2_(&m, &n, a, &lda, tau, work, &info);

  free(work);
  rb_info = INT2NUM(info);
  return rb_ary_new3(3, rb_tau, rb_info, rb_a);
}

void
init_lapack_sgeqr2(VALUE mLapack){
  rb_define_module_function(mLapack, "sgeqr2", rb_sgeqr2, -1);
}

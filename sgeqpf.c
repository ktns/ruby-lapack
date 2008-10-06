#include "rb_lapack.h"

static VALUE
rb_sgeqpf(int argc, VALUE *argv, VALUE self){
  VALUE rb_m;
  integer m; 
  VALUE rb_a;
  real *a; 
  VALUE rb_jpvt;
  integer *jpvt; 
  VALUE rb_tau;
  real *tau; 
  VALUE rb_info;
  integer info; 
  VALUE rb_a_out__;
  real *a_out__;
  VALUE rb_jpvt_out__;
  integer *jpvt_out__;
  real *work;

  integer lda;
  integer n;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  tau, info, a, jpvt = NumRu::Lapack.sgeqpf( m, a, jpvt)\n    or\n  NumRu::Lapack.sgeqpf  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE SGEQPF( M, N, A, LDA, JPVT, TAU, WORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  This routine is deprecated and has been replaced by routine SGEQP3.\n*\n*  SGEQPF computes a QR factorization with column pivoting of a\n*  real M-by-N matrix A: A*P = Q*R.\n*\n\n*  Arguments\n*  =========\n*\n*  M       (input) INTEGER\n*          The number of rows of the matrix A. M >= 0.\n*\n*  N       (input) INTEGER\n*          The number of columns of the matrix A. N >= 0\n*\n*  A       (input/output) REAL array, dimension (LDA,N)\n*          On entry, the M-by-N matrix A.\n*          On exit, the upper triangle of the array contains the\n*          min(M,N)-by-N upper triangular matrix R; the elements\n*          below the diagonal, together with the array TAU,\n*          represent the orthogonal matrix Q as a product of\n*          min(m,n) elementary reflectors.\n*\n*  LDA     (input) INTEGER\n*          The leading dimension of the array A. LDA >= max(1,M).\n*\n*  JPVT    (input/output) INTEGER array, dimension (N)\n*          On entry, if JPVT(i) .ne. 0, the i-th column of A is permuted\n*          to the front of A*P (a leading column); if JPVT(i) = 0,\n*          the i-th column of A is a free column.\n*          On exit, if JPVT(i) = k, then the i-th column of A*P\n*          was the k-th column of A.\n*\n*  TAU     (output) REAL array, dimension (min(M,N))\n*          The scalar factors of the elementary reflectors.\n*\n*  WORK    (workspace) REAL array, dimension (3*N)\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit\n*          < 0:  if INFO = -i, the i-th argument had an illegal value\n*\n\n*  Further Details\n*  ===============\n*\n*  The matrix Q is represented as a product of elementary reflectors\n*\n*     Q = H(1) H(2) . . . H(n)\n*\n*  Each H(i) has the form\n*\n*     H = I - tau * v * v'\n*\n*  where tau is a real scalar, and v is a real vector with\n*  v(1:i-1) = 0 and v(i) = 1; v(i+1:m) is stored on exit in A(i+1:m,i).\n*\n*  The matrix P is represented in jpvt as follows: If\n*     jpvt(j) = i\n*  then the jth column of P is the ith canonical unit vector.\n*\n*  Partial column norm updating strategy modified by\n*    Z. Drmac and Z. Bujanovic, Dept. of Mathematics,\n*    University of Zagreb, Croatia.\n*    June 2006.\n*  For more details see LAPACK Working Note 176.\n*\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 3)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 3)", argc);
  rb_m = argv[0];
  rb_a = argv[1];
  rb_jpvt = argv[2];

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
  if (!NA_IsNArray(rb_jpvt))
    rb_raise(rb_eArgError, "jpvt (3th argument) must be NArray");
  if (NA_RANK(rb_jpvt) != 1)
    rb_raise(rb_eArgError, "rank of jpvt (3th argument) must be %d", 1);
  if (NA_SHAPE0(rb_jpvt) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of jpvt must be the same as shape 1 of a");
  if (NA_TYPE(rb_jpvt) != NA_LINT)
    rb_jpvt = na_change_type(rb_jpvt, NA_LINT);
  jpvt = NA_PTR_TYPE(rb_jpvt, integer*);
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
  {
    int shape[1];
    shape[0] = n;
    rb_jpvt_out__ = na_make_object(NA_LINT, 1, shape, cNArray);
  }
  jpvt_out__ = NA_PTR_TYPE(rb_jpvt_out__, integer*);
  MEMCPY(jpvt_out__, jpvt, integer, NA_TOTAL(rb_jpvt));
  rb_jpvt = rb_jpvt_out__;
  jpvt = jpvt_out__;
  work = ALLOC_N(real, (3*n));

  sgeqpf_(&m, &n, a, &lda, jpvt, tau, work, &info);

  free(work);
  rb_info = INT2NUM(info);
  return rb_ary_new3(4, rb_tau, rb_info, rb_a, rb_jpvt);
}

void
init_lapack_sgeqpf(VALUE mLapack){
  rb_define_module_function(mLapack, "sgeqpf", rb_sgeqpf, -1);
}

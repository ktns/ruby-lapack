#include "rb_lapack.h"

extern VOID zgelq2_(integer *m, integer *n, doublecomplex *a, integer *lda, doublecomplex *tau, doublecomplex *work, integer *info);

static VALUE sHelp, sUsage;

static VALUE
rblapack_zgelq2(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_a;
  doublecomplex *a; 
  VALUE rblapack_tau;
  doublecomplex *tau; 
  VALUE rblapack_info;
  integer info; 
  VALUE rblapack_a_out__;
  doublecomplex *a_out__;
  doublecomplex *work;

  integer lda;
  integer n;
  integer m;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  tau, info, a = NumRu::Lapack.zgelq2( a, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE ZGELQ2( M, N, A, LDA, TAU, WORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  ZGELQ2 computes an LQ factorization of a complex m by n matrix A:\n*  A = L * Q.\n*\n\n*  Arguments\n*  =========\n*\n*  M       (input) INTEGER\n*          The number of rows of the matrix A.  M >= 0.\n*\n*  N       (input) INTEGER\n*          The number of columns of the matrix A.  N >= 0.\n*\n*  A       (input/output) COMPLEX*16 array, dimension (LDA,N)\n*          On entry, the m by n matrix A.\n*          On exit, the elements on and below the diagonal of the array\n*          contain the m by min(m,n) lower trapezoidal matrix L (L is\n*          lower triangular if m <= n); the elements above the diagonal,\n*          with the array TAU, represent the unitary matrix Q as a\n*          product of elementary reflectors (see Further Details).\n*\n*  LDA     (input) INTEGER\n*          The leading dimension of the array A.  LDA >= max(1,M).\n*\n*  TAU     (output) COMPLEX*16 array, dimension (min(M,N))\n*          The scalar factors of the elementary reflectors (see Further\n*          Details).\n*\n*  WORK    (workspace) COMPLEX*16 array, dimension (M)\n*\n*  INFO    (output) INTEGER\n*          = 0: successful exit\n*          < 0: if INFO = -i, the i-th argument had an illegal value\n*\n\n*  Further Details\n*  ===============\n*\n*  The matrix Q is represented as a product of elementary reflectors\n*\n*     Q = H(k)' . . . H(2)' H(1)', where k = min(m,n).\n*\n*  Each H(i) has the form\n*\n*     H(i) = I - tau * v * v'\n*\n*  where tau is a complex scalar, and v is a complex vector with\n*  v(1:i-1) = 0 and v(i) = 1; conjg(v(i+1:n)) is stored on exit in\n*  A(i,i+1:n), and tau in TAU(i).\n*\n*  =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  tau, info, a = NumRu::Lapack.zgelq2( a, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 1)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 1)", argc);
  rblapack_a = argv[0];
  if (rb_options != Qnil) {
  }

  if (!NA_IsNArray(rblapack_a))
    rb_raise(rb_eArgError, "a (1th argument) must be NArray");
  if (NA_RANK(rblapack_a) != 2)
    rb_raise(rb_eArgError, "rank of a (1th argument) must be %d", 2);
  n = NA_SHAPE1(rblapack_a);
  lda = NA_SHAPE0(rblapack_a);
  if (NA_TYPE(rblapack_a) != NA_DCOMPLEX)
    rblapack_a = na_change_type(rblapack_a, NA_DCOMPLEX);
  a = NA_PTR_TYPE(rblapack_a, doublecomplex*);
  m = lda;
  {
    int shape[1];
    shape[0] = MIN(m,n);
    rblapack_tau = na_make_object(NA_DCOMPLEX, 1, shape, cNArray);
  }
  tau = NA_PTR_TYPE(rblapack_tau, doublecomplex*);
  {
    int shape[2];
    shape[0] = lda;
    shape[1] = n;
    rblapack_a_out__ = na_make_object(NA_DCOMPLEX, 2, shape, cNArray);
  }
  a_out__ = NA_PTR_TYPE(rblapack_a_out__, doublecomplex*);
  MEMCPY(a_out__, a, doublecomplex, NA_TOTAL(rblapack_a));
  rblapack_a = rblapack_a_out__;
  a = a_out__;
  work = ALLOC_N(doublecomplex, (m));

  zgelq2_(&m, &n, a, &lda, tau, work, &info);

  free(work);
  rblapack_info = INT2NUM(info);
  return rb_ary_new3(3, rblapack_tau, rblapack_info, rblapack_a);
}

void
init_lapack_zgelq2(VALUE mLapack){
  rb_define_module_function(mLapack, "zgelq2", rblapack_zgelq2, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}

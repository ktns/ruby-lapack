#include "rb_lapack.h"

extern VOID sgeqlf_(integer *m, integer *n, real *a, integer *lda, real *tau, real *work, integer *lwork, integer *info);

static VALUE sHelp, sUsage;

static VALUE
rblapack_sgeqlf(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_m;
  integer m; 
  VALUE rblapack_a;
  real *a; 
  VALUE rblapack_lwork;
  integer lwork; 
  VALUE rblapack_tau;
  real *tau; 
  VALUE rblapack_work;
  real *work; 
  VALUE rblapack_info;
  integer info; 
  VALUE rblapack_a_out__;
  real *a_out__;

  integer lda;
  integer n;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  tau, work, info, a = NumRu::Lapack.sgeqlf( m, a, lwork, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE SGEQLF( M, N, A, LDA, TAU, WORK, LWORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  SGEQLF computes a QL factorization of a real M-by-N matrix A:\n*  A = Q * L.\n*\n\n*  Arguments\n*  =========\n*\n*  M       (input) INTEGER\n*          The number of rows of the matrix A.  M >= 0.\n*\n*  N       (input) INTEGER\n*          The number of columns of the matrix A.  N >= 0.\n*\n*  A       (input/output) REAL array, dimension (LDA,N)\n*          On entry, the M-by-N matrix A.\n*          On exit,\n*          if m >= n, the lower triangle of the subarray\n*          A(m-n+1:m,1:n) contains the N-by-N lower triangular matrix L;\n*          if m <= n, the elements on and below the (n-m)-th\n*          superdiagonal contain the M-by-N lower trapezoidal matrix L;\n*          the remaining elements, with the array TAU, represent the\n*          orthogonal matrix Q as a product of elementary reflectors\n*          (see Further Details).\n*\n*  LDA     (input) INTEGER\n*          The leading dimension of the array A.  LDA >= max(1,M).\n*\n*  TAU     (output) REAL array, dimension (min(M,N))\n*          The scalar factors of the elementary reflectors (see Further\n*          Details).\n*\n*  WORK    (workspace/output) REAL array, dimension (MAX(1,LWORK))\n*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.\n*\n*  LWORK   (input) INTEGER\n*          The dimension of the array WORK.  LWORK >= max(1,N).\n*          For optimum performance LWORK >= N*NB, where NB is the\n*          optimal blocksize.\n*\n*          If LWORK = -1, then a workspace query is assumed; the routine\n*          only calculates the optimal size of the WORK array, returns\n*          this value as the first entry of the WORK array, and no error\n*          message related to LWORK is issued by XERBLA.\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit\n*          < 0:  if INFO = -i, the i-th argument had an illegal value\n*\n\n*  Further Details\n*  ===============\n*\n*  The matrix Q is represented as a product of elementary reflectors\n*\n*     Q = H(k) . . . H(2) H(1), where k = min(m,n).\n*\n*  Each H(i) has the form\n*\n*     H(i) = I - tau * v * v'\n*\n*  where tau is a real scalar, and v is a real vector with\n*  v(m-k+i+1:m) = 0 and v(m-k+i) = 1; v(1:m-k+i-1) is stored on exit in\n*  A(1:m-k+i-1,n-k+i), and tau in TAU(i).\n*\n*  =====================================================================\n*\n*     .. Local Scalars ..\n      LOGICAL            LQUERY\n      INTEGER            I, IB, IINFO, IWS, K, KI, KK, LDWORK, LWKOPT,\n     $                   MU, NB, NBMIN, NU, NX\n*     ..\n*     .. External Subroutines ..\n      EXTERNAL           SGEQL2, SLARFB, SLARFT, XERBLA\n*     ..\n*     .. Intrinsic Functions ..\n      INTRINSIC          MAX, MIN\n*     ..\n*     .. External Functions ..\n      INTEGER            ILAENV\n      EXTERNAL           ILAENV\n*     ..\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  tau, work, info, a = NumRu::Lapack.sgeqlf( m, a, lwork, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 3)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 3)", argc);
  rblapack_m = argv[0];
  rblapack_a = argv[1];
  rblapack_lwork = argv[2];
  if (rb_options != Qnil) {
  }

  if (!NA_IsNArray(rblapack_a))
    rb_raise(rb_eArgError, "a (2th argument) must be NArray");
  if (NA_RANK(rblapack_a) != 2)
    rb_raise(rb_eArgError, "rank of a (2th argument) must be %d", 2);
  n = NA_SHAPE1(rblapack_a);
  lda = NA_SHAPE0(rblapack_a);
  if (NA_TYPE(rblapack_a) != NA_SFLOAT)
    rblapack_a = na_change_type(rblapack_a, NA_SFLOAT);
  a = NA_PTR_TYPE(rblapack_a, real*);
  m = NUM2INT(rblapack_m);
  lwork = NUM2INT(rblapack_lwork);
  {
    int shape[1];
    shape[0] = MIN(m,n);
    rblapack_tau = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  tau = NA_PTR_TYPE(rblapack_tau, real*);
  {
    int shape[1];
    shape[0] = MAX(1,lwork);
    rblapack_work = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  work = NA_PTR_TYPE(rblapack_work, real*);
  {
    int shape[2];
    shape[0] = lda;
    shape[1] = n;
    rblapack_a_out__ = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  a_out__ = NA_PTR_TYPE(rblapack_a_out__, real*);
  MEMCPY(a_out__, a, real, NA_TOTAL(rblapack_a));
  rblapack_a = rblapack_a_out__;
  a = a_out__;

  sgeqlf_(&m, &n, a, &lda, tau, work, &lwork, &info);

  rblapack_info = INT2NUM(info);
  return rb_ary_new3(4, rblapack_tau, rblapack_work, rblapack_info, rblapack_a);
}

void
init_lapack_sgeqlf(VALUE mLapack){
  rb_define_module_function(mLapack, "sgeqlf", rblapack_sgeqlf, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}

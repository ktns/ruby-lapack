#include "rb_lapack.h"

extern VOID zgerqf_(integer *m, integer *n, doublecomplex *a, integer *lda, doublecomplex *tau, doublecomplex *work, integer *lwork, integer *info);

static VALUE sHelp, sUsage;

static VALUE
rblapack_zgerqf(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_m;
  integer m; 
  VALUE rblapack_a;
  doublecomplex *a; 
  VALUE rblapack_lwork;
  integer lwork; 
  VALUE rblapack_tau;
  doublecomplex *tau; 
  VALUE rblapack_work;
  doublecomplex *work; 
  VALUE rblapack_info;
  integer info; 
  VALUE rblapack_a_out__;
  doublecomplex *a_out__;

  integer lda;
  integer n;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  tau, work, info, a = NumRu::Lapack.zgerqf( m, a, lwork, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE ZGERQF( M, N, A, LDA, TAU, WORK, LWORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  ZGERQF computes an RQ factorization of a complex M-by-N matrix A:\n*  A = R * Q.\n*\n\n*  Arguments\n*  =========\n*\n*  M       (input) INTEGER\n*          The number of rows of the matrix A.  M >= 0.\n*\n*  N       (input) INTEGER\n*          The number of columns of the matrix A.  N >= 0.\n*\n*  A       (input/output) COMPLEX*16 array, dimension (LDA,N)\n*          On entry, the M-by-N matrix A.\n*          On exit,\n*          if m <= n, the upper triangle of the subarray\n*          A(1:m,n-m+1:n) contains the M-by-M upper triangular matrix R;\n*          if m >= n, the elements on and above the (m-n)-th subdiagonal\n*          contain the M-by-N upper trapezoidal matrix R;\n*          the remaining elements, with the array TAU, represent the\n*          unitary matrix Q as a product of min(m,n) elementary\n*          reflectors (see Further Details).\n*\n*  LDA     (input) INTEGER\n*          The leading dimension of the array A.  LDA >= max(1,M).\n*\n*  TAU     (output) COMPLEX*16 array, dimension (min(M,N))\n*          The scalar factors of the elementary reflectors (see Further\n*          Details).\n*\n*  WORK    (workspace/output) COMPLEX*16 array, dimension (MAX(1,LWORK))\n*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.\n*\n*  LWORK   (input) INTEGER\n*          The dimension of the array WORK.  LWORK >= max(1,M).\n*          For optimum performance LWORK >= M*NB, where NB is\n*          the optimal blocksize.\n*\n*          If LWORK = -1, then a workspace query is assumed; the routine\n*          only calculates the optimal size of the WORK array, returns\n*          this value as the first entry of the WORK array, and no error\n*          message related to LWORK is issued by XERBLA.\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit\n*          < 0:  if INFO = -i, the i-th argument had an illegal value\n*\n\n*  Further Details\n*  ===============\n*\n*  The matrix Q is represented as a product of elementary reflectors\n*\n*     Q = H(1)' H(2)' . . . H(k)', where k = min(m,n).\n*\n*  Each H(i) has the form\n*\n*     H(i) = I - tau * v * v'\n*\n*  where tau is a complex scalar, and v is a complex vector with\n*  v(n-k+i+1:n) = 0 and v(n-k+i) = 1; conjg(v(1:n-k+i-1)) is stored on\n*  exit in A(m-k+i,1:n-k+i-1), and tau in TAU(i).\n*\n*  =====================================================================\n*\n*     .. Local Scalars ..\n      LOGICAL            LQUERY\n      INTEGER            I, IB, IINFO, IWS, K, KI, KK, LDWORK, LWKOPT,\n     $                   MU, NB, NBMIN, NU, NX\n*     ..\n*     .. External Subroutines ..\n      EXTERNAL           XERBLA, ZGERQ2, ZLARFB, ZLARFT\n*     ..\n*     .. Intrinsic Functions ..\n      INTRINSIC          MAX, MIN\n*     ..\n*     .. External Functions ..\n      INTEGER            ILAENV\n      EXTERNAL           ILAENV\n*     ..\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  tau, work, info, a = NumRu::Lapack.zgerqf( m, a, lwork, [:usage => usage, :help => help])\n");
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
  if (NA_TYPE(rblapack_a) != NA_DCOMPLEX)
    rblapack_a = na_change_type(rblapack_a, NA_DCOMPLEX);
  a = NA_PTR_TYPE(rblapack_a, doublecomplex*);
  m = NUM2INT(rblapack_m);
  lwork = NUM2INT(rblapack_lwork);
  {
    int shape[1];
    shape[0] = MIN(m,n);
    rblapack_tau = na_make_object(NA_DCOMPLEX, 1, shape, cNArray);
  }
  tau = NA_PTR_TYPE(rblapack_tau, doublecomplex*);
  {
    int shape[1];
    shape[0] = MAX(1,lwork);
    rblapack_work = na_make_object(NA_DCOMPLEX, 1, shape, cNArray);
  }
  work = NA_PTR_TYPE(rblapack_work, doublecomplex*);
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

  zgerqf_(&m, &n, a, &lda, tau, work, &lwork, &info);

  rblapack_info = INT2NUM(info);
  return rb_ary_new3(4, rblapack_tau, rblapack_work, rblapack_info, rblapack_a);
}

void
init_lapack_zgerqf(VALUE mLapack){
  rb_define_module_function(mLapack, "zgerqf", rblapack_zgerqf, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}

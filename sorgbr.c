#include "rb_lapack.h"

extern VOID sorgbr_(char *vect, integer *m, integer *n, integer *k, real *a, integer *lda, real *tau, real *work, integer *lwork, integer *info);

static VALUE
rb_sorgbr(int argc, VALUE *argv, VALUE self){
  VALUE rb_vect;
  char vect; 
  VALUE rb_m;
  integer m; 
  VALUE rb_k;
  integer k; 
  VALUE rb_a;
  real *a; 
  VALUE rb_tau;
  real *tau; 
  VALUE rb_lwork;
  integer lwork; 
  VALUE rb_work;
  real *work; 
  VALUE rb_info;
  integer info; 
  VALUE rb_a_out__;
  real *a_out__;

  integer lda;
  integer n;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  work, info, a = NumRu::Lapack.sorgbr( vect, m, k, a, tau, lwork)\n    or\n  NumRu::Lapack.sorgbr  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE SORGBR( VECT, M, N, K, A, LDA, TAU, WORK, LWORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  SORGBR generates one of the real orthogonal matrices Q or P**T\n*  determined by SGEBRD when reducing a real matrix A to bidiagonal\n*  form: A = Q * B * P**T.  Q and P**T are defined as products of\n*  elementary reflectors H(i) or G(i) respectively.\n*\n*  If VECT = 'Q', A is assumed to have been an M-by-K matrix, and Q\n*  is of order M:\n*  if m >= k, Q = H(1) H(2) . . . H(k) and SORGBR returns the first n\n*  columns of Q, where m >= n >= k;\n*  if m < k, Q = H(1) H(2) . . . H(m-1) and SORGBR returns Q as an\n*  M-by-M matrix.\n*\n*  If VECT = 'P', A is assumed to have been a K-by-N matrix, and P**T\n*  is of order N:\n*  if k < n, P**T = G(k) . . . G(2) G(1) and SORGBR returns the first m\n*  rows of P**T, where n >= m >= k;\n*  if k >= n, P**T = G(n-1) . . . G(2) G(1) and SORGBR returns P**T as\n*  an N-by-N matrix.\n*\n\n*  Arguments\n*  =========\n*\n*  VECT    (input) CHARACTER*1\n*          Specifies whether the matrix Q or the matrix P**T is\n*          required, as defined in the transformation applied by SGEBRD:\n*          = 'Q':  generate Q;\n*          = 'P':  generate P**T.\n*\n*  M       (input) INTEGER\n*          The number of rows of the matrix Q or P**T to be returned.\n*          M >= 0.\n*\n*  N       (input) INTEGER\n*          The number of columns of the matrix Q or P**T to be returned.\n*          N >= 0.\n*          If VECT = 'Q', M >= N >= min(M,K);\n*          if VECT = 'P', N >= M >= min(N,K).\n*\n*  K       (input) INTEGER\n*          If VECT = 'Q', the number of columns in the original M-by-K\n*          matrix reduced by SGEBRD.\n*          If VECT = 'P', the number of rows in the original K-by-N\n*          matrix reduced by SGEBRD.\n*          K >= 0.\n*\n*  A       (input/output) REAL array, dimension (LDA,N)\n*          On entry, the vectors which define the elementary reflectors,\n*          as returned by SGEBRD.\n*          On exit, the M-by-N matrix Q or P**T.\n*\n*  LDA     (input) INTEGER\n*          The leading dimension of the array A. LDA >= max(1,M).\n*\n*  TAU     (input) REAL array, dimension\n*                                (min(M,K)) if VECT = 'Q'\n*                                (min(N,K)) if VECT = 'P'\n*          TAU(i) must contain the scalar factor of the elementary\n*          reflector H(i) or G(i), which determines Q or P**T, as\n*          returned by SGEBRD in its array argument TAUQ or TAUP.\n*\n*  WORK    (workspace/output) REAL array, dimension (MAX(1,LWORK))\n*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.\n*\n*  LWORK   (input) INTEGER\n*          The dimension of the array WORK. LWORK >= max(1,min(M,N)).\n*          For optimum performance LWORK >= min(M,N)*NB, where NB\n*          is the optimal blocksize.\n*\n*          If LWORK = -1, then a workspace query is assumed; the routine\n*          only calculates the optimal size of the WORK array, returns\n*          this value as the first entry of the WORK array, and no error\n*          message related to LWORK is issued by XERBLA.\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit\n*          < 0:  if INFO = -i, the i-th argument had an illegal value\n*\n\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 6)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 6)", argc);
  rb_vect = argv[0];
  rb_m = argv[1];
  rb_k = argv[2];
  rb_a = argv[3];
  rb_tau = argv[4];
  rb_lwork = argv[5];

  k = NUM2INT(rb_k);
  if (!NA_IsNArray(rb_a))
    rb_raise(rb_eArgError, "a (4th argument) must be NArray");
  if (NA_RANK(rb_a) != 2)
    rb_raise(rb_eArgError, "rank of a (4th argument) must be %d", 2);
  n = NA_SHAPE1(rb_a);
  lda = NA_SHAPE0(rb_a);
  if (NA_TYPE(rb_a) != NA_SFLOAT)
    rb_a = na_change_type(rb_a, NA_SFLOAT);
  a = NA_PTR_TYPE(rb_a, real*);
  m = NUM2INT(rb_m);
  lwork = NUM2INT(rb_lwork);
  vect = StringValueCStr(rb_vect)[0];
  if (!NA_IsNArray(rb_tau))
    rb_raise(rb_eArgError, "tau (5th argument) must be NArray");
  if (NA_RANK(rb_tau) != 1)
    rb_raise(rb_eArgError, "rank of tau (5th argument) must be %d", 1);
  if (NA_SHAPE0(rb_tau) != (MIN(m,k)))
    rb_raise(rb_eRuntimeError, "shape 0 of tau must be %d", MIN(m,k));
  if (NA_TYPE(rb_tau) != NA_SFLOAT)
    rb_tau = na_change_type(rb_tau, NA_SFLOAT);
  tau = NA_PTR_TYPE(rb_tau, real*);
  {
    int shape[1];
    shape[0] = MAX(1,lwork);
    rb_work = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  work = NA_PTR_TYPE(rb_work, real*);
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

  sorgbr_(&vect, &m, &n, &k, a, &lda, tau, work, &lwork, &info);

  rb_info = INT2NUM(info);
  return rb_ary_new3(3, rb_work, rb_info, rb_a);
}

void
init_lapack_sorgbr(VALUE mLapack){
  rb_define_module_function(mLapack, "sorgbr", rb_sorgbr, -1);
}

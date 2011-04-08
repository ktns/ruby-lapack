#include "rb_lapack.h"

extern VOID cungbr_(char *vect, integer *m, integer *n, integer *k, complex *a, integer *lda, complex *tau, complex *work, integer *lwork, integer *info);

static VALUE sHelp, sUsage;

static VALUE
rblapack_cungbr(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_vect;
  char vect; 
  VALUE rblapack_m;
  integer m; 
  VALUE rblapack_k;
  integer k; 
  VALUE rblapack_a;
  complex *a; 
  VALUE rblapack_tau;
  complex *tau; 
  VALUE rblapack_lwork;
  integer lwork; 
  VALUE rblapack_work;
  complex *work; 
  VALUE rblapack_info;
  integer info; 
  VALUE rblapack_a_out__;
  complex *a_out__;

  integer lda;
  integer n;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  work, info, a = NumRu::Lapack.cungbr( vect, m, k, a, tau, lwork, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE CUNGBR( VECT, M, N, K, A, LDA, TAU, WORK, LWORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  CUNGBR generates one of the complex unitary matrices Q or P**H\n*  determined by CGEBRD when reducing a complex matrix A to bidiagonal\n*  form: A = Q * B * P**H.  Q and P**H are defined as products of\n*  elementary reflectors H(i) or G(i) respectively.\n*\n*  If VECT = 'Q', A is assumed to have been an M-by-K matrix, and Q\n*  is of order M:\n*  if m >= k, Q = H(1) H(2) . . . H(k) and CUNGBR returns the first n\n*  columns of Q, where m >= n >= k;\n*  if m < k, Q = H(1) H(2) . . . H(m-1) and CUNGBR returns Q as an\n*  M-by-M matrix.\n*\n*  If VECT = 'P', A is assumed to have been a K-by-N matrix, and P**H\n*  is of order N:\n*  if k < n, P**H = G(k) . . . G(2) G(1) and CUNGBR returns the first m\n*  rows of P**H, where n >= m >= k;\n*  if k >= n, P**H = G(n-1) . . . G(2) G(1) and CUNGBR returns P**H as\n*  an N-by-N matrix.\n*\n\n*  Arguments\n*  =========\n*\n*  VECT    (input) CHARACTER*1\n*          Specifies whether the matrix Q or the matrix P**H is\n*          required, as defined in the transformation applied by CGEBRD:\n*          = 'Q':  generate Q;\n*          = 'P':  generate P**H.\n*\n*  M       (input) INTEGER\n*          The number of rows of the matrix Q or P**H to be returned.\n*          M >= 0.\n*\n*  N       (input) INTEGER\n*          The number of columns of the matrix Q or P**H to be returned.\n*          N >= 0.\n*          If VECT = 'Q', M >= N >= min(M,K);\n*          if VECT = 'P', N >= M >= min(N,K).\n*\n*  K       (input) INTEGER\n*          If VECT = 'Q', the number of columns in the original M-by-K\n*          matrix reduced by CGEBRD.\n*          If VECT = 'P', the number of rows in the original K-by-N\n*          matrix reduced by CGEBRD.\n*          K >= 0.\n*\n*  A       (input/output) COMPLEX array, dimension (LDA,N)\n*          On entry, the vectors which define the elementary reflectors,\n*          as returned by CGEBRD.\n*          On exit, the M-by-N matrix Q or P**H.\n*\n*  LDA     (input) INTEGER\n*          The leading dimension of the array A. LDA >= M.\n*\n*  TAU     (input) COMPLEX array, dimension\n*                                (min(M,K)) if VECT = 'Q'\n*                                (min(N,K)) if VECT = 'P'\n*          TAU(i) must contain the scalar factor of the elementary\n*          reflector H(i) or G(i), which determines Q or P**H, as\n*          returned by CGEBRD in its array argument TAUQ or TAUP.\n*\n*  WORK    (workspace/output) COMPLEX array, dimension (MAX(1,LWORK))\n*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.\n*\n*  LWORK   (input) INTEGER\n*          The dimension of the array WORK. LWORK >= max(1,min(M,N)).\n*          For optimum performance LWORK >= min(M,N)*NB, where NB\n*          is the optimal blocksize.\n*\n*          If LWORK = -1, then a workspace query is assumed; the routine\n*          only calculates the optimal size of the WORK array, returns\n*          this value as the first entry of the WORK array, and no error\n*          message related to LWORK is issued by XERBLA.\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit\n*          < 0:  if INFO = -i, the i-th argument had an illegal value\n*\n\n*  =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  work, info, a = NumRu::Lapack.cungbr( vect, m, k, a, tau, lwork, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 6)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 6)", argc);
  rblapack_vect = argv[0];
  rblapack_m = argv[1];
  rblapack_k = argv[2];
  rblapack_a = argv[3];
  rblapack_tau = argv[4];
  rblapack_lwork = argv[5];
  if (rb_options != Qnil) {
  }

  k = NUM2INT(rblapack_k);
  if (!NA_IsNArray(rblapack_a))
    rb_raise(rb_eArgError, "a (4th argument) must be NArray");
  if (NA_RANK(rblapack_a) != 2)
    rb_raise(rb_eArgError, "rank of a (4th argument) must be %d", 2);
  n = NA_SHAPE1(rblapack_a);
  lda = NA_SHAPE0(rblapack_a);
  if (NA_TYPE(rblapack_a) != NA_SCOMPLEX)
    rblapack_a = na_change_type(rblapack_a, NA_SCOMPLEX);
  a = NA_PTR_TYPE(rblapack_a, complex*);
  m = NUM2INT(rblapack_m);
  lwork = NUM2INT(rblapack_lwork);
  vect = StringValueCStr(rblapack_vect)[0];
  if (!NA_IsNArray(rblapack_tau))
    rb_raise(rb_eArgError, "tau (5th argument) must be NArray");
  if (NA_RANK(rblapack_tau) != 1)
    rb_raise(rb_eArgError, "rank of tau (5th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_tau) != (MIN(m,k)))
    rb_raise(rb_eRuntimeError, "shape 0 of tau must be %d", MIN(m,k));
  if (NA_TYPE(rblapack_tau) != NA_SCOMPLEX)
    rblapack_tau = na_change_type(rblapack_tau, NA_SCOMPLEX);
  tau = NA_PTR_TYPE(rblapack_tau, complex*);
  {
    int shape[1];
    shape[0] = MAX(1,lwork);
    rblapack_work = na_make_object(NA_SCOMPLEX, 1, shape, cNArray);
  }
  work = NA_PTR_TYPE(rblapack_work, complex*);
  {
    int shape[2];
    shape[0] = lda;
    shape[1] = n;
    rblapack_a_out__ = na_make_object(NA_SCOMPLEX, 2, shape, cNArray);
  }
  a_out__ = NA_PTR_TYPE(rblapack_a_out__, complex*);
  MEMCPY(a_out__, a, complex, NA_TOTAL(rblapack_a));
  rblapack_a = rblapack_a_out__;
  a = a_out__;

  cungbr_(&vect, &m, &n, &k, a, &lda, tau, work, &lwork, &info);

  rblapack_info = INT2NUM(info);
  return rb_ary_new3(3, rblapack_work, rblapack_info, rblapack_a);
}

void
init_lapack_cungbr(VALUE mLapack){
  rb_define_module_function(mLapack, "cungbr", rblapack_cungbr, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}

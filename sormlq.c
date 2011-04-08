#include "rb_lapack.h"

extern VOID sormlq_(char *side, char *trans, integer *m, integer *n, integer *k, real *a, integer *lda, real *tau, real *c, integer *ldc, real *work, integer *lwork, integer *info);

static VALUE sHelp, sUsage;

static VALUE
rblapack_sormlq(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_side;
  char side; 
  VALUE rblapack_trans;
  char trans; 
  VALUE rblapack_a;
  real *a; 
  VALUE rblapack_tau;
  real *tau; 
  VALUE rblapack_c;
  real *c; 
  VALUE rblapack_lwork;
  integer lwork; 
  VALUE rblapack_work;
  real *work; 
  VALUE rblapack_info;
  integer info; 
  VALUE rblapack_c_out__;
  real *c_out__;

  integer lda;
  integer m;
  integer k;
  integer ldc;
  integer n;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  work, info, c = NumRu::Lapack.sormlq( side, trans, a, tau, c, lwork, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE SORMLQ( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, WORK, LWORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  SORMLQ overwrites the general real M-by-N matrix C with\n*\n*                  SIDE = 'L'     SIDE = 'R'\n*  TRANS = 'N':      Q * C          C * Q\n*  TRANS = 'T':      Q**T * C       C * Q**T\n*\n*  where Q is a real orthogonal matrix defined as the product of k\n*  elementary reflectors\n*\n*        Q = H(k) . . . H(2) H(1)\n*\n*  as returned by SGELQF. Q is of order M if SIDE = 'L' and of order N\n*  if SIDE = 'R'.\n*\n\n*  Arguments\n*  =========\n*\n*  SIDE    (input) CHARACTER*1\n*          = 'L': apply Q or Q**T from the Left;\n*          = 'R': apply Q or Q**T from the Right.\n*\n*  TRANS   (input) CHARACTER*1\n*          = 'N':  No transpose, apply Q;\n*          = 'T':  Transpose, apply Q**T.\n*\n*  M       (input) INTEGER\n*          The number of rows of the matrix C. M >= 0.\n*\n*  N       (input) INTEGER\n*          The number of columns of the matrix C. N >= 0.\n*\n*  K       (input) INTEGER\n*          The number of elementary reflectors whose product defines\n*          the matrix Q.\n*          If SIDE = 'L', M >= K >= 0;\n*          if SIDE = 'R', N >= K >= 0.\n*\n*  A       (input) REAL array, dimension\n*                               (LDA,M) if SIDE = 'L',\n*                               (LDA,N) if SIDE = 'R'\n*          The i-th row must contain the vector which defines the\n*          elementary reflector H(i), for i = 1,2,...,k, as returned by\n*          SGELQF in the first k rows of its array argument A.\n*          A is modified by the routine but restored on exit.\n*\n*  LDA     (input) INTEGER\n*          The leading dimension of the array A. LDA >= max(1,K).\n*\n*  TAU     (input) REAL array, dimension (K)\n*          TAU(i) must contain the scalar factor of the elementary\n*          reflector H(i), as returned by SGELQF.\n*\n*  C       (input/output) REAL array, dimension (LDC,N)\n*          On entry, the M-by-N matrix C.\n*          On exit, C is overwritten by Q*C or Q**T*C or C*Q**T or C*Q.\n*\n*  LDC     (input) INTEGER\n*          The leading dimension of the array C. LDC >= max(1,M).\n*\n*  WORK    (workspace/output) REAL array, dimension (MAX(1,LWORK))\n*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.\n*\n*  LWORK   (input) INTEGER\n*          The dimension of the array WORK.\n*          If SIDE = 'L', LWORK >= max(1,N);\n*          if SIDE = 'R', LWORK >= max(1,M).\n*          For optimum performance LWORK >= N*NB if SIDE = 'L', and\n*          LWORK >= M*NB if SIDE = 'R', where NB is the optimal\n*          blocksize.\n*\n*          If LWORK = -1, then a workspace query is assumed; the routine\n*          only calculates the optimal size of the WORK array, returns\n*          this value as the first entry of the WORK array, and no error\n*          message related to LWORK is issued by XERBLA.\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit\n*          < 0:  if INFO = -i, the i-th argument had an illegal value\n*\n\n*  =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  work, info, c = NumRu::Lapack.sormlq( side, trans, a, tau, c, lwork, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 6)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 6)", argc);
  rblapack_side = argv[0];
  rblapack_trans = argv[1];
  rblapack_a = argv[2];
  rblapack_tau = argv[3];
  rblapack_c = argv[4];
  rblapack_lwork = argv[5];
  if (rb_options != Qnil) {
  }

  if (!NA_IsNArray(rblapack_a))
    rb_raise(rb_eArgError, "a (3th argument) must be NArray");
  if (NA_RANK(rblapack_a) != 2)
    rb_raise(rb_eArgError, "rank of a (3th argument) must be %d", 2);
  m = NA_SHAPE1(rblapack_a);
  lda = NA_SHAPE0(rblapack_a);
  if (NA_TYPE(rblapack_a) != NA_SFLOAT)
    rblapack_a = na_change_type(rblapack_a, NA_SFLOAT);
  a = NA_PTR_TYPE(rblapack_a, real*);
  side = StringValueCStr(rblapack_side)[0];
  if (!NA_IsNArray(rblapack_c))
    rb_raise(rb_eArgError, "c (5th argument) must be NArray");
  if (NA_RANK(rblapack_c) != 2)
    rb_raise(rb_eArgError, "rank of c (5th argument) must be %d", 2);
  n = NA_SHAPE1(rblapack_c);
  ldc = NA_SHAPE0(rblapack_c);
  if (NA_TYPE(rblapack_c) != NA_SFLOAT)
    rblapack_c = na_change_type(rblapack_c, NA_SFLOAT);
  c = NA_PTR_TYPE(rblapack_c, real*);
  lwork = NUM2INT(rblapack_lwork);
  if (!NA_IsNArray(rblapack_tau))
    rb_raise(rb_eArgError, "tau (4th argument) must be NArray");
  if (NA_RANK(rblapack_tau) != 1)
    rb_raise(rb_eArgError, "rank of tau (4th argument) must be %d", 1);
  k = NA_SHAPE0(rblapack_tau);
  if (NA_TYPE(rblapack_tau) != NA_SFLOAT)
    rblapack_tau = na_change_type(rblapack_tau, NA_SFLOAT);
  tau = NA_PTR_TYPE(rblapack_tau, real*);
  trans = StringValueCStr(rblapack_trans)[0];
  {
    int shape[1];
    shape[0] = MAX(1,lwork);
    rblapack_work = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  work = NA_PTR_TYPE(rblapack_work, real*);
  {
    int shape[2];
    shape[0] = ldc;
    shape[1] = n;
    rblapack_c_out__ = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  c_out__ = NA_PTR_TYPE(rblapack_c_out__, real*);
  MEMCPY(c_out__, c, real, NA_TOTAL(rblapack_c));
  rblapack_c = rblapack_c_out__;
  c = c_out__;

  sormlq_(&side, &trans, &m, &n, &k, a, &lda, tau, c, &ldc, work, &lwork, &info);

  rblapack_info = INT2NUM(info);
  return rb_ary_new3(3, rblapack_work, rblapack_info, rblapack_c);
}

void
init_lapack_sormlq(VALUE mLapack){
  rb_define_module_function(mLapack, "sormlq", rblapack_sormlq, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}

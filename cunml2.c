#include "rb_lapack.h"

extern VOID cunml2_(char *side, char *trans, integer *m, integer *n, integer *k, complex *a, integer *lda, complex *tau, complex *c, integer *ldc, complex *work, integer *info);

static VALUE sHelp, sUsage;

static VALUE
rblapack_cunml2(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_side;
  char side; 
  VALUE rblapack_trans;
  char trans; 
  VALUE rblapack_a;
  complex *a; 
  VALUE rblapack_tau;
  complex *tau; 
  VALUE rblapack_c;
  complex *c; 
  VALUE rblapack_info;
  integer info; 
  VALUE rblapack_c_out__;
  complex *c_out__;
  complex *work;

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
      printf("%s\n", "USAGE:\n  info, c = NumRu::Lapack.cunml2( side, trans, a, tau, c, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE CUNML2( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, WORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  CUNML2 overwrites the general complex m-by-n matrix C with\n*\n*        Q * C  if SIDE = 'L' and TRANS = 'N', or\n*\n*        Q'* C  if SIDE = 'L' and TRANS = 'C', or\n*\n*        C * Q  if SIDE = 'R' and TRANS = 'N', or\n*\n*        C * Q' if SIDE = 'R' and TRANS = 'C',\n*\n*  where Q is a complex unitary matrix defined as the product of k\n*  elementary reflectors\n*\n*        Q = H(k)' . . . H(2)' H(1)'\n*\n*  as returned by CGELQF. Q is of order m if SIDE = 'L' and of order n\n*  if SIDE = 'R'.\n*\n\n*  Arguments\n*  =========\n*\n*  SIDE    (input) CHARACTER*1\n*          = 'L': apply Q or Q' from the Left\n*          = 'R': apply Q or Q' from the Right\n*\n*  TRANS   (input) CHARACTER*1\n*          = 'N': apply Q  (No transpose)\n*          = 'C': apply Q' (Conjugate transpose)\n*\n*  M       (input) INTEGER\n*          The number of rows of the matrix C. M >= 0.\n*\n*  N       (input) INTEGER\n*          The number of columns of the matrix C. N >= 0.\n*\n*  K       (input) INTEGER\n*          The number of elementary reflectors whose product defines\n*          the matrix Q.\n*          If SIDE = 'L', M >= K >= 0;\n*          if SIDE = 'R', N >= K >= 0.\n*\n*  A       (input) COMPLEX array, dimension\n*                               (LDA,M) if SIDE = 'L',\n*                               (LDA,N) if SIDE = 'R'\n*          The i-th row must contain the vector which defines the\n*          elementary reflector H(i), for i = 1,2,...,k, as returned by\n*          CGELQF in the first k rows of its array argument A.\n*          A is modified by the routine but restored on exit.\n*\n*  LDA     (input) INTEGER\n*          The leading dimension of the array A. LDA >= max(1,K).\n*\n*  TAU     (input) COMPLEX array, dimension (K)\n*          TAU(i) must contain the scalar factor of the elementary\n*          reflector H(i), as returned by CGELQF.\n*\n*  C       (input/output) COMPLEX array, dimension (LDC,N)\n*          On entry, the m-by-n matrix C.\n*          On exit, C is overwritten by Q*C or Q'*C or C*Q' or C*Q.\n*\n*  LDC     (input) INTEGER\n*          The leading dimension of the array C. LDC >= max(1,M).\n*\n*  WORK    (workspace) COMPLEX array, dimension\n*                                   (N) if SIDE = 'L',\n*                                   (M) if SIDE = 'R'\n*\n*  INFO    (output) INTEGER\n*          = 0: successful exit\n*          < 0: if INFO = -i, the i-th argument had an illegal value\n*\n\n*  =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  info, c = NumRu::Lapack.cunml2( side, trans, a, tau, c, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 5)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 5)", argc);
  rblapack_side = argv[0];
  rblapack_trans = argv[1];
  rblapack_a = argv[2];
  rblapack_tau = argv[3];
  rblapack_c = argv[4];
  if (rb_options != Qnil) {
  }

  if (!NA_IsNArray(rblapack_a))
    rb_raise(rb_eArgError, "a (3th argument) must be NArray");
  if (NA_RANK(rblapack_a) != 2)
    rb_raise(rb_eArgError, "rank of a (3th argument) must be %d", 2);
  m = NA_SHAPE1(rblapack_a);
  lda = NA_SHAPE0(rblapack_a);
  if (NA_TYPE(rblapack_a) != NA_SCOMPLEX)
    rblapack_a = na_change_type(rblapack_a, NA_SCOMPLEX);
  a = NA_PTR_TYPE(rblapack_a, complex*);
  side = StringValueCStr(rblapack_side)[0];
  if (!NA_IsNArray(rblapack_c))
    rb_raise(rb_eArgError, "c (5th argument) must be NArray");
  if (NA_RANK(rblapack_c) != 2)
    rb_raise(rb_eArgError, "rank of c (5th argument) must be %d", 2);
  n = NA_SHAPE1(rblapack_c);
  ldc = NA_SHAPE0(rblapack_c);
  if (NA_TYPE(rblapack_c) != NA_SCOMPLEX)
    rblapack_c = na_change_type(rblapack_c, NA_SCOMPLEX);
  c = NA_PTR_TYPE(rblapack_c, complex*);
  if (!NA_IsNArray(rblapack_tau))
    rb_raise(rb_eArgError, "tau (4th argument) must be NArray");
  if (NA_RANK(rblapack_tau) != 1)
    rb_raise(rb_eArgError, "rank of tau (4th argument) must be %d", 1);
  k = NA_SHAPE0(rblapack_tau);
  if (NA_TYPE(rblapack_tau) != NA_SCOMPLEX)
    rblapack_tau = na_change_type(rblapack_tau, NA_SCOMPLEX);
  tau = NA_PTR_TYPE(rblapack_tau, complex*);
  trans = StringValueCStr(rblapack_trans)[0];
  {
    int shape[2];
    shape[0] = ldc;
    shape[1] = n;
    rblapack_c_out__ = na_make_object(NA_SCOMPLEX, 2, shape, cNArray);
  }
  c_out__ = NA_PTR_TYPE(rblapack_c_out__, complex*);
  MEMCPY(c_out__, c, complex, NA_TOTAL(rblapack_c));
  rblapack_c = rblapack_c_out__;
  c = c_out__;
  work = ALLOC_N(complex, (lsame_(&side,"L") ? n : lsame_(&side,"R") ? m : 0));

  cunml2_(&side, &trans, &m, &n, &k, a, &lda, tau, c, &ldc, work, &info);

  free(work);
  rblapack_info = INT2NUM(info);
  return rb_ary_new3(2, rblapack_info, rblapack_c);
}

void
init_lapack_cunml2(VALUE mLapack){
  rb_define_module_function(mLapack, "cunml2", rblapack_cunml2, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}

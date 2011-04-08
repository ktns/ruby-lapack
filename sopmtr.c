#include "rb_lapack.h"

extern VOID sopmtr_(char *side, char *uplo, char *trans, integer *m, integer *n, real *ap, real *tau, real *c, integer *ldc, real *work, integer *info);

static VALUE sHelp, sUsage;

static VALUE
rblapack_sopmtr(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_side;
  char side; 
  VALUE rblapack_uplo;
  char uplo; 
  VALUE rblapack_trans;
  char trans; 
  VALUE rblapack_m;
  integer m; 
  VALUE rblapack_ap;
  real *ap; 
  VALUE rblapack_tau;
  real *tau; 
  VALUE rblapack_c;
  real *c; 
  VALUE rblapack_info;
  integer info; 
  VALUE rblapack_c_out__;
  real *c_out__;
  real *work;

  integer ldc;
  integer n;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  info, c = NumRu::Lapack.sopmtr( side, uplo, trans, m, ap, tau, c, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE SOPMTR( SIDE, UPLO, TRANS, M, N, AP, TAU, C, LDC, WORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  SOPMTR overwrites the general real M-by-N matrix C with\n*\n*                  SIDE = 'L'     SIDE = 'R'\n*  TRANS = 'N':      Q * C          C * Q\n*  TRANS = 'T':      Q**T * C       C * Q**T\n*\n*  where Q is a real orthogonal matrix of order nq, with nq = m if\n*  SIDE = 'L' and nq = n if SIDE = 'R'. Q is defined as the product of\n*  nq-1 elementary reflectors, as returned by SSPTRD using packed\n*  storage:\n*\n*  if UPLO = 'U', Q = H(nq-1) . . . H(2) H(1);\n*\n*  if UPLO = 'L', Q = H(1) H(2) . . . H(nq-1).\n*\n\n*  Arguments\n*  =========\n*\n*  SIDE    (input) CHARACTER*1\n*          = 'L': apply Q or Q**T from the Left;\n*          = 'R': apply Q or Q**T from the Right.\n*\n*  UPLO    (input) CHARACTER*1\n*          = 'U': Upper triangular packed storage used in previous\n*                 call to SSPTRD;\n*          = 'L': Lower triangular packed storage used in previous\n*                 call to SSPTRD.\n*\n*  TRANS   (input) CHARACTER*1\n*          = 'N':  No transpose, apply Q;\n*          = 'T':  Transpose, apply Q**T.\n*\n*  M       (input) INTEGER\n*          The number of rows of the matrix C. M >= 0.\n*\n*  N       (input) INTEGER\n*          The number of columns of the matrix C. N >= 0.\n*\n*  AP      (input) REAL array, dimension\n*                               (M*(M+1)/2) if SIDE = 'L'\n*                               (N*(N+1)/2) if SIDE = 'R'\n*          The vectors which define the elementary reflectors, as\n*          returned by SSPTRD.  AP is modified by the routine but\n*          restored on exit.\n*\n*  TAU     (input) REAL array, dimension (M-1) if SIDE = 'L'\n*                                     or (N-1) if SIDE = 'R'\n*          TAU(i) must contain the scalar factor of the elementary\n*          reflector H(i), as returned by SSPTRD.\n*\n*  C       (input/output) REAL array, dimension (LDC,N)\n*          On entry, the M-by-N matrix C.\n*          On exit, C is overwritten by Q*C or Q**T*C or C*Q**T or C*Q.\n*\n*  LDC     (input) INTEGER\n*          The leading dimension of the array C. LDC >= max(1,M).\n*\n*  WORK    (workspace) REAL array, dimension\n*                                   (N) if SIDE = 'L'\n*                                   (M) if SIDE = 'R'\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit\n*          < 0:  if INFO = -i, the i-th argument had an illegal value\n*\n\n*  =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  info, c = NumRu::Lapack.sopmtr( side, uplo, trans, m, ap, tau, c, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 7)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 7)", argc);
  rblapack_side = argv[0];
  rblapack_uplo = argv[1];
  rblapack_trans = argv[2];
  rblapack_m = argv[3];
  rblapack_ap = argv[4];
  rblapack_tau = argv[5];
  rblapack_c = argv[6];
  if (rb_options != Qnil) {
  }

  side = StringValueCStr(rblapack_side)[0];
  m = NUM2INT(rblapack_m);
  if (!NA_IsNArray(rblapack_c))
    rb_raise(rb_eArgError, "c (7th argument) must be NArray");
  if (NA_RANK(rblapack_c) != 2)
    rb_raise(rb_eArgError, "rank of c (7th argument) must be %d", 2);
  n = NA_SHAPE1(rblapack_c);
  ldc = NA_SHAPE0(rblapack_c);
  if (NA_TYPE(rblapack_c) != NA_SFLOAT)
    rblapack_c = na_change_type(rblapack_c, NA_SFLOAT);
  c = NA_PTR_TYPE(rblapack_c, real*);
  trans = StringValueCStr(rblapack_trans)[0];
  uplo = StringValueCStr(rblapack_uplo)[0];
  if (!NA_IsNArray(rblapack_tau))
    rb_raise(rb_eArgError, "tau (6th argument) must be NArray");
  if (NA_RANK(rblapack_tau) != 1)
    rb_raise(rb_eArgError, "rank of tau (6th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_tau) != (m-1))
    rb_raise(rb_eRuntimeError, "shape 0 of tau must be %d", m-1);
  if (NA_TYPE(rblapack_tau) != NA_SFLOAT)
    rblapack_tau = na_change_type(rblapack_tau, NA_SFLOAT);
  tau = NA_PTR_TYPE(rblapack_tau, real*);
  if (!NA_IsNArray(rblapack_ap))
    rb_raise(rb_eArgError, "ap (5th argument) must be NArray");
  if (NA_RANK(rblapack_ap) != 1)
    rb_raise(rb_eArgError, "rank of ap (5th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_ap) != (m*(m+1)/2))
    rb_raise(rb_eRuntimeError, "shape 0 of ap must be %d", m*(m+1)/2);
  if (NA_TYPE(rblapack_ap) != NA_SFLOAT)
    rblapack_ap = na_change_type(rblapack_ap, NA_SFLOAT);
  ap = NA_PTR_TYPE(rblapack_ap, real*);
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
  work = ALLOC_N(real, (lsame_(&side,"L") ? n : lsame_(&side,"R") ? m : 0));

  sopmtr_(&side, &uplo, &trans, &m, &n, ap, tau, c, &ldc, work, &info);

  free(work);
  rblapack_info = INT2NUM(info);
  return rb_ary_new3(2, rblapack_info, rblapack_c);
}

void
init_lapack_sopmtr(VALUE mLapack){
  rb_define_module_function(mLapack, "sopmtr", rblapack_sopmtr, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}

#include "rb_lapack.h"

extern VOID zlarfb_(char *side, char *trans, char *direct, char *storev, integer *m, integer *n, integer *k, doublecomplex *v, integer *ldv, doublecomplex *t, integer *ldt, doublecomplex *c, integer *ldc, doublecomplex *work, integer *ldwork);

static VALUE sHelp, sUsage;

static VALUE
rblapack_zlarfb(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_side;
  char side; 
  VALUE rblapack_trans;
  char trans; 
  VALUE rblapack_direct;
  char direct; 
  VALUE rblapack_storev;
  char storev; 
  VALUE rblapack_m;
  integer m; 
  VALUE rblapack_v;
  doublecomplex *v; 
  VALUE rblapack_t;
  doublecomplex *t; 
  VALUE rblapack_c;
  doublecomplex *c; 
  VALUE rblapack_c_out__;
  doublecomplex *c_out__;
  doublecomplex *work;

  integer ldv;
  integer k;
  integer ldt;
  integer ldc;
  integer n;
  integer ldwork;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  c = NumRu::Lapack.zlarfb( side, trans, direct, storev, m, v, t, c, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE ZLARFB( SIDE, TRANS, DIRECT, STOREV, M, N, K, V, LDV, T, LDT, C, LDC, WORK, LDWORK )\n\n*  Purpose\n*  =======\n*\n*  ZLARFB applies a complex block reflector H or its transpose H' to a\n*  complex M-by-N matrix C, from either the left or the right.\n*\n\n*  Arguments\n*  =========\n*\n*  SIDE    (input) CHARACTER*1\n*          = 'L': apply H or H' from the Left\n*          = 'R': apply H or H' from the Right\n*\n*  TRANS   (input) CHARACTER*1\n*          = 'N': apply H (No transpose)\n*          = 'C': apply H' (Conjugate transpose)\n*\n*  DIRECT  (input) CHARACTER*1\n*          Indicates how H is formed from a product of elementary\n*          reflectors\n*          = 'F': H = H(1) H(2) . . . H(k) (Forward)\n*          = 'B': H = H(k) . . . H(2) H(1) (Backward)\n*\n*  STOREV  (input) CHARACTER*1\n*          Indicates how the vectors which define the elementary\n*          reflectors are stored:\n*          = 'C': Columnwise\n*          = 'R': Rowwise\n*\n*  M       (input) INTEGER\n*          The number of rows of the matrix C.\n*\n*  N       (input) INTEGER\n*          The number of columns of the matrix C.\n*\n*  K       (input) INTEGER\n*          The order of the matrix T (= the number of elementary\n*          reflectors whose product defines the block reflector).\n*\n*  V       (input) COMPLEX*16 array, dimension\n*                                (LDV,K) if STOREV = 'C'\n*                                (LDV,M) if STOREV = 'R' and SIDE = 'L'\n*                                (LDV,N) if STOREV = 'R' and SIDE = 'R'\n*          The matrix V. See further details.\n*\n*  LDV     (input) INTEGER\n*          The leading dimension of the array V.\n*          If STOREV = 'C' and SIDE = 'L', LDV >= max(1,M);\n*          if STOREV = 'C' and SIDE = 'R', LDV >= max(1,N);\n*          if STOREV = 'R', LDV >= K.\n*\n*  T       (input) COMPLEX*16 array, dimension (LDT,K)\n*          The triangular K-by-K matrix T in the representation of the\n*          block reflector.\n*\n*  LDT     (input) INTEGER\n*          The leading dimension of the array T. LDT >= K.\n*\n*  C       (input/output) COMPLEX*16 array, dimension (LDC,N)\n*          On entry, the M-by-N matrix C.\n*          On exit, C is overwritten by H*C or H'*C or C*H or C*H'.\n*\n*  LDC     (input) INTEGER\n*          The leading dimension of the array C. LDC >= max(1,M).\n*\n*  WORK    (workspace) COMPLEX*16 array, dimension (LDWORK,K)\n*\n*  LDWORK  (input) INTEGER\n*          The leading dimension of the array WORK.\n*          If SIDE = 'L', LDWORK >= max(1,N);\n*          if SIDE = 'R', LDWORK >= max(1,M).\n*\n\n*  =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  c = NumRu::Lapack.zlarfb( side, trans, direct, storev, m, v, t, c, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 8)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 8)", argc);
  rblapack_side = argv[0];
  rblapack_trans = argv[1];
  rblapack_direct = argv[2];
  rblapack_storev = argv[3];
  rblapack_m = argv[4];
  rblapack_v = argv[5];
  rblapack_t = argv[6];
  rblapack_c = argv[7];
  if (rb_options != Qnil) {
  }

  trans = StringValueCStr(rblapack_trans)[0];
  if (!NA_IsNArray(rblapack_v))
    rb_raise(rb_eArgError, "v (6th argument) must be NArray");
  if (NA_RANK(rblapack_v) != 2)
    rb_raise(rb_eArgError, "rank of v (6th argument) must be %d", 2);
  k = NA_SHAPE1(rblapack_v);
  ldv = NA_SHAPE0(rblapack_v);
  if (NA_TYPE(rblapack_v) != NA_DCOMPLEX)
    rblapack_v = na_change_type(rblapack_v, NA_DCOMPLEX);
  v = NA_PTR_TYPE(rblapack_v, doublecomplex*);
  direct = StringValueCStr(rblapack_direct)[0];
  side = StringValueCStr(rblapack_side)[0];
  storev = StringValueCStr(rblapack_storev)[0];
  if (!NA_IsNArray(rblapack_c))
    rb_raise(rb_eArgError, "c (8th argument) must be NArray");
  if (NA_RANK(rblapack_c) != 2)
    rb_raise(rb_eArgError, "rank of c (8th argument) must be %d", 2);
  n = NA_SHAPE1(rblapack_c);
  ldc = NA_SHAPE0(rblapack_c);
  if (NA_TYPE(rblapack_c) != NA_DCOMPLEX)
    rblapack_c = na_change_type(rblapack_c, NA_DCOMPLEX);
  c = NA_PTR_TYPE(rblapack_c, doublecomplex*);
  if (!NA_IsNArray(rblapack_t))
    rb_raise(rb_eArgError, "t (7th argument) must be NArray");
  if (NA_RANK(rblapack_t) != 2)
    rb_raise(rb_eArgError, "rank of t (7th argument) must be %d", 2);
  if (NA_SHAPE1(rblapack_t) != k)
    rb_raise(rb_eRuntimeError, "shape 1 of t must be the same as shape 1 of v");
  ldt = NA_SHAPE0(rblapack_t);
  if (NA_TYPE(rblapack_t) != NA_DCOMPLEX)
    rblapack_t = na_change_type(rblapack_t, NA_DCOMPLEX);
  t = NA_PTR_TYPE(rblapack_t, doublecomplex*);
  m = NUM2INT(rblapack_m);
  ldwork = max(1,n) ? side = 'l' : max(1,m) ? side = 'r' : 0;
  {
    int shape[2];
    shape[0] = ldc;
    shape[1] = n;
    rblapack_c_out__ = na_make_object(NA_DCOMPLEX, 2, shape, cNArray);
  }
  c_out__ = NA_PTR_TYPE(rblapack_c_out__, doublecomplex*);
  MEMCPY(c_out__, c, doublecomplex, NA_TOTAL(rblapack_c));
  rblapack_c = rblapack_c_out__;
  c = c_out__;
  work = ALLOC_N(doublecomplex, (ldwork)*(k));

  zlarfb_(&side, &trans, &direct, &storev, &m, &n, &k, v, &ldv, t, &ldt, c, &ldc, work, &ldwork);

  free(work);
  return rblapack_c;
}

void
init_lapack_zlarfb(VALUE mLapack){
  rb_define_module_function(mLapack, "zlarfb", rblapack_zlarfb, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}

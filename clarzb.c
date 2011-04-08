#include "rb_lapack.h"

extern VOID clarzb_(char *side, char *trans, char *direct, char *storev, integer *m, integer *n, integer *k, integer *l, complex *v, integer *ldv, complex *t, integer *ldt, complex *c, integer *ldc, complex *work, integer *ldwork);

static VALUE sHelp, sUsage;

static VALUE
rblapack_clarzb(int argc, VALUE *argv, VALUE self){
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
  VALUE rblapack_l;
  integer l; 
  VALUE rblapack_v;
  complex *v; 
  VALUE rblapack_t;
  complex *t; 
  VALUE rblapack_c;
  complex *c; 
  VALUE rblapack_c_out__;
  complex *c_out__;
  complex *work;

  integer ldv;
  integer nv;
  integer ldt;
  integer k;
  integer ldc;
  integer n;
  integer ldwork;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  c = NumRu::Lapack.clarzb( side, trans, direct, storev, m, l, v, t, c, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE CLARZB( SIDE, TRANS, DIRECT, STOREV, M, N, K, L, V, LDV, T, LDT, C, LDC, WORK, LDWORK )\n\n*  Purpose\n*  =======\n*\n*  CLARZB applies a complex block reflector H or its transpose H**H\n*  to a complex distributed M-by-N  C from the left or the right.\n*\n*  Currently, only STOREV = 'R' and DIRECT = 'B' are supported.\n*\n\n*  Arguments\n*  =========\n*\n*  SIDE    (input) CHARACTER*1\n*          = 'L': apply H or H' from the Left\n*          = 'R': apply H or H' from the Right\n*\n*  TRANS   (input) CHARACTER*1\n*          = 'N': apply H (No transpose)\n*          = 'C': apply H' (Conjugate transpose)\n*\n*  DIRECT  (input) CHARACTER*1\n*          Indicates how H is formed from a product of elementary\n*          reflectors\n*          = 'F': H = H(1) H(2) . . . H(k) (Forward, not supported yet)\n*          = 'B': H = H(k) . . . H(2) H(1) (Backward)\n*\n*  STOREV  (input) CHARACTER*1\n*          Indicates how the vectors which define the elementary\n*          reflectors are stored:\n*          = 'C': Columnwise                        (not supported yet)\n*          = 'R': Rowwise\n*\n*  M       (input) INTEGER\n*          The number of rows of the matrix C.\n*\n*  N       (input) INTEGER\n*          The number of columns of the matrix C.\n*\n*  K       (input) INTEGER\n*          The order of the matrix T (= the number of elementary\n*          reflectors whose product defines the block reflector).\n*\n*  L       (input) INTEGER\n*          The number of columns of the matrix V containing the\n*          meaningful part of the Householder reflectors.\n*          If SIDE = 'L', M >= L >= 0, if SIDE = 'R', N >= L >= 0.\n*\n*  V       (input) COMPLEX array, dimension (LDV,NV).\n*          If STOREV = 'C', NV = K; if STOREV = 'R', NV = L.\n*\n*  LDV     (input) INTEGER\n*          The leading dimension of the array V.\n*          If STOREV = 'C', LDV >= L; if STOREV = 'R', LDV >= K.\n*\n*  T       (input) COMPLEX array, dimension (LDT,K)\n*          The triangular K-by-K matrix T in the representation of the\n*          block reflector.\n*\n*  LDT     (input) INTEGER\n*          The leading dimension of the array T. LDT >= K.\n*\n*  C       (input/output) COMPLEX array, dimension (LDC,N)\n*          On entry, the M-by-N matrix C.\n*          On exit, C is overwritten by H*C or H'*C or C*H or C*H'.\n*\n*  LDC     (input) INTEGER\n*          The leading dimension of the array C. LDC >= max(1,M).\n*\n*  WORK    (workspace) COMPLEX array, dimension (LDWORK,K)\n*\n*  LDWORK  (input) INTEGER\n*          The leading dimension of the array WORK.\n*          If SIDE = 'L', LDWORK >= max(1,N);\n*          if SIDE = 'R', LDWORK >= max(1,M).\n*\n\n*  Further Details\n*  ===============\n*\n*  Based on contributions by\n*    A. Petitet, Computer Science Dept., Univ. of Tenn., Knoxville, USA\n*\n*  =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  c = NumRu::Lapack.clarzb( side, trans, direct, storev, m, l, v, t, c, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 9)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 9)", argc);
  rblapack_side = argv[0];
  rblapack_trans = argv[1];
  rblapack_direct = argv[2];
  rblapack_storev = argv[3];
  rblapack_m = argv[4];
  rblapack_l = argv[5];
  rblapack_v = argv[6];
  rblapack_t = argv[7];
  rblapack_c = argv[8];
  if (rb_options != Qnil) {
  }

  trans = StringValueCStr(rblapack_trans)[0];
  if (!NA_IsNArray(rblapack_v))
    rb_raise(rb_eArgError, "v (7th argument) must be NArray");
  if (NA_RANK(rblapack_v) != 2)
    rb_raise(rb_eArgError, "rank of v (7th argument) must be %d", 2);
  nv = NA_SHAPE1(rblapack_v);
  ldv = NA_SHAPE0(rblapack_v);
  if (NA_TYPE(rblapack_v) != NA_SCOMPLEX)
    rblapack_v = na_change_type(rblapack_v, NA_SCOMPLEX);
  v = NA_PTR_TYPE(rblapack_v, complex*);
  direct = StringValueCStr(rblapack_direct)[0];
  l = NUM2INT(rblapack_l);
  side = StringValueCStr(rblapack_side)[0];
  storev = StringValueCStr(rblapack_storev)[0];
  if (!NA_IsNArray(rblapack_c))
    rb_raise(rb_eArgError, "c (9th argument) must be NArray");
  if (NA_RANK(rblapack_c) != 2)
    rb_raise(rb_eArgError, "rank of c (9th argument) must be %d", 2);
  n = NA_SHAPE1(rblapack_c);
  ldc = NA_SHAPE0(rblapack_c);
  if (NA_TYPE(rblapack_c) != NA_SCOMPLEX)
    rblapack_c = na_change_type(rblapack_c, NA_SCOMPLEX);
  c = NA_PTR_TYPE(rblapack_c, complex*);
  if (!NA_IsNArray(rblapack_t))
    rb_raise(rb_eArgError, "t (8th argument) must be NArray");
  if (NA_RANK(rblapack_t) != 2)
    rb_raise(rb_eArgError, "rank of t (8th argument) must be %d", 2);
  k = NA_SHAPE1(rblapack_t);
  ldt = NA_SHAPE0(rblapack_t);
  if (NA_TYPE(rblapack_t) != NA_SCOMPLEX)
    rblapack_t = na_change_type(rblapack_t, NA_SCOMPLEX);
  t = NA_PTR_TYPE(rblapack_t, complex*);
  m = NUM2INT(rblapack_m);
  ldwork = max(1,n) ? side = 'l' : max(1,m) ? side = 'r' : 0;
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
  work = ALLOC_N(complex, (ldwork)*(k));

  clarzb_(&side, &trans, &direct, &storev, &m, &n, &k, &l, v, &ldv, t, &ldt, c, &ldc, work, &ldwork);

  free(work);
  return rblapack_c;
}

void
init_lapack_clarzb(VALUE mLapack){
  rb_define_module_function(mLapack, "clarzb", rblapack_clarzb, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}

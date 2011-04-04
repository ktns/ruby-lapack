#include "rb_lapack.h"

extern VOID dlarfb_(char *side, char *trans, char *direct, char *storev, integer *m, integer *n, integer *k, doublereal *v, integer *ldv, doublereal *t, integer *ldt, doublereal *c, integer *ldc, doublereal *work, integer *ldwork);

static VALUE
rb_dlarfb(int argc, VALUE *argv, VALUE self){
  VALUE rb_side;
  char side; 
  VALUE rb_trans;
  char trans; 
  VALUE rb_direct;
  char direct; 
  VALUE rb_storev;
  char storev; 
  VALUE rb_m;
  integer m; 
  VALUE rb_v;
  doublereal *v; 
  VALUE rb_t;
  doublereal *t; 
  VALUE rb_c;
  doublereal *c; 
  VALUE rb_c_out__;
  doublereal *c_out__;
  doublereal *work;

  integer ldv;
  integer k;
  integer ldt;
  integer ldc;
  integer n;
  integer ldwork;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  c = NumRu::Lapack.dlarfb( side, trans, direct, storev, m, v, t, c)\n    or\n  NumRu::Lapack.dlarfb  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE DLARFB( SIDE, TRANS, DIRECT, STOREV, M, N, K, V, LDV, T, LDT, C, LDC, WORK, LDWORK )\n\n*  Purpose\n*  =======\n*\n*  DLARFB applies a real block reflector H or its transpose H' to a\n*  real m by n matrix C, from either the left or the right.\n*\n\n*  Arguments\n*  =========\n*\n*  SIDE    (input) CHARACTER*1\n*          = 'L': apply H or H' from the Left\n*          = 'R': apply H or H' from the Right\n*\n*  TRANS   (input) CHARACTER*1\n*          = 'N': apply H (No transpose)\n*          = 'T': apply H' (Transpose)\n*\n*  DIRECT  (input) CHARACTER*1\n*          Indicates how H is formed from a product of elementary\n*          reflectors\n*          = 'F': H = H(1) H(2) . . . H(k) (Forward)\n*          = 'B': H = H(k) . . . H(2) H(1) (Backward)\n*\n*  STOREV  (input) CHARACTER*1\n*          Indicates how the vectors which define the elementary\n*          reflectors are stored:\n*          = 'C': Columnwise\n*          = 'R': Rowwise\n*\n*  M       (input) INTEGER\n*          The number of rows of the matrix C.\n*\n*  N       (input) INTEGER\n*          The number of columns of the matrix C.\n*\n*  K       (input) INTEGER\n*          The order of the matrix T (= the number of elementary\n*          reflectors whose product defines the block reflector).\n*\n*  V       (input) DOUBLE PRECISION array, dimension\n*                                (LDV,K) if STOREV = 'C'\n*                                (LDV,M) if STOREV = 'R' and SIDE = 'L'\n*                                (LDV,N) if STOREV = 'R' and SIDE = 'R'\n*          The matrix V. See further details.\n*\n*  LDV     (input) INTEGER\n*          The leading dimension of the array V.\n*          If STOREV = 'C' and SIDE = 'L', LDV >= max(1,M);\n*          if STOREV = 'C' and SIDE = 'R', LDV >= max(1,N);\n*          if STOREV = 'R', LDV >= K.\n*\n*  T       (input) DOUBLE PRECISION array, dimension (LDT,K)\n*          The triangular k by k matrix T in the representation of the\n*          block reflector.\n*\n*  LDT     (input) INTEGER\n*          The leading dimension of the array T. LDT >= K.\n*\n*  C       (input/output) DOUBLE PRECISION array, dimension (LDC,N)\n*          On entry, the m by n matrix C.\n*          On exit, C is overwritten by H*C or H'*C or C*H or C*H'.\n*\n*  LDC     (input) INTEGER\n*          The leading dimension of the array C. LDA >= max(1,M).\n*\n*  WORK    (workspace) DOUBLE PRECISION array, dimension (LDWORK,K)\n*\n*  LDWORK  (input) INTEGER\n*          The leading dimension of the array WORK.\n*          If SIDE = 'L', LDWORK >= max(1,N);\n*          if SIDE = 'R', LDWORK >= max(1,M).\n*\n\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 8)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 8)", argc);
  rb_side = argv[0];
  rb_trans = argv[1];
  rb_direct = argv[2];
  rb_storev = argv[3];
  rb_m = argv[4];
  rb_v = argv[5];
  rb_t = argv[6];
  rb_c = argv[7];

  trans = StringValueCStr(rb_trans)[0];
  if (!NA_IsNArray(rb_v))
    rb_raise(rb_eArgError, "v (6th argument) must be NArray");
  if (NA_RANK(rb_v) != 2)
    rb_raise(rb_eArgError, "rank of v (6th argument) must be %d", 2);
  k = NA_SHAPE1(rb_v);
  ldv = NA_SHAPE0(rb_v);
  if (NA_TYPE(rb_v) != NA_DFLOAT)
    rb_v = na_change_type(rb_v, NA_DFLOAT);
  v = NA_PTR_TYPE(rb_v, doublereal*);
  direct = StringValueCStr(rb_direct)[0];
  side = StringValueCStr(rb_side)[0];
  storev = StringValueCStr(rb_storev)[0];
  if (!NA_IsNArray(rb_c))
    rb_raise(rb_eArgError, "c (8th argument) must be NArray");
  if (NA_RANK(rb_c) != 2)
    rb_raise(rb_eArgError, "rank of c (8th argument) must be %d", 2);
  n = NA_SHAPE1(rb_c);
  ldc = NA_SHAPE0(rb_c);
  if (NA_TYPE(rb_c) != NA_DFLOAT)
    rb_c = na_change_type(rb_c, NA_DFLOAT);
  c = NA_PTR_TYPE(rb_c, doublereal*);
  if (!NA_IsNArray(rb_t))
    rb_raise(rb_eArgError, "t (7th argument) must be NArray");
  if (NA_RANK(rb_t) != 2)
    rb_raise(rb_eArgError, "rank of t (7th argument) must be %d", 2);
  if (NA_SHAPE1(rb_t) != k)
    rb_raise(rb_eRuntimeError, "shape 1 of t must be the same as shape 1 of v");
  ldt = NA_SHAPE0(rb_t);
  if (NA_TYPE(rb_t) != NA_DFLOAT)
    rb_t = na_change_type(rb_t, NA_DFLOAT);
  t = NA_PTR_TYPE(rb_t, doublereal*);
  m = NUM2INT(rb_m);
  ldwork = max(1,n) ? side = 'l' : max(1,m) ? side = 'r' : 0;
  {
    int shape[2];
    shape[0] = ldc;
    shape[1] = n;
    rb_c_out__ = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  c_out__ = NA_PTR_TYPE(rb_c_out__, doublereal*);
  MEMCPY(c_out__, c, doublereal, NA_TOTAL(rb_c));
  rb_c = rb_c_out__;
  c = c_out__;
  work = ALLOC_N(doublereal, (ldwork)*(k));

  dlarfb_(&side, &trans, &direct, &storev, &m, &n, &k, v, &ldv, t, &ldt, c, &ldc, work, &ldwork);

  free(work);
  return rb_c;
}

void
init_lapack_dlarfb(VALUE mLapack){
  rb_define_module_function(mLapack, "dlarfb", rb_dlarfb, -1);
}

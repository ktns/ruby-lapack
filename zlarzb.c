#include "rb_lapack.h"

static VALUE
rb_zlarzb(int argc, VALUE *argv, VALUE self){
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
  VALUE rb_l;
  integer l; 
  VALUE rb_v;
  doublecomplex *v; 
  VALUE rb_t;
  doublecomplex *t; 
  VALUE rb_c;
  doublecomplex *c; 
  VALUE rb_c_out__;
  doublecomplex *c_out__;
  doublecomplex *work;

  integer ldv;
  integer nv;
  integer ldt;
  integer k;
  integer ldc;
  integer n;
  integer ldwork;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  c = NumRu::Lapack.zlarzb( side, trans, direct, storev, m, l, v, t, c)\n    or\n  NumRu::Lapack.zlarzb  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE ZLARZB( SIDE, TRANS, DIRECT, STOREV, M, N, K, L, V, LDV, T, LDT, C, LDC, WORK, LDWORK )\n\n*  Purpose\n*  =======\n*\n*  ZLARZB applies a complex block reflector H or its transpose H**H\n*  to a complex distributed M-by-N  C from the left or the right.\n*\n*  Currently, only STOREV = 'R' and DIRECT = 'B' are supported.\n*\n\n*  Arguments\n*  =========\n*\n*  SIDE    (input) CHARACTER*1\n*          = 'L': apply H or H' from the Left\n*          = 'R': apply H or H' from the Right\n*\n*  TRANS   (input) CHARACTER*1\n*          = 'N': apply H (No transpose)\n*          = 'C': apply H' (Conjugate transpose)\n*\n*  DIRECT  (input) CHARACTER*1\n*          Indicates how H is formed from a product of elementary\n*          reflectors\n*          = 'F': H = H(1) H(2) . . . H(k) (Forward, not supported yet)\n*          = 'B': H = H(k) . . . H(2) H(1) (Backward)\n*\n*  STOREV  (input) CHARACTER*1\n*          Indicates how the vectors which define the elementary\n*          reflectors are stored:\n*          = 'C': Columnwise                        (not supported yet)\n*          = 'R': Rowwise\n*\n*  M       (input) INTEGER\n*          The number of rows of the matrix C.\n*\n*  N       (input) INTEGER\n*          The number of columns of the matrix C.\n*\n*  K       (input) INTEGER\n*          The order of the matrix T (= the number of elementary\n*          reflectors whose product defines the block reflector).\n*\n*  L       (input) INTEGER\n*          The number of columns of the matrix V containing the\n*          meaningful part of the Householder reflectors.\n*          If SIDE = 'L', M >= L >= 0, if SIDE = 'R', N >= L >= 0.\n*\n*  V       (input) COMPLEX*16 array, dimension (LDV,NV).\n*          If STOREV = 'C', NV = K; if STOREV = 'R', NV = L.\n*\n*  LDV     (input) INTEGER\n*          The leading dimension of the array V.\n*          If STOREV = 'C', LDV >= L; if STOREV = 'R', LDV >= K.\n*\n*  T       (input) COMPLEX*16 array, dimension (LDT,K)\n*          The triangular K-by-K matrix T in the representation of the\n*          block reflector.\n*\n*  LDT     (input) INTEGER\n*          The leading dimension of the array T. LDT >= K.\n*\n*  C       (input/output) COMPLEX*16 array, dimension (LDC,N)\n*          On entry, the M-by-N matrix C.\n*          On exit, C is overwritten by H*C or H'*C or C*H or C*H'.\n*\n*  LDC     (input) INTEGER\n*          The leading dimension of the array C. LDC >= max(1,M).\n*\n*  WORK    (workspace) COMPLEX*16 array, dimension (LDWORK,K)\n*\n*  LDWORK  (input) INTEGER\n*          The leading dimension of the array WORK.\n*          If SIDE = 'L', LDWORK >= max(1,N);\n*          if SIDE = 'R', LDWORK >= max(1,M).\n*\n\n*  Further Details\n*  ===============\n*\n*  Based on contributions by\n*    A. Petitet, Computer Science Dept., Univ. of Tenn., Knoxville, USA\n*\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 9)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 9)", argc);
  rb_side = argv[0];
  rb_trans = argv[1];
  rb_direct = argv[2];
  rb_storev = argv[3];
  rb_m = argv[4];
  rb_l = argv[5];
  rb_v = argv[6];
  rb_t = argv[7];
  rb_c = argv[8];

  side = StringValueCStr(rb_side)[0];
  trans = StringValueCStr(rb_trans)[0];
  direct = StringValueCStr(rb_direct)[0];
  storev = StringValueCStr(rb_storev)[0];
  m = NUM2INT(rb_m);
  l = NUM2INT(rb_l);
  if (!NA_IsNArray(rb_v))
    rb_raise(rb_eArgError, "v (7th argument) must be NArray");
  if (NA_RANK(rb_v) != 2)
    rb_raise(rb_eArgError, "rank of v (7th argument) must be %d", 2);
  ldv = NA_SHAPE0(rb_v);
  nv = NA_SHAPE1(rb_v);
  if (NA_TYPE(rb_v) != NA_DCOMPLEX)
    rb_v = na_change_type(rb_v, NA_DCOMPLEX);
  v = NA_PTR_TYPE(rb_v, doublecomplex*);
  if (!NA_IsNArray(rb_t))
    rb_raise(rb_eArgError, "t (8th argument) must be NArray");
  if (NA_RANK(rb_t) != 2)
    rb_raise(rb_eArgError, "rank of t (8th argument) must be %d", 2);
  ldt = NA_SHAPE0(rb_t);
  k = NA_SHAPE1(rb_t);
  if (NA_TYPE(rb_t) != NA_DCOMPLEX)
    rb_t = na_change_type(rb_t, NA_DCOMPLEX);
  t = NA_PTR_TYPE(rb_t, doublecomplex*);
  if (!NA_IsNArray(rb_c))
    rb_raise(rb_eArgError, "c (9th argument) must be NArray");
  if (NA_RANK(rb_c) != 2)
    rb_raise(rb_eArgError, "rank of c (9th argument) must be %d", 2);
  ldc = NA_SHAPE0(rb_c);
  n = NA_SHAPE1(rb_c);
  if (NA_TYPE(rb_c) != NA_DCOMPLEX)
    rb_c = na_change_type(rb_c, NA_DCOMPLEX);
  c = NA_PTR_TYPE(rb_c, doublecomplex*);
  {
    int shape[2];
    shape[0] = ldc;
    shape[1] = n;
    rb_c_out__ = na_make_object(NA_DCOMPLEX, 2, shape, cNArray);
  }
  c_out__ = NA_PTR_TYPE(rb_c_out__, doublecomplex*);
  MEMCPY(c_out__, c, doublecomplex, NA_TOTAL(rb_c));
  rb_c = rb_c_out__;
  c = c_out__;
  ldwork = max(1,n) ? side = 'l' : max(1,m) ? side = 'r' : 0;
  work = ALLOC_N(doublecomplex, (ldwork)*(k));

  zlarzb_(&side, &trans, &direct, &storev, &m, &n, &k, &l, v, &ldv, t, &ldt, c, &ldc, work, &ldwork);

  free(work);
  return rb_c;
}

void
init_lapack_zlarzb(VALUE mLapack){
  rb_define_module_function(mLapack, "zlarzb", rb_zlarzb, -1);
}

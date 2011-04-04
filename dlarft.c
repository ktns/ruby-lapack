#include "rb_lapack.h"

extern VOID dlarft_(char *direct, char *storev, integer *n, integer *k, doublereal *v, integer *ldv, doublereal *tau, doublereal *t, integer *ldt);

static VALUE
rb_dlarft(int argc, VALUE *argv, VALUE self){
  VALUE rb_direct;
  char direct; 
  VALUE rb_storev;
  char storev; 
  VALUE rb_n;
  integer n; 
  VALUE rb_v;
  doublereal *v; 
  VALUE rb_tau;
  doublereal *tau; 
  VALUE rb_t;
  doublereal *t; 
  VALUE rb_v_out__;
  doublereal *v_out__;

  integer ldv;
  integer k;
  integer ldt;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  t, v = NumRu::Lapack.dlarft( direct, storev, n, v, tau)\n    or\n  NumRu::Lapack.dlarft  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE DLARFT( DIRECT, STOREV, N, K, V, LDV, TAU, T, LDT )\n\n*  Purpose\n*  =======\n*\n*  DLARFT forms the triangular factor T of a real block reflector H\n*  of order n, which is defined as a product of k elementary reflectors.\n*\n*  If DIRECT = 'F', H = H(1) H(2) . . . H(k) and T is upper triangular;\n*\n*  If DIRECT = 'B', H = H(k) . . . H(2) H(1) and T is lower triangular.\n*\n*  If STOREV = 'C', the vector which defines the elementary reflector\n*  H(i) is stored in the i-th column of the array V, and\n*\n*     H  =  I - V * T * V'\n*\n*  If STOREV = 'R', the vector which defines the elementary reflector\n*  H(i) is stored in the i-th row of the array V, and\n*\n*     H  =  I - V' * T * V\n*\n\n*  Arguments\n*  =========\n*\n*  DIRECT  (input) CHARACTER*1\n*          Specifies the order in which the elementary reflectors are\n*          multiplied to form the block reflector:\n*          = 'F': H = H(1) H(2) . . . H(k) (Forward)\n*          = 'B': H = H(k) . . . H(2) H(1) (Backward)\n*\n*  STOREV  (input) CHARACTER*1\n*          Specifies how the vectors which define the elementary\n*          reflectors are stored (see also Further Details):\n*          = 'C': columnwise\n*          = 'R': rowwise\n*\n*  N       (input) INTEGER\n*          The order of the block reflector H. N >= 0.\n*\n*  K       (input) INTEGER\n*          The order of the triangular factor T (= the number of\n*          elementary reflectors). K >= 1.\n*\n*  V       (input/output) DOUBLE PRECISION array, dimension\n*                               (LDV,K) if STOREV = 'C'\n*                               (LDV,N) if STOREV = 'R'\n*          The matrix V. See further details.\n*\n*  LDV     (input) INTEGER\n*          The leading dimension of the array V.\n*          If STOREV = 'C', LDV >= max(1,N); if STOREV = 'R', LDV >= K.\n*\n*  TAU     (input) DOUBLE PRECISION array, dimension (K)\n*          TAU(i) must contain the scalar factor of the elementary\n*          reflector H(i).\n*\n*  T       (output) DOUBLE PRECISION array, dimension (LDT,K)\n*          The k by k triangular factor T of the block reflector.\n*          If DIRECT = 'F', T is upper triangular; if DIRECT = 'B', T is\n*          lower triangular. The rest of the array is not used.\n*\n*  LDT     (input) INTEGER\n*          The leading dimension of the array T. LDT >= K.\n*\n\n*  Further Details\n*  ===============\n*\n*  The shape of the matrix V and the storage of the vectors which define\n*  the H(i) is best illustrated by the following example with n = 5 and\n*  k = 3. The elements equal to 1 are not stored; the corresponding\n*  array elements are modified but restored on exit. The rest of the\n*  array is not used.\n*\n*  DIRECT = 'F' and STOREV = 'C':         DIRECT = 'F' and STOREV = 'R':\n*\n*               V = (  1       )                 V = (  1 v1 v1 v1 v1 )\n*                   ( v1  1    )                     (     1 v2 v2 v2 )\n*                   ( v1 v2  1 )                     (        1 v3 v3 )\n*                   ( v1 v2 v3 )\n*                   ( v1 v2 v3 )\n*\n*  DIRECT = 'B' and STOREV = 'C':         DIRECT = 'B' and STOREV = 'R':\n*\n*               V = ( v1 v2 v3 )                 V = ( v1 v1  1       )\n*                   ( v1 v2 v3 )                     ( v2 v2 v2  1    )\n*                   (  1 v2 v3 )                     ( v3 v3 v3 v3  1 )\n*                   (     1 v3 )\n*                   (        1 )\n*\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 5)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 5)", argc);
  rb_direct = argv[0];
  rb_storev = argv[1];
  rb_n = argv[2];
  rb_v = argv[3];
  rb_tau = argv[4];

  storev = StringValueCStr(rb_storev)[0];
  direct = StringValueCStr(rb_direct)[0];
  n = NUM2INT(rb_n);
  if (!NA_IsNArray(rb_tau))
    rb_raise(rb_eArgError, "tau (5th argument) must be NArray");
  if (NA_RANK(rb_tau) != 1)
    rb_raise(rb_eArgError, "rank of tau (5th argument) must be %d", 1);
  k = NA_SHAPE0(rb_tau);
  if (NA_TYPE(rb_tau) != NA_DFLOAT)
    rb_tau = na_change_type(rb_tau, NA_DFLOAT);
  tau = NA_PTR_TYPE(rb_tau, doublereal*);
  ldt = k;
  if (!NA_IsNArray(rb_v))
    rb_raise(rb_eArgError, "v (4th argument) must be NArray");
  if (NA_RANK(rb_v) != 2)
    rb_raise(rb_eArgError, "rank of v (4th argument) must be %d", 2);
  if (NA_SHAPE1(rb_v) != (lsame_(&storev,"C") ? k : lsame_(&storev,"R") ? n : 0))
    rb_raise(rb_eRuntimeError, "shape 1 of v must be %d", lsame_(&storev,"C") ? k : lsame_(&storev,"R") ? n : 0);
  ldv = NA_SHAPE0(rb_v);
  if (NA_TYPE(rb_v) != NA_DFLOAT)
    rb_v = na_change_type(rb_v, NA_DFLOAT);
  v = NA_PTR_TYPE(rb_v, doublereal*);
  {
    int shape[2];
    shape[0] = ldt;
    shape[1] = k;
    rb_t = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  t = NA_PTR_TYPE(rb_t, doublereal*);
  {
    int shape[2];
    shape[0] = ldv;
    shape[1] = lsame_(&storev,"C") ? k : lsame_(&storev,"R") ? n : 0;
    rb_v_out__ = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  v_out__ = NA_PTR_TYPE(rb_v_out__, doublereal*);
  MEMCPY(v_out__, v, doublereal, NA_TOTAL(rb_v));
  rb_v = rb_v_out__;
  v = v_out__;

  dlarft_(&direct, &storev, &n, &k, v, &ldv, tau, t, &ldt);

  return rb_ary_new3(2, rb_t, rb_v);
}

void
init_lapack_dlarft(VALUE mLapack){
  rb_define_module_function(mLapack, "dlarft", rb_dlarft, -1);
}

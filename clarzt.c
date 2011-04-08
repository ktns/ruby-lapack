#include "rb_lapack.h"

extern VOID clarzt_(char *direct, char *storev, integer *n, integer *k, complex *v, integer *ldv, complex *tau, complex *t, integer *ldt);

static VALUE sHelp, sUsage;

static VALUE
rblapack_clarzt(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_direct;
  char direct; 
  VALUE rblapack_storev;
  char storev; 
  VALUE rblapack_n;
  integer n; 
  VALUE rblapack_v;
  complex *v; 
  VALUE rblapack_tau;
  complex *tau; 
  VALUE rblapack_t;
  complex *t; 
  VALUE rblapack_v_out__;
  complex *v_out__;

  integer ldv;
  integer k;
  integer ldt;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  t, v = NumRu::Lapack.clarzt( direct, storev, n, v, tau, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE CLARZT( DIRECT, STOREV, N, K, V, LDV, TAU, T, LDT )\n\n*  Purpose\n*  =======\n*\n*  CLARZT forms the triangular factor T of a complex block reflector\n*  H of order > n, which is defined as a product of k elementary\n*  reflectors.\n*\n*  If DIRECT = 'F', H = H(1) H(2) . . . H(k) and T is upper triangular;\n*\n*  If DIRECT = 'B', H = H(k) . . . H(2) H(1) and T is lower triangular.\n*\n*  If STOREV = 'C', the vector which defines the elementary reflector\n*  H(i) is stored in the i-th column of the array V, and\n*\n*     H  =  I - V * T * V'\n*\n*  If STOREV = 'R', the vector which defines the elementary reflector\n*  H(i) is stored in the i-th row of the array V, and\n*\n*     H  =  I - V' * T * V\n*\n*  Currently, only STOREV = 'R' and DIRECT = 'B' are supported.\n*\n\n*  Arguments\n*  =========\n*\n*  DIRECT  (input) CHARACTER*1\n*          Specifies the order in which the elementary reflectors are\n*          multiplied to form the block reflector:\n*          = 'F': H = H(1) H(2) . . . H(k) (Forward, not supported yet)\n*          = 'B': H = H(k) . . . H(2) H(1) (Backward)\n*\n*  STOREV  (input) CHARACTER*1\n*          Specifies how the vectors which define the elementary\n*          reflectors are stored (see also Further Details):\n*          = 'C': columnwise                        (not supported yet)\n*          = 'R': rowwise\n*\n*  N       (input) INTEGER\n*          The order of the block reflector H. N >= 0.\n*\n*  K       (input) INTEGER\n*          The order of the triangular factor T (= the number of\n*          elementary reflectors). K >= 1.\n*\n*  V       (input/output) COMPLEX array, dimension\n*                               (LDV,K) if STOREV = 'C'\n*                               (LDV,N) if STOREV = 'R'\n*          The matrix V. See further details.\n*\n*  LDV     (input) INTEGER\n*          The leading dimension of the array V.\n*          If STOREV = 'C', LDV >= max(1,N); if STOREV = 'R', LDV >= K.\n*\n*  TAU     (input) COMPLEX array, dimension (K)\n*          TAU(i) must contain the scalar factor of the elementary\n*          reflector H(i).\n*\n*  T       (output) COMPLEX array, dimension (LDT,K)\n*          The k by k triangular factor T of the block reflector.\n*          If DIRECT = 'F', T is upper triangular; if DIRECT = 'B', T is\n*          lower triangular. The rest of the array is not used.\n*\n*  LDT     (input) INTEGER\n*          The leading dimension of the array T. LDT >= K.\n*\n\n*  Further Details\n*  ===============\n*\n*  Based on contributions by\n*    A. Petitet, Computer Science Dept., Univ. of Tenn., Knoxville, USA\n*\n*  The shape of the matrix V and the storage of the vectors which define\n*  the H(i) is best illustrated by the following example with n = 5 and\n*  k = 3. The elements equal to 1 are not stored; the corresponding\n*  array elements are modified but restored on exit. The rest of the\n*  array is not used.\n*\n*  DIRECT = 'F' and STOREV = 'C':         DIRECT = 'F' and STOREV = 'R':\n*\n*                                              ______V_____\n*         ( v1 v2 v3 )                        /            \\\n*         ( v1 v2 v3 )                      ( v1 v1 v1 v1 v1 . . . . 1 )\n*     V = ( v1 v2 v3 )                      ( v2 v2 v2 v2 v2 . . . 1   )\n*         ( v1 v2 v3 )                      ( v3 v3 v3 v3 v3 . . 1     )\n*         ( v1 v2 v3 )\n*            .  .  .\n*            .  .  .\n*            1  .  .\n*               1  .\n*                  1\n*\n*  DIRECT = 'B' and STOREV = 'C':         DIRECT = 'B' and STOREV = 'R':\n*\n*                                                        ______V_____\n*            1                                          /            \\\n*            .  1                           ( 1 . . . . v1 v1 v1 v1 v1 )\n*            .  .  1                        ( . 1 . . . v2 v2 v2 v2 v2 )\n*            .  .  .                        ( . . 1 . . v3 v3 v3 v3 v3 )\n*            .  .  .\n*         ( v1 v2 v3 )\n*         ( v1 v2 v3 )\n*     V = ( v1 v2 v3 )\n*         ( v1 v2 v3 )\n*         ( v1 v2 v3 )\n*\n*  =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  t, v = NumRu::Lapack.clarzt( direct, storev, n, v, tau, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 5)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 5)", argc);
  rblapack_direct = argv[0];
  rblapack_storev = argv[1];
  rblapack_n = argv[2];
  rblapack_v = argv[3];
  rblapack_tau = argv[4];
  if (rb_options != Qnil) {
  }

  storev = StringValueCStr(rblapack_storev)[0];
  direct = StringValueCStr(rblapack_direct)[0];
  n = NUM2INT(rblapack_n);
  if (!NA_IsNArray(rblapack_tau))
    rb_raise(rb_eArgError, "tau (5th argument) must be NArray");
  if (NA_RANK(rblapack_tau) != 1)
    rb_raise(rb_eArgError, "rank of tau (5th argument) must be %d", 1);
  k = NA_SHAPE0(rblapack_tau);
  if (NA_TYPE(rblapack_tau) != NA_SCOMPLEX)
    rblapack_tau = na_change_type(rblapack_tau, NA_SCOMPLEX);
  tau = NA_PTR_TYPE(rblapack_tau, complex*);
  ldt = k;
  if (!NA_IsNArray(rblapack_v))
    rb_raise(rb_eArgError, "v (4th argument) must be NArray");
  if (NA_RANK(rblapack_v) != 2)
    rb_raise(rb_eArgError, "rank of v (4th argument) must be %d", 2);
  if (NA_SHAPE1(rblapack_v) != (lsame_(&storev,"C") ? k : lsame_(&storev,"R") ? n : 0))
    rb_raise(rb_eRuntimeError, "shape 1 of v must be %d", lsame_(&storev,"C") ? k : lsame_(&storev,"R") ? n : 0);
  ldv = NA_SHAPE0(rblapack_v);
  if (NA_TYPE(rblapack_v) != NA_SCOMPLEX)
    rblapack_v = na_change_type(rblapack_v, NA_SCOMPLEX);
  v = NA_PTR_TYPE(rblapack_v, complex*);
  {
    int shape[2];
    shape[0] = ldt;
    shape[1] = k;
    rblapack_t = na_make_object(NA_SCOMPLEX, 2, shape, cNArray);
  }
  t = NA_PTR_TYPE(rblapack_t, complex*);
  {
    int shape[2];
    shape[0] = ldv;
    shape[1] = lsame_(&storev,"C") ? k : lsame_(&storev,"R") ? n : 0;
    rblapack_v_out__ = na_make_object(NA_SCOMPLEX, 2, shape, cNArray);
  }
  v_out__ = NA_PTR_TYPE(rblapack_v_out__, complex*);
  MEMCPY(v_out__, v, complex, NA_TOTAL(rblapack_v));
  rblapack_v = rblapack_v_out__;
  v = v_out__;

  clarzt_(&direct, &storev, &n, &k, v, &ldv, tau, t, &ldt);

  return rb_ary_new3(2, rblapack_t, rblapack_v);
}

void
init_lapack_clarzt(VALUE mLapack){
  rb_define_module_function(mLapack, "clarzt", rblapack_clarzt, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}

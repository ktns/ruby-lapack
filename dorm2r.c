#include "rb_lapack.h"

extern VOID dorm2r_(char *side, char *trans, integer *m, integer *n, integer *k, doublereal *a, integer *lda, doublereal *tau, doublereal *c, integer *ldc, doublereal *work, integer *info);

static VALUE sHelp, sUsage;

static VALUE
rblapack_dorm2r(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_side;
  char side; 
  VALUE rblapack_trans;
  char trans; 
  VALUE rblapack_m;
  integer m; 
  VALUE rblapack_a;
  doublereal *a; 
  VALUE rblapack_tau;
  doublereal *tau; 
  VALUE rblapack_c;
  doublereal *c; 
  VALUE rblapack_info;
  integer info; 
  VALUE rblapack_c_out__;
  doublereal *c_out__;
  doublereal *work;

  integer lda;
  integer k;
  integer ldc;
  integer n;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  info, c = NumRu::Lapack.dorm2r( side, trans, m, a, tau, c, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE DORM2R( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, WORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  DORM2R overwrites the general real m by n matrix C with\n*\n*        Q * C  if SIDE = 'L' and TRANS = 'N', or\n*\n*        Q'* C  if SIDE = 'L' and TRANS = 'T', or\n*\n*        C * Q  if SIDE = 'R' and TRANS = 'N', or\n*\n*        C * Q' if SIDE = 'R' and TRANS = 'T',\n*\n*  where Q is a real orthogonal matrix defined as the product of k\n*  elementary reflectors\n*\n*        Q = H(1) H(2) . . . H(k)\n*\n*  as returned by DGEQRF. Q is of order m if SIDE = 'L' and of order n\n*  if SIDE = 'R'.\n*\n\n*  Arguments\n*  =========\n*\n*  SIDE    (input) CHARACTER*1\n*          = 'L': apply Q or Q' from the Left\n*          = 'R': apply Q or Q' from the Right\n*\n*  TRANS   (input) CHARACTER*1\n*          = 'N': apply Q  (No transpose)\n*          = 'T': apply Q' (Transpose)\n*\n*  M       (input) INTEGER\n*          The number of rows of the matrix C. M >= 0.\n*\n*  N       (input) INTEGER\n*          The number of columns of the matrix C. N >= 0.\n*\n*  K       (input) INTEGER\n*          The number of elementary reflectors whose product defines\n*          the matrix Q.\n*          If SIDE = 'L', M >= K >= 0;\n*          if SIDE = 'R', N >= K >= 0.\n*\n*  A       (input) DOUBLE PRECISION array, dimension (LDA,K)\n*          The i-th column must contain the vector which defines the\n*          elementary reflector H(i), for i = 1,2,...,k, as returned by\n*          DGEQRF in the first k columns of its array argument A.\n*          A is modified by the routine but restored on exit.\n*\n*  LDA     (input) INTEGER\n*          The leading dimension of the array A.\n*          If SIDE = 'L', LDA >= max(1,M);\n*          if SIDE = 'R', LDA >= max(1,N).\n*\n*  TAU     (input) DOUBLE PRECISION array, dimension (K)\n*          TAU(i) must contain the scalar factor of the elementary\n*          reflector H(i), as returned by DGEQRF.\n*\n*  C       (input/output) DOUBLE PRECISION array, dimension (LDC,N)\n*          On entry, the m by n matrix C.\n*          On exit, C is overwritten by Q*C or Q'*C or C*Q' or C*Q.\n*\n*  LDC     (input) INTEGER\n*          The leading dimension of the array C. LDC >= max(1,M).\n*\n*  WORK    (workspace) DOUBLE PRECISION array, dimension\n*                                   (N) if SIDE = 'L',\n*                                   (M) if SIDE = 'R'\n*\n*  INFO    (output) INTEGER\n*          = 0: successful exit\n*          < 0: if INFO = -i, the i-th argument had an illegal value\n*\n\n*  =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  info, c = NumRu::Lapack.dorm2r( side, trans, m, a, tau, c, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 6)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 6)", argc);
  rblapack_side = argv[0];
  rblapack_trans = argv[1];
  rblapack_m = argv[2];
  rblapack_a = argv[3];
  rblapack_tau = argv[4];
  rblapack_c = argv[5];
  if (rb_options != Qnil) {
  }

  if (!NA_IsNArray(rblapack_a))
    rb_raise(rb_eArgError, "a (4th argument) must be NArray");
  if (NA_RANK(rblapack_a) != 2)
    rb_raise(rb_eArgError, "rank of a (4th argument) must be %d", 2);
  k = NA_SHAPE1(rblapack_a);
  lda = NA_SHAPE0(rblapack_a);
  if (NA_TYPE(rblapack_a) != NA_DFLOAT)
    rblapack_a = na_change_type(rblapack_a, NA_DFLOAT);
  a = NA_PTR_TYPE(rblapack_a, doublereal*);
  side = StringValueCStr(rblapack_side)[0];
  m = NUM2INT(rblapack_m);
  if (!NA_IsNArray(rblapack_c))
    rb_raise(rb_eArgError, "c (6th argument) must be NArray");
  if (NA_RANK(rblapack_c) != 2)
    rb_raise(rb_eArgError, "rank of c (6th argument) must be %d", 2);
  n = NA_SHAPE1(rblapack_c);
  ldc = NA_SHAPE0(rblapack_c);
  if (NA_TYPE(rblapack_c) != NA_DFLOAT)
    rblapack_c = na_change_type(rblapack_c, NA_DFLOAT);
  c = NA_PTR_TYPE(rblapack_c, doublereal*);
  if (!NA_IsNArray(rblapack_tau))
    rb_raise(rb_eArgError, "tau (5th argument) must be NArray");
  if (NA_RANK(rblapack_tau) != 1)
    rb_raise(rb_eArgError, "rank of tau (5th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_tau) != k)
    rb_raise(rb_eRuntimeError, "shape 0 of tau must be the same as shape 1 of a");
  if (NA_TYPE(rblapack_tau) != NA_DFLOAT)
    rblapack_tau = na_change_type(rblapack_tau, NA_DFLOAT);
  tau = NA_PTR_TYPE(rblapack_tau, doublereal*);
  trans = StringValueCStr(rblapack_trans)[0];
  {
    int shape[2];
    shape[0] = ldc;
    shape[1] = n;
    rblapack_c_out__ = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  c_out__ = NA_PTR_TYPE(rblapack_c_out__, doublereal*);
  MEMCPY(c_out__, c, doublereal, NA_TOTAL(rblapack_c));
  rblapack_c = rblapack_c_out__;
  c = c_out__;
  work = ALLOC_N(doublereal, (lsame_(&side,"L") ? n : lsame_(&side,"R") ? m : 0));

  dorm2r_(&side, &trans, &m, &n, &k, a, &lda, tau, c, &ldc, work, &info);

  free(work);
  rblapack_info = INT2NUM(info);
  return rb_ary_new3(2, rblapack_info, rblapack_c);
}

void
init_lapack_dorm2r(VALUE mLapack){
  rb_define_module_function(mLapack, "dorm2r", rblapack_dorm2r, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}

#include "rb_lapack.h"

extern VOID zlarz_(char *side, integer *m, integer *n, integer *l, doublecomplex *v, integer *incv, doublecomplex *tau, doublecomplex *c, integer *ldc, doublecomplex *work);

static VALUE sHelp, sUsage;

static VALUE
rblapack_zlarz(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_side;
  char side; 
  VALUE rblapack_m;
  integer m; 
  VALUE rblapack_l;
  integer l; 
  VALUE rblapack_v;
  doublecomplex *v; 
  VALUE rblapack_incv;
  integer incv; 
  VALUE rblapack_tau;
  doublecomplex tau; 
  VALUE rblapack_c;
  doublecomplex *c; 
  VALUE rblapack_c_out__;
  doublecomplex *c_out__;
  doublecomplex *work;

  integer ldc;
  integer n;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  c = NumRu::Lapack.zlarz( side, m, l, v, incv, tau, c, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE ZLARZ( SIDE, M, N, L, V, INCV, TAU, C, LDC, WORK )\n\n*  Purpose\n*  =======\n*\n*  ZLARZ applies a complex elementary reflector H to a complex\n*  M-by-N matrix C, from either the left or the right. H is represented\n*  in the form\n*\n*        H = I - tau * v * v'\n*\n*  where tau is a complex scalar and v is a complex vector.\n*\n*  If tau = 0, then H is taken to be the unit matrix.\n*\n*  To apply H' (the conjugate transpose of H), supply conjg(tau) instead\n*  tau.\n*\n*  H is a product of k elementary reflectors as returned by ZTZRZF.\n*\n\n*  Arguments\n*  =========\n*\n*  SIDE    (input) CHARACTER*1\n*          = 'L': form  H * C\n*          = 'R': form  C * H\n*\n*  M       (input) INTEGER\n*          The number of rows of the matrix C.\n*\n*  N       (input) INTEGER\n*          The number of columns of the matrix C.\n*\n*  L       (input) INTEGER\n*          The number of entries of the vector V containing\n*          the meaningful part of the Householder vectors.\n*          If SIDE = 'L', M >= L >= 0, if SIDE = 'R', N >= L >= 0.\n*\n*  V       (input) COMPLEX*16 array, dimension (1+(L-1)*abs(INCV))\n*          The vector v in the representation of H as returned by\n*          ZTZRZF. V is not used if TAU = 0.\n*\n*  INCV    (input) INTEGER\n*          The increment between elements of v. INCV <> 0.\n*\n*  TAU     (input) COMPLEX*16\n*          The value tau in the representation of H.\n*\n*  C       (input/output) COMPLEX*16 array, dimension (LDC,N)\n*          On entry, the M-by-N matrix C.\n*          On exit, C is overwritten by the matrix H * C if SIDE = 'L',\n*          or C * H if SIDE = 'R'.\n*\n*  LDC     (input) INTEGER\n*          The leading dimension of the array C. LDC >= max(1,M).\n*\n*  WORK    (workspace) COMPLEX*16 array, dimension\n*                         (N) if SIDE = 'L'\n*                      or (M) if SIDE = 'R'\n*\n\n*  Further Details\n*  ===============\n*\n*  Based on contributions by\n*    A. Petitet, Computer Science Dept., Univ. of Tenn., Knoxville, USA\n*\n*  =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  c = NumRu::Lapack.zlarz( side, m, l, v, incv, tau, c, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 7)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 7)", argc);
  rblapack_side = argv[0];
  rblapack_m = argv[1];
  rblapack_l = argv[2];
  rblapack_v = argv[3];
  rblapack_incv = argv[4];
  rblapack_tau = argv[5];
  rblapack_c = argv[6];
  if (rb_options != Qnil) {
  }

  tau.r = NUM2DBL(rb_funcall(rblapack_tau, rb_intern("real"), 0));
  tau.i = NUM2DBL(rb_funcall(rblapack_tau, rb_intern("imag"), 0));
  l = NUM2INT(rblapack_l);
  side = StringValueCStr(rblapack_side)[0];
  incv = NUM2INT(rblapack_incv);
  if (!NA_IsNArray(rblapack_c))
    rb_raise(rb_eArgError, "c (7th argument) must be NArray");
  if (NA_RANK(rblapack_c) != 2)
    rb_raise(rb_eArgError, "rank of c (7th argument) must be %d", 2);
  n = NA_SHAPE1(rblapack_c);
  ldc = NA_SHAPE0(rblapack_c);
  if (NA_TYPE(rblapack_c) != NA_DCOMPLEX)
    rblapack_c = na_change_type(rblapack_c, NA_DCOMPLEX);
  c = NA_PTR_TYPE(rblapack_c, doublecomplex*);
  m = NUM2INT(rblapack_m);
  if (!NA_IsNArray(rblapack_v))
    rb_raise(rb_eArgError, "v (4th argument) must be NArray");
  if (NA_RANK(rblapack_v) != 1)
    rb_raise(rb_eArgError, "rank of v (4th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_v) != (1+(l-1)*abs(incv)))
    rb_raise(rb_eRuntimeError, "shape 0 of v must be %d", 1+(l-1)*abs(incv));
  if (NA_TYPE(rblapack_v) != NA_DCOMPLEX)
    rblapack_v = na_change_type(rblapack_v, NA_DCOMPLEX);
  v = NA_PTR_TYPE(rblapack_v, doublecomplex*);
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
  work = ALLOC_N(doublecomplex, (lsame_(&side,"L") ? n : lsame_(&side,"R") ? m : 0));

  zlarz_(&side, &m, &n, &l, v, &incv, &tau, c, &ldc, work);

  free(work);
  return rblapack_c;
}

void
init_lapack_zlarz(VALUE mLapack){
  rb_define_module_function(mLapack, "zlarz", rblapack_zlarz, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}

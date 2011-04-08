#include "rb_lapack.h"

extern VOID dlarfx_(char *side, integer *m, integer *n, doublereal *v, doublereal *tau, doublereal *c, integer *ldc, doublereal *work);

static VALUE sHelp, sUsage;

static VALUE
rblapack_dlarfx(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_side;
  char side; 
  VALUE rblapack_v;
  doublereal *v; 
  VALUE rblapack_tau;
  doublereal tau; 
  VALUE rblapack_c;
  doublereal *c; 
  VALUE rblapack_c_out__;
  doublereal *c_out__;
  doublereal *work;

  integer m;
  integer ldc;
  integer n;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  c = NumRu::Lapack.dlarfx( side, v, tau, c, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE DLARFX( SIDE, M, N, V, TAU, C, LDC, WORK )\n\n*  Purpose\n*  =======\n*\n*  DLARFX applies a real elementary reflector H to a real m by n\n*  matrix C, from either the left or the right. H is represented in the\n*  form\n*\n*        H = I - tau * v * v'\n*\n*  where tau is a real scalar and v is a real vector.\n*\n*  If tau = 0, then H is taken to be the unit matrix\n*\n*  This version uses inline code if H has order < 11.\n*\n\n*  Arguments\n*  =========\n*\n*  SIDE    (input) CHARACTER*1\n*          = 'L': form  H * C\n*          = 'R': form  C * H\n*\n*  M       (input) INTEGER\n*          The number of rows of the matrix C.\n*\n*  N       (input) INTEGER\n*          The number of columns of the matrix C.\n*\n*  V       (input) DOUBLE PRECISION array, dimension (M) if SIDE = 'L'\n*                                     or (N) if SIDE = 'R'\n*          The vector v in the representation of H.\n*\n*  TAU     (input) DOUBLE PRECISION\n*          The value tau in the representation of H.\n*\n*  C       (input/output) DOUBLE PRECISION array, dimension (LDC,N)\n*          On entry, the m by n matrix C.\n*          On exit, C is overwritten by the matrix H * C if SIDE = 'L',\n*          or C * H if SIDE = 'R'.\n*\n*  LDC     (input) INTEGER\n*          The leading dimension of the array C. LDA >= (1,M).\n*\n*  WORK    (workspace) DOUBLE PRECISION array, dimension\n*                      (N) if SIDE = 'L'\n*                      or (M) if SIDE = 'R'\n*          WORK is not referenced if H has order < 11.\n*\n\n*  =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  c = NumRu::Lapack.dlarfx( side, v, tau, c, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 4)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 4)", argc);
  rblapack_side = argv[0];
  rblapack_v = argv[1];
  rblapack_tau = argv[2];
  rblapack_c = argv[3];
  if (rb_options != Qnil) {
  }

  if (!NA_IsNArray(rblapack_v))
    rb_raise(rb_eArgError, "v (2th argument) must be NArray");
  if (NA_RANK(rblapack_v) != 1)
    rb_raise(rb_eArgError, "rank of v (2th argument) must be %d", 1);
  m = NA_SHAPE0(rblapack_v);
  if (NA_TYPE(rblapack_v) != NA_DFLOAT)
    rblapack_v = na_change_type(rblapack_v, NA_DFLOAT);
  v = NA_PTR_TYPE(rblapack_v, doublereal*);
  side = StringValueCStr(rblapack_side)[0];
  if (!NA_IsNArray(rblapack_c))
    rb_raise(rb_eArgError, "c (4th argument) must be NArray");
  if (NA_RANK(rblapack_c) != 2)
    rb_raise(rb_eArgError, "rank of c (4th argument) must be %d", 2);
  n = NA_SHAPE1(rblapack_c);
  ldc = NA_SHAPE0(rblapack_c);
  if (NA_TYPE(rblapack_c) != NA_DFLOAT)
    rblapack_c = na_change_type(rblapack_c, NA_DFLOAT);
  c = NA_PTR_TYPE(rblapack_c, doublereal*);
  tau = NUM2DBL(rblapack_tau);
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

  dlarfx_(&side, &m, &n, v, &tau, c, &ldc, work);

  free(work);
  return rblapack_c;
}

void
init_lapack_dlarfx(VALUE mLapack){
  rb_define_module_function(mLapack, "dlarfx", rblapack_dlarfx, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}

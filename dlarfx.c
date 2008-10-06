#include "rb_lapack.h"

static VALUE
rb_dlarfx(int argc, VALUE *argv, VALUE self){
  VALUE rb_side;
  char side; 
  VALUE rb_v;
  doublereal *v; 
  VALUE rb_tau;
  doublereal tau; 
  VALUE rb_c;
  doublereal *c; 
  VALUE rb_c_out__;
  doublereal *c_out__;
  doublereal *work;

  integer m;
  integer ldc;
  integer n;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  c = NumRu::Lapack.dlarfx( side, v, tau, c)\n    or\n  NumRu::Lapack.dlarfx  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE DLARFX( SIDE, M, N, V, TAU, C, LDC, WORK )\n\n*  Purpose\n*  =======\n*\n*  DLARFX applies a real elementary reflector H to a real m by n\n*  matrix C, from either the left or the right. H is represented in the\n*  form\n*\n*        H = I - tau * v * v'\n*\n*  where tau is a real scalar and v is a real vector.\n*\n*  If tau = 0, then H is taken to be the unit matrix\n*\n*  This version uses inline code if H has order < 11.\n*\n\n*  Arguments\n*  =========\n*\n*  SIDE    (input) CHARACTER*1\n*          = 'L': form  H * C\n*          = 'R': form  C * H\n*\n*  M       (input) INTEGER\n*          The number of rows of the matrix C.\n*\n*  N       (input) INTEGER\n*          The number of columns of the matrix C.\n*\n*  V       (input) DOUBLE PRECISION array, dimension (M) if SIDE = 'L'\n*                                     or (N) if SIDE = 'R'\n*          The vector v in the representation of H.\n*\n*  TAU     (input) DOUBLE PRECISION\n*          The value tau in the representation of H.\n*\n*  C       (input/output) DOUBLE PRECISION array, dimension (LDC,N)\n*          On entry, the m by n matrix C.\n*          On exit, C is overwritten by the matrix H * C if SIDE = 'L',\n*          or C * H if SIDE = 'R'.\n*\n*  LDC     (input) INTEGER\n*          The leading dimension of the array C. LDA >= (1,M).\n*\n*  WORK    (workspace) DOUBLE PRECISION array, dimension\n*                      (N) if SIDE = 'L'\n*                      or (M) if SIDE = 'R'\n*          WORK is not referenced if H has order < 11.\n*\n\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 4)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 4)", argc);
  rb_side = argv[0];
  rb_v = argv[1];
  rb_tau = argv[2];
  rb_c = argv[3];

  side = StringValueCStr(rb_side)[0];
  tau = NUM2DBL(rb_tau);
  if (!NA_IsNArray(rb_v))
    rb_raise(rb_eArgError, "v (2th argument) must be NArray");
  if (NA_RANK(rb_v) != 1)
    rb_raise(rb_eArgError, "rank of v (2th argument) must be %d", 1);
  m = NA_SHAPE0(rb_v);
  if (NA_TYPE(rb_v) != NA_DFLOAT)
    rb_v = na_change_type(rb_v, NA_DFLOAT);
  v = NA_PTR_TYPE(rb_v, doublereal*);
  if (!NA_IsNArray(rb_c))
    rb_raise(rb_eArgError, "c (4th argument) must be NArray");
  if (NA_RANK(rb_c) != 2)
    rb_raise(rb_eArgError, "rank of c (4th argument) must be %d", 2);
  ldc = NA_SHAPE0(rb_c);
  n = NA_SHAPE1(rb_c);
  if (NA_TYPE(rb_c) != NA_DFLOAT)
    rb_c = na_change_type(rb_c, NA_DFLOAT);
  c = NA_PTR_TYPE(rb_c, doublereal*);
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
  work = ALLOC_N(doublereal, (lsame_(&side,"L") ? n : lsame_(&side,"R") ? m : 0));

  dlarfx_(&side, &m, &n, v, &tau, c, &ldc, work);

  free(work);
  return rb_c;
}

void
init_lapack_dlarfx(VALUE mLapack){
  rb_define_module_function(mLapack, "dlarfx", rb_dlarfx, -1);
}

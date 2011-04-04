#include "rb_lapack.h"

extern VOID slarf_(char *side, integer *m, integer *n, real *v, integer *incv, real *tau, real *c, integer *ldc, real *work);

static VALUE
rb_slarf(int argc, VALUE *argv, VALUE self){
  VALUE rb_side;
  char side; 
  VALUE rb_m;
  integer m; 
  VALUE rb_v;
  real *v; 
  VALUE rb_incv;
  integer incv; 
  VALUE rb_tau;
  real tau; 
  VALUE rb_c;
  real *c; 
  VALUE rb_c_out__;
  real *c_out__;
  real *work;

  integer ldc;
  integer n;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  c = NumRu::Lapack.slarf( side, m, v, incv, tau, c)\n    or\n  NumRu::Lapack.slarf  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE SLARF( SIDE, M, N, V, INCV, TAU, C, LDC, WORK )\n\n*  Purpose\n*  =======\n*\n*  SLARF applies a real elementary reflector H to a real m by n matrix\n*  C, from either the left or the right. H is represented in the form\n*\n*        H = I - tau * v * v'\n*\n*  where tau is a real scalar and v is a real vector.\n*\n*  If tau = 0, then H is taken to be the unit matrix.\n*\n\n*  Arguments\n*  =========\n*\n*  SIDE    (input) CHARACTER*1\n*          = 'L': form  H * C\n*          = 'R': form  C * H\n*\n*  M       (input) INTEGER\n*          The number of rows of the matrix C.\n*\n*  N       (input) INTEGER\n*          The number of columns of the matrix C.\n*\n*  V       (input) REAL array, dimension\n*                     (1 + (M-1)*abs(INCV)) if SIDE = 'L'\n*                  or (1 + (N-1)*abs(INCV)) if SIDE = 'R'\n*          The vector v in the representation of H. V is not used if\n*          TAU = 0.\n*\n*  INCV    (input) INTEGER\n*          The increment between elements of v. INCV <> 0.\n*\n*  TAU     (input) REAL\n*          The value tau in the representation of H.\n*\n*  C       (input/output) REAL array, dimension (LDC,N)\n*          On entry, the m by n matrix C.\n*          On exit, C is overwritten by the matrix H * C if SIDE = 'L',\n*          or C * H if SIDE = 'R'.\n*\n*  LDC     (input) INTEGER\n*          The leading dimension of the array C. LDC >= max(1,M).\n*\n*  WORK    (workspace) REAL array, dimension\n*                         (N) if SIDE = 'L'\n*                      or (M) if SIDE = 'R'\n*\n\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 6)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 6)", argc);
  rb_side = argv[0];
  rb_m = argv[1];
  rb_v = argv[2];
  rb_incv = argv[3];
  rb_tau = argv[4];
  rb_c = argv[5];

  tau = (real)NUM2DBL(rb_tau);
  side = StringValueCStr(rb_side)[0];
  m = NUM2INT(rb_m);
  incv = NUM2INT(rb_incv);
  if (!NA_IsNArray(rb_c))
    rb_raise(rb_eArgError, "c (6th argument) must be NArray");
  if (NA_RANK(rb_c) != 2)
    rb_raise(rb_eArgError, "rank of c (6th argument) must be %d", 2);
  n = NA_SHAPE1(rb_c);
  ldc = NA_SHAPE0(rb_c);
  if (NA_TYPE(rb_c) != NA_SFLOAT)
    rb_c = na_change_type(rb_c, NA_SFLOAT);
  c = NA_PTR_TYPE(rb_c, real*);
  if (!NA_IsNArray(rb_v))
    rb_raise(rb_eArgError, "v (3th argument) must be NArray");
  if (NA_RANK(rb_v) != 1)
    rb_raise(rb_eArgError, "rank of v (3th argument) must be %d", 1);
  if (NA_SHAPE0(rb_v) != (1 + (m-1)*abs(incv)))
    rb_raise(rb_eRuntimeError, "shape 0 of v must be %d", 1 + (m-1)*abs(incv));
  if (NA_TYPE(rb_v) != NA_SFLOAT)
    rb_v = na_change_type(rb_v, NA_SFLOAT);
  v = NA_PTR_TYPE(rb_v, real*);
  {
    int shape[2];
    shape[0] = ldc;
    shape[1] = n;
    rb_c_out__ = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  c_out__ = NA_PTR_TYPE(rb_c_out__, real*);
  MEMCPY(c_out__, c, real, NA_TOTAL(rb_c));
  rb_c = rb_c_out__;
  c = c_out__;
  work = ALLOC_N(real, (lsame_(&side,"L") ? n : lsame_(&side,"R") ? m : 0));

  slarf_(&side, &m, &n, v, &incv, &tau, c, &ldc, work);

  free(work);
  return rb_c;
}

void
init_lapack_slarf(VALUE mLapack){
  rb_define_module_function(mLapack, "slarf", rb_slarf, -1);
}

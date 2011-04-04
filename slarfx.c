#include "rb_lapack.h"

extern VOID slarfx_(char *side, integer *m, integer *n, real *v, real *tau, real *c, integer *ldc, real *work);

static VALUE
rb_slarfx(int argc, VALUE *argv, VALUE self){
  VALUE rb_side;
  char side; 
  VALUE rb_v;
  real *v; 
  VALUE rb_tau;
  real tau; 
  VALUE rb_c;
  real *c; 
  VALUE rb_c_out__;
  real *c_out__;
  real *work;

  integer m;
  integer ldc;
  integer n;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  c = NumRu::Lapack.slarfx( side, v, tau, c)\n    or\n  NumRu::Lapack.slarfx  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 4)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 4)", argc);
  rb_side = argv[0];
  rb_v = argv[1];
  rb_tau = argv[2];
  rb_c = argv[3];

  if (!NA_IsNArray(rb_v))
    rb_raise(rb_eArgError, "v (2th argument) must be NArray");
  if (NA_RANK(rb_v) != 1)
    rb_raise(rb_eArgError, "rank of v (2th argument) must be %d", 1);
  m = NA_SHAPE0(rb_v);
  if (NA_TYPE(rb_v) != NA_SFLOAT)
    rb_v = na_change_type(rb_v, NA_SFLOAT);
  v = NA_PTR_TYPE(rb_v, real*);
  side = StringValueCStr(rb_side)[0];
  if (!NA_IsNArray(rb_c))
    rb_raise(rb_eArgError, "c (4th argument) must be NArray");
  if (NA_RANK(rb_c) != 2)
    rb_raise(rb_eArgError, "rank of c (4th argument) must be %d", 2);
  n = NA_SHAPE1(rb_c);
  ldc = NA_SHAPE0(rb_c);
  if (NA_TYPE(rb_c) != NA_SFLOAT)
    rb_c = na_change_type(rb_c, NA_SFLOAT);
  c = NA_PTR_TYPE(rb_c, real*);
  tau = (real)NUM2DBL(rb_tau);
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

  slarfx_(&side, &m, &n, v, &tau, c, &ldc, work);

  free(work);
  return rb_c;
}

void
init_lapack_slarfx(VALUE mLapack){
  rb_define_module_function(mLapack, "slarfx", rb_slarfx, -1);
}

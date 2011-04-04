#include "rb_lapack.h"

extern VOID dlarfx_(char *side, integer *m, integer *n, doublereal *v, doublereal *tau, doublereal *c, integer *ldc, doublereal *work);

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
    printf("%s\n", "USAGE:\n  c = NumRu::Lapack.dlarfx( side, v, tau, c)\n    or\n  NumRu::Lapack.dlarfx  # print help\n\n\nFORTRAN MANUAL\n\n");
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
  if (NA_TYPE(rb_v) != NA_DFLOAT)
    rb_v = na_change_type(rb_v, NA_DFLOAT);
  v = NA_PTR_TYPE(rb_v, doublereal*);
  side = StringValueCStr(rb_side)[0];
  if (!NA_IsNArray(rb_c))
    rb_raise(rb_eArgError, "c (4th argument) must be NArray");
  if (NA_RANK(rb_c) != 2)
    rb_raise(rb_eArgError, "rank of c (4th argument) must be %d", 2);
  n = NA_SHAPE1(rb_c);
  ldc = NA_SHAPE0(rb_c);
  if (NA_TYPE(rb_c) != NA_DFLOAT)
    rb_c = na_change_type(rb_c, NA_DFLOAT);
  c = NA_PTR_TYPE(rb_c, doublereal*);
  tau = NUM2DBL(rb_tau);
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

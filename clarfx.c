#include "rb_lapack.h"

extern VOID clarfx_(char *side, integer *m, integer *n, complex *v, complex *tau, complex *c, integer *ldc, complex *work);

static VALUE
rb_clarfx(int argc, VALUE *argv, VALUE self){
  VALUE rb_side;
  char side; 
  VALUE rb_v;
  complex *v; 
  VALUE rb_tau;
  complex tau; 
  VALUE rb_c;
  complex *c; 
  VALUE rb_c_out__;
  complex *c_out__;
  complex *work;

  integer m;
  integer ldc;
  integer n;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  c = NumRu::Lapack.clarfx( side, v, tau, c)\n    or\n  NumRu::Lapack.clarfx  # print help\n\n\nFORTRAN MANUAL\n\n");
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
  if (NA_TYPE(rb_v) != NA_SCOMPLEX)
    rb_v = na_change_type(rb_v, NA_SCOMPLEX);
  v = NA_PTR_TYPE(rb_v, complex*);
  side = StringValueCStr(rb_side)[0];
  if (!NA_IsNArray(rb_c))
    rb_raise(rb_eArgError, "c (4th argument) must be NArray");
  if (NA_RANK(rb_c) != 2)
    rb_raise(rb_eArgError, "rank of c (4th argument) must be %d", 2);
  n = NA_SHAPE1(rb_c);
  ldc = NA_SHAPE0(rb_c);
  if (NA_TYPE(rb_c) != NA_SCOMPLEX)
    rb_c = na_change_type(rb_c, NA_SCOMPLEX);
  c = NA_PTR_TYPE(rb_c, complex*);
  tau.r = (real)NUM2DBL(rb_funcall(rb_tau, rb_intern("real"), 0));
  tau.i = (real)NUM2DBL(rb_funcall(rb_tau, rb_intern("imag"), 0));
  {
    int shape[2];
    shape[0] = ldc;
    shape[1] = n;
    rb_c_out__ = na_make_object(NA_SCOMPLEX, 2, shape, cNArray);
  }
  c_out__ = NA_PTR_TYPE(rb_c_out__, complex*);
  MEMCPY(c_out__, c, complex, NA_TOTAL(rb_c));
  rb_c = rb_c_out__;
  c = c_out__;
  work = ALLOC_N(complex, (lsame_(&side,"L") ? n : lsame_(&side,"R") ? m : 0));

  clarfx_(&side, &m, &n, v, &tau, c, &ldc, work);

  free(work);
  return rb_c;
}

void
init_lapack_clarfx(VALUE mLapack){
  rb_define_module_function(mLapack, "clarfx", rb_clarfx, -1);
}

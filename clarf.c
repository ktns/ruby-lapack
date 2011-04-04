#include "rb_lapack.h"

extern VOID clarf_(char *side, integer *m, integer *n, complex *v, integer *incv, complex *tau, complex *c, integer *ldc, complex *work);

static VALUE
rb_clarf(int argc, VALUE *argv, VALUE self){
  VALUE rb_side;
  char side; 
  VALUE rb_m;
  integer m; 
  VALUE rb_v;
  complex *v; 
  VALUE rb_incv;
  integer incv; 
  VALUE rb_tau;
  complex tau; 
  VALUE rb_c;
  complex *c; 
  VALUE rb_c_out__;
  complex *c_out__;
  complex *work;

  integer ldc;
  integer n;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  c = NumRu::Lapack.clarf( side, m, v, incv, tau, c)\n    or\n  NumRu::Lapack.clarf  # print help\n\n\nFORTRAN MANUAL\n\n");
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

  tau.r = (real)NUM2DBL(rb_funcall(rb_tau, rb_intern("real"), 0));
  tau.i = (real)NUM2DBL(rb_funcall(rb_tau, rb_intern("imag"), 0));
  side = StringValueCStr(rb_side)[0];
  m = NUM2INT(rb_m);
  incv = NUM2INT(rb_incv);
  if (!NA_IsNArray(rb_c))
    rb_raise(rb_eArgError, "c (6th argument) must be NArray");
  if (NA_RANK(rb_c) != 2)
    rb_raise(rb_eArgError, "rank of c (6th argument) must be %d", 2);
  n = NA_SHAPE1(rb_c);
  ldc = NA_SHAPE0(rb_c);
  if (NA_TYPE(rb_c) != NA_SCOMPLEX)
    rb_c = na_change_type(rb_c, NA_SCOMPLEX);
  c = NA_PTR_TYPE(rb_c, complex*);
  if (!NA_IsNArray(rb_v))
    rb_raise(rb_eArgError, "v (3th argument) must be NArray");
  if (NA_RANK(rb_v) != 1)
    rb_raise(rb_eArgError, "rank of v (3th argument) must be %d", 1);
  if (NA_SHAPE0(rb_v) != (1 + (m-1)*abs(incv)))
    rb_raise(rb_eRuntimeError, "shape 0 of v must be %d", 1 + (m-1)*abs(incv));
  if (NA_TYPE(rb_v) != NA_SCOMPLEX)
    rb_v = na_change_type(rb_v, NA_SCOMPLEX);
  v = NA_PTR_TYPE(rb_v, complex*);
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

  clarf_(&side, &m, &n, v, &incv, &tau, c, &ldc, work);

  free(work);
  return rb_c;
}

void
init_lapack_clarf(VALUE mLapack){
  rb_define_module_function(mLapack, "clarf", rb_clarf, -1);
}

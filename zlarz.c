#include "rb_lapack.h"

extern VOID zlarz_(char *side, integer *m, integer *n, integer *l, doublecomplex *v, integer *incv, doublecomplex *tau, doublecomplex *c, integer *ldc, doublecomplex *work);

static VALUE
rb_zlarz(int argc, VALUE *argv, VALUE self){
  VALUE rb_side;
  char side; 
  VALUE rb_m;
  integer m; 
  VALUE rb_l;
  integer l; 
  VALUE rb_v;
  doublecomplex *v; 
  VALUE rb_incv;
  integer incv; 
  VALUE rb_tau;
  doublecomplex tau; 
  VALUE rb_c;
  doublecomplex *c; 
  VALUE rb_c_out__;
  doublecomplex *c_out__;
  doublecomplex *work;

  integer ldc;
  integer n;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  c = NumRu::Lapack.zlarz( side, m, l, v, incv, tau, c)\n    or\n  NumRu::Lapack.zlarz  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 7)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 7)", argc);
  rb_side = argv[0];
  rb_m = argv[1];
  rb_l = argv[2];
  rb_v = argv[3];
  rb_incv = argv[4];
  rb_tau = argv[5];
  rb_c = argv[6];

  tau.r = NUM2DBL(rb_funcall(rb_tau, rb_intern("real"), 0));
  tau.i = NUM2DBL(rb_funcall(rb_tau, rb_intern("imag"), 0));
  l = NUM2INT(rb_l);
  side = StringValueCStr(rb_side)[0];
  incv = NUM2INT(rb_incv);
  if (!NA_IsNArray(rb_c))
    rb_raise(rb_eArgError, "c (7th argument) must be NArray");
  if (NA_RANK(rb_c) != 2)
    rb_raise(rb_eArgError, "rank of c (7th argument) must be %d", 2);
  n = NA_SHAPE1(rb_c);
  ldc = NA_SHAPE0(rb_c);
  if (NA_TYPE(rb_c) != NA_DCOMPLEX)
    rb_c = na_change_type(rb_c, NA_DCOMPLEX);
  c = NA_PTR_TYPE(rb_c, doublecomplex*);
  m = NUM2INT(rb_m);
  if (!NA_IsNArray(rb_v))
    rb_raise(rb_eArgError, "v (4th argument) must be NArray");
  if (NA_RANK(rb_v) != 1)
    rb_raise(rb_eArgError, "rank of v (4th argument) must be %d", 1);
  if (NA_SHAPE0(rb_v) != (1+(l-1)*abs(incv)))
    rb_raise(rb_eRuntimeError, "shape 0 of v must be %d", 1+(l-1)*abs(incv));
  if (NA_TYPE(rb_v) != NA_DCOMPLEX)
    rb_v = na_change_type(rb_v, NA_DCOMPLEX);
  v = NA_PTR_TYPE(rb_v, doublecomplex*);
  {
    int shape[2];
    shape[0] = ldc;
    shape[1] = n;
    rb_c_out__ = na_make_object(NA_DCOMPLEX, 2, shape, cNArray);
  }
  c_out__ = NA_PTR_TYPE(rb_c_out__, doublecomplex*);
  MEMCPY(c_out__, c, doublecomplex, NA_TOTAL(rb_c));
  rb_c = rb_c_out__;
  c = c_out__;
  work = ALLOC_N(doublecomplex, (lsame_(&side,"L") ? n : lsame_(&side,"R") ? m : 0));

  zlarz_(&side, &m, &n, &l, v, &incv, &tau, c, &ldc, work);

  free(work);
  return rb_c;
}

void
init_lapack_zlarz(VALUE mLapack){
  rb_define_module_function(mLapack, "zlarz", rb_zlarz, -1);
}

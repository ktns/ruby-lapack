#include "rb_lapack.h"

extern VOID slarzt_(char *direct, char *storev, integer *n, integer *k, real *v, integer *ldv, real *tau, real *t, integer *ldt);

static VALUE
rb_slarzt(int argc, VALUE *argv, VALUE self){
  VALUE rb_direct;
  char direct; 
  VALUE rb_storev;
  char storev; 
  VALUE rb_n;
  integer n; 
  VALUE rb_v;
  real *v; 
  VALUE rb_tau;
  real *tau; 
  VALUE rb_t;
  real *t; 
  VALUE rb_v_out__;
  real *v_out__;

  integer ldv;
  integer k;
  integer ldt;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  t, v = NumRu::Lapack.slarzt( direct, storev, n, v, tau)\n    or\n  NumRu::Lapack.slarzt  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 5)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 5)", argc);
  rb_direct = argv[0];
  rb_storev = argv[1];
  rb_n = argv[2];
  rb_v = argv[3];
  rb_tau = argv[4];

  storev = StringValueCStr(rb_storev)[0];
  direct = StringValueCStr(rb_direct)[0];
  n = NUM2INT(rb_n);
  if (!NA_IsNArray(rb_tau))
    rb_raise(rb_eArgError, "tau (5th argument) must be NArray");
  if (NA_RANK(rb_tau) != 1)
    rb_raise(rb_eArgError, "rank of tau (5th argument) must be %d", 1);
  k = NA_SHAPE0(rb_tau);
  if (NA_TYPE(rb_tau) != NA_SFLOAT)
    rb_tau = na_change_type(rb_tau, NA_SFLOAT);
  tau = NA_PTR_TYPE(rb_tau, real*);
  ldt = k;
  if (!NA_IsNArray(rb_v))
    rb_raise(rb_eArgError, "v (4th argument) must be NArray");
  if (NA_RANK(rb_v) != 2)
    rb_raise(rb_eArgError, "rank of v (4th argument) must be %d", 2);
  if (NA_SHAPE1(rb_v) != (lsame_(&storev,"C") ? k : lsame_(&storev,"R") ? n : 0))
    rb_raise(rb_eRuntimeError, "shape 1 of v must be %d", lsame_(&storev,"C") ? k : lsame_(&storev,"R") ? n : 0);
  ldv = NA_SHAPE0(rb_v);
  if (NA_TYPE(rb_v) != NA_SFLOAT)
    rb_v = na_change_type(rb_v, NA_SFLOAT);
  v = NA_PTR_TYPE(rb_v, real*);
  {
    int shape[2];
    shape[0] = ldt;
    shape[1] = k;
    rb_t = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  t = NA_PTR_TYPE(rb_t, real*);
  {
    int shape[2];
    shape[0] = ldv;
    shape[1] = lsame_(&storev,"C") ? k : lsame_(&storev,"R") ? n : 0;
    rb_v_out__ = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  v_out__ = NA_PTR_TYPE(rb_v_out__, real*);
  MEMCPY(v_out__, v, real, NA_TOTAL(rb_v));
  rb_v = rb_v_out__;
  v = v_out__;

  slarzt_(&direct, &storev, &n, &k, v, &ldv, tau, t, &ldt);

  return rb_ary_new3(2, rb_t, rb_v);
}

void
init_lapack_slarzt(VALUE mLapack){
  rb_define_module_function(mLapack, "slarzt", rb_slarzt, -1);
}

#include "rb_lapack.h"

extern VOID slarfb_(char *side, char *trans, char *direct, char *storev, integer *m, integer *n, integer *k, real *v, integer *ldv, real *t, integer *ldt, real *c, integer *ldc, real *work, integer *ldwork);

static VALUE
rb_slarfb(int argc, VALUE *argv, VALUE self){
  VALUE rb_side;
  char side; 
  VALUE rb_trans;
  char trans; 
  VALUE rb_direct;
  char direct; 
  VALUE rb_storev;
  char storev; 
  VALUE rb_m;
  integer m; 
  VALUE rb_v;
  real *v; 
  VALUE rb_t;
  real *t; 
  VALUE rb_c;
  real *c; 
  VALUE rb_c_out__;
  real *c_out__;
  real *work;

  integer ldv;
  integer k;
  integer ldt;
  integer ldc;
  integer n;
  integer ldwork;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  c = NumRu::Lapack.slarfb( side, trans, direct, storev, m, v, t, c)\n    or\n  NumRu::Lapack.slarfb  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 8)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 8)", argc);
  rb_side = argv[0];
  rb_trans = argv[1];
  rb_direct = argv[2];
  rb_storev = argv[3];
  rb_m = argv[4];
  rb_v = argv[5];
  rb_t = argv[6];
  rb_c = argv[7];

  trans = StringValueCStr(rb_trans)[0];
  if (!NA_IsNArray(rb_v))
    rb_raise(rb_eArgError, "v (6th argument) must be NArray");
  if (NA_RANK(rb_v) != 2)
    rb_raise(rb_eArgError, "rank of v (6th argument) must be %d", 2);
  k = NA_SHAPE1(rb_v);
  ldv = NA_SHAPE0(rb_v);
  if (NA_TYPE(rb_v) != NA_SFLOAT)
    rb_v = na_change_type(rb_v, NA_SFLOAT);
  v = NA_PTR_TYPE(rb_v, real*);
  direct = StringValueCStr(rb_direct)[0];
  side = StringValueCStr(rb_side)[0];
  storev = StringValueCStr(rb_storev)[0];
  if (!NA_IsNArray(rb_c))
    rb_raise(rb_eArgError, "c (8th argument) must be NArray");
  if (NA_RANK(rb_c) != 2)
    rb_raise(rb_eArgError, "rank of c (8th argument) must be %d", 2);
  n = NA_SHAPE1(rb_c);
  ldc = NA_SHAPE0(rb_c);
  if (NA_TYPE(rb_c) != NA_SFLOAT)
    rb_c = na_change_type(rb_c, NA_SFLOAT);
  c = NA_PTR_TYPE(rb_c, real*);
  if (!NA_IsNArray(rb_t))
    rb_raise(rb_eArgError, "t (7th argument) must be NArray");
  if (NA_RANK(rb_t) != 2)
    rb_raise(rb_eArgError, "rank of t (7th argument) must be %d", 2);
  if (NA_SHAPE1(rb_t) != k)
    rb_raise(rb_eRuntimeError, "shape 1 of t must be the same as shape 1 of v");
  ldt = NA_SHAPE0(rb_t);
  if (NA_TYPE(rb_t) != NA_SFLOAT)
    rb_t = na_change_type(rb_t, NA_SFLOAT);
  t = NA_PTR_TYPE(rb_t, real*);
  m = NUM2INT(rb_m);
  ldwork = max(1,n) ? side = 'l' : max(1,m) ? side = 'r' : 0;
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
  work = ALLOC_N(real, (ldwork)*(k));

  slarfb_(&side, &trans, &direct, &storev, &m, &n, &k, v, &ldv, t, &ldt, c, &ldc, work, &ldwork);

  free(work);
  return rb_c;
}

void
init_lapack_slarfb(VALUE mLapack){
  rb_define_module_function(mLapack, "slarfb", rb_slarfb, -1);
}

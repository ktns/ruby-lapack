#include "rb_lapack.h"

extern VOID slarzb_(char *side, char *trans, char *direct, char *storev, integer *m, integer *n, integer *k, integer *l, real *v, integer *ldv, real *t, integer *ldt, real *c, integer *ldc, real *work, integer *ldwork);

static VALUE
rb_slarzb(int argc, VALUE *argv, VALUE self){
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
  VALUE rb_l;
  integer l; 
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
  integer nv;
  integer ldt;
  integer k;
  integer ldc;
  integer n;
  integer ldwork;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  c = NumRu::Lapack.slarzb( side, trans, direct, storev, m, l, v, t, c)\n    or\n  NumRu::Lapack.slarzb  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 9)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 9)", argc);
  rb_side = argv[0];
  rb_trans = argv[1];
  rb_direct = argv[2];
  rb_storev = argv[3];
  rb_m = argv[4];
  rb_l = argv[5];
  rb_v = argv[6];
  rb_t = argv[7];
  rb_c = argv[8];

  trans = StringValueCStr(rb_trans)[0];
  if (!NA_IsNArray(rb_v))
    rb_raise(rb_eArgError, "v (7th argument) must be NArray");
  if (NA_RANK(rb_v) != 2)
    rb_raise(rb_eArgError, "rank of v (7th argument) must be %d", 2);
  nv = NA_SHAPE1(rb_v);
  ldv = NA_SHAPE0(rb_v);
  if (NA_TYPE(rb_v) != NA_SFLOAT)
    rb_v = na_change_type(rb_v, NA_SFLOAT);
  v = NA_PTR_TYPE(rb_v, real*);
  direct = StringValueCStr(rb_direct)[0];
  l = NUM2INT(rb_l);
  side = StringValueCStr(rb_side)[0];
  storev = StringValueCStr(rb_storev)[0];
  if (!NA_IsNArray(rb_c))
    rb_raise(rb_eArgError, "c (9th argument) must be NArray");
  if (NA_RANK(rb_c) != 2)
    rb_raise(rb_eArgError, "rank of c (9th argument) must be %d", 2);
  n = NA_SHAPE1(rb_c);
  ldc = NA_SHAPE0(rb_c);
  if (NA_TYPE(rb_c) != NA_SFLOAT)
    rb_c = na_change_type(rb_c, NA_SFLOAT);
  c = NA_PTR_TYPE(rb_c, real*);
  if (!NA_IsNArray(rb_t))
    rb_raise(rb_eArgError, "t (8th argument) must be NArray");
  if (NA_RANK(rb_t) != 2)
    rb_raise(rb_eArgError, "rank of t (8th argument) must be %d", 2);
  k = NA_SHAPE1(rb_t);
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

  slarzb_(&side, &trans, &direct, &storev, &m, &n, &k, &l, v, &ldv, t, &ldt, c, &ldc, work, &ldwork);

  free(work);
  return rb_c;
}

void
init_lapack_slarzb(VALUE mLapack){
  rb_define_module_function(mLapack, "slarzb", rb_slarzb, -1);
}

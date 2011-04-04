#include "rb_lapack.h"

extern VOID slaed9_(integer *k, integer *kstart, integer *kstop, integer *n, real *d, real *q, integer *ldq, real *rho, real *dlamda, real *w, real *s, integer *lds, integer *info);

static VALUE
rb_slaed9(int argc, VALUE *argv, VALUE self){
  VALUE rb_kstart;
  integer kstart; 
  VALUE rb_kstop;
  integer kstop; 
  VALUE rb_n;
  integer n; 
  VALUE rb_rho;
  real rho; 
  VALUE rb_dlamda;
  real *dlamda; 
  VALUE rb_w;
  real *w; 
  VALUE rb_d;
  real *d; 
  VALUE rb_s;
  real *s; 
  VALUE rb_info;
  integer info; 
  real *q;

  integer k;
  integer lds;
  integer ldq;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  d, s, info = NumRu::Lapack.slaed9( kstart, kstop, n, rho, dlamda, w)\n    or\n  NumRu::Lapack.slaed9  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 6)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 6)", argc);
  rb_kstart = argv[0];
  rb_kstop = argv[1];
  rb_n = argv[2];
  rb_rho = argv[3];
  rb_dlamda = argv[4];
  rb_w = argv[5];

  if (!NA_IsNArray(rb_w))
    rb_raise(rb_eArgError, "w (6th argument) must be NArray");
  if (NA_RANK(rb_w) != 1)
    rb_raise(rb_eArgError, "rank of w (6th argument) must be %d", 1);
  k = NA_SHAPE0(rb_w);
  if (NA_TYPE(rb_w) != NA_SFLOAT)
    rb_w = na_change_type(rb_w, NA_SFLOAT);
  w = NA_PTR_TYPE(rb_w, real*);
  kstart = NUM2INT(rb_kstart);
  if (!NA_IsNArray(rb_dlamda))
    rb_raise(rb_eArgError, "dlamda (5th argument) must be NArray");
  if (NA_RANK(rb_dlamda) != 1)
    rb_raise(rb_eArgError, "rank of dlamda (5th argument) must be %d", 1);
  if (NA_SHAPE0(rb_dlamda) != k)
    rb_raise(rb_eRuntimeError, "shape 0 of dlamda must be the same as shape 0 of w");
  if (NA_TYPE(rb_dlamda) != NA_SFLOAT)
    rb_dlamda = na_change_type(rb_dlamda, NA_SFLOAT);
  dlamda = NA_PTR_TYPE(rb_dlamda, real*);
  rho = (real)NUM2DBL(rb_rho);
  kstop = NUM2INT(rb_kstop);
  n = NUM2INT(rb_n);
  lds = MAX( 1, k );
  ldq = MAX( 1, n );
  {
    int shape[1];
    shape[0] = MAX(1,n);
    rb_d = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  d = NA_PTR_TYPE(rb_d, real*);
  {
    int shape[2];
    shape[0] = lds;
    shape[1] = k;
    rb_s = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  s = NA_PTR_TYPE(rb_s, real*);
  q = ALLOC_N(real, (ldq)*(MAX(1,n)));

  slaed9_(&k, &kstart, &kstop, &n, d, q, &ldq, &rho, dlamda, w, s, &lds, &info);

  free(q);
  rb_info = INT2NUM(info);
  return rb_ary_new3(3, rb_d, rb_s, rb_info);
}

void
init_lapack_slaed9(VALUE mLapack){
  rb_define_module_function(mLapack, "slaed9", rb_slaed9, -1);
}

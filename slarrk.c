#include "rb_lapack.h"

extern VOID slarrk_(integer *n, integer *iw, real *gl, real *gu, real *d, real *e2, real *pivmin, real *reltol, real *w, real *werr, integer *info);

static VALUE
rb_slarrk(int argc, VALUE *argv, VALUE self){
  VALUE rb_iw;
  integer iw; 
  VALUE rb_gl;
  real gl; 
  VALUE rb_gu;
  real gu; 
  VALUE rb_d;
  real *d; 
  VALUE rb_e2;
  real *e2; 
  VALUE rb_pivmin;
  real pivmin; 
  VALUE rb_reltol;
  real reltol; 
  VALUE rb_w;
  real w; 
  VALUE rb_werr;
  real werr; 
  VALUE rb_info;
  integer info; 

  integer n;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  w, werr, info = NumRu::Lapack.slarrk( iw, gl, gu, d, e2, pivmin, reltol)\n    or\n  NumRu::Lapack.slarrk  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 7)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 7)", argc);
  rb_iw = argv[0];
  rb_gl = argv[1];
  rb_gu = argv[2];
  rb_d = argv[3];
  rb_e2 = argv[4];
  rb_pivmin = argv[5];
  rb_reltol = argv[6];

  pivmin = (real)NUM2DBL(rb_pivmin);
  gu = (real)NUM2DBL(rb_gu);
  iw = NUM2INT(rb_iw);
  gl = (real)NUM2DBL(rb_gl);
  if (!NA_IsNArray(rb_d))
    rb_raise(rb_eArgError, "d (4th argument) must be NArray");
  if (NA_RANK(rb_d) != 1)
    rb_raise(rb_eArgError, "rank of d (4th argument) must be %d", 1);
  n = NA_SHAPE0(rb_d);
  if (NA_TYPE(rb_d) != NA_SFLOAT)
    rb_d = na_change_type(rb_d, NA_SFLOAT);
  d = NA_PTR_TYPE(rb_d, real*);
  reltol = (real)NUM2DBL(rb_reltol);
  if (!NA_IsNArray(rb_e2))
    rb_raise(rb_eArgError, "e2 (5th argument) must be NArray");
  if (NA_RANK(rb_e2) != 1)
    rb_raise(rb_eArgError, "rank of e2 (5th argument) must be %d", 1);
  if (NA_SHAPE0(rb_e2) != (n-1))
    rb_raise(rb_eRuntimeError, "shape 0 of e2 must be %d", n-1);
  if (NA_TYPE(rb_e2) != NA_SFLOAT)
    rb_e2 = na_change_type(rb_e2, NA_SFLOAT);
  e2 = NA_PTR_TYPE(rb_e2, real*);

  slarrk_(&n, &iw, &gl, &gu, d, e2, &pivmin, &reltol, &w, &werr, &info);

  rb_w = rb_float_new((double)w);
  rb_werr = rb_float_new((double)werr);
  rb_info = INT2NUM(info);
  return rb_ary_new3(3, rb_w, rb_werr, rb_info);
}

void
init_lapack_slarrk(VALUE mLapack){
  rb_define_module_function(mLapack, "slarrk", rb_slarrk, -1);
}

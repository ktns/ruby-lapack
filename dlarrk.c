#include "rb_lapack.h"

extern VOID dlarrk_(integer *n, integer *iw, doublereal *gl, doublereal *gu, doublereal *d, doublereal *e2, doublereal *pivmin, doublereal *reltol, doublereal *w, doublereal *werr, integer *info);

static VALUE
rb_dlarrk(int argc, VALUE *argv, VALUE self){
  VALUE rb_iw;
  integer iw; 
  VALUE rb_gl;
  doublereal gl; 
  VALUE rb_gu;
  doublereal gu; 
  VALUE rb_d;
  doublereal *d; 
  VALUE rb_e2;
  doublereal *e2; 
  VALUE rb_pivmin;
  doublereal pivmin; 
  VALUE rb_reltol;
  doublereal reltol; 
  VALUE rb_w;
  doublereal w; 
  VALUE rb_werr;
  doublereal werr; 
  VALUE rb_info;
  integer info; 

  integer n;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  w, werr, info = NumRu::Lapack.dlarrk( iw, gl, gu, d, e2, pivmin, reltol)\n    or\n  NumRu::Lapack.dlarrk  # print help\n\n\nFORTRAN MANUAL\n\n");
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

  pivmin = NUM2DBL(rb_pivmin);
  gu = NUM2DBL(rb_gu);
  iw = NUM2INT(rb_iw);
  gl = NUM2DBL(rb_gl);
  if (!NA_IsNArray(rb_d))
    rb_raise(rb_eArgError, "d (4th argument) must be NArray");
  if (NA_RANK(rb_d) != 1)
    rb_raise(rb_eArgError, "rank of d (4th argument) must be %d", 1);
  n = NA_SHAPE0(rb_d);
  if (NA_TYPE(rb_d) != NA_DFLOAT)
    rb_d = na_change_type(rb_d, NA_DFLOAT);
  d = NA_PTR_TYPE(rb_d, doublereal*);
  reltol = NUM2DBL(rb_reltol);
  if (!NA_IsNArray(rb_e2))
    rb_raise(rb_eArgError, "e2 (5th argument) must be NArray");
  if (NA_RANK(rb_e2) != 1)
    rb_raise(rb_eArgError, "rank of e2 (5th argument) must be %d", 1);
  if (NA_SHAPE0(rb_e2) != (n-1))
    rb_raise(rb_eRuntimeError, "shape 0 of e2 must be %d", n-1);
  if (NA_TYPE(rb_e2) != NA_DFLOAT)
    rb_e2 = na_change_type(rb_e2, NA_DFLOAT);
  e2 = NA_PTR_TYPE(rb_e2, doublereal*);

  dlarrk_(&n, &iw, &gl, &gu, d, e2, &pivmin, &reltol, &w, &werr, &info);

  rb_w = rb_float_new((double)w);
  rb_werr = rb_float_new((double)werr);
  rb_info = INT2NUM(info);
  return rb_ary_new3(3, rb_w, rb_werr, rb_info);
}

void
init_lapack_dlarrk(VALUE mLapack){
  rb_define_module_function(mLapack, "dlarrk", rb_dlarrk, -1);
}

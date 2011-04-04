#include "rb_lapack.h"

extern VOID slarre_(char *range, integer *n, real *vl, real *vu, integer *il, integer *iu, real *d, real *e, real *e2, real *rtol1, real *rtol2, real *spltol, integer *nsplit, integer *isplit, integer *m, real *w, real *werr, real *wgap, integer *iblock, integer *indexw, real *gers, real *pivmin, real *work, integer *iwork, integer *info);

static VALUE
rb_slarre(int argc, VALUE *argv, VALUE self){
  VALUE rb_range;
  char range; 
  VALUE rb_vl;
  real vl; 
  VALUE rb_vu;
  real vu; 
  VALUE rb_il;
  integer il; 
  VALUE rb_iu;
  integer iu; 
  VALUE rb_d;
  real *d; 
  VALUE rb_e;
  real *e; 
  VALUE rb_e2;
  real *e2; 
  VALUE rb_rtol1;
  real rtol1; 
  VALUE rb_rtol2;
  real rtol2; 
  VALUE rb_spltol;
  real spltol; 
  VALUE rb_nsplit;
  integer nsplit; 
  VALUE rb_isplit;
  integer *isplit; 
  VALUE rb_m;
  integer m; 
  VALUE rb_w;
  real *w; 
  VALUE rb_werr;
  real *werr; 
  VALUE rb_wgap;
  real *wgap; 
  VALUE rb_iblock;
  integer *iblock; 
  VALUE rb_indexw;
  integer *indexw; 
  VALUE rb_gers;
  real *gers; 
  VALUE rb_pivmin;
  real pivmin; 
  VALUE rb_info;
  integer info; 
  VALUE rb_d_out__;
  real *d_out__;
  VALUE rb_e_out__;
  real *e_out__;
  VALUE rb_e2_out__;
  real *e2_out__;
  real *work;
  integer *iwork;

  integer n;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  nsplit, isplit, m, w, werr, wgap, iblock, indexw, gers, pivmin, info, vl, vu, d, e, e2 = NumRu::Lapack.slarre( range, vl, vu, il, iu, d, e, e2, rtol1, rtol2, spltol)\n    or\n  NumRu::Lapack.slarre  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 11)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 11)", argc);
  rb_range = argv[0];
  rb_vl = argv[1];
  rb_vu = argv[2];
  rb_il = argv[3];
  rb_iu = argv[4];
  rb_d = argv[5];
  rb_e = argv[6];
  rb_e2 = argv[7];
  rb_rtol1 = argv[8];
  rb_rtol2 = argv[9];
  rb_spltol = argv[10];

  vl = (real)NUM2DBL(rb_vl);
  iu = NUM2INT(rb_iu);
  if (!NA_IsNArray(rb_e2))
    rb_raise(rb_eArgError, "e2 (8th argument) must be NArray");
  if (NA_RANK(rb_e2) != 1)
    rb_raise(rb_eArgError, "rank of e2 (8th argument) must be %d", 1);
  n = NA_SHAPE0(rb_e2);
  if (NA_TYPE(rb_e2) != NA_SFLOAT)
    rb_e2 = na_change_type(rb_e2, NA_SFLOAT);
  e2 = NA_PTR_TYPE(rb_e2, real*);
  il = NUM2INT(rb_il);
  spltol = (real)NUM2DBL(rb_spltol);
  if (!NA_IsNArray(rb_d))
    rb_raise(rb_eArgError, "d (6th argument) must be NArray");
  if (NA_RANK(rb_d) != 1)
    rb_raise(rb_eArgError, "rank of d (6th argument) must be %d", 1);
  if (NA_SHAPE0(rb_d) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of d must be the same as shape 0 of e2");
  if (NA_TYPE(rb_d) != NA_SFLOAT)
    rb_d = na_change_type(rb_d, NA_SFLOAT);
  d = NA_PTR_TYPE(rb_d, real*);
  if (!NA_IsNArray(rb_e))
    rb_raise(rb_eArgError, "e (7th argument) must be NArray");
  if (NA_RANK(rb_e) != 1)
    rb_raise(rb_eArgError, "rank of e (7th argument) must be %d", 1);
  if (NA_SHAPE0(rb_e) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of e must be the same as shape 0 of e2");
  if (NA_TYPE(rb_e) != NA_SFLOAT)
    rb_e = na_change_type(rb_e, NA_SFLOAT);
  e = NA_PTR_TYPE(rb_e, real*);
  range = StringValueCStr(rb_range)[0];
  rtol1 = (real)NUM2DBL(rb_rtol1);
  vu = (real)NUM2DBL(rb_vu);
  rtol2 = (real)NUM2DBL(rb_rtol2);
  {
    int shape[1];
    shape[0] = n;
    rb_isplit = na_make_object(NA_LINT, 1, shape, cNArray);
  }
  isplit = NA_PTR_TYPE(rb_isplit, integer*);
  {
    int shape[1];
    shape[0] = n;
    rb_w = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  w = NA_PTR_TYPE(rb_w, real*);
  {
    int shape[1];
    shape[0] = n;
    rb_werr = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  werr = NA_PTR_TYPE(rb_werr, real*);
  {
    int shape[1];
    shape[0] = n;
    rb_wgap = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  wgap = NA_PTR_TYPE(rb_wgap, real*);
  {
    int shape[1];
    shape[0] = n;
    rb_iblock = na_make_object(NA_LINT, 1, shape, cNArray);
  }
  iblock = NA_PTR_TYPE(rb_iblock, integer*);
  {
    int shape[1];
    shape[0] = n;
    rb_indexw = na_make_object(NA_LINT, 1, shape, cNArray);
  }
  indexw = NA_PTR_TYPE(rb_indexw, integer*);
  {
    int shape[1];
    shape[0] = 2*n;
    rb_gers = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  gers = NA_PTR_TYPE(rb_gers, real*);
  {
    int shape[1];
    shape[0] = n;
    rb_d_out__ = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  d_out__ = NA_PTR_TYPE(rb_d_out__, real*);
  MEMCPY(d_out__, d, real, NA_TOTAL(rb_d));
  rb_d = rb_d_out__;
  d = d_out__;
  {
    int shape[1];
    shape[0] = n;
    rb_e_out__ = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  e_out__ = NA_PTR_TYPE(rb_e_out__, real*);
  MEMCPY(e_out__, e, real, NA_TOTAL(rb_e));
  rb_e = rb_e_out__;
  e = e_out__;
  {
    int shape[1];
    shape[0] = n;
    rb_e2_out__ = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  e2_out__ = NA_PTR_TYPE(rb_e2_out__, real*);
  MEMCPY(e2_out__, e2, real, NA_TOTAL(rb_e2));
  rb_e2 = rb_e2_out__;
  e2 = e2_out__;
  work = ALLOC_N(real, (6*n));
  iwork = ALLOC_N(integer, (5*n));

  slarre_(&range, &n, &vl, &vu, &il, &iu, d, e, e2, &rtol1, &rtol2, &spltol, &nsplit, isplit, &m, w, werr, wgap, iblock, indexw, gers, &pivmin, work, iwork, &info);

  free(work);
  free(iwork);
  rb_nsplit = INT2NUM(nsplit);
  rb_m = INT2NUM(m);
  rb_pivmin = rb_float_new((double)pivmin);
  rb_info = INT2NUM(info);
  rb_vl = rb_float_new((double)vl);
  rb_vu = rb_float_new((double)vu);
  return rb_ary_new3(16, rb_nsplit, rb_isplit, rb_m, rb_w, rb_werr, rb_wgap, rb_iblock, rb_indexw, rb_gers, rb_pivmin, rb_info, rb_vl, rb_vu, rb_d, rb_e, rb_e2);
}

void
init_lapack_slarre(VALUE mLapack){
  rb_define_module_function(mLapack, "slarre", rb_slarre, -1);
}

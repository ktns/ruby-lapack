#include "rb_lapack.h"

extern VOID dlarre_(char *range, integer *n, doublereal *vl, doublereal *vu, integer *il, integer *iu, doublereal *d, doublereal *e, doublereal *e2, doublereal *rtol1, doublereal *rtol2, doublereal *spltol, integer *nsplit, integer *isplit, integer *m, doublereal *w, doublereal *werr, doublereal *wgap, integer *iblock, integer *indexw, doublereal *gers, doublereal *pivmin, doublereal *work, integer *iwork, integer *info);

static VALUE
rb_dlarre(int argc, VALUE *argv, VALUE self){
  VALUE rb_range;
  char range; 
  VALUE rb_vl;
  doublereal vl; 
  VALUE rb_vu;
  doublereal vu; 
  VALUE rb_il;
  integer il; 
  VALUE rb_iu;
  integer iu; 
  VALUE rb_d;
  doublereal *d; 
  VALUE rb_e;
  doublereal *e; 
  VALUE rb_e2;
  doublereal *e2; 
  VALUE rb_rtol1;
  doublereal rtol1; 
  VALUE rb_rtol2;
  doublereal rtol2; 
  VALUE rb_spltol;
  doublereal spltol; 
  VALUE rb_nsplit;
  integer nsplit; 
  VALUE rb_isplit;
  integer *isplit; 
  VALUE rb_m;
  integer m; 
  VALUE rb_w;
  doublereal *w; 
  VALUE rb_werr;
  doublereal *werr; 
  VALUE rb_wgap;
  doublereal *wgap; 
  VALUE rb_iblock;
  integer *iblock; 
  VALUE rb_indexw;
  integer *indexw; 
  VALUE rb_gers;
  doublereal *gers; 
  VALUE rb_pivmin;
  doublereal pivmin; 
  VALUE rb_info;
  integer info; 
  VALUE rb_d_out__;
  doublereal *d_out__;
  VALUE rb_e_out__;
  doublereal *e_out__;
  VALUE rb_e2_out__;
  doublereal *e2_out__;
  doublereal *work;
  integer *iwork;

  integer n;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  nsplit, isplit, m, w, werr, wgap, iblock, indexw, gers, pivmin, info, vl, vu, d, e, e2 = NumRu::Lapack.dlarre( range, vl, vu, il, iu, d, e, e2, rtol1, rtol2, spltol)\n    or\n  NumRu::Lapack.dlarre  # print help\n\n\nFORTRAN MANUAL\n\n");
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

  vl = NUM2DBL(rb_vl);
  iu = NUM2INT(rb_iu);
  if (!NA_IsNArray(rb_e2))
    rb_raise(rb_eArgError, "e2 (8th argument) must be NArray");
  if (NA_RANK(rb_e2) != 1)
    rb_raise(rb_eArgError, "rank of e2 (8th argument) must be %d", 1);
  n = NA_SHAPE0(rb_e2);
  if (NA_TYPE(rb_e2) != NA_DFLOAT)
    rb_e2 = na_change_type(rb_e2, NA_DFLOAT);
  e2 = NA_PTR_TYPE(rb_e2, doublereal*);
  il = NUM2INT(rb_il);
  spltol = NUM2DBL(rb_spltol);
  if (!NA_IsNArray(rb_d))
    rb_raise(rb_eArgError, "d (6th argument) must be NArray");
  if (NA_RANK(rb_d) != 1)
    rb_raise(rb_eArgError, "rank of d (6th argument) must be %d", 1);
  if (NA_SHAPE0(rb_d) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of d must be the same as shape 0 of e2");
  if (NA_TYPE(rb_d) != NA_DFLOAT)
    rb_d = na_change_type(rb_d, NA_DFLOAT);
  d = NA_PTR_TYPE(rb_d, doublereal*);
  if (!NA_IsNArray(rb_e))
    rb_raise(rb_eArgError, "e (7th argument) must be NArray");
  if (NA_RANK(rb_e) != 1)
    rb_raise(rb_eArgError, "rank of e (7th argument) must be %d", 1);
  if (NA_SHAPE0(rb_e) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of e must be the same as shape 0 of e2");
  if (NA_TYPE(rb_e) != NA_DFLOAT)
    rb_e = na_change_type(rb_e, NA_DFLOAT);
  e = NA_PTR_TYPE(rb_e, doublereal*);
  range = StringValueCStr(rb_range)[0];
  rtol1 = NUM2DBL(rb_rtol1);
  vu = NUM2DBL(rb_vu);
  rtol2 = NUM2DBL(rb_rtol2);
  {
    int shape[1];
    shape[0] = n;
    rb_isplit = na_make_object(NA_LINT, 1, shape, cNArray);
  }
  isplit = NA_PTR_TYPE(rb_isplit, integer*);
  {
    int shape[1];
    shape[0] = n;
    rb_w = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  w = NA_PTR_TYPE(rb_w, doublereal*);
  {
    int shape[1];
    shape[0] = n;
    rb_werr = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  werr = NA_PTR_TYPE(rb_werr, doublereal*);
  {
    int shape[1];
    shape[0] = n;
    rb_wgap = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  wgap = NA_PTR_TYPE(rb_wgap, doublereal*);
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
    rb_gers = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  gers = NA_PTR_TYPE(rb_gers, doublereal*);
  {
    int shape[1];
    shape[0] = n;
    rb_d_out__ = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  d_out__ = NA_PTR_TYPE(rb_d_out__, doublereal*);
  MEMCPY(d_out__, d, doublereal, NA_TOTAL(rb_d));
  rb_d = rb_d_out__;
  d = d_out__;
  {
    int shape[1];
    shape[0] = n;
    rb_e_out__ = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  e_out__ = NA_PTR_TYPE(rb_e_out__, doublereal*);
  MEMCPY(e_out__, e, doublereal, NA_TOTAL(rb_e));
  rb_e = rb_e_out__;
  e = e_out__;
  {
    int shape[1];
    shape[0] = n;
    rb_e2_out__ = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  e2_out__ = NA_PTR_TYPE(rb_e2_out__, doublereal*);
  MEMCPY(e2_out__, e2, doublereal, NA_TOTAL(rb_e2));
  rb_e2 = rb_e2_out__;
  e2 = e2_out__;
  work = ALLOC_N(doublereal, (6*n));
  iwork = ALLOC_N(integer, (5*n));

  dlarre_(&range, &n, &vl, &vu, &il, &iu, d, e, e2, &rtol1, &rtol2, &spltol, &nsplit, isplit, &m, w, werr, wgap, iblock, indexw, gers, &pivmin, work, iwork, &info);

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
init_lapack_dlarre(VALUE mLapack){
  rb_define_module_function(mLapack, "dlarre", rb_dlarre, -1);
}

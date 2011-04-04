#include "rb_lapack.h"

extern VOID zlarrv_(integer *n, doublereal *vl, doublereal *vu, doublereal *d, doublereal *l, doublereal *pivmin, integer *isplit, integer *m, integer *dol, integer *dou, doublereal *minrgp, doublereal *rtol1, doublereal *rtol2, doublereal *w, doublereal *werr, doublereal *wgap, integer *iblock, integer *indexw, doublereal *gers, doublecomplex *z, integer *ldz, integer *isuppz, doublereal *work, integer *iwork, integer *info);

static VALUE
rb_zlarrv(int argc, VALUE *argv, VALUE self){
  VALUE rb_vl;
  doublereal vl; 
  VALUE rb_vu;
  doublereal vu; 
  VALUE rb_d;
  doublereal *d; 
  VALUE rb_l;
  doublereal *l; 
  VALUE rb_pivmin;
  doublereal pivmin; 
  VALUE rb_isplit;
  integer *isplit; 
  VALUE rb_m;
  integer m; 
  VALUE rb_dol;
  integer dol; 
  VALUE rb_dou;
  integer dou; 
  VALUE rb_minrgp;
  doublereal minrgp; 
  VALUE rb_rtol1;
  doublereal rtol1; 
  VALUE rb_rtol2;
  doublereal rtol2; 
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
  VALUE rb_z;
  doublecomplex *z; 
  VALUE rb_isuppz;
  integer *isuppz; 
  VALUE rb_info;
  integer info; 
  VALUE rb_d_out__;
  doublereal *d_out__;
  VALUE rb_l_out__;
  doublereal *l_out__;
  VALUE rb_w_out__;
  doublereal *w_out__;
  VALUE rb_werr_out__;
  doublereal *werr_out__;
  VALUE rb_wgap_out__;
  doublereal *wgap_out__;
  doublereal *work;
  integer *iwork;

  integer n;
  integer ldz;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  z, isuppz, info, d, l, w, werr, wgap = NumRu::Lapack.zlarrv( vl, vu, d, l, pivmin, isplit, m, dol, dou, minrgp, rtol1, rtol2, w, werr, wgap, iblock, indexw, gers)\n    or\n  NumRu::Lapack.zlarrv  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 18)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 18)", argc);
  rb_vl = argv[0];
  rb_vu = argv[1];
  rb_d = argv[2];
  rb_l = argv[3];
  rb_pivmin = argv[4];
  rb_isplit = argv[5];
  rb_m = argv[6];
  rb_dol = argv[7];
  rb_dou = argv[8];
  rb_minrgp = argv[9];
  rb_rtol1 = argv[10];
  rb_rtol2 = argv[11];
  rb_w = argv[12];
  rb_werr = argv[13];
  rb_wgap = argv[14];
  rb_iblock = argv[15];
  rb_indexw = argv[16];
  rb_gers = argv[17];

  vl = NUM2DBL(rb_vl);
  if (!NA_IsNArray(rb_w))
    rb_raise(rb_eArgError, "w (13th argument) must be NArray");
  if (NA_RANK(rb_w) != 1)
    rb_raise(rb_eArgError, "rank of w (13th argument) must be %d", 1);
  n = NA_SHAPE0(rb_w);
  if (NA_TYPE(rb_w) != NA_DFLOAT)
    rb_w = na_change_type(rb_w, NA_DFLOAT);
  w = NA_PTR_TYPE(rb_w, doublereal*);
  dol = NUM2INT(rb_dol);
  if (!NA_IsNArray(rb_l))
    rb_raise(rb_eArgError, "l (4th argument) must be NArray");
  if (NA_RANK(rb_l) != 1)
    rb_raise(rb_eArgError, "rank of l (4th argument) must be %d", 1);
  if (NA_SHAPE0(rb_l) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of l must be the same as shape 0 of w");
  if (NA_TYPE(rb_l) != NA_DFLOAT)
    rb_l = na_change_type(rb_l, NA_DFLOAT);
  l = NA_PTR_TYPE(rb_l, doublereal*);
  pivmin = NUM2DBL(rb_pivmin);
  dou = NUM2INT(rb_dou);
  if (!NA_IsNArray(rb_wgap))
    rb_raise(rb_eArgError, "wgap (15th argument) must be NArray");
  if (NA_RANK(rb_wgap) != 1)
    rb_raise(rb_eArgError, "rank of wgap (15th argument) must be %d", 1);
  if (NA_SHAPE0(rb_wgap) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of wgap must be the same as shape 0 of w");
  if (NA_TYPE(rb_wgap) != NA_DFLOAT)
    rb_wgap = na_change_type(rb_wgap, NA_DFLOAT);
  wgap = NA_PTR_TYPE(rb_wgap, doublereal*);
  m = NUM2INT(rb_m);
  minrgp = NUM2DBL(rb_minrgp);
  rtol2 = NUM2DBL(rb_rtol2);
  if (!NA_IsNArray(rb_isplit))
    rb_raise(rb_eArgError, "isplit (6th argument) must be NArray");
  if (NA_RANK(rb_isplit) != 1)
    rb_raise(rb_eArgError, "rank of isplit (6th argument) must be %d", 1);
  if (NA_SHAPE0(rb_isplit) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of isplit must be the same as shape 0 of w");
  if (NA_TYPE(rb_isplit) != NA_LINT)
    rb_isplit = na_change_type(rb_isplit, NA_LINT);
  isplit = NA_PTR_TYPE(rb_isplit, integer*);
  if (!NA_IsNArray(rb_indexw))
    rb_raise(rb_eArgError, "indexw (17th argument) must be NArray");
  if (NA_RANK(rb_indexw) != 1)
    rb_raise(rb_eArgError, "rank of indexw (17th argument) must be %d", 1);
  if (NA_SHAPE0(rb_indexw) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of indexw must be the same as shape 0 of w");
  if (NA_TYPE(rb_indexw) != NA_LINT)
    rb_indexw = na_change_type(rb_indexw, NA_LINT);
  indexw = NA_PTR_TYPE(rb_indexw, integer*);
  if (!NA_IsNArray(rb_werr))
    rb_raise(rb_eArgError, "werr (14th argument) must be NArray");
  if (NA_RANK(rb_werr) != 1)
    rb_raise(rb_eArgError, "rank of werr (14th argument) must be %d", 1);
  if (NA_SHAPE0(rb_werr) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of werr must be the same as shape 0 of w");
  if (NA_TYPE(rb_werr) != NA_DFLOAT)
    rb_werr = na_change_type(rb_werr, NA_DFLOAT);
  werr = NA_PTR_TYPE(rb_werr, doublereal*);
  if (!NA_IsNArray(rb_iblock))
    rb_raise(rb_eArgError, "iblock (16th argument) must be NArray");
  if (NA_RANK(rb_iblock) != 1)
    rb_raise(rb_eArgError, "rank of iblock (16th argument) must be %d", 1);
  if (NA_SHAPE0(rb_iblock) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of iblock must be the same as shape 0 of w");
  if (NA_TYPE(rb_iblock) != NA_LINT)
    rb_iblock = na_change_type(rb_iblock, NA_LINT);
  iblock = NA_PTR_TYPE(rb_iblock, integer*);
  rtol1 = NUM2DBL(rb_rtol1);
  vu = NUM2DBL(rb_vu);
  if (!NA_IsNArray(rb_d))
    rb_raise(rb_eArgError, "d (3th argument) must be NArray");
  if (NA_RANK(rb_d) != 1)
    rb_raise(rb_eArgError, "rank of d (3th argument) must be %d", 1);
  if (NA_SHAPE0(rb_d) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of d must be the same as shape 0 of w");
  if (NA_TYPE(rb_d) != NA_DFLOAT)
    rb_d = na_change_type(rb_d, NA_DFLOAT);
  d = NA_PTR_TYPE(rb_d, doublereal*);
  if (!NA_IsNArray(rb_gers))
    rb_raise(rb_eArgError, "gers (18th argument) must be NArray");
  if (NA_RANK(rb_gers) != 1)
    rb_raise(rb_eArgError, "rank of gers (18th argument) must be %d", 1);
  if (NA_SHAPE0(rb_gers) != (2*n))
    rb_raise(rb_eRuntimeError, "shape 0 of gers must be %d", 2*n);
  if (NA_TYPE(rb_gers) != NA_DFLOAT)
    rb_gers = na_change_type(rb_gers, NA_DFLOAT);
  gers = NA_PTR_TYPE(rb_gers, doublereal*);
  ldz = n;
  {
    int shape[2];
    shape[0] = ldz;
    shape[1] = MAX(1,m);
    rb_z = na_make_object(NA_DCOMPLEX, 2, shape, cNArray);
  }
  z = NA_PTR_TYPE(rb_z, doublecomplex*);
  {
    int shape[1];
    shape[0] = 2*MAX(1,m);
    rb_isuppz = na_make_object(NA_LINT, 1, shape, cNArray);
  }
  isuppz = NA_PTR_TYPE(rb_isuppz, integer*);
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
    rb_l_out__ = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  l_out__ = NA_PTR_TYPE(rb_l_out__, doublereal*);
  MEMCPY(l_out__, l, doublereal, NA_TOTAL(rb_l));
  rb_l = rb_l_out__;
  l = l_out__;
  {
    int shape[1];
    shape[0] = n;
    rb_w_out__ = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  w_out__ = NA_PTR_TYPE(rb_w_out__, doublereal*);
  MEMCPY(w_out__, w, doublereal, NA_TOTAL(rb_w));
  rb_w = rb_w_out__;
  w = w_out__;
  {
    int shape[1];
    shape[0] = n;
    rb_werr_out__ = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  werr_out__ = NA_PTR_TYPE(rb_werr_out__, doublereal*);
  MEMCPY(werr_out__, werr, doublereal, NA_TOTAL(rb_werr));
  rb_werr = rb_werr_out__;
  werr = werr_out__;
  {
    int shape[1];
    shape[0] = n;
    rb_wgap_out__ = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  wgap_out__ = NA_PTR_TYPE(rb_wgap_out__, doublereal*);
  MEMCPY(wgap_out__, wgap, doublereal, NA_TOTAL(rb_wgap));
  rb_wgap = rb_wgap_out__;
  wgap = wgap_out__;
  work = ALLOC_N(doublereal, (12*n));
  iwork = ALLOC_N(integer, (7*n));

  zlarrv_(&n, &vl, &vu, d, l, &pivmin, isplit, &m, &dol, &dou, &minrgp, &rtol1, &rtol2, w, werr, wgap, iblock, indexw, gers, z, &ldz, isuppz, work, iwork, &info);

  free(work);
  free(iwork);
  rb_info = INT2NUM(info);
  return rb_ary_new3(8, rb_z, rb_isuppz, rb_info, rb_d, rb_l, rb_w, rb_werr, rb_wgap);
}

void
init_lapack_zlarrv(VALUE mLapack){
  rb_define_module_function(mLapack, "zlarrv", rb_zlarrv, -1);
}

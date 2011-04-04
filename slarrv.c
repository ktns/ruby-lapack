#include "rb_lapack.h"

extern VOID slarrv_(integer *n, real *vl, real *vu, real *d, real *l, real *pivmin, integer *isplit, integer *m, integer *dol, integer *dou, real *minrgp, real *rtol1, real *rtol2, real *w, real *werr, real *wgap, integer *iblock, integer *indexw, real *gers, real *z, integer *ldz, integer *isuppz, real *work, integer *iwork, integer *info);

static VALUE
rb_slarrv(int argc, VALUE *argv, VALUE self){
  VALUE rb_vl;
  real vl; 
  VALUE rb_vu;
  real vu; 
  VALUE rb_d;
  real *d; 
  VALUE rb_l;
  real *l; 
  VALUE rb_pivmin;
  real pivmin; 
  VALUE rb_isplit;
  integer *isplit; 
  VALUE rb_m;
  integer m; 
  VALUE rb_dol;
  integer dol; 
  VALUE rb_dou;
  integer dou; 
  VALUE rb_minrgp;
  real minrgp; 
  VALUE rb_rtol1;
  real rtol1; 
  VALUE rb_rtol2;
  real rtol2; 
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
  VALUE rb_z;
  real *z; 
  VALUE rb_isuppz;
  integer *isuppz; 
  VALUE rb_info;
  integer info; 
  VALUE rb_d_out__;
  real *d_out__;
  VALUE rb_l_out__;
  real *l_out__;
  VALUE rb_w_out__;
  real *w_out__;
  VALUE rb_werr_out__;
  real *werr_out__;
  VALUE rb_wgap_out__;
  real *wgap_out__;
  real *work;
  integer *iwork;

  integer n;
  integer ldz;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  z, isuppz, info, d, l, w, werr, wgap = NumRu::Lapack.slarrv( vl, vu, d, l, pivmin, isplit, m, dol, dou, minrgp, rtol1, rtol2, w, werr, wgap, iblock, indexw, gers)\n    or\n  NumRu::Lapack.slarrv  # print help\n\n\nFORTRAN MANUAL\n\n");
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

  vl = (real)NUM2DBL(rb_vl);
  if (!NA_IsNArray(rb_w))
    rb_raise(rb_eArgError, "w (13th argument) must be NArray");
  if (NA_RANK(rb_w) != 1)
    rb_raise(rb_eArgError, "rank of w (13th argument) must be %d", 1);
  n = NA_SHAPE0(rb_w);
  if (NA_TYPE(rb_w) != NA_SFLOAT)
    rb_w = na_change_type(rb_w, NA_SFLOAT);
  w = NA_PTR_TYPE(rb_w, real*);
  dol = NUM2INT(rb_dol);
  if (!NA_IsNArray(rb_l))
    rb_raise(rb_eArgError, "l (4th argument) must be NArray");
  if (NA_RANK(rb_l) != 1)
    rb_raise(rb_eArgError, "rank of l (4th argument) must be %d", 1);
  if (NA_SHAPE0(rb_l) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of l must be the same as shape 0 of w");
  if (NA_TYPE(rb_l) != NA_SFLOAT)
    rb_l = na_change_type(rb_l, NA_SFLOAT);
  l = NA_PTR_TYPE(rb_l, real*);
  pivmin = (real)NUM2DBL(rb_pivmin);
  dou = NUM2INT(rb_dou);
  if (!NA_IsNArray(rb_wgap))
    rb_raise(rb_eArgError, "wgap (15th argument) must be NArray");
  if (NA_RANK(rb_wgap) != 1)
    rb_raise(rb_eArgError, "rank of wgap (15th argument) must be %d", 1);
  if (NA_SHAPE0(rb_wgap) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of wgap must be the same as shape 0 of w");
  if (NA_TYPE(rb_wgap) != NA_SFLOAT)
    rb_wgap = na_change_type(rb_wgap, NA_SFLOAT);
  wgap = NA_PTR_TYPE(rb_wgap, real*);
  m = NUM2INT(rb_m);
  minrgp = (real)NUM2DBL(rb_minrgp);
  rtol2 = (real)NUM2DBL(rb_rtol2);
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
  if (NA_TYPE(rb_werr) != NA_SFLOAT)
    rb_werr = na_change_type(rb_werr, NA_SFLOAT);
  werr = NA_PTR_TYPE(rb_werr, real*);
  if (!NA_IsNArray(rb_iblock))
    rb_raise(rb_eArgError, "iblock (16th argument) must be NArray");
  if (NA_RANK(rb_iblock) != 1)
    rb_raise(rb_eArgError, "rank of iblock (16th argument) must be %d", 1);
  if (NA_SHAPE0(rb_iblock) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of iblock must be the same as shape 0 of w");
  if (NA_TYPE(rb_iblock) != NA_LINT)
    rb_iblock = na_change_type(rb_iblock, NA_LINT);
  iblock = NA_PTR_TYPE(rb_iblock, integer*);
  rtol1 = (real)NUM2DBL(rb_rtol1);
  vu = (real)NUM2DBL(rb_vu);
  if (!NA_IsNArray(rb_d))
    rb_raise(rb_eArgError, "d (3th argument) must be NArray");
  if (NA_RANK(rb_d) != 1)
    rb_raise(rb_eArgError, "rank of d (3th argument) must be %d", 1);
  if (NA_SHAPE0(rb_d) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of d must be the same as shape 0 of w");
  if (NA_TYPE(rb_d) != NA_SFLOAT)
    rb_d = na_change_type(rb_d, NA_SFLOAT);
  d = NA_PTR_TYPE(rb_d, real*);
  if (!NA_IsNArray(rb_gers))
    rb_raise(rb_eArgError, "gers (18th argument) must be NArray");
  if (NA_RANK(rb_gers) != 1)
    rb_raise(rb_eArgError, "rank of gers (18th argument) must be %d", 1);
  if (NA_SHAPE0(rb_gers) != (2*n))
    rb_raise(rb_eRuntimeError, "shape 0 of gers must be %d", 2*n);
  if (NA_TYPE(rb_gers) != NA_SFLOAT)
    rb_gers = na_change_type(rb_gers, NA_SFLOAT);
  gers = NA_PTR_TYPE(rb_gers, real*);
  ldz = n;
  {
    int shape[2];
    shape[0] = ldz;
    shape[1] = MAX(1,m);
    rb_z = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  z = NA_PTR_TYPE(rb_z, real*);
  {
    int shape[1];
    shape[0] = 2*MAX(1,m);
    rb_isuppz = na_make_object(NA_LINT, 1, shape, cNArray);
  }
  isuppz = NA_PTR_TYPE(rb_isuppz, integer*);
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
    rb_l_out__ = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  l_out__ = NA_PTR_TYPE(rb_l_out__, real*);
  MEMCPY(l_out__, l, real, NA_TOTAL(rb_l));
  rb_l = rb_l_out__;
  l = l_out__;
  {
    int shape[1];
    shape[0] = n;
    rb_w_out__ = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  w_out__ = NA_PTR_TYPE(rb_w_out__, real*);
  MEMCPY(w_out__, w, real, NA_TOTAL(rb_w));
  rb_w = rb_w_out__;
  w = w_out__;
  {
    int shape[1];
    shape[0] = n;
    rb_werr_out__ = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  werr_out__ = NA_PTR_TYPE(rb_werr_out__, real*);
  MEMCPY(werr_out__, werr, real, NA_TOTAL(rb_werr));
  rb_werr = rb_werr_out__;
  werr = werr_out__;
  {
    int shape[1];
    shape[0] = n;
    rb_wgap_out__ = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  wgap_out__ = NA_PTR_TYPE(rb_wgap_out__, real*);
  MEMCPY(wgap_out__, wgap, real, NA_TOTAL(rb_wgap));
  rb_wgap = rb_wgap_out__;
  wgap = wgap_out__;
  work = ALLOC_N(real, (12*n));
  iwork = ALLOC_N(integer, (7*n));

  slarrv_(&n, &vl, &vu, d, l, &pivmin, isplit, &m, &dol, &dou, &minrgp, &rtol1, &rtol2, w, werr, wgap, iblock, indexw, gers, z, &ldz, isuppz, work, iwork, &info);

  free(work);
  free(iwork);
  rb_info = INT2NUM(info);
  return rb_ary_new3(8, rb_z, rb_isuppz, rb_info, rb_d, rb_l, rb_w, rb_werr, rb_wgap);
}

void
init_lapack_slarrv(VALUE mLapack){
  rb_define_module_function(mLapack, "slarrv", rb_slarrv, -1);
}

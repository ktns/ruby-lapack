#include "rb_lapack.h"

extern VOID dlarrb_(integer *n, doublereal *d, doublereal *lld, integer *ifirst, integer *ilast, doublereal *rtol1, doublereal *rtol2, integer *offset, doublereal *w, doublereal *wgap, doublereal *werr, doublereal *work, integer *iwork, doublereal *pivmin, doublereal *spdiam, integer *twist, integer *info);

static VALUE
rb_dlarrb(int argc, VALUE *argv, VALUE self){
  VALUE rb_d;
  doublereal *d; 
  VALUE rb_lld;
  doublereal *lld; 
  VALUE rb_ifirst;
  integer ifirst; 
  VALUE rb_ilast;
  integer ilast; 
  VALUE rb_rtol1;
  doublereal rtol1; 
  VALUE rb_rtol2;
  doublereal rtol2; 
  VALUE rb_offset;
  integer offset; 
  VALUE rb_w;
  doublereal *w; 
  VALUE rb_wgap;
  doublereal *wgap; 
  VALUE rb_werr;
  doublereal *werr; 
  VALUE rb_pivmin;
  doublereal pivmin; 
  VALUE rb_spdiam;
  doublereal spdiam; 
  VALUE rb_twist;
  integer twist; 
  VALUE rb_info;
  integer info; 
  VALUE rb_w_out__;
  doublereal *w_out__;
  VALUE rb_wgap_out__;
  doublereal *wgap_out__;
  VALUE rb_werr_out__;
  doublereal *werr_out__;
  doublereal *work;
  integer *iwork;

  integer n;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  info, w, wgap, werr = NumRu::Lapack.dlarrb( d, lld, ifirst, ilast, rtol1, rtol2, offset, w, wgap, werr, pivmin, spdiam, twist)\n    or\n  NumRu::Lapack.dlarrb  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 13)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 13)", argc);
  rb_d = argv[0];
  rb_lld = argv[1];
  rb_ifirst = argv[2];
  rb_ilast = argv[3];
  rb_rtol1 = argv[4];
  rb_rtol2 = argv[5];
  rb_offset = argv[6];
  rb_w = argv[7];
  rb_wgap = argv[8];
  rb_werr = argv[9];
  rb_pivmin = argv[10];
  rb_spdiam = argv[11];
  rb_twist = argv[12];

  ilast = NUM2INT(rb_ilast);
  if (!NA_IsNArray(rb_w))
    rb_raise(rb_eArgError, "w (8th argument) must be NArray");
  if (NA_RANK(rb_w) != 1)
    rb_raise(rb_eArgError, "rank of w (8th argument) must be %d", 1);
  n = NA_SHAPE0(rb_w);
  if (NA_TYPE(rb_w) != NA_DFLOAT)
    rb_w = na_change_type(rb_w, NA_DFLOAT);
  w = NA_PTR_TYPE(rb_w, doublereal*);
  rtol2 = NUM2DBL(rb_rtol2);
  spdiam = NUM2DBL(rb_spdiam);
  offset = NUM2INT(rb_offset);
  if (!NA_IsNArray(rb_d))
    rb_raise(rb_eArgError, "d (1th argument) must be NArray");
  if (NA_RANK(rb_d) != 1)
    rb_raise(rb_eArgError, "rank of d (1th argument) must be %d", 1);
  if (NA_SHAPE0(rb_d) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of d must be the same as shape 0 of w");
  if (NA_TYPE(rb_d) != NA_DFLOAT)
    rb_d = na_change_type(rb_d, NA_DFLOAT);
  d = NA_PTR_TYPE(rb_d, doublereal*);
  pivmin = NUM2DBL(rb_pivmin);
  twist = NUM2INT(rb_twist);
  rtol1 = NUM2DBL(rb_rtol1);
  ifirst = NUM2INT(rb_ifirst);
  if (!NA_IsNArray(rb_werr))
    rb_raise(rb_eArgError, "werr (10th argument) must be NArray");
  if (NA_RANK(rb_werr) != 1)
    rb_raise(rb_eArgError, "rank of werr (10th argument) must be %d", 1);
  if (NA_SHAPE0(rb_werr) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of werr must be the same as shape 0 of w");
  if (NA_TYPE(rb_werr) != NA_DFLOAT)
    rb_werr = na_change_type(rb_werr, NA_DFLOAT);
  werr = NA_PTR_TYPE(rb_werr, doublereal*);
  if (!NA_IsNArray(rb_lld))
    rb_raise(rb_eArgError, "lld (2th argument) must be NArray");
  if (NA_RANK(rb_lld) != 1)
    rb_raise(rb_eArgError, "rank of lld (2th argument) must be %d", 1);
  if (NA_SHAPE0(rb_lld) != (n-1))
    rb_raise(rb_eRuntimeError, "shape 0 of lld must be %d", n-1);
  if (NA_TYPE(rb_lld) != NA_DFLOAT)
    rb_lld = na_change_type(rb_lld, NA_DFLOAT);
  lld = NA_PTR_TYPE(rb_lld, doublereal*);
  if (!NA_IsNArray(rb_wgap))
    rb_raise(rb_eArgError, "wgap (9th argument) must be NArray");
  if (NA_RANK(rb_wgap) != 1)
    rb_raise(rb_eArgError, "rank of wgap (9th argument) must be %d", 1);
  if (NA_SHAPE0(rb_wgap) != (n-1))
    rb_raise(rb_eRuntimeError, "shape 0 of wgap must be %d", n-1);
  if (NA_TYPE(rb_wgap) != NA_DFLOAT)
    rb_wgap = na_change_type(rb_wgap, NA_DFLOAT);
  wgap = NA_PTR_TYPE(rb_wgap, doublereal*);
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
    shape[0] = n-1;
    rb_wgap_out__ = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  wgap_out__ = NA_PTR_TYPE(rb_wgap_out__, doublereal*);
  MEMCPY(wgap_out__, wgap, doublereal, NA_TOTAL(rb_wgap));
  rb_wgap = rb_wgap_out__;
  wgap = wgap_out__;
  {
    int shape[1];
    shape[0] = n;
    rb_werr_out__ = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  werr_out__ = NA_PTR_TYPE(rb_werr_out__, doublereal*);
  MEMCPY(werr_out__, werr, doublereal, NA_TOTAL(rb_werr));
  rb_werr = rb_werr_out__;
  werr = werr_out__;
  work = ALLOC_N(doublereal, (2*n));
  iwork = ALLOC_N(integer, (2*n));

  dlarrb_(&n, d, lld, &ifirst, &ilast, &rtol1, &rtol2, &offset, w, wgap, werr, work, iwork, &pivmin, &spdiam, &twist, &info);

  free(work);
  free(iwork);
  rb_info = INT2NUM(info);
  return rb_ary_new3(4, rb_info, rb_w, rb_wgap, rb_werr);
}

void
init_lapack_dlarrb(VALUE mLapack){
  rb_define_module_function(mLapack, "dlarrb", rb_dlarrb, -1);
}

#include "rb_lapack.h"

extern VOID dlaqr5_(logical *wantt, logical *wantz, integer *kacc22, integer *n, integer *ktop, integer *kbot, integer *nshfts, doublereal *sr, doublereal *si, doublereal *h, integer *ldh, integer *iloz, integer *ihiz, doublereal *z, integer *ldz, doublereal *v, integer *ldv, doublereal *u, integer *ldu, integer *nv, doublereal *wv, integer *ldwv, integer *nh, doublereal *wh, integer *ldwh);

static VALUE
rb_dlaqr5(int argc, VALUE *argv, VALUE self){
  VALUE rb_wantt;
  logical wantt; 
  VALUE rb_wantz;
  logical wantz; 
  VALUE rb_kacc22;
  integer kacc22; 
  VALUE rb_ktop;
  integer ktop; 
  VALUE rb_kbot;
  integer kbot; 
  VALUE rb_sr;
  doublereal *sr; 
  VALUE rb_si;
  doublereal *si; 
  VALUE rb_h;
  doublereal *h; 
  VALUE rb_iloz;
  integer iloz; 
  VALUE rb_ihiz;
  integer ihiz; 
  VALUE rb_z;
  doublereal *z; 
  VALUE rb_ldz;
  integer ldz; 
  VALUE rb_nv;
  integer nv; 
  VALUE rb_nh;
  integer nh; 
  VALUE rb_sr_out__;
  doublereal *sr_out__;
  VALUE rb_si_out__;
  doublereal *si_out__;
  VALUE rb_h_out__;
  doublereal *h_out__;
  VALUE rb_z_out__;
  doublereal *z_out__;
  doublereal *v;
  doublereal *u;
  doublereal *wv;
  doublereal *wh;

  integer nshfts;
  integer ldh;
  integer n;
  integer ldv;
  integer ldu;
  integer ldwv;
  integer ldwh;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  sr, si, h, z = NumRu::Lapack.dlaqr5( wantt, wantz, kacc22, ktop, kbot, sr, si, h, iloz, ihiz, z, ldz, nv, nh)\n    or\n  NumRu::Lapack.dlaqr5  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 14)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 14)", argc);
  rb_wantt = argv[0];
  rb_wantz = argv[1];
  rb_kacc22 = argv[2];
  rb_ktop = argv[3];
  rb_kbot = argv[4];
  rb_sr = argv[5];
  rb_si = argv[6];
  rb_h = argv[7];
  rb_iloz = argv[8];
  rb_ihiz = argv[9];
  rb_z = argv[10];
  rb_ldz = argv[11];
  rb_nv = argv[12];
  rb_nh = argv[13];

  if (!NA_IsNArray(rb_si))
    rb_raise(rb_eArgError, "si (7th argument) must be NArray");
  if (NA_RANK(rb_si) != 1)
    rb_raise(rb_eArgError, "rank of si (7th argument) must be %d", 1);
  nshfts = NA_SHAPE0(rb_si);
  if (NA_TYPE(rb_si) != NA_DFLOAT)
    rb_si = na_change_type(rb_si, NA_DFLOAT);
  si = NA_PTR_TYPE(rb_si, doublereal*);
  kacc22 = NUM2INT(rb_kacc22);
  ktop = NUM2INT(rb_ktop);
  wantz = (rb_wantz == Qtrue);
  if (!NA_IsNArray(rb_sr))
    rb_raise(rb_eArgError, "sr (6th argument) must be NArray");
  if (NA_RANK(rb_sr) != 1)
    rb_raise(rb_eArgError, "rank of sr (6th argument) must be %d", 1);
  if (NA_SHAPE0(rb_sr) != nshfts)
    rb_raise(rb_eRuntimeError, "shape 0 of sr must be the same as shape 0 of si");
  if (NA_TYPE(rb_sr) != NA_DFLOAT)
    rb_sr = na_change_type(rb_sr, NA_DFLOAT);
  sr = NA_PTR_TYPE(rb_sr, doublereal*);
  kbot = NUM2INT(rb_kbot);
  nh = NUM2INT(rb_nh);
  if (!NA_IsNArray(rb_h))
    rb_raise(rb_eArgError, "h (8th argument) must be NArray");
  if (NA_RANK(rb_h) != 2)
    rb_raise(rb_eArgError, "rank of h (8th argument) must be %d", 2);
  n = NA_SHAPE1(rb_h);
  ldh = NA_SHAPE0(rb_h);
  if (NA_TYPE(rb_h) != NA_DFLOAT)
    rb_h = na_change_type(rb_h, NA_DFLOAT);
  h = NA_PTR_TYPE(rb_h, doublereal*);
  nv = NUM2INT(rb_nv);
  ihiz = NUM2INT(rb_ihiz);
  wantt = (rb_wantt == Qtrue);
  ldv = 3;
  iloz = NUM2INT(rb_iloz);
  ldz = n;
  if (!NA_IsNArray(rb_z))
    rb_raise(rb_eArgError, "z (11th argument) must be NArray");
  if (NA_RANK(rb_z) != 2)
    rb_raise(rb_eArgError, "rank of z (11th argument) must be %d", 2);
  if (NA_SHAPE1(rb_z) != (wantz ? ihiz : 0))
    rb_raise(rb_eRuntimeError, "shape 1 of z must be %d", wantz ? ihiz : 0);
  if (NA_SHAPE0(rb_z) != (wantz ? ldz : 0))
    rb_raise(rb_eRuntimeError, "shape 0 of z must be %d", wantz ? ldz : 0);
  if (NA_TYPE(rb_z) != NA_DFLOAT)
    rb_z = na_change_type(rb_z, NA_DFLOAT);
  z = NA_PTR_TYPE(rb_z, doublereal*);
  ldwv = nv;
  ldu = 3*nshfts-3;
  ldwh = 3*nshfts-3;
  {
    int shape[1];
    shape[0] = nshfts;
    rb_sr_out__ = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  sr_out__ = NA_PTR_TYPE(rb_sr_out__, doublereal*);
  MEMCPY(sr_out__, sr, doublereal, NA_TOTAL(rb_sr));
  rb_sr = rb_sr_out__;
  sr = sr_out__;
  {
    int shape[1];
    shape[0] = nshfts;
    rb_si_out__ = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  si_out__ = NA_PTR_TYPE(rb_si_out__, doublereal*);
  MEMCPY(si_out__, si, doublereal, NA_TOTAL(rb_si));
  rb_si = rb_si_out__;
  si = si_out__;
  {
    int shape[2];
    shape[0] = ldh;
    shape[1] = n;
    rb_h_out__ = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  h_out__ = NA_PTR_TYPE(rb_h_out__, doublereal*);
  MEMCPY(h_out__, h, doublereal, NA_TOTAL(rb_h));
  rb_h = rb_h_out__;
  h = h_out__;
  {
    int shape[2];
    shape[0] = wantz ? ldz : 0;
    shape[1] = wantz ? ihiz : 0;
    rb_z_out__ = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  z_out__ = NA_PTR_TYPE(rb_z_out__, doublereal*);
  MEMCPY(z_out__, z, doublereal, NA_TOTAL(rb_z));
  rb_z = rb_z_out__;
  z = z_out__;
  v = ALLOC_N(doublereal, (ldv)*(nshfts/2));
  u = ALLOC_N(doublereal, (ldu)*(3*nshfts-3));
  wv = ALLOC_N(doublereal, (ldwv)*(3*nshfts-3));
  wh = ALLOC_N(doublereal, (ldwh)*(MAX(1,nh)));

  dlaqr5_(&wantt, &wantz, &kacc22, &n, &ktop, &kbot, &nshfts, sr, si, h, &ldh, &iloz, &ihiz, z, &ldz, v, &ldv, u, &ldu, &nv, wv, &ldwv, &nh, wh, &ldwh);

  free(v);
  free(u);
  free(wv);
  free(wh);
  return rb_ary_new3(4, rb_sr, rb_si, rb_h, rb_z);
}

void
init_lapack_dlaqr5(VALUE mLapack){
  rb_define_module_function(mLapack, "dlaqr5", rb_dlaqr5, -1);
}

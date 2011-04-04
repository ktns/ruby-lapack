#include "rb_lapack.h"

extern VOID slaqr5_(logical *wantt, logical *wantz, integer *kacc22, integer *n, integer *ktop, integer *kbot, integer *nshfts, real *sr, real *si, real *h, integer *ldh, integer *iloz, integer *ihiz, real *z, integer *ldz, real *v, integer *ldv, real *u, integer *ldu, integer *nv, real *wv, integer *ldwv, integer *nh, real *wh, integer *ldwh);

static VALUE
rb_slaqr5(int argc, VALUE *argv, VALUE self){
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
  real *sr; 
  VALUE rb_si;
  real *si; 
  VALUE rb_h;
  real *h; 
  VALUE rb_iloz;
  integer iloz; 
  VALUE rb_ihiz;
  integer ihiz; 
  VALUE rb_z;
  real *z; 
  VALUE rb_ldz;
  integer ldz; 
  VALUE rb_nv;
  integer nv; 
  VALUE rb_nh;
  integer nh; 
  VALUE rb_sr_out__;
  real *sr_out__;
  VALUE rb_si_out__;
  real *si_out__;
  VALUE rb_h_out__;
  real *h_out__;
  VALUE rb_z_out__;
  real *z_out__;
  real *v;
  real *u;
  real *wv;
  real *wh;

  integer nshfts;
  integer ldh;
  integer n;
  integer ldv;
  integer ldu;
  integer ldwv;
  integer ldwh;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  sr, si, h, z = NumRu::Lapack.slaqr5( wantt, wantz, kacc22, ktop, kbot, sr, si, h, iloz, ihiz, z, ldz, nv, nh)\n    or\n  NumRu::Lapack.slaqr5  # print help\n\n\nFORTRAN MANUAL\n\n");
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
  if (NA_TYPE(rb_si) != NA_SFLOAT)
    rb_si = na_change_type(rb_si, NA_SFLOAT);
  si = NA_PTR_TYPE(rb_si, real*);
  kacc22 = NUM2INT(rb_kacc22);
  ktop = NUM2INT(rb_ktop);
  wantz = (rb_wantz == Qtrue);
  if (!NA_IsNArray(rb_sr))
    rb_raise(rb_eArgError, "sr (6th argument) must be NArray");
  if (NA_RANK(rb_sr) != 1)
    rb_raise(rb_eArgError, "rank of sr (6th argument) must be %d", 1);
  if (NA_SHAPE0(rb_sr) != nshfts)
    rb_raise(rb_eRuntimeError, "shape 0 of sr must be the same as shape 0 of si");
  if (NA_TYPE(rb_sr) != NA_SFLOAT)
    rb_sr = na_change_type(rb_sr, NA_SFLOAT);
  sr = NA_PTR_TYPE(rb_sr, real*);
  kbot = NUM2INT(rb_kbot);
  nh = NUM2INT(rb_nh);
  if (!NA_IsNArray(rb_h))
    rb_raise(rb_eArgError, "h (8th argument) must be NArray");
  if (NA_RANK(rb_h) != 2)
    rb_raise(rb_eArgError, "rank of h (8th argument) must be %d", 2);
  n = NA_SHAPE1(rb_h);
  ldh = NA_SHAPE0(rb_h);
  if (NA_TYPE(rb_h) != NA_SFLOAT)
    rb_h = na_change_type(rb_h, NA_SFLOAT);
  h = NA_PTR_TYPE(rb_h, real*);
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
  if (NA_TYPE(rb_z) != NA_SFLOAT)
    rb_z = na_change_type(rb_z, NA_SFLOAT);
  z = NA_PTR_TYPE(rb_z, real*);
  ldwv = nv;
  ldu = 3*nshfts-3;
  ldwh = 3*nshfts-3;
  {
    int shape[1];
    shape[0] = nshfts;
    rb_sr_out__ = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  sr_out__ = NA_PTR_TYPE(rb_sr_out__, real*);
  MEMCPY(sr_out__, sr, real, NA_TOTAL(rb_sr));
  rb_sr = rb_sr_out__;
  sr = sr_out__;
  {
    int shape[1];
    shape[0] = nshfts;
    rb_si_out__ = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  si_out__ = NA_PTR_TYPE(rb_si_out__, real*);
  MEMCPY(si_out__, si, real, NA_TOTAL(rb_si));
  rb_si = rb_si_out__;
  si = si_out__;
  {
    int shape[2];
    shape[0] = ldh;
    shape[1] = n;
    rb_h_out__ = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  h_out__ = NA_PTR_TYPE(rb_h_out__, real*);
  MEMCPY(h_out__, h, real, NA_TOTAL(rb_h));
  rb_h = rb_h_out__;
  h = h_out__;
  {
    int shape[2];
    shape[0] = wantz ? ldz : 0;
    shape[1] = wantz ? ihiz : 0;
    rb_z_out__ = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  z_out__ = NA_PTR_TYPE(rb_z_out__, real*);
  MEMCPY(z_out__, z, real, NA_TOTAL(rb_z));
  rb_z = rb_z_out__;
  z = z_out__;
  v = ALLOC_N(real, (ldv)*(nshfts/2));
  u = ALLOC_N(real, (ldu)*(3*nshfts-3));
  wv = ALLOC_N(real, (ldwv)*(3*nshfts-3));
  wh = ALLOC_N(real, (ldwh)*(MAX(1,nh)));

  slaqr5_(&wantt, &wantz, &kacc22, &n, &ktop, &kbot, &nshfts, sr, si, h, &ldh, &iloz, &ihiz, z, &ldz, v, &ldv, u, &ldu, &nv, wv, &ldwv, &nh, wh, &ldwh);

  free(v);
  free(u);
  free(wv);
  free(wh);
  return rb_ary_new3(4, rb_sr, rb_si, rb_h, rb_z);
}

void
init_lapack_slaqr5(VALUE mLapack){
  rb_define_module_function(mLapack, "slaqr5", rb_slaqr5, -1);
}

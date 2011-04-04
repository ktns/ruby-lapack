#include "rb_lapack.h"

extern VOID zlaqr5_(logical *wantt, logical *wantz, integer *kacc22, integer *n, integer *ktop, integer *kbot, integer *nshfts, doublecomplex *s, doublecomplex *h, integer *ldh, integer *iloz, integer *ihiz, doublecomplex *z, integer *ldz, doublecomplex *v, integer *ldv, doublecomplex *u, integer *ldu, integer *nv, doublecomplex *wv, integer *ldwv, integer *nh, doublecomplex *wh, integer *ldwh);

static VALUE
rb_zlaqr5(int argc, VALUE *argv, VALUE self){
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
  VALUE rb_s;
  doublecomplex *s; 
  VALUE rb_h;
  doublecomplex *h; 
  VALUE rb_iloz;
  integer iloz; 
  VALUE rb_ihiz;
  integer ihiz; 
  VALUE rb_z;
  doublecomplex *z; 
  VALUE rb_ldz;
  integer ldz; 
  VALUE rb_nv;
  integer nv; 
  VALUE rb_nh;
  integer nh; 
  VALUE rb_s_out__;
  doublecomplex *s_out__;
  VALUE rb_h_out__;
  doublecomplex *h_out__;
  VALUE rb_z_out__;
  doublecomplex *z_out__;
  doublecomplex *v;
  doublecomplex *u;
  doublecomplex *wv;
  doublecomplex *wh;

  integer nshfts;
  integer ldh;
  integer n;
  integer ldv;
  integer ldu;
  integer ldwv;
  integer ldwh;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  s, h, z = NumRu::Lapack.zlaqr5( wantt, wantz, kacc22, ktop, kbot, s, h, iloz, ihiz, z, ldz, nv, nh)\n    or\n  NumRu::Lapack.zlaqr5  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 13)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 13)", argc);
  rb_wantt = argv[0];
  rb_wantz = argv[1];
  rb_kacc22 = argv[2];
  rb_ktop = argv[3];
  rb_kbot = argv[4];
  rb_s = argv[5];
  rb_h = argv[6];
  rb_iloz = argv[7];
  rb_ihiz = argv[8];
  rb_z = argv[9];
  rb_ldz = argv[10];
  rb_nv = argv[11];
  rb_nh = argv[12];

  kacc22 = NUM2INT(rb_kacc22);
  ktop = NUM2INT(rb_ktop);
  wantz = (rb_wantz == Qtrue);
  kbot = NUM2INT(rb_kbot);
  ldz = NUM2INT(rb_ldz);
  nh = NUM2INT(rb_nh);
  if (!NA_IsNArray(rb_s))
    rb_raise(rb_eArgError, "s (6th argument) must be NArray");
  if (NA_RANK(rb_s) != 1)
    rb_raise(rb_eArgError, "rank of s (6th argument) must be %d", 1);
  nshfts = NA_SHAPE0(rb_s);
  if (NA_TYPE(rb_s) != NA_DCOMPLEX)
    rb_s = na_change_type(rb_s, NA_DCOMPLEX);
  s = NA_PTR_TYPE(rb_s, doublecomplex*);
  if (!NA_IsNArray(rb_h))
    rb_raise(rb_eArgError, "h (7th argument) must be NArray");
  if (NA_RANK(rb_h) != 2)
    rb_raise(rb_eArgError, "rank of h (7th argument) must be %d", 2);
  n = NA_SHAPE1(rb_h);
  ldh = NA_SHAPE0(rb_h);
  if (NA_TYPE(rb_h) != NA_DCOMPLEX)
    rb_h = na_change_type(rb_h, NA_DCOMPLEX);
  h = NA_PTR_TYPE(rb_h, doublecomplex*);
  ldv = 3;
  wantt = (rb_wantt == Qtrue);
  nv = NUM2INT(rb_nv);
  ihiz = NUM2INT(rb_ihiz);
  iloz = NUM2INT(rb_iloz);
  ldu = 3*nshfts-3;
  ldwv = nv;
  if (!NA_IsNArray(rb_z))
    rb_raise(rb_eArgError, "z (10th argument) must be NArray");
  if (NA_RANK(rb_z) != 2)
    rb_raise(rb_eArgError, "rank of z (10th argument) must be %d", 2);
  if (NA_SHAPE1(rb_z) != (wantz ? ihiz : 0))
    rb_raise(rb_eRuntimeError, "shape 1 of z must be %d", wantz ? ihiz : 0);
  if (NA_SHAPE0(rb_z) != (wantz ? ldz : 0))
    rb_raise(rb_eRuntimeError, "shape 0 of z must be %d", wantz ? ldz : 0);
  if (NA_TYPE(rb_z) != NA_DCOMPLEX)
    rb_z = na_change_type(rb_z, NA_DCOMPLEX);
  z = NA_PTR_TYPE(rb_z, doublecomplex*);
  ldwh = 3*nshfts-3;
  {
    int shape[1];
    shape[0] = nshfts;
    rb_s_out__ = na_make_object(NA_DCOMPLEX, 1, shape, cNArray);
  }
  s_out__ = NA_PTR_TYPE(rb_s_out__, doublecomplex*);
  MEMCPY(s_out__, s, doublecomplex, NA_TOTAL(rb_s));
  rb_s = rb_s_out__;
  s = s_out__;
  {
    int shape[2];
    shape[0] = ldh;
    shape[1] = n;
    rb_h_out__ = na_make_object(NA_DCOMPLEX, 2, shape, cNArray);
  }
  h_out__ = NA_PTR_TYPE(rb_h_out__, doublecomplex*);
  MEMCPY(h_out__, h, doublecomplex, NA_TOTAL(rb_h));
  rb_h = rb_h_out__;
  h = h_out__;
  {
    int shape[2];
    shape[0] = wantz ? ldz : 0;
    shape[1] = wantz ? ihiz : 0;
    rb_z_out__ = na_make_object(NA_DCOMPLEX, 2, shape, cNArray);
  }
  z_out__ = NA_PTR_TYPE(rb_z_out__, doublecomplex*);
  MEMCPY(z_out__, z, doublecomplex, NA_TOTAL(rb_z));
  rb_z = rb_z_out__;
  z = z_out__;
  v = ALLOC_N(doublecomplex, (ldv)*(nshfts/2));
  u = ALLOC_N(doublecomplex, (ldu)*(3*nshfts-3));
  wv = ALLOC_N(doublecomplex, (ldwv)*(3*nshfts-3));
  wh = ALLOC_N(doublecomplex, (ldwh)*(MAX(1,nh)));

  zlaqr5_(&wantt, &wantz, &kacc22, &n, &ktop, &kbot, &nshfts, s, h, &ldh, &iloz, &ihiz, z, &ldz, v, &ldv, u, &ldu, &nv, wv, &ldwv, &nh, wh, &ldwh);

  free(v);
  free(u);
  free(wv);
  free(wh);
  return rb_ary_new3(3, rb_s, rb_h, rb_z);
}

void
init_lapack_zlaqr5(VALUE mLapack){
  rb_define_module_function(mLapack, "zlaqr5", rb_zlaqr5, -1);
}

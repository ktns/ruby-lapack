#include "rb_lapack.h"

extern VOID dlaqr3_(logical *wantt, logical *wantz, integer *n, integer *ktop, integer *kbot, integer *nw, doublereal *h, integer *ldh, integer *iloz, integer *ihiz, doublereal *z, integer *ldz, integer *ns, integer *nd, doublereal *sr, doublereal *si, doublereal *v, integer *ldv, integer *nh, doublereal *t, integer *ldt, integer *nv, doublereal *wv, integer *ldwv, doublereal *work, integer *lwork);

static VALUE
rb_dlaqr3(int argc, VALUE *argv, VALUE self){
  VALUE rb_wantt;
  logical wantt; 
  VALUE rb_wantz;
  logical wantz; 
  VALUE rb_ktop;
  integer ktop; 
  VALUE rb_kbot;
  integer kbot; 
  VALUE rb_nw;
  integer nw; 
  VALUE rb_h;
  doublereal *h; 
  VALUE rb_iloz;
  integer iloz; 
  VALUE rb_ihiz;
  integer ihiz; 
  VALUE rb_z;
  doublereal *z; 
  VALUE rb_nh;
  integer nh; 
  VALUE rb_nv;
  integer nv; 
  VALUE rb_lwork;
  integer lwork; 
  VALUE rb_ns;
  integer ns; 
  VALUE rb_nd;
  integer nd; 
  VALUE rb_sr;
  doublereal *sr; 
  VALUE rb_si;
  doublereal *si; 
  VALUE rb_h_out__;
  doublereal *h_out__;
  VALUE rb_z_out__;
  doublereal *z_out__;
  doublereal *v;
  doublereal *t;
  doublereal *wv;
  doublereal *work;

  integer ldh;
  integer n;
  integer ldz;
  integer ldv;
  integer ldt;
  integer ldwv;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  ns, nd, sr, si, h, z = NumRu::Lapack.dlaqr3( wantt, wantz, ktop, kbot, nw, h, iloz, ihiz, z, nh, nv, lwork)\n    or\n  NumRu::Lapack.dlaqr3  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 12)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 12)", argc);
  rb_wantt = argv[0];
  rb_wantz = argv[1];
  rb_ktop = argv[2];
  rb_kbot = argv[3];
  rb_nw = argv[4];
  rb_h = argv[5];
  rb_iloz = argv[6];
  rb_ihiz = argv[7];
  rb_z = argv[8];
  rb_nh = argv[9];
  rb_nv = argv[10];
  rb_lwork = argv[11];

  ktop = NUM2INT(rb_ktop);
  wantz = (rb_wantz == Qtrue);
  nh = NUM2INT(rb_nh);
  if (!NA_IsNArray(rb_z))
    rb_raise(rb_eArgError, "z (9th argument) must be NArray");
  if (NA_RANK(rb_z) != 2)
    rb_raise(rb_eArgError, "rank of z (9th argument) must be %d", 2);
  n = NA_SHAPE1(rb_z);
  ldz = NA_SHAPE0(rb_z);
  if (NA_TYPE(rb_z) != NA_DFLOAT)
    rb_z = na_change_type(rb_z, NA_DFLOAT);
  z = NA_PTR_TYPE(rb_z, doublereal*);
  kbot = NUM2INT(rb_kbot);
  lwork = NUM2INT(rb_lwork);
  iloz = NUM2INT(rb_iloz);
  nv = NUM2INT(rb_nv);
  wantt = (rb_wantt == Qtrue);
  ihiz = NUM2INT(rb_ihiz);
  if (!NA_IsNArray(rb_h))
    rb_raise(rb_eArgError, "h (6th argument) must be NArray");
  if (NA_RANK(rb_h) != 2)
    rb_raise(rb_eArgError, "rank of h (6th argument) must be %d", 2);
  if (NA_SHAPE1(rb_h) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of h must be the same as shape 1 of z");
  ldh = NA_SHAPE0(rb_h);
  if (NA_TYPE(rb_h) != NA_DFLOAT)
    rb_h = na_change_type(rb_h, NA_DFLOAT);
  h = NA_PTR_TYPE(rb_h, doublereal*);
  nw = NUM2INT(rb_nw);
  ldv = nw;
  ldwv = nw;
  ldt = nw;
  {
    int shape[1];
    shape[0] = MAX(1,kbot);
    rb_sr = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  sr = NA_PTR_TYPE(rb_sr, doublereal*);
  {
    int shape[1];
    shape[0] = MAX(1,kbot);
    rb_si = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  si = NA_PTR_TYPE(rb_si, doublereal*);
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
    shape[0] = ldz;
    shape[1] = n;
    rb_z_out__ = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  z_out__ = NA_PTR_TYPE(rb_z_out__, doublereal*);
  MEMCPY(z_out__, z, doublereal, NA_TOTAL(rb_z));
  rb_z = rb_z_out__;
  z = z_out__;
  v = ALLOC_N(doublereal, (ldv)*(MAX(1,nw)));
  t = ALLOC_N(doublereal, (ldt)*(MAX(1,nw)));
  wv = ALLOC_N(doublereal, (ldwv)*(MAX(1,nw)));
  work = ALLOC_N(doublereal, (MAX(1,lwork)));

  dlaqr3_(&wantt, &wantz, &n, &ktop, &kbot, &nw, h, &ldh, &iloz, &ihiz, z, &ldz, &ns, &nd, sr, si, v, &ldv, &nh, t, &ldt, &nv, wv, &ldwv, work, &lwork);

  free(v);
  free(t);
  free(wv);
  free(work);
  rb_ns = INT2NUM(ns);
  rb_nd = INT2NUM(nd);
  return rb_ary_new3(6, rb_ns, rb_nd, rb_sr, rb_si, rb_h, rb_z);
}

void
init_lapack_dlaqr3(VALUE mLapack){
  rb_define_module_function(mLapack, "dlaqr3", rb_dlaqr3, -1);
}

#include "rb_lapack.h"

extern VOID slaqr2_(logical *wantt, logical *wantz, integer *n, integer *ktop, integer *kbot, integer *nw, real *h, integer *ldh, integer *iloz, integer *ihiz, real *z, integer *ldz, integer *ns, integer *nd, real *sr, real *si, real *v, integer *ldv, integer *nh, real *t, integer *ldt, integer *nv, real *wv, integer *ldwv, real *work, integer *lwork);

static VALUE
rb_slaqr2(int argc, VALUE *argv, VALUE self){
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
  real *h; 
  VALUE rb_iloz;
  integer iloz; 
  VALUE rb_ihiz;
  integer ihiz; 
  VALUE rb_z;
  real *z; 
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
  real *sr; 
  VALUE rb_si;
  real *si; 
  VALUE rb_h_out__;
  real *h_out__;
  VALUE rb_z_out__;
  real *z_out__;
  real *v;
  real *t;
  real *wv;
  real *work;

  integer ldh;
  integer n;
  integer ldz;
  integer ldv;
  integer ldt;
  integer ldwv;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  ns, nd, sr, si, h, z = NumRu::Lapack.slaqr2( wantt, wantz, ktop, kbot, nw, h, iloz, ihiz, z, nh, nv, lwork)\n    or\n  NumRu::Lapack.slaqr2  # print help\n\n\nFORTRAN MANUAL\n\n");
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
  if (NA_TYPE(rb_z) != NA_SFLOAT)
    rb_z = na_change_type(rb_z, NA_SFLOAT);
  z = NA_PTR_TYPE(rb_z, real*);
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
  if (NA_TYPE(rb_h) != NA_SFLOAT)
    rb_h = na_change_type(rb_h, NA_SFLOAT);
  h = NA_PTR_TYPE(rb_h, real*);
  nw = NUM2INT(rb_nw);
  ldv = nw;
  ldwv = nw;
  ldt = nw;
  {
    int shape[1];
    shape[0] = MAX(1,kbot);
    rb_sr = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  sr = NA_PTR_TYPE(rb_sr, real*);
  {
    int shape[1];
    shape[0] = MAX(1,kbot);
    rb_si = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  si = NA_PTR_TYPE(rb_si, real*);
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
    shape[0] = ldz;
    shape[1] = n;
    rb_z_out__ = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  z_out__ = NA_PTR_TYPE(rb_z_out__, real*);
  MEMCPY(z_out__, z, real, NA_TOTAL(rb_z));
  rb_z = rb_z_out__;
  z = z_out__;
  v = ALLOC_N(real, (ldv)*(MAX(1,nw)));
  t = ALLOC_N(real, (ldt)*(MAX(1,nw)));
  wv = ALLOC_N(real, (ldwv)*(MAX(1,nw)));
  work = ALLOC_N(real, (MAX(1,lwork)));

  slaqr2_(&wantt, &wantz, &n, &ktop, &kbot, &nw, h, &ldh, &iloz, &ihiz, z, &ldz, &ns, &nd, sr, si, v, &ldv, &nh, t, &ldt, &nv, wv, &ldwv, work, &lwork);

  free(v);
  free(t);
  free(wv);
  free(work);
  rb_ns = INT2NUM(ns);
  rb_nd = INT2NUM(nd);
  return rb_ary_new3(6, rb_ns, rb_nd, rb_sr, rb_si, rb_h, rb_z);
}

void
init_lapack_slaqr2(VALUE mLapack){
  rb_define_module_function(mLapack, "slaqr2", rb_slaqr2, -1);
}

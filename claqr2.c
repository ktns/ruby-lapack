#include "rb_lapack.h"

extern VOID claqr2_(logical *wantt, logical *wantz, integer *n, integer *ktop, integer *kbot, integer *nw, complex *h, integer *ldh, integer *iloz, integer *ihiz, complex *z, integer *ldz, integer *ns, integer *nd, complex *sh, complex *v, integer *ldv, integer *nh, complex *t, integer *ldt, integer *nv, complex *wv, integer *ldwv, complex *work, integer *lwork);

static VALUE
rb_claqr2(int argc, VALUE *argv, VALUE self){
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
  complex *h; 
  VALUE rb_iloz;
  integer iloz; 
  VALUE rb_ihiz;
  integer ihiz; 
  VALUE rb_z;
  complex *z; 
  VALUE rb_nh;
  integer nh; 
  VALUE rb_ldt;
  integer ldt; 
  VALUE rb_nv;
  integer nv; 
  VALUE rb_ldwv;
  integer ldwv; 
  VALUE rb_lwork;
  integer lwork; 
  VALUE rb_ns;
  integer ns; 
  VALUE rb_nd;
  integer nd; 
  VALUE rb_sh;
  complex *sh; 
  VALUE rb_h_out__;
  complex *h_out__;
  VALUE rb_z_out__;
  complex *z_out__;
  complex *v;
  complex *t;
  complex *wv;
  complex *work;

  integer ldh;
  integer n;
  integer ldz;
  integer ldv;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  ns, nd, sh, h, z = NumRu::Lapack.claqr2( wantt, wantz, ktop, kbot, nw, h, iloz, ihiz, z, nh, ldt, nv, ldwv, lwork)\n    or\n  NumRu::Lapack.claqr2  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 14)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 14)", argc);
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
  rb_ldt = argv[10];
  rb_nv = argv[11];
  rb_ldwv = argv[12];
  rb_lwork = argv[13];

  ktop = NUM2INT(rb_ktop);
  wantz = (rb_wantz == Qtrue);
  nh = NUM2INT(rb_nh);
  if (!NA_IsNArray(rb_z))
    rb_raise(rb_eArgError, "z (9th argument) must be NArray");
  if (NA_RANK(rb_z) != 2)
    rb_raise(rb_eArgError, "rank of z (9th argument) must be %d", 2);
  n = NA_SHAPE1(rb_z);
  ldz = NA_SHAPE0(rb_z);
  if (NA_TYPE(rb_z) != NA_SCOMPLEX)
    rb_z = na_change_type(rb_z, NA_SCOMPLEX);
  z = NA_PTR_TYPE(rb_z, complex*);
  lwork = NUM2INT(rb_lwork);
  kbot = NUM2INT(rb_kbot);
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
  if (NA_TYPE(rb_h) != NA_SCOMPLEX)
    rb_h = na_change_type(rb_h, NA_SCOMPLEX);
  h = NA_PTR_TYPE(rb_h, complex*);
  nw = NUM2INT(rb_nw);
  ldv = nw;
  ldwv = nw;
  ldt = nw;
  {
    int shape[1];
    shape[0] = MAX(1,kbot);
    rb_sh = na_make_object(NA_SCOMPLEX, 1, shape, cNArray);
  }
  sh = NA_PTR_TYPE(rb_sh, complex*);
  {
    int shape[2];
    shape[0] = ldh;
    shape[1] = n;
    rb_h_out__ = na_make_object(NA_SCOMPLEX, 2, shape, cNArray);
  }
  h_out__ = NA_PTR_TYPE(rb_h_out__, complex*);
  MEMCPY(h_out__, h, complex, NA_TOTAL(rb_h));
  rb_h = rb_h_out__;
  h = h_out__;
  {
    int shape[2];
    shape[0] = ldz;
    shape[1] = n;
    rb_z_out__ = na_make_object(NA_SCOMPLEX, 2, shape, cNArray);
  }
  z_out__ = NA_PTR_TYPE(rb_z_out__, complex*);
  MEMCPY(z_out__, z, complex, NA_TOTAL(rb_z));
  rb_z = rb_z_out__;
  z = z_out__;
  v = ALLOC_N(complex, (ldv)*(MAX(1,nw)));
  t = ALLOC_N(complex, (ldv)*(MAX(1,nw)));
  wv = ALLOC_N(complex, (ldv)*(MAX(1,nw)));
  work = ALLOC_N(complex, (MAX(1,lwork)));

  claqr2_(&wantt, &wantz, &n, &ktop, &kbot, &nw, h, &ldh, &iloz, &ihiz, z, &ldz, &ns, &nd, sh, v, &ldv, &nh, t, &ldt, &nv, wv, &ldwv, work, &lwork);

  free(v);
  free(t);
  free(wv);
  free(work);
  rb_ns = INT2NUM(ns);
  rb_nd = INT2NUM(nd);
  return rb_ary_new3(5, rb_ns, rb_nd, rb_sh, rb_h, rb_z);
}

void
init_lapack_claqr2(VALUE mLapack){
  rb_define_module_function(mLapack, "claqr2", rb_claqr2, -1);
}

#include "rb_lapack.h"

extern VOID dlahqr_(logical *wantt, logical *wantz, integer *n, integer *ilo, integer *ihi, doublereal *h, integer *ldh, doublereal *wr, doublereal *wi, integer *iloz, integer *ihiz, doublereal *z, integer *ldz, integer *info);

static VALUE
rb_dlahqr(int argc, VALUE *argv, VALUE self){
  VALUE rb_wantt;
  logical wantt; 
  VALUE rb_wantz;
  logical wantz; 
  VALUE rb_ilo;
  integer ilo; 
  VALUE rb_ihi;
  integer ihi; 
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
  VALUE rb_wr;
  doublereal *wr; 
  VALUE rb_wi;
  doublereal *wi; 
  VALUE rb_info;
  integer info; 
  VALUE rb_h_out__;
  doublereal *h_out__;
  VALUE rb_z_out__;
  doublereal *z_out__;

  integer ldh;
  integer n;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  wr, wi, info, h, z = NumRu::Lapack.dlahqr( wantt, wantz, ilo, ihi, h, iloz, ihiz, z, ldz)\n    or\n  NumRu::Lapack.dlahqr  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 9)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 9)", argc);
  rb_wantt = argv[0];
  rb_wantz = argv[1];
  rb_ilo = argv[2];
  rb_ihi = argv[3];
  rb_h = argv[4];
  rb_iloz = argv[5];
  rb_ihiz = argv[6];
  rb_z = argv[7];
  rb_ldz = argv[8];

  ilo = NUM2INT(rb_ilo);
  wantz = (rb_wantz == Qtrue);
  ldz = NUM2INT(rb_ldz);
  ihi = NUM2INT(rb_ihi);
  if (!NA_IsNArray(rb_h))
    rb_raise(rb_eArgError, "h (5th argument) must be NArray");
  if (NA_RANK(rb_h) != 2)
    rb_raise(rb_eArgError, "rank of h (5th argument) must be %d", 2);
  n = NA_SHAPE1(rb_h);
  ldh = NA_SHAPE0(rb_h);
  if (NA_TYPE(rb_h) != NA_DFLOAT)
    rb_h = na_change_type(rb_h, NA_DFLOAT);
  h = NA_PTR_TYPE(rb_h, doublereal*);
  wantt = (rb_wantt == Qtrue);
  ihiz = NUM2INT(rb_ihiz);
  iloz = NUM2INT(rb_iloz);
  if (!NA_IsNArray(rb_z))
    rb_raise(rb_eArgError, "z (8th argument) must be NArray");
  if (NA_RANK(rb_z) != 2)
    rb_raise(rb_eArgError, "rank of z (8th argument) must be %d", 2);
  if (NA_SHAPE1(rb_z) != (wantz ? n : 0))
    rb_raise(rb_eRuntimeError, "shape 1 of z must be %d", wantz ? n : 0);
  if (NA_SHAPE0(rb_z) != (wantz ? ldz : 0))
    rb_raise(rb_eRuntimeError, "shape 0 of z must be %d", wantz ? ldz : 0);
  if (NA_TYPE(rb_z) != NA_DFLOAT)
    rb_z = na_change_type(rb_z, NA_DFLOAT);
  z = NA_PTR_TYPE(rb_z, doublereal*);
  {
    int shape[1];
    shape[0] = n;
    rb_wr = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  wr = NA_PTR_TYPE(rb_wr, doublereal*);
  {
    int shape[1];
    shape[0] = n;
    rb_wi = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  wi = NA_PTR_TYPE(rb_wi, doublereal*);
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
    shape[1] = wantz ? n : 0;
    rb_z_out__ = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  z_out__ = NA_PTR_TYPE(rb_z_out__, doublereal*);
  MEMCPY(z_out__, z, doublereal, NA_TOTAL(rb_z));
  rb_z = rb_z_out__;
  z = z_out__;

  dlahqr_(&wantt, &wantz, &n, &ilo, &ihi, h, &ldh, wr, wi, &iloz, &ihiz, z, &ldz, &info);

  rb_info = INT2NUM(info);
  return rb_ary_new3(5, rb_wr, rb_wi, rb_info, rb_h, rb_z);
}

void
init_lapack_dlahqr(VALUE mLapack){
  rb_define_module_function(mLapack, "dlahqr", rb_dlahqr, -1);
}

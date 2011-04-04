#include "rb_lapack.h"

extern VOID claqr0_(logical *wantt, logical *wantz, integer *n, integer *ilo, integer *ihi, complex *h, integer *ldh, complex *w, integer *iloz, integer *ihiz, complex *z, integer *ldz, complex *work, integer *lwork, integer *info);

static VALUE
rb_claqr0(int argc, VALUE *argv, VALUE self){
  VALUE rb_wantt;
  logical wantt; 
  VALUE rb_wantz;
  logical wantz; 
  VALUE rb_ilo;
  integer ilo; 
  VALUE rb_h;
  complex *h; 
  VALUE rb_iloz;
  integer iloz; 
  VALUE rb_ihiz;
  integer ihiz; 
  VALUE rb_z;
  complex *z; 
  VALUE rb_lwork;
  integer lwork; 
  VALUE rb_w;
  complex *w; 
  VALUE rb_work;
  complex *work; 
  VALUE rb_info;
  integer info; 
  VALUE rb_h_out__;
  complex *h_out__;
  VALUE rb_z_out__;
  complex *z_out__;

  integer ldh;
  integer n;
  integer ldz;
  integer ihi;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  w, work, info, h, z = NumRu::Lapack.claqr0( wantt, wantz, ilo, h, iloz, ihiz, z, lwork)\n    or\n  NumRu::Lapack.claqr0  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 8)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 8)", argc);
  rb_wantt = argv[0];
  rb_wantz = argv[1];
  rb_ilo = argv[2];
  rb_h = argv[3];
  rb_iloz = argv[4];
  rb_ihiz = argv[5];
  rb_z = argv[6];
  rb_lwork = argv[7];

  ilo = NUM2INT(rb_ilo);
  wantz = (rb_wantz == Qtrue);
  if (!NA_IsNArray(rb_z))
    rb_raise(rb_eArgError, "z (7th argument) must be NArray");
  if (NA_RANK(rb_z) != 2)
    rb_raise(rb_eArgError, "rank of z (7th argument) must be %d", 2);
  ihi = NA_SHAPE1(rb_z);
  ldz = NA_SHAPE0(rb_z);
  if (NA_TYPE(rb_z) != NA_SCOMPLEX)
    rb_z = na_change_type(rb_z, NA_SCOMPLEX);
  z = NA_PTR_TYPE(rb_z, complex*);
  lwork = NUM2INT(rb_lwork);
  iloz = NUM2INT(rb_iloz);
  wantt = (rb_wantt == Qtrue);
  ihiz = NUM2INT(rb_ihiz);
  if (!NA_IsNArray(rb_h))
    rb_raise(rb_eArgError, "h (4th argument) must be NArray");
  if (NA_RANK(rb_h) != 2)
    rb_raise(rb_eArgError, "rank of h (4th argument) must be %d", 2);
  n = NA_SHAPE1(rb_h);
  ldh = NA_SHAPE0(rb_h);
  if (NA_TYPE(rb_h) != NA_SCOMPLEX)
    rb_h = na_change_type(rb_h, NA_SCOMPLEX);
  h = NA_PTR_TYPE(rb_h, complex*);
  {
    int shape[1];
    shape[0] = n;
    rb_w = na_make_object(NA_SCOMPLEX, 1, shape, cNArray);
  }
  w = NA_PTR_TYPE(rb_w, complex*);
  {
    int shape[1];
    shape[0] = MAX(1,lwork);
    rb_work = na_make_object(NA_SCOMPLEX, 1, shape, cNArray);
  }
  work = NA_PTR_TYPE(rb_work, complex*);
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
    shape[1] = ihi;
    rb_z_out__ = na_make_object(NA_SCOMPLEX, 2, shape, cNArray);
  }
  z_out__ = NA_PTR_TYPE(rb_z_out__, complex*);
  MEMCPY(z_out__, z, complex, NA_TOTAL(rb_z));
  rb_z = rb_z_out__;
  z = z_out__;

  claqr0_(&wantt, &wantz, &n, &ilo, &ihi, h, &ldh, w, &iloz, &ihiz, z, &ldz, work, &lwork, &info);

  rb_info = INT2NUM(info);
  return rb_ary_new3(5, rb_w, rb_work, rb_info, rb_h, rb_z);
}

void
init_lapack_claqr0(VALUE mLapack){
  rb_define_module_function(mLapack, "claqr0", rb_claqr0, -1);
}

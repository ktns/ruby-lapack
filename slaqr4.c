#include "rb_lapack.h"

extern VOID slaqr4_(logical *wantt, logical *wantz, integer *n, integer *ilo, integer *ihi, real *h, integer *ldh, real *wr, real *wi, integer *iloz, integer *ihiz, real *z, integer *ldz, real *work, integer *lwork, integer *info);

static VALUE
rb_slaqr4(int argc, VALUE *argv, VALUE self){
  VALUE rb_wantt;
  logical wantt; 
  VALUE rb_wantz;
  logical wantz; 
  VALUE rb_ilo;
  integer ilo; 
  VALUE rb_h;
  real *h; 
  VALUE rb_iloz;
  integer iloz; 
  VALUE rb_ihiz;
  integer ihiz; 
  VALUE rb_z;
  real *z; 
  VALUE rb_lwork;
  integer lwork; 
  VALUE rb_wr;
  real *wr; 
  VALUE rb_wi;
  real *wi; 
  VALUE rb_work;
  real *work; 
  VALUE rb_info;
  integer info; 
  VALUE rb_h_out__;
  real *h_out__;
  VALUE rb_z_out__;
  real *z_out__;

  integer ldh;
  integer n;
  integer ldz;
  integer ihi;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  wr, wi, work, info, h, z = NumRu::Lapack.slaqr4( wantt, wantz, ilo, h, iloz, ihiz, z, lwork)\n    or\n  NumRu::Lapack.slaqr4  # print help\n\n\nFORTRAN MANUAL\n\n");
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
  if (NA_TYPE(rb_z) != NA_SFLOAT)
    rb_z = na_change_type(rb_z, NA_SFLOAT);
  z = NA_PTR_TYPE(rb_z, real*);
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
  if (NA_TYPE(rb_h) != NA_SFLOAT)
    rb_h = na_change_type(rb_h, NA_SFLOAT);
  h = NA_PTR_TYPE(rb_h, real*);
  {
    int shape[1];
    shape[0] = ihi;
    rb_wr = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  wr = NA_PTR_TYPE(rb_wr, real*);
  {
    int shape[1];
    shape[0] = ihi;
    rb_wi = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  wi = NA_PTR_TYPE(rb_wi, real*);
  {
    int shape[1];
    shape[0] = MAX(1,lwork);
    rb_work = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  work = NA_PTR_TYPE(rb_work, real*);
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
    shape[1] = ihi;
    rb_z_out__ = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  z_out__ = NA_PTR_TYPE(rb_z_out__, real*);
  MEMCPY(z_out__, z, real, NA_TOTAL(rb_z));
  rb_z = rb_z_out__;
  z = z_out__;

  slaqr4_(&wantt, &wantz, &n, &ilo, &ihi, h, &ldh, wr, wi, &iloz, &ihiz, z, &ldz, work, &lwork, &info);

  rb_info = INT2NUM(info);
  return rb_ary_new3(6, rb_wr, rb_wi, rb_work, rb_info, rb_h, rb_z);
}

void
init_lapack_slaqr4(VALUE mLapack){
  rb_define_module_function(mLapack, "slaqr4", rb_slaqr4, -1);
}

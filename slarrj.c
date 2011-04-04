#include "rb_lapack.h"

extern VOID slarrj_(integer *n, real *d, real *e2, integer *ifirst, integer *ilast, real *rtol, integer *offset, real *w, real *werr, real *work, integer *iwork, real *pivmin, real *spdiam, integer *info);

static VALUE
rb_slarrj(int argc, VALUE *argv, VALUE self){
  VALUE rb_d;
  real *d; 
  VALUE rb_e2;
  real *e2; 
  VALUE rb_ifirst;
  integer ifirst; 
  VALUE rb_ilast;
  integer ilast; 
  VALUE rb_rtol;
  real rtol; 
  VALUE rb_offset;
  integer offset; 
  VALUE rb_w;
  real *w; 
  VALUE rb_werr;
  real *werr; 
  VALUE rb_pivmin;
  real pivmin; 
  VALUE rb_spdiam;
  real spdiam; 
  VALUE rb_info;
  integer info; 
  VALUE rb_w_out__;
  real *w_out__;
  VALUE rb_werr_out__;
  real *werr_out__;
  real *work;
  integer *iwork;

  integer n;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  info, w, werr = NumRu::Lapack.slarrj( d, e2, ifirst, ilast, rtol, offset, w, werr, pivmin, spdiam)\n    or\n  NumRu::Lapack.slarrj  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 10)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 10)", argc);
  rb_d = argv[0];
  rb_e2 = argv[1];
  rb_ifirst = argv[2];
  rb_ilast = argv[3];
  rb_rtol = argv[4];
  rb_offset = argv[5];
  rb_w = argv[6];
  rb_werr = argv[7];
  rb_pivmin = argv[8];
  rb_spdiam = argv[9];

  ilast = NUM2INT(rb_ilast);
  if (!NA_IsNArray(rb_w))
    rb_raise(rb_eArgError, "w (7th argument) must be NArray");
  if (NA_RANK(rb_w) != 1)
    rb_raise(rb_eArgError, "rank of w (7th argument) must be %d", 1);
  n = NA_SHAPE0(rb_w);
  if (NA_TYPE(rb_w) != NA_SFLOAT)
    rb_w = na_change_type(rb_w, NA_SFLOAT);
  w = NA_PTR_TYPE(rb_w, real*);
  rtol = (real)NUM2DBL(rb_rtol);
  offset = NUM2INT(rb_offset);
  spdiam = (real)NUM2DBL(rb_spdiam);
  pivmin = (real)NUM2DBL(rb_pivmin);
  if (!NA_IsNArray(rb_werr))
    rb_raise(rb_eArgError, "werr (8th argument) must be NArray");
  if (NA_RANK(rb_werr) != 1)
    rb_raise(rb_eArgError, "rank of werr (8th argument) must be %d", 1);
  if (NA_SHAPE0(rb_werr) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of werr must be the same as shape 0 of w");
  if (NA_TYPE(rb_werr) != NA_SFLOAT)
    rb_werr = na_change_type(rb_werr, NA_SFLOAT);
  werr = NA_PTR_TYPE(rb_werr, real*);
  ifirst = NUM2INT(rb_ifirst);
  if (!NA_IsNArray(rb_d))
    rb_raise(rb_eArgError, "d (1th argument) must be NArray");
  if (NA_RANK(rb_d) != 1)
    rb_raise(rb_eArgError, "rank of d (1th argument) must be %d", 1);
  if (NA_SHAPE0(rb_d) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of d must be the same as shape 0 of w");
  if (NA_TYPE(rb_d) != NA_SFLOAT)
    rb_d = na_change_type(rb_d, NA_SFLOAT);
  d = NA_PTR_TYPE(rb_d, real*);
  if (!NA_IsNArray(rb_e2))
    rb_raise(rb_eArgError, "e2 (2th argument) must be NArray");
  if (NA_RANK(rb_e2) != 1)
    rb_raise(rb_eArgError, "rank of e2 (2th argument) must be %d", 1);
  if (NA_SHAPE0(rb_e2) != (n-1))
    rb_raise(rb_eRuntimeError, "shape 0 of e2 must be %d", n-1);
  if (NA_TYPE(rb_e2) != NA_SFLOAT)
    rb_e2 = na_change_type(rb_e2, NA_SFLOAT);
  e2 = NA_PTR_TYPE(rb_e2, real*);
  {
    int shape[1];
    shape[0] = n;
    rb_w_out__ = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  w_out__ = NA_PTR_TYPE(rb_w_out__, real*);
  MEMCPY(w_out__, w, real, NA_TOTAL(rb_w));
  rb_w = rb_w_out__;
  w = w_out__;
  {
    int shape[1];
    shape[0] = n;
    rb_werr_out__ = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  werr_out__ = NA_PTR_TYPE(rb_werr_out__, real*);
  MEMCPY(werr_out__, werr, real, NA_TOTAL(rb_werr));
  rb_werr = rb_werr_out__;
  werr = werr_out__;
  work = ALLOC_N(real, (2*n));
  iwork = ALLOC_N(integer, (2*n));

  slarrj_(&n, d, e2, &ifirst, &ilast, &rtol, &offset, w, werr, work, iwork, &pivmin, &spdiam, &info);

  free(work);
  free(iwork);
  rb_info = INT2NUM(info);
  return rb_ary_new3(3, rb_info, rb_w, rb_werr);
}

void
init_lapack_slarrj(VALUE mLapack){
  rb_define_module_function(mLapack, "slarrj", rb_slarrj, -1);
}

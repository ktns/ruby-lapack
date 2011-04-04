#include "rb_lapack.h"

extern VOID dlarrj_(integer *n, doublereal *d, doublereal *e2, integer *ifirst, integer *ilast, doublereal *rtol, integer *offset, doublereal *w, doublereal *werr, doublereal *work, integer *iwork, doublereal *pivmin, doublereal *spdiam, integer *info);

static VALUE
rb_dlarrj(int argc, VALUE *argv, VALUE self){
  VALUE rb_d;
  doublereal *d; 
  VALUE rb_e2;
  doublereal *e2; 
  VALUE rb_ifirst;
  integer ifirst; 
  VALUE rb_ilast;
  integer ilast; 
  VALUE rb_rtol;
  doublereal rtol; 
  VALUE rb_offset;
  integer offset; 
  VALUE rb_w;
  doublereal *w; 
  VALUE rb_werr;
  doublereal *werr; 
  VALUE rb_pivmin;
  doublereal pivmin; 
  VALUE rb_spdiam;
  doublereal spdiam; 
  VALUE rb_info;
  integer info; 
  VALUE rb_w_out__;
  doublereal *w_out__;
  VALUE rb_werr_out__;
  doublereal *werr_out__;
  doublereal *work;
  integer *iwork;

  integer n;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  info, w, werr = NumRu::Lapack.dlarrj( d, e2, ifirst, ilast, rtol, offset, w, werr, pivmin, spdiam)\n    or\n  NumRu::Lapack.dlarrj  # print help\n\n\nFORTRAN MANUAL\n\n");
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
  if (NA_TYPE(rb_w) != NA_DFLOAT)
    rb_w = na_change_type(rb_w, NA_DFLOAT);
  w = NA_PTR_TYPE(rb_w, doublereal*);
  rtol = NUM2DBL(rb_rtol);
  offset = NUM2INT(rb_offset);
  spdiam = NUM2DBL(rb_spdiam);
  pivmin = NUM2DBL(rb_pivmin);
  if (!NA_IsNArray(rb_werr))
    rb_raise(rb_eArgError, "werr (8th argument) must be NArray");
  if (NA_RANK(rb_werr) != 1)
    rb_raise(rb_eArgError, "rank of werr (8th argument) must be %d", 1);
  if (NA_SHAPE0(rb_werr) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of werr must be the same as shape 0 of w");
  if (NA_TYPE(rb_werr) != NA_DFLOAT)
    rb_werr = na_change_type(rb_werr, NA_DFLOAT);
  werr = NA_PTR_TYPE(rb_werr, doublereal*);
  ifirst = NUM2INT(rb_ifirst);
  if (!NA_IsNArray(rb_d))
    rb_raise(rb_eArgError, "d (1th argument) must be NArray");
  if (NA_RANK(rb_d) != 1)
    rb_raise(rb_eArgError, "rank of d (1th argument) must be %d", 1);
  if (NA_SHAPE0(rb_d) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of d must be the same as shape 0 of w");
  if (NA_TYPE(rb_d) != NA_DFLOAT)
    rb_d = na_change_type(rb_d, NA_DFLOAT);
  d = NA_PTR_TYPE(rb_d, doublereal*);
  if (!NA_IsNArray(rb_e2))
    rb_raise(rb_eArgError, "e2 (2th argument) must be NArray");
  if (NA_RANK(rb_e2) != 1)
    rb_raise(rb_eArgError, "rank of e2 (2th argument) must be %d", 1);
  if (NA_SHAPE0(rb_e2) != (n-1))
    rb_raise(rb_eRuntimeError, "shape 0 of e2 must be %d", n-1);
  if (NA_TYPE(rb_e2) != NA_DFLOAT)
    rb_e2 = na_change_type(rb_e2, NA_DFLOAT);
  e2 = NA_PTR_TYPE(rb_e2, doublereal*);
  {
    int shape[1];
    shape[0] = n;
    rb_w_out__ = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  w_out__ = NA_PTR_TYPE(rb_w_out__, doublereal*);
  MEMCPY(w_out__, w, doublereal, NA_TOTAL(rb_w));
  rb_w = rb_w_out__;
  w = w_out__;
  {
    int shape[1];
    shape[0] = n;
    rb_werr_out__ = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  werr_out__ = NA_PTR_TYPE(rb_werr_out__, doublereal*);
  MEMCPY(werr_out__, werr, doublereal, NA_TOTAL(rb_werr));
  rb_werr = rb_werr_out__;
  werr = werr_out__;
  work = ALLOC_N(doublereal, (2*n));
  iwork = ALLOC_N(integer, (2*n));

  dlarrj_(&n, d, e2, &ifirst, &ilast, &rtol, &offset, w, werr, work, iwork, &pivmin, &spdiam, &info);

  free(work);
  free(iwork);
  rb_info = INT2NUM(info);
  return rb_ary_new3(3, rb_info, rb_w, rb_werr);
}

void
init_lapack_dlarrj(VALUE mLapack){
  rb_define_module_function(mLapack, "dlarrj", rb_dlarrj, -1);
}

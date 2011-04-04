#include "rb_lapack.h"

extern VOID dlarrd_(char *range, char *order, integer *n, doublereal *vl, doublereal *vu, integer *il, integer *iu, doublereal *gers, doublereal *reltol, doublereal *d, doublereal *e, doublereal *e2, doublereal *pivmin, integer *nsplit, integer *isplit, integer *m, doublereal *w, doublereal *werr, doublereal *wl, doublereal *wu, integer *iblock, integer *indexw, doublereal *work, integer *iwork, integer *info);

static VALUE
rb_dlarrd(int argc, VALUE *argv, VALUE self){
  VALUE rb_range;
  char range; 
  VALUE rb_order;
  char order; 
  VALUE rb_vl;
  doublereal vl; 
  VALUE rb_vu;
  doublereal vu; 
  VALUE rb_il;
  integer il; 
  VALUE rb_iu;
  integer iu; 
  VALUE rb_gers;
  doublereal *gers; 
  VALUE rb_reltol;
  doublereal reltol; 
  VALUE rb_d;
  doublereal *d; 
  VALUE rb_e;
  doublereal *e; 
  VALUE rb_e2;
  doublereal *e2; 
  VALUE rb_pivmin;
  doublereal pivmin; 
  VALUE rb_nsplit;
  integer nsplit; 
  VALUE rb_isplit;
  integer *isplit; 
  VALUE rb_m;
  integer m; 
  VALUE rb_w;
  doublereal *w; 
  VALUE rb_werr;
  doublereal *werr; 
  VALUE rb_wl;
  doublereal wl; 
  VALUE rb_wu;
  doublereal wu; 
  VALUE rb_iblock;
  integer *iblock; 
  VALUE rb_indexw;
  integer *indexw; 
  VALUE rb_info;
  integer info; 
  doublereal *work;
  integer *iwork;

  integer n;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  m, w, werr, wl, wu, iblock, indexw, info = NumRu::Lapack.dlarrd( range, order, vl, vu, il, iu, gers, reltol, d, e, e2, pivmin, nsplit, isplit)\n    or\n  NumRu::Lapack.dlarrd  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 14)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 14)", argc);
  rb_range = argv[0];
  rb_order = argv[1];
  rb_vl = argv[2];
  rb_vu = argv[3];
  rb_il = argv[4];
  rb_iu = argv[5];
  rb_gers = argv[6];
  rb_reltol = argv[7];
  rb_d = argv[8];
  rb_e = argv[9];
  rb_e2 = argv[10];
  rb_pivmin = argv[11];
  rb_nsplit = argv[12];
  rb_isplit = argv[13];

  vl = NUM2DBL(rb_vl);
  iu = NUM2INT(rb_iu);
  pivmin = NUM2DBL(rb_pivmin);
  vu = NUM2DBL(rb_vu);
  nsplit = NUM2INT(rb_nsplit);
  il = NUM2INT(rb_il);
  if (!NA_IsNArray(rb_d))
    rb_raise(rb_eArgError, "d (9th argument) must be NArray");
  if (NA_RANK(rb_d) != 1)
    rb_raise(rb_eArgError, "rank of d (9th argument) must be %d", 1);
  n = NA_SHAPE0(rb_d);
  if (NA_TYPE(rb_d) != NA_DFLOAT)
    rb_d = na_change_type(rb_d, NA_DFLOAT);
  d = NA_PTR_TYPE(rb_d, doublereal*);
  if (!NA_IsNArray(rb_isplit))
    rb_raise(rb_eArgError, "isplit (14th argument) must be NArray");
  if (NA_RANK(rb_isplit) != 1)
    rb_raise(rb_eArgError, "rank of isplit (14th argument) must be %d", 1);
  if (NA_SHAPE0(rb_isplit) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of isplit must be the same as shape 0 of d");
  if (NA_TYPE(rb_isplit) != NA_LINT)
    rb_isplit = na_change_type(rb_isplit, NA_LINT);
  isplit = NA_PTR_TYPE(rb_isplit, integer*);
  range = StringValueCStr(rb_range)[0];
  order = StringValueCStr(rb_order)[0];
  reltol = NUM2DBL(rb_reltol);
  if (!NA_IsNArray(rb_e))
    rb_raise(rb_eArgError, "e (10th argument) must be NArray");
  if (NA_RANK(rb_e) != 1)
    rb_raise(rb_eArgError, "rank of e (10th argument) must be %d", 1);
  if (NA_SHAPE0(rb_e) != (n-1))
    rb_raise(rb_eRuntimeError, "shape 0 of e must be %d", n-1);
  if (NA_TYPE(rb_e) != NA_DFLOAT)
    rb_e = na_change_type(rb_e, NA_DFLOAT);
  e = NA_PTR_TYPE(rb_e, doublereal*);
  if (!NA_IsNArray(rb_gers))
    rb_raise(rb_eArgError, "gers (7th argument) must be NArray");
  if (NA_RANK(rb_gers) != 1)
    rb_raise(rb_eArgError, "rank of gers (7th argument) must be %d", 1);
  if (NA_SHAPE0(rb_gers) != (2*n))
    rb_raise(rb_eRuntimeError, "shape 0 of gers must be %d", 2*n);
  if (NA_TYPE(rb_gers) != NA_DFLOAT)
    rb_gers = na_change_type(rb_gers, NA_DFLOAT);
  gers = NA_PTR_TYPE(rb_gers, doublereal*);
  if (!NA_IsNArray(rb_e2))
    rb_raise(rb_eArgError, "e2 (11th argument) must be NArray");
  if (NA_RANK(rb_e2) != 1)
    rb_raise(rb_eArgError, "rank of e2 (11th argument) must be %d", 1);
  if (NA_SHAPE0(rb_e2) != (n-1))
    rb_raise(rb_eRuntimeError, "shape 0 of e2 must be %d", n-1);
  if (NA_TYPE(rb_e2) != NA_DFLOAT)
    rb_e2 = na_change_type(rb_e2, NA_DFLOAT);
  e2 = NA_PTR_TYPE(rb_e2, doublereal*);
  {
    int shape[1];
    shape[0] = n;
    rb_w = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  w = NA_PTR_TYPE(rb_w, doublereal*);
  {
    int shape[1];
    shape[0] = n;
    rb_werr = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  werr = NA_PTR_TYPE(rb_werr, doublereal*);
  {
    int shape[1];
    shape[0] = n;
    rb_iblock = na_make_object(NA_LINT, 1, shape, cNArray);
  }
  iblock = NA_PTR_TYPE(rb_iblock, integer*);
  {
    int shape[1];
    shape[0] = n;
    rb_indexw = na_make_object(NA_LINT, 1, shape, cNArray);
  }
  indexw = NA_PTR_TYPE(rb_indexw, integer*);
  work = ALLOC_N(doublereal, (4*n));
  iwork = ALLOC_N(integer, (3*n));

  dlarrd_(&range, &order, &n, &vl, &vu, &il, &iu, gers, &reltol, d, e, e2, &pivmin, &nsplit, isplit, &m, w, werr, &wl, &wu, iblock, indexw, work, iwork, &info);

  free(work);
  free(iwork);
  rb_m = INT2NUM(m);
  rb_wl = rb_float_new((double)wl);
  rb_wu = rb_float_new((double)wu);
  rb_info = INT2NUM(info);
  return rb_ary_new3(8, rb_m, rb_w, rb_werr, rb_wl, rb_wu, rb_iblock, rb_indexw, rb_info);
}

void
init_lapack_dlarrd(VALUE mLapack){
  rb_define_module_function(mLapack, "dlarrd", rb_dlarrd, -1);
}

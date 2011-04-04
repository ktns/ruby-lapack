#include "rb_lapack.h"

extern VOID zhsein_(char *side, char *eigsrc, char *initv, logical *select, integer *n, doublecomplex *h, integer *ldh, doublecomplex *w, doublecomplex *vl, integer *ldvl, doublecomplex *vr, integer *ldvr, integer *mm, integer *m, doublecomplex *work, doublereal *rwork, integer *ifaill, integer *ifailr, integer *info);

static VALUE
rb_zhsein(int argc, VALUE *argv, VALUE self){
  VALUE rb_side;
  char side; 
  VALUE rb_eigsrc;
  char eigsrc; 
  VALUE rb_initv;
  char initv; 
  VALUE rb_select;
  logical *select; 
  VALUE rb_h;
  doublecomplex *h; 
  VALUE rb_w;
  doublecomplex *w; 
  VALUE rb_vl;
  doublecomplex *vl; 
  VALUE rb_vr;
  doublecomplex *vr; 
  VALUE rb_m;
  integer m; 
  VALUE rb_ifaill;
  integer *ifaill; 
  VALUE rb_ifailr;
  integer *ifailr; 
  VALUE rb_info;
  integer info; 
  VALUE rb_w_out__;
  doublecomplex *w_out__;
  VALUE rb_vl_out__;
  doublecomplex *vl_out__;
  VALUE rb_vr_out__;
  doublecomplex *vr_out__;
  doublecomplex *work;
  doublereal *rwork;

  integer n;
  integer ldh;
  integer ldvl;
  integer mm;
  integer ldvr;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  m, ifaill, ifailr, info, w, vl, vr = NumRu::Lapack.zhsein( side, eigsrc, initv, select, h, w, vl, vr)\n    or\n  NumRu::Lapack.zhsein  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 8)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 8)", argc);
  rb_side = argv[0];
  rb_eigsrc = argv[1];
  rb_initv = argv[2];
  rb_select = argv[3];
  rb_h = argv[4];
  rb_w = argv[5];
  rb_vl = argv[6];
  rb_vr = argv[7];

  if (!NA_IsNArray(rb_vl))
    rb_raise(rb_eArgError, "vl (7th argument) must be NArray");
  if (NA_RANK(rb_vl) != 2)
    rb_raise(rb_eArgError, "rank of vl (7th argument) must be %d", 2);
  mm = NA_SHAPE1(rb_vl);
  ldvl = NA_SHAPE0(rb_vl);
  if (NA_TYPE(rb_vl) != NA_DCOMPLEX)
    rb_vl = na_change_type(rb_vl, NA_DCOMPLEX);
  vl = NA_PTR_TYPE(rb_vl, doublecomplex*);
  if (!NA_IsNArray(rb_w))
    rb_raise(rb_eArgError, "w (6th argument) must be NArray");
  if (NA_RANK(rb_w) != 1)
    rb_raise(rb_eArgError, "rank of w (6th argument) must be %d", 1);
  n = NA_SHAPE0(rb_w);
  if (NA_TYPE(rb_w) != NA_DCOMPLEX)
    rb_w = na_change_type(rb_w, NA_DCOMPLEX);
  w = NA_PTR_TYPE(rb_w, doublecomplex*);
  side = StringValueCStr(rb_side)[0];
  eigsrc = StringValueCStr(rb_eigsrc)[0];
  if (!NA_IsNArray(rb_vr))
    rb_raise(rb_eArgError, "vr (8th argument) must be NArray");
  if (NA_RANK(rb_vr) != 2)
    rb_raise(rb_eArgError, "rank of vr (8th argument) must be %d", 2);
  if (NA_SHAPE1(rb_vr) != mm)
    rb_raise(rb_eRuntimeError, "shape 1 of vr must be the same as shape 1 of vl");
  ldvr = NA_SHAPE0(rb_vr);
  if (NA_TYPE(rb_vr) != NA_DCOMPLEX)
    rb_vr = na_change_type(rb_vr, NA_DCOMPLEX);
  vr = NA_PTR_TYPE(rb_vr, doublecomplex*);
  initv = StringValueCStr(rb_initv)[0];
  if (!NA_IsNArray(rb_h))
    rb_raise(rb_eArgError, "h (5th argument) must be NArray");
  if (NA_RANK(rb_h) != 2)
    rb_raise(rb_eArgError, "rank of h (5th argument) must be %d", 2);
  if (NA_SHAPE1(rb_h) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of h must be the same as shape 0 of w");
  ldh = NA_SHAPE0(rb_h);
  if (NA_TYPE(rb_h) != NA_DCOMPLEX)
    rb_h = na_change_type(rb_h, NA_DCOMPLEX);
  h = NA_PTR_TYPE(rb_h, doublecomplex*);
  if (!NA_IsNArray(rb_select))
    rb_raise(rb_eArgError, "select (4th argument) must be NArray");
  if (NA_RANK(rb_select) != 1)
    rb_raise(rb_eArgError, "rank of select (4th argument) must be %d", 1);
  if (NA_SHAPE0(rb_select) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of select must be the same as shape 0 of w");
  if (NA_TYPE(rb_select) != NA_LINT)
    rb_select = na_change_type(rb_select, NA_LINT);
  select = NA_PTR_TYPE(rb_select, logical*);
  {
    int shape[1];
    shape[0] = mm;
    rb_ifaill = na_make_object(NA_LINT, 1, shape, cNArray);
  }
  ifaill = NA_PTR_TYPE(rb_ifaill, integer*);
  {
    int shape[1];
    shape[0] = mm;
    rb_ifailr = na_make_object(NA_LINT, 1, shape, cNArray);
  }
  ifailr = NA_PTR_TYPE(rb_ifailr, integer*);
  {
    int shape[1];
    shape[0] = n;
    rb_w_out__ = na_make_object(NA_DCOMPLEX, 1, shape, cNArray);
  }
  w_out__ = NA_PTR_TYPE(rb_w_out__, doublecomplex*);
  MEMCPY(w_out__, w, doublecomplex, NA_TOTAL(rb_w));
  rb_w = rb_w_out__;
  w = w_out__;
  {
    int shape[2];
    shape[0] = ldvl;
    shape[1] = mm;
    rb_vl_out__ = na_make_object(NA_DCOMPLEX, 2, shape, cNArray);
  }
  vl_out__ = NA_PTR_TYPE(rb_vl_out__, doublecomplex*);
  MEMCPY(vl_out__, vl, doublecomplex, NA_TOTAL(rb_vl));
  rb_vl = rb_vl_out__;
  vl = vl_out__;
  {
    int shape[2];
    shape[0] = ldvr;
    shape[1] = mm;
    rb_vr_out__ = na_make_object(NA_DCOMPLEX, 2, shape, cNArray);
  }
  vr_out__ = NA_PTR_TYPE(rb_vr_out__, doublecomplex*);
  MEMCPY(vr_out__, vr, doublecomplex, NA_TOTAL(rb_vr));
  rb_vr = rb_vr_out__;
  vr = vr_out__;
  work = ALLOC_N(doublecomplex, (n*n));
  rwork = ALLOC_N(doublereal, (n));

  zhsein_(&side, &eigsrc, &initv, select, &n, h, &ldh, w, vl, &ldvl, vr, &ldvr, &mm, &m, work, rwork, ifaill, ifailr, &info);

  free(work);
  free(rwork);
  rb_m = INT2NUM(m);
  rb_info = INT2NUM(info);
  return rb_ary_new3(7, rb_m, rb_ifaill, rb_ifailr, rb_info, rb_w, rb_vl, rb_vr);
}

void
init_lapack_zhsein(VALUE mLapack){
  rb_define_module_function(mLapack, "zhsein", rb_zhsein, -1);
}

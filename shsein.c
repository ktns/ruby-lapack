#include "rb_lapack.h"

extern VOID shsein_(char *side, char *eigsrc, char *initv, logical *select, integer *n, real *h, integer *ldh, real *wr, real *wi, real *vl, integer *ldvl, real *vr, integer *ldvr, integer *mm, integer *m, real *work, integer *ifaill, integer *ifailr, integer *info);

static VALUE
rb_shsein(int argc, VALUE *argv, VALUE self){
  VALUE rb_side;
  char side; 
  VALUE rb_eigsrc;
  char eigsrc; 
  VALUE rb_initv;
  char initv; 
  VALUE rb_select;
  logical *select; 
  VALUE rb_h;
  real *h; 
  VALUE rb_wr;
  real *wr; 
  VALUE rb_wi;
  real *wi; 
  VALUE rb_vl;
  real *vl; 
  VALUE rb_vr;
  real *vr; 
  VALUE rb_m;
  integer m; 
  VALUE rb_ifaill;
  integer *ifaill; 
  VALUE rb_ifailr;
  integer *ifailr; 
  VALUE rb_info;
  integer info; 
  VALUE rb_select_out__;
  logical *select_out__;
  VALUE rb_wr_out__;
  real *wr_out__;
  VALUE rb_vl_out__;
  real *vl_out__;
  VALUE rb_vr_out__;
  real *vr_out__;
  real *work;

  integer n;
  integer ldh;
  integer ldvl;
  integer mm;
  integer ldvr;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  m, ifaill, ifailr, info, select, wr, vl, vr = NumRu::Lapack.shsein( side, eigsrc, initv, select, h, wr, wi, vl, vr)\n    or\n  NumRu::Lapack.shsein  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 9)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 9)", argc);
  rb_side = argv[0];
  rb_eigsrc = argv[1];
  rb_initv = argv[2];
  rb_select = argv[3];
  rb_h = argv[4];
  rb_wr = argv[5];
  rb_wi = argv[6];
  rb_vl = argv[7];
  rb_vr = argv[8];

  if (!NA_IsNArray(rb_vl))
    rb_raise(rb_eArgError, "vl (8th argument) must be NArray");
  if (NA_RANK(rb_vl) != 2)
    rb_raise(rb_eArgError, "rank of vl (8th argument) must be %d", 2);
  mm = NA_SHAPE1(rb_vl);
  ldvl = NA_SHAPE0(rb_vl);
  if (NA_TYPE(rb_vl) != NA_SFLOAT)
    rb_vl = na_change_type(rb_vl, NA_SFLOAT);
  vl = NA_PTR_TYPE(rb_vl, real*);
  side = StringValueCStr(rb_side)[0];
  eigsrc = StringValueCStr(rb_eigsrc)[0];
  if (!NA_IsNArray(rb_wr))
    rb_raise(rb_eArgError, "wr (6th argument) must be NArray");
  if (NA_RANK(rb_wr) != 1)
    rb_raise(rb_eArgError, "rank of wr (6th argument) must be %d", 1);
  n = NA_SHAPE0(rb_wr);
  if (NA_TYPE(rb_wr) != NA_SFLOAT)
    rb_wr = na_change_type(rb_wr, NA_SFLOAT);
  wr = NA_PTR_TYPE(rb_wr, real*);
  if (!NA_IsNArray(rb_vr))
    rb_raise(rb_eArgError, "vr (9th argument) must be NArray");
  if (NA_RANK(rb_vr) != 2)
    rb_raise(rb_eArgError, "rank of vr (9th argument) must be %d", 2);
  if (NA_SHAPE1(rb_vr) != mm)
    rb_raise(rb_eRuntimeError, "shape 1 of vr must be the same as shape 1 of vl");
  ldvr = NA_SHAPE0(rb_vr);
  if (NA_TYPE(rb_vr) != NA_SFLOAT)
    rb_vr = na_change_type(rb_vr, NA_SFLOAT);
  vr = NA_PTR_TYPE(rb_vr, real*);
  initv = StringValueCStr(rb_initv)[0];
  if (!NA_IsNArray(rb_wi))
    rb_raise(rb_eArgError, "wi (7th argument) must be NArray");
  if (NA_RANK(rb_wi) != 1)
    rb_raise(rb_eArgError, "rank of wi (7th argument) must be %d", 1);
  if (NA_SHAPE0(rb_wi) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of wi must be the same as shape 0 of wr");
  if (NA_TYPE(rb_wi) != NA_SFLOAT)
    rb_wi = na_change_type(rb_wi, NA_SFLOAT);
  wi = NA_PTR_TYPE(rb_wi, real*);
  if (!NA_IsNArray(rb_h))
    rb_raise(rb_eArgError, "h (5th argument) must be NArray");
  if (NA_RANK(rb_h) != 2)
    rb_raise(rb_eArgError, "rank of h (5th argument) must be %d", 2);
  if (NA_SHAPE1(rb_h) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of h must be the same as shape 0 of wr");
  ldh = NA_SHAPE0(rb_h);
  if (NA_TYPE(rb_h) != NA_SFLOAT)
    rb_h = na_change_type(rb_h, NA_SFLOAT);
  h = NA_PTR_TYPE(rb_h, real*);
  if (!NA_IsNArray(rb_select))
    rb_raise(rb_eArgError, "select (4th argument) must be NArray");
  if (NA_RANK(rb_select) != 1)
    rb_raise(rb_eArgError, "rank of select (4th argument) must be %d", 1);
  if (NA_SHAPE0(rb_select) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of select must be the same as shape 0 of wr");
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
    rb_select_out__ = na_make_object(NA_LINT, 1, shape, cNArray);
  }
  select_out__ = NA_PTR_TYPE(rb_select_out__, logical*);
  MEMCPY(select_out__, select, logical, NA_TOTAL(rb_select));
  rb_select = rb_select_out__;
  select = select_out__;
  {
    int shape[1];
    shape[0] = n;
    rb_wr_out__ = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  wr_out__ = NA_PTR_TYPE(rb_wr_out__, real*);
  MEMCPY(wr_out__, wr, real, NA_TOTAL(rb_wr));
  rb_wr = rb_wr_out__;
  wr = wr_out__;
  {
    int shape[2];
    shape[0] = ldvl;
    shape[1] = mm;
    rb_vl_out__ = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  vl_out__ = NA_PTR_TYPE(rb_vl_out__, real*);
  MEMCPY(vl_out__, vl, real, NA_TOTAL(rb_vl));
  rb_vl = rb_vl_out__;
  vl = vl_out__;
  {
    int shape[2];
    shape[0] = ldvr;
    shape[1] = mm;
    rb_vr_out__ = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  vr_out__ = NA_PTR_TYPE(rb_vr_out__, real*);
  MEMCPY(vr_out__, vr, real, NA_TOTAL(rb_vr));
  rb_vr = rb_vr_out__;
  vr = vr_out__;
  work = ALLOC_N(real, ((n+2)*n));

  shsein_(&side, &eigsrc, &initv, select, &n, h, &ldh, wr, wi, vl, &ldvl, vr, &ldvr, &mm, &m, work, ifaill, ifailr, &info);

  free(work);
  rb_m = INT2NUM(m);
  rb_info = INT2NUM(info);
  return rb_ary_new3(8, rb_m, rb_ifaill, rb_ifailr, rb_info, rb_select, rb_wr, rb_vl, rb_vr);
}

void
init_lapack_shsein(VALUE mLapack){
  rb_define_module_function(mLapack, "shsein", rb_shsein, -1);
}

#include "rb_lapack.h"

extern VOID ztrevc_(char *side, char *howmny, logical *select, integer *n, doublecomplex *t, integer *ldt, doublecomplex *vl, integer *ldvl, doublecomplex *vr, integer *ldvr, integer *mm, integer *m, doublecomplex *work, doublereal *rwork, integer *info);

static VALUE
rb_ztrevc(int argc, VALUE *argv, VALUE self){
  VALUE rb_side;
  char side; 
  VALUE rb_howmny;
  char howmny; 
  VALUE rb_select;
  logical *select; 
  VALUE rb_t;
  doublecomplex *t; 
  VALUE rb_vl;
  doublecomplex *vl; 
  VALUE rb_vr;
  doublecomplex *vr; 
  VALUE rb_m;
  integer m; 
  VALUE rb_info;
  integer info; 
  VALUE rb_t_out__;
  doublecomplex *t_out__;
  VALUE rb_vl_out__;
  doublecomplex *vl_out__;
  VALUE rb_vr_out__;
  doublecomplex *vr_out__;
  doublecomplex *work;
  doublereal *rwork;

  integer n;
  integer ldt;
  integer ldvl;
  integer mm;
  integer ldvr;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  m, info, t, vl, vr = NumRu::Lapack.ztrevc( side, howmny, select, t, vl, vr)\n    or\n  NumRu::Lapack.ztrevc  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 6)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 6)", argc);
  rb_side = argv[0];
  rb_howmny = argv[1];
  rb_select = argv[2];
  rb_t = argv[3];
  rb_vl = argv[4];
  rb_vr = argv[5];

  if (!NA_IsNArray(rb_vl))
    rb_raise(rb_eArgError, "vl (5th argument) must be NArray");
  if (NA_RANK(rb_vl) != 2)
    rb_raise(rb_eArgError, "rank of vl (5th argument) must be %d", 2);
  mm = NA_SHAPE1(rb_vl);
  ldvl = NA_SHAPE0(rb_vl);
  if (NA_TYPE(rb_vl) != NA_DCOMPLEX)
    rb_vl = na_change_type(rb_vl, NA_DCOMPLEX);
  vl = NA_PTR_TYPE(rb_vl, doublecomplex*);
  side = StringValueCStr(rb_side)[0];
  if (!NA_IsNArray(rb_vr))
    rb_raise(rb_eArgError, "vr (6th argument) must be NArray");
  if (NA_RANK(rb_vr) != 2)
    rb_raise(rb_eArgError, "rank of vr (6th argument) must be %d", 2);
  if (NA_SHAPE1(rb_vr) != mm)
    rb_raise(rb_eRuntimeError, "shape 1 of vr must be the same as shape 1 of vl");
  ldvr = NA_SHAPE0(rb_vr);
  if (NA_TYPE(rb_vr) != NA_DCOMPLEX)
    rb_vr = na_change_type(rb_vr, NA_DCOMPLEX);
  vr = NA_PTR_TYPE(rb_vr, doublecomplex*);
  if (!NA_IsNArray(rb_select))
    rb_raise(rb_eArgError, "select (3th argument) must be NArray");
  if (NA_RANK(rb_select) != 1)
    rb_raise(rb_eArgError, "rank of select (3th argument) must be %d", 1);
  n = NA_SHAPE0(rb_select);
  if (NA_TYPE(rb_select) != NA_LINT)
    rb_select = na_change_type(rb_select, NA_LINT);
  select = NA_PTR_TYPE(rb_select, logical*);
  if (!NA_IsNArray(rb_t))
    rb_raise(rb_eArgError, "t (4th argument) must be NArray");
  if (NA_RANK(rb_t) != 2)
    rb_raise(rb_eArgError, "rank of t (4th argument) must be %d", 2);
  if (NA_SHAPE1(rb_t) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of t must be the same as shape 0 of select");
  ldt = NA_SHAPE0(rb_t);
  if (NA_TYPE(rb_t) != NA_DCOMPLEX)
    rb_t = na_change_type(rb_t, NA_DCOMPLEX);
  t = NA_PTR_TYPE(rb_t, doublecomplex*);
  howmny = StringValueCStr(rb_howmny)[0];
  {
    int shape[2];
    shape[0] = ldt;
    shape[1] = n;
    rb_t_out__ = na_make_object(NA_DCOMPLEX, 2, shape, cNArray);
  }
  t_out__ = NA_PTR_TYPE(rb_t_out__, doublecomplex*);
  MEMCPY(t_out__, t, doublecomplex, NA_TOTAL(rb_t));
  rb_t = rb_t_out__;
  t = t_out__;
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
  work = ALLOC_N(doublecomplex, (2*n));
  rwork = ALLOC_N(doublereal, (n));

  ztrevc_(&side, &howmny, select, &n, t, &ldt, vl, &ldvl, vr, &ldvr, &mm, &m, work, rwork, &info);

  free(work);
  free(rwork);
  rb_m = INT2NUM(m);
  rb_info = INT2NUM(info);
  return rb_ary_new3(5, rb_m, rb_info, rb_t, rb_vl, rb_vr);
}

void
init_lapack_ztrevc(VALUE mLapack){
  rb_define_module_function(mLapack, "ztrevc", rb_ztrevc, -1);
}

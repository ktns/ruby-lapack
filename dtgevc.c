#include "rb_lapack.h"

extern VOID dtgevc_(char *side, char *howmny, logical *select, integer *n, doublereal *s, integer *lds, doublereal *p, integer *ldp, doublereal *vl, integer *ldvl, doublereal *vr, integer *ldvr, integer *mm, integer *m, doublereal *work, integer *info);

static VALUE
rb_dtgevc(int argc, VALUE *argv, VALUE self){
  VALUE rb_side;
  char side; 
  VALUE rb_howmny;
  char howmny; 
  VALUE rb_select;
  logical *select; 
  VALUE rb_s;
  doublereal *s; 
  VALUE rb_p;
  doublereal *p; 
  VALUE rb_vl;
  doublereal *vl; 
  VALUE rb_vr;
  doublereal *vr; 
  VALUE rb_m;
  integer m; 
  VALUE rb_info;
  integer info; 
  VALUE rb_vl_out__;
  doublereal *vl_out__;
  VALUE rb_vr_out__;
  doublereal *vr_out__;
  doublereal *work;

  integer n;
  integer lds;
  integer ldp;
  integer ldvl;
  integer mm;
  integer ldvr;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  m, info, vl, vr = NumRu::Lapack.dtgevc( side, howmny, select, s, p, vl, vr)\n    or\n  NumRu::Lapack.dtgevc  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 7)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 7)", argc);
  rb_side = argv[0];
  rb_howmny = argv[1];
  rb_select = argv[2];
  rb_s = argv[3];
  rb_p = argv[4];
  rb_vl = argv[5];
  rb_vr = argv[6];

  if (!NA_IsNArray(rb_vl))
    rb_raise(rb_eArgError, "vl (6th argument) must be NArray");
  if (NA_RANK(rb_vl) != 2)
    rb_raise(rb_eArgError, "rank of vl (6th argument) must be %d", 2);
  mm = NA_SHAPE1(rb_vl);
  ldvl = NA_SHAPE0(rb_vl);
  if (NA_TYPE(rb_vl) != NA_DFLOAT)
    rb_vl = na_change_type(rb_vl, NA_DFLOAT);
  vl = NA_PTR_TYPE(rb_vl, doublereal*);
  side = StringValueCStr(rb_side)[0];
  if (!NA_IsNArray(rb_p))
    rb_raise(rb_eArgError, "p (5th argument) must be NArray");
  if (NA_RANK(rb_p) != 2)
    rb_raise(rb_eArgError, "rank of p (5th argument) must be %d", 2);
  n = NA_SHAPE1(rb_p);
  ldp = NA_SHAPE0(rb_p);
  if (NA_TYPE(rb_p) != NA_DFLOAT)
    rb_p = na_change_type(rb_p, NA_DFLOAT);
  p = NA_PTR_TYPE(rb_p, doublereal*);
  if (!NA_IsNArray(rb_vr))
    rb_raise(rb_eArgError, "vr (7th argument) must be NArray");
  if (NA_RANK(rb_vr) != 2)
    rb_raise(rb_eArgError, "rank of vr (7th argument) must be %d", 2);
  if (NA_SHAPE1(rb_vr) != mm)
    rb_raise(rb_eRuntimeError, "shape 1 of vr must be the same as shape 1 of vl");
  ldvr = NA_SHAPE0(rb_vr);
  if (NA_TYPE(rb_vr) != NA_DFLOAT)
    rb_vr = na_change_type(rb_vr, NA_DFLOAT);
  vr = NA_PTR_TYPE(rb_vr, doublereal*);
  if (!NA_IsNArray(rb_s))
    rb_raise(rb_eArgError, "s (4th argument) must be NArray");
  if (NA_RANK(rb_s) != 2)
    rb_raise(rb_eArgError, "rank of s (4th argument) must be %d", 2);
  if (NA_SHAPE1(rb_s) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of s must be the same as shape 1 of p");
  lds = NA_SHAPE0(rb_s);
  if (NA_TYPE(rb_s) != NA_DFLOAT)
    rb_s = na_change_type(rb_s, NA_DFLOAT);
  s = NA_PTR_TYPE(rb_s, doublereal*);
  if (!NA_IsNArray(rb_select))
    rb_raise(rb_eArgError, "select (3th argument) must be NArray");
  if (NA_RANK(rb_select) != 1)
    rb_raise(rb_eArgError, "rank of select (3th argument) must be %d", 1);
  if (NA_SHAPE0(rb_select) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of select must be the same as shape 1 of p");
  if (NA_TYPE(rb_select) != NA_LINT)
    rb_select = na_change_type(rb_select, NA_LINT);
  select = NA_PTR_TYPE(rb_select, logical*);
  howmny = StringValueCStr(rb_howmny)[0];
  {
    int shape[2];
    shape[0] = ldvl;
    shape[1] = mm;
    rb_vl_out__ = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  vl_out__ = NA_PTR_TYPE(rb_vl_out__, doublereal*);
  MEMCPY(vl_out__, vl, doublereal, NA_TOTAL(rb_vl));
  rb_vl = rb_vl_out__;
  vl = vl_out__;
  {
    int shape[2];
    shape[0] = ldvr;
    shape[1] = mm;
    rb_vr_out__ = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  vr_out__ = NA_PTR_TYPE(rb_vr_out__, doublereal*);
  MEMCPY(vr_out__, vr, doublereal, NA_TOTAL(rb_vr));
  rb_vr = rb_vr_out__;
  vr = vr_out__;
  work = ALLOC_N(doublereal, (6*n));

  dtgevc_(&side, &howmny, select, &n, s, &lds, p, &ldp, vl, &ldvl, vr, &ldvr, &mm, &m, work, &info);

  free(work);
  rb_m = INT2NUM(m);
  rb_info = INT2NUM(info);
  return rb_ary_new3(4, rb_m, rb_info, rb_vl, rb_vr);
}

void
init_lapack_dtgevc(VALUE mLapack){
  rb_define_module_function(mLapack, "dtgevc", rb_dtgevc, -1);
}

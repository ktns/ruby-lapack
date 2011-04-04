#include "rb_lapack.h"

extern VOID dlarrc_(char *jobt, integer *n, doublereal *vl, doublereal *vu, doublereal *d, doublereal *e, doublereal *pivmin, integer *eigcnt, integer *lcnt, integer *rcnt, integer *info);

static VALUE
rb_dlarrc(int argc, VALUE *argv, VALUE self){
  VALUE rb_jobt;
  char jobt; 
  VALUE rb_vl;
  doublereal vl; 
  VALUE rb_vu;
  doublereal vu; 
  VALUE rb_d;
  doublereal *d; 
  VALUE rb_e;
  doublereal *e; 
  VALUE rb_pivmin;
  doublereal pivmin; 
  VALUE rb_eigcnt;
  integer eigcnt; 
  VALUE rb_lcnt;
  integer lcnt; 
  VALUE rb_rcnt;
  integer rcnt; 
  VALUE rb_info;
  integer info; 

  integer n;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  eigcnt, lcnt, rcnt, info = NumRu::Lapack.dlarrc( jobt, vl, vu, d, e, pivmin)\n    or\n  NumRu::Lapack.dlarrc  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 6)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 6)", argc);
  rb_jobt = argv[0];
  rb_vl = argv[1];
  rb_vu = argv[2];
  rb_d = argv[3];
  rb_e = argv[4];
  rb_pivmin = argv[5];

  vl = NUM2DBL(rb_vl);
  jobt = StringValueCStr(rb_jobt)[0];
  if (!NA_IsNArray(rb_d))
    rb_raise(rb_eArgError, "d (4th argument) must be NArray");
  if (NA_RANK(rb_d) != 1)
    rb_raise(rb_eArgError, "rank of d (4th argument) must be %d", 1);
  n = NA_SHAPE0(rb_d);
  if (NA_TYPE(rb_d) != NA_DFLOAT)
    rb_d = na_change_type(rb_d, NA_DFLOAT);
  d = NA_PTR_TYPE(rb_d, doublereal*);
  if (!NA_IsNArray(rb_e))
    rb_raise(rb_eArgError, "e (5th argument) must be NArray");
  if (NA_RANK(rb_e) != 1)
    rb_raise(rb_eArgError, "rank of e (5th argument) must be %d", 1);
  if (NA_SHAPE0(rb_e) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of e must be the same as shape 0 of d");
  if (NA_TYPE(rb_e) != NA_DFLOAT)
    rb_e = na_change_type(rb_e, NA_DFLOAT);
  e = NA_PTR_TYPE(rb_e, doublereal*);
  vu = NUM2DBL(rb_vu);
  pivmin = NUM2DBL(rb_pivmin);

  dlarrc_(&jobt, &n, &vl, &vu, d, e, &pivmin, &eigcnt, &lcnt, &rcnt, &info);

  rb_eigcnt = INT2NUM(eigcnt);
  rb_lcnt = INT2NUM(lcnt);
  rb_rcnt = INT2NUM(rcnt);
  rb_info = INT2NUM(info);
  return rb_ary_new3(4, rb_eigcnt, rb_lcnt, rb_rcnt, rb_info);
}

void
init_lapack_dlarrc(VALUE mLapack){
  rb_define_module_function(mLapack, "dlarrc", rb_dlarrc, -1);
}

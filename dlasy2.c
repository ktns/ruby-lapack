#include "rb_lapack.h"

extern VOID dlasy2_(logical *ltranl, logical *ltranr, integer *isgn, integer *n1, integer *n2, doublereal *tl, integer *ldtl, doublereal *tr, integer *ldtr, doublereal *b, integer *ldb, doublereal *scale, doublereal *x, integer *ldx, doublereal *xnorm, integer *info);

static VALUE
rb_dlasy2(int argc, VALUE *argv, VALUE self){
  VALUE rb_ltranl;
  logical ltranl; 
  VALUE rb_ltranr;
  logical ltranr; 
  VALUE rb_isgn;
  integer isgn; 
  VALUE rb_n1;
  integer n1; 
  VALUE rb_n2;
  integer n2; 
  VALUE rb_tl;
  doublereal *tl; 
  VALUE rb_tr;
  doublereal *tr; 
  VALUE rb_b;
  doublereal *b; 
  VALUE rb_scale;
  doublereal scale; 
  VALUE rb_x;
  doublereal *x; 
  VALUE rb_xnorm;
  doublereal xnorm; 
  VALUE rb_info;
  integer info; 

  integer ldtl;
  integer ldtr;
  integer ldb;
  integer ldx;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  scale, x, xnorm, info = NumRu::Lapack.dlasy2( ltranl, ltranr, isgn, n1, n2, tl, tr, b)\n    or\n  NumRu::Lapack.dlasy2  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 8)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 8)", argc);
  rb_ltranl = argv[0];
  rb_ltranr = argv[1];
  rb_isgn = argv[2];
  rb_n1 = argv[3];
  rb_n2 = argv[4];
  rb_tl = argv[5];
  rb_tr = argv[6];
  rb_b = argv[7];

  ltranl = (rb_ltranl == Qtrue);
  if (!NA_IsNArray(rb_b))
    rb_raise(rb_eArgError, "b (8th argument) must be NArray");
  if (NA_RANK(rb_b) != 2)
    rb_raise(rb_eArgError, "rank of b (8th argument) must be %d", 2);
  if (NA_SHAPE1(rb_b) != (2))
    rb_raise(rb_eRuntimeError, "shape 1 of b must be %d", 2);
  ldb = NA_SHAPE0(rb_b);
  if (NA_TYPE(rb_b) != NA_DFLOAT)
    rb_b = na_change_type(rb_b, NA_DFLOAT);
  b = NA_PTR_TYPE(rb_b, doublereal*);
  if (!NA_IsNArray(rb_tl))
    rb_raise(rb_eArgError, "tl (6th argument) must be NArray");
  if (NA_RANK(rb_tl) != 2)
    rb_raise(rb_eArgError, "rank of tl (6th argument) must be %d", 2);
  if (NA_SHAPE1(rb_tl) != (2))
    rb_raise(rb_eRuntimeError, "shape 1 of tl must be %d", 2);
  ldtl = NA_SHAPE0(rb_tl);
  if (NA_TYPE(rb_tl) != NA_DFLOAT)
    rb_tl = na_change_type(rb_tl, NA_DFLOAT);
  tl = NA_PTR_TYPE(rb_tl, doublereal*);
  n1 = NUM2INT(rb_n1);
  isgn = NUM2INT(rb_isgn);
  ltranr = (rb_ltranr == Qtrue);
  if (!NA_IsNArray(rb_tr))
    rb_raise(rb_eArgError, "tr (7th argument) must be NArray");
  if (NA_RANK(rb_tr) != 2)
    rb_raise(rb_eArgError, "rank of tr (7th argument) must be %d", 2);
  if (NA_SHAPE1(rb_tr) != (2))
    rb_raise(rb_eRuntimeError, "shape 1 of tr must be %d", 2);
  ldtr = NA_SHAPE0(rb_tr);
  if (NA_TYPE(rb_tr) != NA_DFLOAT)
    rb_tr = na_change_type(rb_tr, NA_DFLOAT);
  tr = NA_PTR_TYPE(rb_tr, doublereal*);
  n2 = NUM2INT(rb_n2);
  ldx = MAX(1,n1);
  {
    int shape[2];
    shape[0] = ldx;
    shape[1] = 2;
    rb_x = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  x = NA_PTR_TYPE(rb_x, doublereal*);

  dlasy2_(&ltranl, &ltranr, &isgn, &n1, &n2, tl, &ldtl, tr, &ldtr, b, &ldb, &scale, x, &ldx, &xnorm, &info);

  rb_scale = rb_float_new((double)scale);
  rb_xnorm = rb_float_new((double)xnorm);
  rb_info = INT2NUM(info);
  return rb_ary_new3(4, rb_scale, rb_x, rb_xnorm, rb_info);
}

void
init_lapack_dlasy2(VALUE mLapack){
  rb_define_module_function(mLapack, "dlasy2", rb_dlasy2, -1);
}

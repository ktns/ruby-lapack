#include "rb_lapack.h"

extern VOID dstebz_(char *range, char *order, integer *n, doublereal *vl, doublereal *vu, integer *il, integer *iu, doublereal *abstol, doublereal *d, doublereal *e, integer *m, integer *nsplit, doublereal *w, integer *iblock, integer *isplit, doublereal *work, integer *iwork, integer *info);

static VALUE
rb_dstebz(int argc, VALUE *argv, VALUE self){
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
  VALUE rb_abstol;
  doublereal abstol; 
  VALUE rb_d;
  doublereal *d; 
  VALUE rb_e;
  doublereal *e; 
  VALUE rb_m;
  integer m; 
  VALUE rb_nsplit;
  integer nsplit; 
  VALUE rb_w;
  doublereal *w; 
  VALUE rb_iblock;
  integer *iblock; 
  VALUE rb_isplit;
  integer *isplit; 
  VALUE rb_info;
  integer info; 
  doublereal *work;
  integer *iwork;

  integer n;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  m, nsplit, w, iblock, isplit, info = NumRu::Lapack.dstebz( range, order, vl, vu, il, iu, abstol, d, e)\n    or\n  NumRu::Lapack.dstebz  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 9)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 9)", argc);
  rb_range = argv[0];
  rb_order = argv[1];
  rb_vl = argv[2];
  rb_vu = argv[3];
  rb_il = argv[4];
  rb_iu = argv[5];
  rb_abstol = argv[6];
  rb_d = argv[7];
  rb_e = argv[8];

  abstol = NUM2DBL(rb_abstol);
  vl = NUM2DBL(rb_vl);
  iu = NUM2INT(rb_iu);
  il = NUM2INT(rb_il);
  range = StringValueCStr(rb_range)[0];
  vu = NUM2DBL(rb_vu);
  order = StringValueCStr(rb_order)[0];
  if (!NA_IsNArray(rb_d))
    rb_raise(rb_eArgError, "d (8th argument) must be NArray");
  if (NA_RANK(rb_d) != 1)
    rb_raise(rb_eArgError, "rank of d (8th argument) must be %d", 1);
  n = NA_SHAPE0(rb_d);
  if (NA_TYPE(rb_d) != NA_DFLOAT)
    rb_d = na_change_type(rb_d, NA_DFLOAT);
  d = NA_PTR_TYPE(rb_d, doublereal*);
  if (!NA_IsNArray(rb_e))
    rb_raise(rb_eArgError, "e (9th argument) must be NArray");
  if (NA_RANK(rb_e) != 1)
    rb_raise(rb_eArgError, "rank of e (9th argument) must be %d", 1);
  if (NA_SHAPE0(rb_e) != (n-1))
    rb_raise(rb_eRuntimeError, "shape 0 of e must be %d", n-1);
  if (NA_TYPE(rb_e) != NA_DFLOAT)
    rb_e = na_change_type(rb_e, NA_DFLOAT);
  e = NA_PTR_TYPE(rb_e, doublereal*);
  {
    int shape[1];
    shape[0] = n;
    rb_w = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  w = NA_PTR_TYPE(rb_w, doublereal*);
  {
    int shape[1];
    shape[0] = n;
    rb_iblock = na_make_object(NA_LINT, 1, shape, cNArray);
  }
  iblock = NA_PTR_TYPE(rb_iblock, integer*);
  {
    int shape[1];
    shape[0] = n;
    rb_isplit = na_make_object(NA_LINT, 1, shape, cNArray);
  }
  isplit = NA_PTR_TYPE(rb_isplit, integer*);
  work = ALLOC_N(doublereal, (4*n));
  iwork = ALLOC_N(integer, (3*n));

  dstebz_(&range, &order, &n, &vl, &vu, &il, &iu, &abstol, d, e, &m, &nsplit, w, iblock, isplit, work, iwork, &info);

  free(work);
  free(iwork);
  rb_m = INT2NUM(m);
  rb_nsplit = INT2NUM(nsplit);
  rb_info = INT2NUM(info);
  return rb_ary_new3(6, rb_m, rb_nsplit, rb_w, rb_iblock, rb_isplit, rb_info);
}

void
init_lapack_dstebz(VALUE mLapack){
  rb_define_module_function(mLapack, "dstebz", rb_dstebz, -1);
}

#include "rb_lapack.h"

extern VOID dlarra_(integer *n, doublereal *d, doublereal *e, doublereal *e2, doublereal *spltol, doublereal *tnrm, integer *nsplit, integer *isplit, integer *info);

static VALUE
rb_dlarra(int argc, VALUE *argv, VALUE self){
  VALUE rb_d;
  doublereal *d; 
  VALUE rb_e;
  doublereal *e; 
  VALUE rb_e2;
  doublereal *e2; 
  VALUE rb_spltol;
  doublereal spltol; 
  VALUE rb_tnrm;
  doublereal tnrm; 
  VALUE rb_nsplit;
  integer nsplit; 
  VALUE rb_isplit;
  integer *isplit; 
  VALUE rb_info;
  integer info; 
  VALUE rb_e_out__;
  doublereal *e_out__;
  VALUE rb_e2_out__;
  doublereal *e2_out__;

  integer n;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  nsplit, isplit, info, e, e2 = NumRu::Lapack.dlarra( d, e, e2, spltol, tnrm)\n    or\n  NumRu::Lapack.dlarra  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 5)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 5)", argc);
  rb_d = argv[0];
  rb_e = argv[1];
  rb_e2 = argv[2];
  rb_spltol = argv[3];
  rb_tnrm = argv[4];

  tnrm = NUM2DBL(rb_tnrm);
  if (!NA_IsNArray(rb_e2))
    rb_raise(rb_eArgError, "e2 (3th argument) must be NArray");
  if (NA_RANK(rb_e2) != 1)
    rb_raise(rb_eArgError, "rank of e2 (3th argument) must be %d", 1);
  n = NA_SHAPE0(rb_e2);
  if (NA_TYPE(rb_e2) != NA_DFLOAT)
    rb_e2 = na_change_type(rb_e2, NA_DFLOAT);
  e2 = NA_PTR_TYPE(rb_e2, doublereal*);
  spltol = NUM2DBL(rb_spltol);
  if (!NA_IsNArray(rb_d))
    rb_raise(rb_eArgError, "d (1th argument) must be NArray");
  if (NA_RANK(rb_d) != 1)
    rb_raise(rb_eArgError, "rank of d (1th argument) must be %d", 1);
  if (NA_SHAPE0(rb_d) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of d must be the same as shape 0 of e2");
  if (NA_TYPE(rb_d) != NA_DFLOAT)
    rb_d = na_change_type(rb_d, NA_DFLOAT);
  d = NA_PTR_TYPE(rb_d, doublereal*);
  if (!NA_IsNArray(rb_e))
    rb_raise(rb_eArgError, "e (2th argument) must be NArray");
  if (NA_RANK(rb_e) != 1)
    rb_raise(rb_eArgError, "rank of e (2th argument) must be %d", 1);
  if (NA_SHAPE0(rb_e) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of e must be the same as shape 0 of e2");
  if (NA_TYPE(rb_e) != NA_DFLOAT)
    rb_e = na_change_type(rb_e, NA_DFLOAT);
  e = NA_PTR_TYPE(rb_e, doublereal*);
  {
    int shape[1];
    shape[0] = n;
    rb_isplit = na_make_object(NA_LINT, 1, shape, cNArray);
  }
  isplit = NA_PTR_TYPE(rb_isplit, integer*);
  {
    int shape[1];
    shape[0] = n;
    rb_e_out__ = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  e_out__ = NA_PTR_TYPE(rb_e_out__, doublereal*);
  MEMCPY(e_out__, e, doublereal, NA_TOTAL(rb_e));
  rb_e = rb_e_out__;
  e = e_out__;
  {
    int shape[1];
    shape[0] = n;
    rb_e2_out__ = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  e2_out__ = NA_PTR_TYPE(rb_e2_out__, doublereal*);
  MEMCPY(e2_out__, e2, doublereal, NA_TOTAL(rb_e2));
  rb_e2 = rb_e2_out__;
  e2 = e2_out__;

  dlarra_(&n, d, e, e2, &spltol, &tnrm, &nsplit, isplit, &info);

  rb_nsplit = INT2NUM(nsplit);
  rb_info = INT2NUM(info);
  return rb_ary_new3(5, rb_nsplit, rb_isplit, rb_info, rb_e, rb_e2);
}

void
init_lapack_dlarra(VALUE mLapack){
  rb_define_module_function(mLapack, "dlarra", rb_dlarra, -1);
}

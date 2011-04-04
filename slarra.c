#include "rb_lapack.h"

extern VOID slarra_(integer *n, real *d, real *e, real *e2, real *spltol, real *tnrm, integer *nsplit, integer *isplit, integer *info);

static VALUE
rb_slarra(int argc, VALUE *argv, VALUE self){
  VALUE rb_d;
  real *d; 
  VALUE rb_e;
  real *e; 
  VALUE rb_e2;
  real *e2; 
  VALUE rb_spltol;
  real spltol; 
  VALUE rb_tnrm;
  real tnrm; 
  VALUE rb_nsplit;
  integer nsplit; 
  VALUE rb_isplit;
  integer *isplit; 
  VALUE rb_info;
  integer info; 
  VALUE rb_e_out__;
  real *e_out__;
  VALUE rb_e2_out__;
  real *e2_out__;

  integer n;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  nsplit, isplit, info, e, e2 = NumRu::Lapack.slarra( d, e, e2, spltol, tnrm)\n    or\n  NumRu::Lapack.slarra  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 5)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 5)", argc);
  rb_d = argv[0];
  rb_e = argv[1];
  rb_e2 = argv[2];
  rb_spltol = argv[3];
  rb_tnrm = argv[4];

  tnrm = (real)NUM2DBL(rb_tnrm);
  if (!NA_IsNArray(rb_e2))
    rb_raise(rb_eArgError, "e2 (3th argument) must be NArray");
  if (NA_RANK(rb_e2) != 1)
    rb_raise(rb_eArgError, "rank of e2 (3th argument) must be %d", 1);
  n = NA_SHAPE0(rb_e2);
  if (NA_TYPE(rb_e2) != NA_SFLOAT)
    rb_e2 = na_change_type(rb_e2, NA_SFLOAT);
  e2 = NA_PTR_TYPE(rb_e2, real*);
  spltol = (real)NUM2DBL(rb_spltol);
  if (!NA_IsNArray(rb_d))
    rb_raise(rb_eArgError, "d (1th argument) must be NArray");
  if (NA_RANK(rb_d) != 1)
    rb_raise(rb_eArgError, "rank of d (1th argument) must be %d", 1);
  if (NA_SHAPE0(rb_d) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of d must be the same as shape 0 of e2");
  if (NA_TYPE(rb_d) != NA_SFLOAT)
    rb_d = na_change_type(rb_d, NA_SFLOAT);
  d = NA_PTR_TYPE(rb_d, real*);
  if (!NA_IsNArray(rb_e))
    rb_raise(rb_eArgError, "e (2th argument) must be NArray");
  if (NA_RANK(rb_e) != 1)
    rb_raise(rb_eArgError, "rank of e (2th argument) must be %d", 1);
  if (NA_SHAPE0(rb_e) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of e must be the same as shape 0 of e2");
  if (NA_TYPE(rb_e) != NA_SFLOAT)
    rb_e = na_change_type(rb_e, NA_SFLOAT);
  e = NA_PTR_TYPE(rb_e, real*);
  {
    int shape[1];
    shape[0] = n;
    rb_isplit = na_make_object(NA_LINT, 1, shape, cNArray);
  }
  isplit = NA_PTR_TYPE(rb_isplit, integer*);
  {
    int shape[1];
    shape[0] = n;
    rb_e_out__ = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  e_out__ = NA_PTR_TYPE(rb_e_out__, real*);
  MEMCPY(e_out__, e, real, NA_TOTAL(rb_e));
  rb_e = rb_e_out__;
  e = e_out__;
  {
    int shape[1];
    shape[0] = n;
    rb_e2_out__ = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  e2_out__ = NA_PTR_TYPE(rb_e2_out__, real*);
  MEMCPY(e2_out__, e2, real, NA_TOTAL(rb_e2));
  rb_e2 = rb_e2_out__;
  e2 = e2_out__;

  slarra_(&n, d, e, e2, &spltol, &tnrm, &nsplit, isplit, &info);

  rb_nsplit = INT2NUM(nsplit);
  rb_info = INT2NUM(info);
  return rb_ary_new3(5, rb_nsplit, rb_isplit, rb_info, rb_e, rb_e2);
}

void
init_lapack_slarra(VALUE mLapack){
  rb_define_module_function(mLapack, "slarra", rb_slarra, -1);
}

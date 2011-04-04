#include "rb_lapack.h"

extern VOID dla_lin_berr_(integer *n, integer *nz, integer *nrhs, doublereal *res, doublereal *ayb, doublereal *berr);

static VALUE
rb_dla_lin_berr(int argc, VALUE *argv, VALUE self){
  VALUE rb_nz;
  integer nz; 
  VALUE rb_res;
  doublereal *res; 
  VALUE rb_ayb;
  doublereal *ayb; 
  VALUE rb_berr;
  doublereal *berr; 

  integer n;
  integer nrhs;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  berr = NumRu::Lapack.dla_lin_berr( nz, res, ayb)\n    or\n  NumRu::Lapack.dla_lin_berr  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 3)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 3)", argc);
  rb_nz = argv[0];
  rb_res = argv[1];
  rb_ayb = argv[2];

  if (!NA_IsNArray(rb_res))
    rb_raise(rb_eArgError, "res (2th argument) must be NArray");
  if (NA_RANK(rb_res) != 2)
    rb_raise(rb_eArgError, "rank of res (2th argument) must be %d", 2);
  nrhs = NA_SHAPE1(rb_res);
  n = NA_SHAPE0(rb_res);
  if (NA_TYPE(rb_res) != NA_DFLOAT)
    rb_res = na_change_type(rb_res, NA_DFLOAT);
  res = NA_PTR_TYPE(rb_res, doublereal*);
  nz = NUM2INT(rb_nz);
  if (!NA_IsNArray(rb_ayb))
    rb_raise(rb_eArgError, "ayb (3th argument) must be NArray");
  if (NA_RANK(rb_ayb) != 2)
    rb_raise(rb_eArgError, "rank of ayb (3th argument) must be %d", 2);
  if (NA_SHAPE1(rb_ayb) != nrhs)
    rb_raise(rb_eRuntimeError, "shape 1 of ayb must be the same as shape 1 of res");
  if (NA_SHAPE0(rb_ayb) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of ayb must be the same as shape 0 of res");
  if (NA_TYPE(rb_ayb) != NA_DFLOAT)
    rb_ayb = na_change_type(rb_ayb, NA_DFLOAT);
  ayb = NA_PTR_TYPE(rb_ayb, doublereal*);
  {
    int shape[1];
    shape[0] = nrhs;
    rb_berr = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  berr = NA_PTR_TYPE(rb_berr, doublereal*);

  dla_lin_berr_(&n, &nz, &nrhs, res, ayb, berr);

  return rb_berr;
}

void
init_lapack_dla_lin_berr(VALUE mLapack){
  rb_define_module_function(mLapack, "dla_lin_berr", rb_dla_lin_berr, -1);
}

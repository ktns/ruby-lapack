#include "rb_lapack.h"

extern VOID dlaed5_(integer *i, doublereal *d, doublereal *z, doublereal *delta, doublereal *rho, doublereal *dlam);

static VALUE
rb_dlaed5(int argc, VALUE *argv, VALUE self){
  VALUE rb_i;
  integer i; 
  VALUE rb_d;
  doublereal *d; 
  VALUE rb_z;
  doublereal *z; 
  VALUE rb_rho;
  doublereal rho; 
  VALUE rb_delta;
  doublereal *delta; 
  VALUE rb_dlam;
  doublereal dlam; 


  if (argc == 0) {
    printf("%s\n", "USAGE:\n  delta, dlam = NumRu::Lapack.dlaed5( i, d, z, rho)\n    or\n  NumRu::Lapack.dlaed5  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 4)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 4)", argc);
  rb_i = argv[0];
  rb_d = argv[1];
  rb_z = argv[2];
  rb_rho = argv[3];

  rho = NUM2DBL(rb_rho);
  if (!NA_IsNArray(rb_z))
    rb_raise(rb_eArgError, "z (3th argument) must be NArray");
  if (NA_RANK(rb_z) != 1)
    rb_raise(rb_eArgError, "rank of z (3th argument) must be %d", 1);
  if (NA_SHAPE0(rb_z) != (2))
    rb_raise(rb_eRuntimeError, "shape 0 of z must be %d", 2);
  if (NA_TYPE(rb_z) != NA_DFLOAT)
    rb_z = na_change_type(rb_z, NA_DFLOAT);
  z = NA_PTR_TYPE(rb_z, doublereal*);
  if (!NA_IsNArray(rb_d))
    rb_raise(rb_eArgError, "d (2th argument) must be NArray");
  if (NA_RANK(rb_d) != 1)
    rb_raise(rb_eArgError, "rank of d (2th argument) must be %d", 1);
  if (NA_SHAPE0(rb_d) != (2))
    rb_raise(rb_eRuntimeError, "shape 0 of d must be %d", 2);
  if (NA_TYPE(rb_d) != NA_DFLOAT)
    rb_d = na_change_type(rb_d, NA_DFLOAT);
  d = NA_PTR_TYPE(rb_d, doublereal*);
  i = NUM2INT(rb_i);
  {
    int shape[1];
    shape[0] = 2;
    rb_delta = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  delta = NA_PTR_TYPE(rb_delta, doublereal*);

  dlaed5_(&i, d, z, delta, &rho, &dlam);

  rb_dlam = rb_float_new((double)dlam);
  return rb_ary_new3(2, rb_delta, rb_dlam);
}

void
init_lapack_dlaed5(VALUE mLapack){
  rb_define_module_function(mLapack, "dlaed5", rb_dlaed5, -1);
}

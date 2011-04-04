#include "rb_lapack.h"

extern VOID dlaed4_(integer *n, integer *i, doublereal *d, doublereal *z, doublereal *delta, doublereal *rho, doublereal *dlam, integer *info);

static VALUE
rb_dlaed4(int argc, VALUE *argv, VALUE self){
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
  VALUE rb_info;
  integer info; 

  integer n;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  delta, dlam, info = NumRu::Lapack.dlaed4( i, d, z, rho)\n    or\n  NumRu::Lapack.dlaed4  # print help\n\n\nFORTRAN MANUAL\n\n");
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
  n = NA_SHAPE0(rb_z);
  if (NA_TYPE(rb_z) != NA_DFLOAT)
    rb_z = na_change_type(rb_z, NA_DFLOAT);
  z = NA_PTR_TYPE(rb_z, doublereal*);
  if (!NA_IsNArray(rb_d))
    rb_raise(rb_eArgError, "d (2th argument) must be NArray");
  if (NA_RANK(rb_d) != 1)
    rb_raise(rb_eArgError, "rank of d (2th argument) must be %d", 1);
  if (NA_SHAPE0(rb_d) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of d must be the same as shape 0 of z");
  if (NA_TYPE(rb_d) != NA_DFLOAT)
    rb_d = na_change_type(rb_d, NA_DFLOAT);
  d = NA_PTR_TYPE(rb_d, doublereal*);
  i = NUM2INT(rb_i);
  {
    int shape[1];
    shape[0] = n;
    rb_delta = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  delta = NA_PTR_TYPE(rb_delta, doublereal*);

  dlaed4_(&n, &i, d, z, delta, &rho, &dlam, &info);

  rb_dlam = rb_float_new((double)dlam);
  rb_info = INT2NUM(info);
  return rb_ary_new3(3, rb_delta, rb_dlam, rb_info);
}

void
init_lapack_dlaed4(VALUE mLapack){
  rb_define_module_function(mLapack, "dlaed4", rb_dlaed4, -1);
}

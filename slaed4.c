#include "rb_lapack.h"

extern VOID slaed4_(integer *n, integer *i, real *d, real *z, real *delta, real *rho, real *dlam, integer *info);

static VALUE
rb_slaed4(int argc, VALUE *argv, VALUE self){
  VALUE rb_i;
  integer i; 
  VALUE rb_d;
  real *d; 
  VALUE rb_z;
  real *z; 
  VALUE rb_rho;
  real rho; 
  VALUE rb_delta;
  real *delta; 
  VALUE rb_dlam;
  real dlam; 
  VALUE rb_info;
  integer info; 

  integer n;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  delta, dlam, info = NumRu::Lapack.slaed4( i, d, z, rho)\n    or\n  NumRu::Lapack.slaed4  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 4)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 4)", argc);
  rb_i = argv[0];
  rb_d = argv[1];
  rb_z = argv[2];
  rb_rho = argv[3];

  rho = (real)NUM2DBL(rb_rho);
  if (!NA_IsNArray(rb_z))
    rb_raise(rb_eArgError, "z (3th argument) must be NArray");
  if (NA_RANK(rb_z) != 1)
    rb_raise(rb_eArgError, "rank of z (3th argument) must be %d", 1);
  n = NA_SHAPE0(rb_z);
  if (NA_TYPE(rb_z) != NA_SFLOAT)
    rb_z = na_change_type(rb_z, NA_SFLOAT);
  z = NA_PTR_TYPE(rb_z, real*);
  if (!NA_IsNArray(rb_d))
    rb_raise(rb_eArgError, "d (2th argument) must be NArray");
  if (NA_RANK(rb_d) != 1)
    rb_raise(rb_eArgError, "rank of d (2th argument) must be %d", 1);
  if (NA_SHAPE0(rb_d) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of d must be the same as shape 0 of z");
  if (NA_TYPE(rb_d) != NA_SFLOAT)
    rb_d = na_change_type(rb_d, NA_SFLOAT);
  d = NA_PTR_TYPE(rb_d, real*);
  i = NUM2INT(rb_i);
  {
    int shape[1];
    shape[0] = n;
    rb_delta = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  delta = NA_PTR_TYPE(rb_delta, real*);

  slaed4_(&n, &i, d, z, delta, &rho, &dlam, &info);

  rb_dlam = rb_float_new((double)dlam);
  rb_info = INT2NUM(info);
  return rb_ary_new3(3, rb_delta, rb_dlam, rb_info);
}

void
init_lapack_slaed4(VALUE mLapack){
  rb_define_module_function(mLapack, "slaed4", rb_slaed4, -1);
}

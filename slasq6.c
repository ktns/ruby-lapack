#include "rb_lapack.h"

extern VOID slasq6_(integer *i0, integer *n0, real *z, integer *pp, real *dmin, real *dmin1, real *dmin2, real *dn, real *dnm1, real *dnm2);

static VALUE
rb_slasq6(int argc, VALUE *argv, VALUE self){
  VALUE rb_i0;
  integer i0; 
  VALUE rb_n0;
  integer n0; 
  VALUE rb_z;
  real *z; 
  VALUE rb_pp;
  integer pp; 
  VALUE rb_dmin;
  real dmin; 
  VALUE rb_dmin1;
  real dmin1; 
  VALUE rb_dmin2;
  real dmin2; 
  VALUE rb_dn;
  real dn; 
  VALUE rb_dnm1;
  real dnm1; 
  VALUE rb_dnm2;
  real dnm2; 


  if (argc == 0) {
    printf("%s\n", "USAGE:\n  dmin, dmin1, dmin2, dn, dnm1, dnm2 = NumRu::Lapack.slasq6( i0, n0, z, pp)\n    or\n  NumRu::Lapack.slasq6  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 4)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 4)", argc);
  rb_i0 = argv[0];
  rb_n0 = argv[1];
  rb_z = argv[2];
  rb_pp = argv[3];

  pp = NUM2INT(rb_pp);
  n0 = NUM2INT(rb_n0);
  i0 = NUM2INT(rb_i0);
  if (!NA_IsNArray(rb_z))
    rb_raise(rb_eArgError, "z (3th argument) must be NArray");
  if (NA_RANK(rb_z) != 1)
    rb_raise(rb_eArgError, "rank of z (3th argument) must be %d", 1);
  if (NA_SHAPE0(rb_z) != (4*n0))
    rb_raise(rb_eRuntimeError, "shape 0 of z must be %d", 4*n0);
  if (NA_TYPE(rb_z) != NA_SFLOAT)
    rb_z = na_change_type(rb_z, NA_SFLOAT);
  z = NA_PTR_TYPE(rb_z, real*);

  slasq6_(&i0, &n0, z, &pp, &dmin, &dmin1, &dmin2, &dn, &dnm1, &dnm2);

  rb_dmin = rb_float_new((double)dmin);
  rb_dmin1 = rb_float_new((double)dmin1);
  rb_dmin2 = rb_float_new((double)dmin2);
  rb_dn = rb_float_new((double)dn);
  rb_dnm1 = rb_float_new((double)dnm1);
  rb_dnm2 = rb_float_new((double)dnm2);
  return rb_ary_new3(6, rb_dmin, rb_dmin1, rb_dmin2, rb_dn, rb_dnm1, rb_dnm2);
}

void
init_lapack_slasq6(VALUE mLapack){
  rb_define_module_function(mLapack, "slasq6", rb_slasq6, -1);
}

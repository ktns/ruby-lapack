#include "rb_lapack.h"

extern VOID dlaed6_(integer *kniter, logical *orgati, doublereal *rho, doublereal *d, doublereal *z, doublereal *finit, doublereal *tau, integer *info);

static VALUE
rb_dlaed6(int argc, VALUE *argv, VALUE self){
  VALUE rb_kniter;
  integer kniter; 
  VALUE rb_orgati;
  logical orgati; 
  VALUE rb_rho;
  doublereal rho; 
  VALUE rb_d;
  doublereal *d; 
  VALUE rb_z;
  doublereal *z; 
  VALUE rb_finit;
  doublereal finit; 
  VALUE rb_tau;
  doublereal tau; 
  VALUE rb_info;
  integer info; 


  if (argc == 0) {
    printf("%s\n", "USAGE:\n  tau, info = NumRu::Lapack.dlaed6( kniter, orgati, rho, d, z, finit)\n    or\n  NumRu::Lapack.dlaed6  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 6)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 6)", argc);
  rb_kniter = argv[0];
  rb_orgati = argv[1];
  rb_rho = argv[2];
  rb_d = argv[3];
  rb_z = argv[4];
  rb_finit = argv[5];

  orgati = (rb_orgati == Qtrue);
  finit = NUM2DBL(rb_finit);
  rho = NUM2DBL(rb_rho);
  if (!NA_IsNArray(rb_z))
    rb_raise(rb_eArgError, "z (5th argument) must be NArray");
  if (NA_RANK(rb_z) != 1)
    rb_raise(rb_eArgError, "rank of z (5th argument) must be %d", 1);
  if (NA_SHAPE0(rb_z) != (3))
    rb_raise(rb_eRuntimeError, "shape 0 of z must be %d", 3);
  if (NA_TYPE(rb_z) != NA_DFLOAT)
    rb_z = na_change_type(rb_z, NA_DFLOAT);
  z = NA_PTR_TYPE(rb_z, doublereal*);
  if (!NA_IsNArray(rb_d))
    rb_raise(rb_eArgError, "d (4th argument) must be NArray");
  if (NA_RANK(rb_d) != 1)
    rb_raise(rb_eArgError, "rank of d (4th argument) must be %d", 1);
  if (NA_SHAPE0(rb_d) != (3))
    rb_raise(rb_eRuntimeError, "shape 0 of d must be %d", 3);
  if (NA_TYPE(rb_d) != NA_DFLOAT)
    rb_d = na_change_type(rb_d, NA_DFLOAT);
  d = NA_PTR_TYPE(rb_d, doublereal*);
  kniter = NUM2INT(rb_kniter);

  dlaed6_(&kniter, &orgati, &rho, d, z, &finit, &tau, &info);

  rb_tau = rb_float_new((double)tau);
  rb_info = INT2NUM(info);
  return rb_ary_new3(2, rb_tau, rb_info);
}

void
init_lapack_dlaed6(VALUE mLapack){
  rb_define_module_function(mLapack, "dlaed6", rb_dlaed6, -1);
}

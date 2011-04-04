#include "rb_lapack.h"

extern doublereal dlapy3_(doublereal *x, doublereal *y, doublereal *z);

static VALUE
rb_dlapy3(int argc, VALUE *argv, VALUE self){
  VALUE rb_x;
  doublereal x; 
  VALUE rb_y;
  doublereal y; 
  VALUE rb_z;
  doublereal z; 
  VALUE rb___out__;
  doublereal __out__; 


  if (argc == 0) {
    printf("%s\n", "USAGE:\n  __out__ = NumRu::Lapack.dlapy3( x, y, z)\n    or\n  NumRu::Lapack.dlapy3  # print help\n\n\nFORTRAN MANUAL\n      DOUBLE PRECISION FUNCTION DLAPY3( X, Y, Z )\n\n*  Purpose\n*  =======\n*\n*  DLAPY3 returns sqrt(x**2+y**2+z**2), taking care not to cause\n*  unnecessary overflow.\n*\n\n*  Arguments\n*  =========\n*\n*  X       (input) DOUBLE PRECISION\n*  Y       (input) DOUBLE PRECISION\n*  Z       (input) DOUBLE PRECISION\n*          X, Y and Z specify the values x, y and z.\n*\n\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 3)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 3)", argc);
  rb_x = argv[0];
  rb_y = argv[1];
  rb_z = argv[2];

  x = NUM2DBL(rb_x);
  y = NUM2DBL(rb_y);
  z = NUM2DBL(rb_z);

  __out__ = dlapy3_(&x, &y, &z);

  rb___out__ = rb_float_new((double)__out__);
  return rb___out__;
}

void
init_lapack_dlapy3(VALUE mLapack){
  rb_define_module_function(mLapack, "dlapy3", rb_dlapy3, -1);
}

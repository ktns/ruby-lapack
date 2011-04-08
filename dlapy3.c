#include "rb_lapack.h"

extern doublereal dlapy3_(doublereal *x, doublereal *y, doublereal *z);

static VALUE sHelp, sUsage;

static VALUE
rblapack_dlapy3(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_x;
  doublereal x; 
  VALUE rblapack_y;
  doublereal y; 
  VALUE rblapack_z;
  doublereal z; 
  VALUE rblapack___out__;
  doublereal __out__; 


  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  __out__ = NumRu::Lapack.dlapy3( x, y, z, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      DOUBLE PRECISION FUNCTION DLAPY3( X, Y, Z )\n\n*  Purpose\n*  =======\n*\n*  DLAPY3 returns sqrt(x**2+y**2+z**2), taking care not to cause\n*  unnecessary overflow.\n*\n\n*  Arguments\n*  =========\n*\n*  X       (input) DOUBLE PRECISION\n*  Y       (input) DOUBLE PRECISION\n*  Z       (input) DOUBLE PRECISION\n*          X, Y and Z specify the values x, y and z.\n*\n\n*  =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  __out__ = NumRu::Lapack.dlapy3( x, y, z, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 3)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 3)", argc);
  rblapack_x = argv[0];
  rblapack_y = argv[1];
  rblapack_z = argv[2];
  if (rb_options != Qnil) {
  }

  x = NUM2DBL(rblapack_x);
  y = NUM2DBL(rblapack_y);
  z = NUM2DBL(rblapack_z);

  __out__ = dlapy3_(&x, &y, &z);

  rblapack___out__ = rb_float_new((double)__out__);
  return rblapack___out__;
}

void
init_lapack_dlapy3(VALUE mLapack){
  rb_define_module_function(mLapack, "dlapy3", rblapack_dlapy3, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}

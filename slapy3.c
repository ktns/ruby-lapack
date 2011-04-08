#include "rb_lapack.h"

extern real slapy3_(real *x, real *y, real *z);

static VALUE sHelp, sUsage;

static VALUE
rblapack_slapy3(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_x;
  real x; 
  VALUE rblapack_y;
  real y; 
  VALUE rblapack_z;
  real z; 
  VALUE rblapack___out__;
  real __out__; 


  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  __out__ = NumRu::Lapack.slapy3( x, y, z, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      REAL             FUNCTION SLAPY3( X, Y, Z )\n\n*  Purpose\n*  =======\n*\n*  SLAPY3 returns sqrt(x**2+y**2+z**2), taking care not to cause\n*  unnecessary overflow.\n*\n\n*  Arguments\n*  =========\n*\n*  X       (input) REAL\n*  Y       (input) REAL\n*  Z       (input) REAL\n*          X, Y and Z specify the values x, y and z.\n*\n\n*  =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  __out__ = NumRu::Lapack.slapy3( x, y, z, [:usage => usage, :help => help])\n");
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

  x = (real)NUM2DBL(rblapack_x);
  y = (real)NUM2DBL(rblapack_y);
  z = (real)NUM2DBL(rblapack_z);

  __out__ = slapy3_(&x, &y, &z);

  rblapack___out__ = rb_float_new((double)__out__);
  return rblapack___out__;
}

void
init_lapack_slapy3(VALUE mLapack){
  rb_define_module_function(mLapack, "slapy3", rblapack_slapy3, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}

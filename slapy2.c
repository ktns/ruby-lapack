#include "rb_lapack.h"

extern real slapy2_(real *x, real *y);

static VALUE sHelp, sUsage;

static VALUE
rblapack_slapy2(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_x;
  real x; 
  VALUE rblapack_y;
  real y; 
  VALUE rblapack___out__;
  real __out__; 


  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  __out__ = NumRu::Lapack.slapy2( x, y, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      REAL             FUNCTION SLAPY2( X, Y )\n\n*  Purpose\n*  =======\n*\n*  SLAPY2 returns sqrt(x**2+y**2), taking care not to cause unnecessary\n*  overflow.\n*\n\n*  Arguments\n*  =========\n*\n*  X       (input) REAL\n*  Y       (input) REAL\n*          X and Y specify the values x and y.\n*\n\n*  =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  __out__ = NumRu::Lapack.slapy2( x, y, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 2)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 2)", argc);
  rblapack_x = argv[0];
  rblapack_y = argv[1];
  if (rb_options != Qnil) {
  }

  x = (real)NUM2DBL(rblapack_x);
  y = (real)NUM2DBL(rblapack_y);

  __out__ = slapy2_(&x, &y);

  rblapack___out__ = rb_float_new((double)__out__);
  return rblapack___out__;
}

void
init_lapack_slapy2(VALUE mLapack){
  rb_define_module_function(mLapack, "slapy2", rblapack_slapy2, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}

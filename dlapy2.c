#include "rb_lapack.h"

extern doublereal dlapy2_(doublereal *x, doublereal *y);

static VALUE sHelp, sUsage;

static VALUE
rblapack_dlapy2(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_x;
  doublereal x; 
  VALUE rblapack_y;
  doublereal y; 
  VALUE rblapack___out__;
  doublereal __out__; 


  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  __out__ = NumRu::Lapack.dlapy2( x, y, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      DOUBLE PRECISION FUNCTION DLAPY2( X, Y )\n\n*  Purpose\n*  =======\n*\n*  DLAPY2 returns sqrt(x**2+y**2), taking care not to cause unnecessary\n*  overflow.\n*\n\n*  Arguments\n*  =========\n*\n*  X       (input) DOUBLE PRECISION\n*  Y       (input) DOUBLE PRECISION\n*          X and Y specify the values x and y.\n*\n\n*  =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  __out__ = NumRu::Lapack.dlapy2( x, y, [:usage => usage, :help => help])\n");
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

  x = NUM2DBL(rblapack_x);
  y = NUM2DBL(rblapack_y);

  __out__ = dlapy2_(&x, &y);

  rblapack___out__ = rb_float_new((double)__out__);
  return rblapack___out__;
}

void
init_lapack_dlapy2(VALUE mLapack){
  rb_define_module_function(mLapack, "dlapy2", rblapack_dlapy2, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}

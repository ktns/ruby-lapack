#include "rb_lapack.h"

extern integer ilauplo_(char *uplo);

static VALUE sHelp, sUsage;

static VALUE
rblapack_ilauplo(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_uplo;
  char uplo; 
  VALUE rblapack___out__;
  integer __out__; 


  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  __out__ = NumRu::Lapack.ilauplo( uplo, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      INTEGER FUNCTION ILAUPLO( UPLO )\n\n*  Purpose\n*  =======\n*\n*  This subroutine translated from a character string specifying a\n*  upper- or lower-triangular matrix to the relevant BLAST-specified\n*  integer constant.\n*\n*  ILAUPLO returns an INTEGER.  If ILAUPLO < 0, then the input is not\n*  a character indicating an upper- or lower-triangular matrix.\n*  Otherwise ILAUPLO returns the constant value corresponding to UPLO.\n*\n\n*  Arguments\n*  =========\n*  UPLO    (input) CHARACTER\n*          = 'U':  A is upper triangular;\n*          = 'L':  A is lower triangular.\n\n*  =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  __out__ = NumRu::Lapack.ilauplo( uplo, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 1)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 1)", argc);
  rblapack_uplo = argv[0];
  if (rb_options != Qnil) {
  }

  uplo = StringValueCStr(rblapack_uplo)[0];

  __out__ = ilauplo_(&uplo);

  rblapack___out__ = INT2NUM(__out__);
  return rblapack___out__;
}

void
init_lapack_ilauplo(VALUE mLapack){
  rb_define_module_function(mLapack, "ilauplo", rblapack_ilauplo, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}

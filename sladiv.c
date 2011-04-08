#include "rb_lapack.h"

extern VOID sladiv_(real *a, real *b, real *c, real *d, real *p, real *q);

static VALUE sHelp, sUsage;

static VALUE
rblapack_sladiv(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_a;
  real a; 
  VALUE rblapack_b;
  real b; 
  VALUE rblapack_c;
  real c; 
  VALUE rblapack_d;
  real d; 
  VALUE rblapack_p;
  real p; 
  VALUE rblapack_q;
  real q; 


  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  p, q = NumRu::Lapack.sladiv( a, b, c, d, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE SLADIV( A, B, C, D, P, Q )\n\n*  Purpose\n*  =======\n*\n*  SLADIV performs complex division in  real arithmetic\n*\n*                        a + i*b\n*             p + i*q = ---------\n*                        c + i*d\n*\n*  The algorithm is due to Robert L. Smith and can be found\n*  in D. Knuth, The art of Computer Programming, Vol.2, p.195\n*\n\n*  Arguments\n*  =========\n*\n*  A       (input) REAL\n*  B       (input) REAL\n*  C       (input) REAL\n*  D       (input) REAL\n*          The scalars a, b, c, and d in the above expression.\n*\n*  P       (output) REAL\n*  Q       (output) REAL\n*          The scalars p and q in the above expression.\n*\n\n*  =====================================================================\n*\n*     .. Local Scalars ..\n      REAL               E, F\n*     ..\n*     .. Intrinsic Functions ..\n      INTRINSIC          ABS\n*     ..\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  p, q = NumRu::Lapack.sladiv( a, b, c, d, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 4)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 4)", argc);
  rblapack_a = argv[0];
  rblapack_b = argv[1];
  rblapack_c = argv[2];
  rblapack_d = argv[3];
  if (rb_options != Qnil) {
  }

  a = (real)NUM2DBL(rblapack_a);
  b = (real)NUM2DBL(rblapack_b);
  c = (real)NUM2DBL(rblapack_c);
  d = (real)NUM2DBL(rblapack_d);

  sladiv_(&a, &b, &c, &d, &p, &q);

  rblapack_p = rb_float_new((double)p);
  rblapack_q = rb_float_new((double)q);
  return rb_ary_new3(2, rblapack_p, rblapack_q);
}

void
init_lapack_sladiv(VALUE mLapack){
  rb_define_module_function(mLapack, "sladiv", rblapack_sladiv, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}

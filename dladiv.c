#include "rb_lapack.h"

extern VOID dladiv_(doublereal *a, doublereal *b, doublereal *c, doublereal *d, doublereal *p, doublereal *q);

static VALUE sHelp, sUsage;

static VALUE
rblapack_dladiv(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_a;
  doublereal a; 
  VALUE rblapack_b;
  doublereal b; 
  VALUE rblapack_c;
  doublereal c; 
  VALUE rblapack_d;
  doublereal d; 
  VALUE rblapack_p;
  doublereal p; 
  VALUE rblapack_q;
  doublereal q; 


  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  p, q = NumRu::Lapack.dladiv( a, b, c, d, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE DLADIV( A, B, C, D, P, Q )\n\n*  Purpose\n*  =======\n*\n*  DLADIV performs complex division in  real arithmetic\n*\n*                        a + i*b\n*             p + i*q = ---------\n*                        c + i*d\n*\n*  The algorithm is due to Robert L. Smith and can be found\n*  in D. Knuth, The art of Computer Programming, Vol.2, p.195\n*\n\n*  Arguments\n*  =========\n*\n*  A       (input) DOUBLE PRECISION\n*  B       (input) DOUBLE PRECISION\n*  C       (input) DOUBLE PRECISION\n*  D       (input) DOUBLE PRECISION\n*          The scalars a, b, c, and d in the above expression.\n*\n*  P       (output) DOUBLE PRECISION\n*  Q       (output) DOUBLE PRECISION\n*          The scalars p and q in the above expression.\n*\n\n*  =====================================================================\n*\n*     .. Local Scalars ..\n      DOUBLE PRECISION   E, F\n*     ..\n*     .. Intrinsic Functions ..\n      INTRINSIC          ABS\n*     ..\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  p, q = NumRu::Lapack.dladiv( a, b, c, d, [:usage => usage, :help => help])\n");
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

  a = NUM2DBL(rblapack_a);
  b = NUM2DBL(rblapack_b);
  c = NUM2DBL(rblapack_c);
  d = NUM2DBL(rblapack_d);

  dladiv_(&a, &b, &c, &d, &p, &q);

  rblapack_p = rb_float_new((double)p);
  rblapack_q = rb_float_new((double)q);
  return rb_ary_new3(2, rblapack_p, rblapack_q);
}

void
init_lapack_dladiv(VALUE mLapack){
  rb_define_module_function(mLapack, "dladiv", rblapack_dladiv, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}

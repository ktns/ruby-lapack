#include "rb_lapack.h"

static VALUE
rb_sladiv(int argc, VALUE *argv, VALUE self){
  VALUE rb_a;
  real a; 
  VALUE rb_b;
  real b; 
  VALUE rb_c;
  real c; 
  VALUE rb_d;
  real d; 
  VALUE rb_p;
  real p; 
  VALUE rb_q;
  real q; 


  if (argc == 0) {
    printf("%s\n", "USAGE:\n  p, q = NumRu::Lapack.sladiv( a, b, c, d)\n    or\n  NumRu::Lapack.sladiv  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE SLADIV( A, B, C, D, P, Q )\n\n*  Purpose\n*  =======\n*\n*  SLADIV performs complex division in  real arithmetic\n*\n*                        a + i*b\n*             p + i*q = ---------\n*                        c + i*d\n*\n*  The algorithm is due to Robert L. Smith and can be found\n*  in D. Knuth, The art of Computer Programming, Vol.2, p.195\n*\n\n*  Arguments\n*  =========\n*\n*  A       (input) REAL\n*  B       (input) REAL\n*  C       (input) REAL\n*  D       (input) REAL\n*          The scalars a, b, c, and d in the above expression.\n*\n*  P       (output) REAL\n*  Q       (output) REAL\n*          The scalars p and q in the above expression.\n*\n\n*  =====================================================================\n*\n*     .. Local Scalars ..\n      REAL               E, F\n*     ..\n*     .. Intrinsic Functions ..\n      INTRINSIC          ABS\n*     ..\n\n");
    return Qnil;
  }
  if (argc != 4)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 4)", argc);
  rb_a = argv[0];
  rb_b = argv[1];
  rb_c = argv[2];
  rb_d = argv[3];

  a = (real)NUM2DBL(rb_a);
  b = (real)NUM2DBL(rb_b);
  c = (real)NUM2DBL(rb_c);
  d = (real)NUM2DBL(rb_d);

  sladiv_(&a, &b, &c, &d, &p, &q);

  rb_p = rb_float_new((double)p);
  rb_q = rb_float_new((double)q);
  return rb_ary_new3(2, rb_p, rb_q);
}

void
init_lapack_sladiv(VALUE mLapack){
  rb_define_module_function(mLapack, "sladiv", rb_sladiv, -1);
}

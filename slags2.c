#include "rb_lapack.h"

extern VOID slags2_(logical *upper, real *a1, real *a2, real *a3, real *b1, real *b2, real *b3, real *csu, real *snu, real *csv, real *snv, real *csq, real *snq);

static VALUE
rb_slags2(int argc, VALUE *argv, VALUE self){
  VALUE rb_upper;
  logical upper; 
  VALUE rb_a1;
  real a1; 
  VALUE rb_a2;
  real a2; 
  VALUE rb_a3;
  real a3; 
  VALUE rb_b1;
  real b1; 
  VALUE rb_b2;
  real b2; 
  VALUE rb_b3;
  real b3; 
  VALUE rb_csu;
  real csu; 
  VALUE rb_snu;
  real snu; 
  VALUE rb_csv;
  real csv; 
  VALUE rb_snv;
  real snv; 
  VALUE rb_csq;
  real csq; 
  VALUE rb_snq;
  real snq; 


  if (argc == 0) {
    printf("%s\n", "USAGE:\n  csu, snu, csv, snv, csq, snq = NumRu::Lapack.slags2( upper, a1, a2, a3, b1, b2, b3)\n    or\n  NumRu::Lapack.slags2  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 7)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 7)", argc);
  rb_upper = argv[0];
  rb_a1 = argv[1];
  rb_a2 = argv[2];
  rb_a3 = argv[3];
  rb_b1 = argv[4];
  rb_b2 = argv[5];
  rb_b3 = argv[6];

  b1 = (real)NUM2DBL(rb_b1);
  upper = (rb_upper == Qtrue);
  b2 = (real)NUM2DBL(rb_b2);
  a1 = (real)NUM2DBL(rb_a1);
  b3 = (real)NUM2DBL(rb_b3);
  a2 = (real)NUM2DBL(rb_a2);
  a3 = (real)NUM2DBL(rb_a3);

  slags2_(&upper, &a1, &a2, &a3, &b1, &b2, &b3, &csu, &snu, &csv, &snv, &csq, &snq);

  rb_csu = rb_float_new((double)csu);
  rb_snu = rb_float_new((double)snu);
  rb_csv = rb_float_new((double)csv);
  rb_snv = rb_float_new((double)snv);
  rb_csq = rb_float_new((double)csq);
  rb_snq = rb_float_new((double)snq);
  return rb_ary_new3(6, rb_csu, rb_snu, rb_csv, rb_snv, rb_csq, rb_snq);
}

void
init_lapack_slags2(VALUE mLapack){
  rb_define_module_function(mLapack, "slags2", rb_slags2, -1);
}

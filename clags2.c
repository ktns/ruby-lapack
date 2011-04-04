#include "rb_lapack.h"

extern VOID clags2_(logical *upper, real *a1, complex *a2, real *a3, real *b1, complex *b2, real *b3, real *csu, complex *snu, real *csv, complex *snv, real *csq, complex *snq);

static VALUE
rb_clags2(int argc, VALUE *argv, VALUE self){
  VALUE rb_upper;
  logical upper; 
  VALUE rb_a1;
  real a1; 
  VALUE rb_a2;
  complex a2; 
  VALUE rb_a3;
  real a3; 
  VALUE rb_b1;
  real b1; 
  VALUE rb_b2;
  complex b2; 
  VALUE rb_b3;
  real b3; 
  VALUE rb_csu;
  real csu; 
  VALUE rb_snu;
  complex snu; 
  VALUE rb_csv;
  real csv; 
  VALUE rb_snv;
  complex snv; 
  VALUE rb_csq;
  real csq; 
  VALUE rb_snq;
  complex snq; 


  if (argc == 0) {
    printf("%s\n", "USAGE:\n  csu, snu, csv, snv, csq, snq = NumRu::Lapack.clags2( upper, a1, a2, a3, b1, b2, b3)\n    or\n  NumRu::Lapack.clags2  # print help\n\n\nFORTRAN MANUAL\n\n");
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
  b2.r = (real)NUM2DBL(rb_funcall(rb_b2, rb_intern("real"), 0));
  b2.i = (real)NUM2DBL(rb_funcall(rb_b2, rb_intern("imag"), 0));
  a1 = (real)NUM2DBL(rb_a1);
  b3 = (real)NUM2DBL(rb_b3);
  a2.r = (real)NUM2DBL(rb_funcall(rb_a2, rb_intern("real"), 0));
  a2.i = (real)NUM2DBL(rb_funcall(rb_a2, rb_intern("imag"), 0));
  a3 = (real)NUM2DBL(rb_a3);

  clags2_(&upper, &a1, &a2, &a3, &b1, &b2, &b3, &csu, &snu, &csv, &snv, &csq, &snq);

  rb_csu = rb_float_new((double)csu);
  rb_snu = rb_funcall(rb_gv_get("Complex"), rb_intern("new"), 2, rb_float_new((double)(snu.r)), rb_float_new((double)(snu.i)));
  rb_csv = rb_float_new((double)csv);
  rb_snv = rb_funcall(rb_gv_get("Complex"), rb_intern("new"), 2, rb_float_new((double)(snv.r)), rb_float_new((double)(snv.i)));
  rb_csq = rb_float_new((double)csq);
  rb_snq = rb_funcall(rb_gv_get("Complex"), rb_intern("new"), 2, rb_float_new((double)(snq.r)), rb_float_new((double)(snq.i)));
  return rb_ary_new3(6, rb_csu, rb_snu, rb_csv, rb_snv, rb_csq, rb_snq);
}

void
init_lapack_clags2(VALUE mLapack){
  rb_define_module_function(mLapack, "clags2", rb_clags2, -1);
}

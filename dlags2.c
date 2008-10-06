#include "rb_lapack.h"

static VALUE
rb_dlags2(int argc, VALUE *argv, VALUE self){
  VALUE rb_upper;
  logical upper; 
  VALUE rb_a1;
  doublereal a1; 
  VALUE rb_a2;
  doublereal a2; 
  VALUE rb_a3;
  doublereal a3; 
  VALUE rb_b1;
  doublereal b1; 
  VALUE rb_b2;
  doublereal b2; 
  VALUE rb_b3;
  doublereal b3; 
  VALUE rb_csu;
  doublereal csu; 
  VALUE rb_snu;
  doublereal snu; 
  VALUE rb_csv;
  doublereal csv; 
  VALUE rb_snv;
  doublereal snv; 
  VALUE rb_csq;
  doublereal csq; 
  VALUE rb_snq;
  doublereal snq; 


  if (argc == 0) {
    printf("%s\n", "USAGE:\n  csu, snu, csv, snv, csq, snq = NumRu::Lapack.dlags2( upper, a1, a2, a3, b1, b2, b3)\n    or\n  NumRu::Lapack.dlags2  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE DLAGS2( UPPER, A1, A2, A3, B1, B2, B3, CSU, SNU, CSV, SNV, CSQ, SNQ )\n\n*  Purpose\n*  =======\n*\n*  DLAGS2 computes 2-by-2 orthogonal matrices U, V and Q, such\n*  that if ( UPPER ) then\n*\n*            U'*A*Q = U'*( A1 A2 )*Q = ( x  0  )\n*                        ( 0  A3 )     ( x  x  )\n*  and\n*            V'*B*Q = V'*( B1 B2 )*Q = ( x  0  )\n*                        ( 0  B3 )     ( x  x  )\n*\n*  or if ( .NOT.UPPER ) then\n*\n*            U'*A*Q = U'*( A1 0  )*Q = ( x  x  )\n*                        ( A2 A3 )     ( 0  x  )\n*  and\n*            V'*B*Q = V'*( B1 0  )*Q = ( x  x  )\n*                        ( B2 B3 )     ( 0  x  )\n*\n*  The rows of the transformed A and B are parallel, where\n*\n*    U = (  CSU  SNU ), V = (  CSV SNV ), Q = (  CSQ   SNQ )\n*        ( -SNU  CSU )      ( -SNV CSV )      ( -SNQ   CSQ )\n*\n*  Z' denotes the transpose of Z.\n*\n*\n\n*  Arguments\n*  =========\n*\n*  UPPER   (input) LOGICAL\n*          = .TRUE.: the input matrices A and B are upper triangular.\n*          = .FALSE.: the input matrices A and B are lower triangular.\n*\n*  A1      (input) DOUBLE PRECISION\n*  A2      (input) DOUBLE PRECISION\n*  A3      (input) DOUBLE PRECISION\n*          On entry, A1, A2 and A3 are elements of the input 2-by-2\n*          upper (lower) triangular matrix A.\n*\n*  B1      (input) DOUBLE PRECISION\n*  B2      (input) DOUBLE PRECISION\n*  B3      (input) DOUBLE PRECISION\n*          On entry, B1, B2 and B3 are elements of the input 2-by-2\n*          upper (lower) triangular matrix B.\n*\n*  CSU     (output) DOUBLE PRECISION\n*  SNU     (output) DOUBLE PRECISION\n*          The desired orthogonal matrix U.\n*\n*  CSV     (output) DOUBLE PRECISION\n*  SNV     (output) DOUBLE PRECISION\n*          The desired orthogonal matrix V.\n*\n*  CSQ     (output) DOUBLE PRECISION\n*  SNQ     (output) DOUBLE PRECISION\n*          The desired orthogonal matrix Q.\n*\n\n*  =====================================================================\n*\n\n");
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

  upper = (rb_upper == Qtrue);
  a1 = NUM2DBL(rb_a1);
  a2 = NUM2DBL(rb_a2);
  a3 = NUM2DBL(rb_a3);
  b1 = NUM2DBL(rb_b1);
  b2 = NUM2DBL(rb_b2);
  b3 = NUM2DBL(rb_b3);

  dlags2_(&upper, &a1, &a2, &a3, &b1, &b2, &b3, &csu, &snu, &csv, &snv, &csq, &snq);

  rb_csu = rb_float_new((double)csu);
  rb_snu = rb_float_new((double)snu);
  rb_csv = rb_float_new((double)csv);
  rb_snv = rb_float_new((double)snv);
  rb_csq = rb_float_new((double)csq);
  rb_snq = rb_float_new((double)snq);
  return rb_ary_new3(6, rb_csu, rb_snu, rb_csv, rb_snv, rb_csq, rb_snq);
}

void
init_lapack_dlags2(VALUE mLapack){
  rb_define_module_function(mLapack, "dlags2", rb_dlags2, -1);
}

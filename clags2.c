#include "rb_lapack.h"

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
    printf("%s\n", "USAGE:\n  csu, snu, csv, snv, csq, snq = NumRu::Lapack.clags2( upper, a1, a2, a3, b1, b2, b3)\n    or\n  NumRu::Lapack.clags2  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE CLAGS2( UPPER, A1, A2, A3, B1, B2, B3, CSU, SNU, CSV, SNV, CSQ, SNQ )\n\n*  Purpose\n*  =======\n*\n*  CLAGS2 computes 2-by-2 unitary matrices U, V and Q, such\n*  that if ( UPPER ) then\n*\n*            U'*A*Q = U'*( A1 A2 )*Q = ( x  0  )\n*                        ( 0  A3 )     ( x  x  )\n*  and\n*            V'*B*Q = V'*( B1 B2 )*Q = ( x  0  )\n*                        ( 0  B3 )     ( x  x  )\n*\n*  or if ( .NOT.UPPER ) then\n*\n*            U'*A*Q = U'*( A1 0  )*Q = ( x  x  )\n*                        ( A2 A3 )     ( 0  x  )\n*  and\n*            V'*B*Q = V'*( B1 0  )*Q = ( x  x  )\n*                        ( B2 B3 )     ( 0  x  )\n*  where\n*\n*    U = (     CSU      SNU ), V = (     CSV     SNV ),\n*        ( -CONJG(SNU)  CSU )      ( -CONJG(SNV) CSV )\n*\n*    Q = (     CSQ      SNQ )\n*        ( -CONJG(SNQ)  CSQ )\n*\n*  Z' denotes the conjugate transpose of Z.\n*\n*  The rows of the transformed A and B are parallel. Moreover, if the\n*  input 2-by-2 matrix A is not zero, then the transformed (1,1) entry\n*  of A is not zero. If the input matrices A and B are both not zero,\n*  then the transformed (2,2) element of B is not zero, except when the\n*  first rows of input A and B are parallel and the second rows are\n*  zero.\n*\n\n*  Arguments\n*  =========\n*\n*  UPPER   (input) LOGICAL\n*          = .TRUE.: the input matrices A and B are upper triangular.\n*          = .FALSE.: the input matrices A and B are lower triangular.\n*\n*  A1      (input) REAL\n*  A2      (input) COMPLEX\n*  A3      (input) REAL\n*          On entry, A1, A2 and A3 are elements of the input 2-by-2\n*          upper (lower) triangular matrix A.\n*\n*  B1      (input) REAL\n*  B2      (input) COMPLEX\n*  B3      (input) REAL\n*          On entry, B1, B2 and B3 are elements of the input 2-by-2\n*          upper (lower) triangular matrix B.\n*\n*  CSU     (output) REAL\n*  SNU     (output) COMPLEX\n*          The desired unitary matrix U.\n*\n*  CSV     (output) REAL\n*  SNV     (output) COMPLEX\n*          The desired unitary matrix V.\n*\n*  CSQ     (output) REAL\n*  SNQ     (output) COMPLEX\n*          The desired unitary matrix Q.\n*\n\n*  =====================================================================\n*\n\n");
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
  a1 = (real)NUM2DBL(rb_a1);
  a2.r = (real)NUM2DBL(rb_funcall(rb_a2, rb_intern("real"), 0));
  a2.i = (real)NUM2DBL(rb_funcall(rb_a2, rb_intern("imag"), 0));
  a3 = (real)NUM2DBL(rb_a3);
  b1 = (real)NUM2DBL(rb_b1);
  b2.r = (real)NUM2DBL(rb_funcall(rb_b2, rb_intern("real"), 0));
  b2.i = (real)NUM2DBL(rb_funcall(rb_b2, rb_intern("imag"), 0));
  b3 = (real)NUM2DBL(rb_b3);

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

#include "rb_lapack.h"

extern VOID dlags2_(logical *upper, doublereal *a1, doublereal *a2, doublereal *a3, doublereal *b1, doublereal *b2, doublereal *b3, doublereal *csu, doublereal *snu, doublereal *csv, doublereal *snv, doublereal *csq, doublereal *snq);

static VALUE sHelp, sUsage;

static VALUE
rblapack_dlags2(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_upper;
  logical upper; 
  VALUE rblapack_a1;
  doublereal a1; 
  VALUE rblapack_a2;
  doublereal a2; 
  VALUE rblapack_a3;
  doublereal a3; 
  VALUE rblapack_b1;
  doublereal b1; 
  VALUE rblapack_b2;
  doublereal b2; 
  VALUE rblapack_b3;
  doublereal b3; 
  VALUE rblapack_csu;
  doublereal csu; 
  VALUE rblapack_snu;
  doublereal snu; 
  VALUE rblapack_csv;
  doublereal csv; 
  VALUE rblapack_snv;
  doublereal snv; 
  VALUE rblapack_csq;
  doublereal csq; 
  VALUE rblapack_snq;
  doublereal snq; 


  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  csu, snu, csv, snv, csq, snq = NumRu::Lapack.dlags2( upper, a1, a2, a3, b1, b2, b3, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE DLAGS2( UPPER, A1, A2, A3, B1, B2, B3, CSU, SNU, CSV, SNV, CSQ, SNQ )\n\n*  Purpose\n*  =======\n*\n*  DLAGS2 computes 2-by-2 orthogonal matrices U, V and Q, such\n*  that if ( UPPER ) then\n*\n*            U'*A*Q = U'*( A1 A2 )*Q = ( x  0  )\n*                        ( 0  A3 )     ( x  x  )\n*  and\n*            V'*B*Q = V'*( B1 B2 )*Q = ( x  0  )\n*                        ( 0  B3 )     ( x  x  )\n*\n*  or if ( .NOT.UPPER ) then\n*\n*            U'*A*Q = U'*( A1 0  )*Q = ( x  x  )\n*                        ( A2 A3 )     ( 0  x  )\n*  and\n*            V'*B*Q = V'*( B1 0  )*Q = ( x  x  )\n*                        ( B2 B3 )     ( 0  x  )\n*\n*  The rows of the transformed A and B are parallel, where\n*\n*    U = (  CSU  SNU ), V = (  CSV SNV ), Q = (  CSQ   SNQ )\n*        ( -SNU  CSU )      ( -SNV CSV )      ( -SNQ   CSQ )\n*\n*  Z' denotes the transpose of Z.\n*\n*\n\n*  Arguments\n*  =========\n*\n*  UPPER   (input) LOGICAL\n*          = .TRUE.: the input matrices A and B are upper triangular.\n*          = .FALSE.: the input matrices A and B are lower triangular.\n*\n*  A1      (input) DOUBLE PRECISION\n*  A2      (input) DOUBLE PRECISION\n*  A3      (input) DOUBLE PRECISION\n*          On entry, A1, A2 and A3 are elements of the input 2-by-2\n*          upper (lower) triangular matrix A.\n*\n*  B1      (input) DOUBLE PRECISION\n*  B2      (input) DOUBLE PRECISION\n*  B3      (input) DOUBLE PRECISION\n*          On entry, B1, B2 and B3 are elements of the input 2-by-2\n*          upper (lower) triangular matrix B.\n*\n*  CSU     (output) DOUBLE PRECISION\n*  SNU     (output) DOUBLE PRECISION\n*          The desired orthogonal matrix U.\n*\n*  CSV     (output) DOUBLE PRECISION\n*  SNV     (output) DOUBLE PRECISION\n*          The desired orthogonal matrix V.\n*\n*  CSQ     (output) DOUBLE PRECISION\n*  SNQ     (output) DOUBLE PRECISION\n*          The desired orthogonal matrix Q.\n*\n\n*  =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  csu, snu, csv, snv, csq, snq = NumRu::Lapack.dlags2( upper, a1, a2, a3, b1, b2, b3, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 7)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 7)", argc);
  rblapack_upper = argv[0];
  rblapack_a1 = argv[1];
  rblapack_a2 = argv[2];
  rblapack_a3 = argv[3];
  rblapack_b1 = argv[4];
  rblapack_b2 = argv[5];
  rblapack_b3 = argv[6];
  if (rb_options != Qnil) {
  }

  b1 = NUM2DBL(rblapack_b1);
  upper = (rblapack_upper == Qtrue);
  b2 = NUM2DBL(rblapack_b2);
  a1 = NUM2DBL(rblapack_a1);
  b3 = NUM2DBL(rblapack_b3);
  a2 = NUM2DBL(rblapack_a2);
  a3 = NUM2DBL(rblapack_a3);

  dlags2_(&upper, &a1, &a2, &a3, &b1, &b2, &b3, &csu, &snu, &csv, &snv, &csq, &snq);

  rblapack_csu = rb_float_new((double)csu);
  rblapack_snu = rb_float_new((double)snu);
  rblapack_csv = rb_float_new((double)csv);
  rblapack_snv = rb_float_new((double)snv);
  rblapack_csq = rb_float_new((double)csq);
  rblapack_snq = rb_float_new((double)snq);
  return rb_ary_new3(6, rblapack_csu, rblapack_snu, rblapack_csv, rblapack_snv, rblapack_csq, rblapack_snq);
}

void
init_lapack_dlags2(VALUE mLapack){
  rb_define_module_function(mLapack, "dlags2", rblapack_dlags2, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}

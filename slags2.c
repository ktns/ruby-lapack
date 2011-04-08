#include "rb_lapack.h"

extern VOID slags2_(logical *upper, real *a1, real *a2, real *a3, real *b1, real *b2, real *b3, real *csu, real *snu, real *csv, real *snv, real *csq, real *snq);

static VALUE sHelp, sUsage;

static VALUE
rblapack_slags2(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_upper;
  logical upper; 
  VALUE rblapack_a1;
  real a1; 
  VALUE rblapack_a2;
  real a2; 
  VALUE rblapack_a3;
  real a3; 
  VALUE rblapack_b1;
  real b1; 
  VALUE rblapack_b2;
  real b2; 
  VALUE rblapack_b3;
  real b3; 
  VALUE rblapack_csu;
  real csu; 
  VALUE rblapack_snu;
  real snu; 
  VALUE rblapack_csv;
  real csv; 
  VALUE rblapack_snv;
  real snv; 
  VALUE rblapack_csq;
  real csq; 
  VALUE rblapack_snq;
  real snq; 


  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  csu, snu, csv, snv, csq, snq = NumRu::Lapack.slags2( upper, a1, a2, a3, b1, b2, b3, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE SLAGS2( UPPER, A1, A2, A3, B1, B2, B3, CSU, SNU, CSV, SNV, CSQ, SNQ )\n\n*  Purpose\n*  =======\n*\n*  SLAGS2 computes 2-by-2 orthogonal matrices U, V and Q, such\n*  that if ( UPPER ) then\n*\n*            U'*A*Q = U'*( A1 A2 )*Q = ( x  0  )\n*                        ( 0  A3 )     ( x  x  )\n*  and\n*            V'*B*Q = V'*( B1 B2 )*Q = ( x  0  )\n*                        ( 0  B3 )     ( x  x  )\n*\n*  or if ( .NOT.UPPER ) then\n*\n*            U'*A*Q = U'*( A1 0  )*Q = ( x  x  )\n*                        ( A2 A3 )     ( 0  x  )\n*  and\n*            V'*B*Q = V'*( B1 0  )*Q = ( x  x  )\n*                        ( B2 B3 )     ( 0  x  )\n*\n*  The rows of the transformed A and B are parallel, where\n*\n*    U = (  CSU  SNU ), V = (  CSV SNV ), Q = (  CSQ   SNQ )\n*        ( -SNU  CSU )      ( -SNV CSV )      ( -SNQ   CSQ )\n*\n*  Z' denotes the transpose of Z.\n*\n*\n\n*  Arguments\n*  =========\n*\n*  UPPER   (input) LOGICAL\n*          = .TRUE.: the input matrices A and B are upper triangular.\n*          = .FALSE.: the input matrices A and B are lower triangular.\n*\n*  A1      (input) REAL\n*  A2      (input) REAL\n*  A3      (input) REAL\n*          On entry, A1, A2 and A3 are elements of the input 2-by-2\n*          upper (lower) triangular matrix A.\n*\n*  B1      (input) REAL\n*  B2      (input) REAL\n*  B3      (input) REAL\n*          On entry, B1, B2 and B3 are elements of the input 2-by-2\n*          upper (lower) triangular matrix B.\n*\n*  CSU     (output) REAL\n*  SNU     (output) REAL\n*          The desired orthogonal matrix U.\n*\n*  CSV     (output) REAL\n*  SNV     (output) REAL\n*          The desired orthogonal matrix V.\n*\n*  CSQ     (output) REAL\n*  SNQ     (output) REAL\n*          The desired orthogonal matrix Q.\n*\n\n*  =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  csu, snu, csv, snv, csq, snq = NumRu::Lapack.slags2( upper, a1, a2, a3, b1, b2, b3, [:usage => usage, :help => help])\n");
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

  b1 = (real)NUM2DBL(rblapack_b1);
  upper = (rblapack_upper == Qtrue);
  b2 = (real)NUM2DBL(rblapack_b2);
  a1 = (real)NUM2DBL(rblapack_a1);
  b3 = (real)NUM2DBL(rblapack_b3);
  a2 = (real)NUM2DBL(rblapack_a2);
  a3 = (real)NUM2DBL(rblapack_a3);

  slags2_(&upper, &a1, &a2, &a3, &b1, &b2, &b3, &csu, &snu, &csv, &snv, &csq, &snq);

  rblapack_csu = rb_float_new((double)csu);
  rblapack_snu = rb_float_new((double)snu);
  rblapack_csv = rb_float_new((double)csv);
  rblapack_snv = rb_float_new((double)snv);
  rblapack_csq = rb_float_new((double)csq);
  rblapack_snq = rb_float_new((double)snq);
  return rb_ary_new3(6, rblapack_csu, rblapack_snu, rblapack_csv, rblapack_snv, rblapack_csq, rblapack_snq);
}

void
init_lapack_slags2(VALUE mLapack){
  rb_define_module_function(mLapack, "slags2", rblapack_slags2, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}

#include "rb_lapack.h"

extern VOID zlags2_(logical *upper, doublereal *a1, doublecomplex *a2, doublereal *a3, doublereal *b1, doublecomplex *b2, doublereal *b3, doublereal *csu, doublecomplex *snu, doublereal *csv, doublecomplex *snv, doublereal *csq, doublecomplex *snq);

static VALUE sHelp, sUsage;

static VALUE
rblapack_zlags2(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_upper;
  logical upper; 
  VALUE rblapack_a1;
  doublereal a1; 
  VALUE rblapack_a2;
  doublecomplex a2; 
  VALUE rblapack_a3;
  doublereal a3; 
  VALUE rblapack_b1;
  doublereal b1; 
  VALUE rblapack_b2;
  doublecomplex b2; 
  VALUE rblapack_b3;
  doublereal b3; 
  VALUE rblapack_csu;
  doublereal csu; 
  VALUE rblapack_snu;
  doublecomplex snu; 
  VALUE rblapack_csv;
  doublereal csv; 
  VALUE rblapack_snv;
  doublecomplex snv; 
  VALUE rblapack_csq;
  doublereal csq; 
  VALUE rblapack_snq;
  doublecomplex snq; 


  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  csu, snu, csv, snv, csq, snq = NumRu::Lapack.zlags2( upper, a1, a2, a3, b1, b2, b3, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE ZLAGS2( UPPER, A1, A2, A3, B1, B2, B3, CSU, SNU, CSV, SNV, CSQ, SNQ )\n\n*  Purpose\n*  =======\n*\n*  ZLAGS2 computes 2-by-2 unitary matrices U, V and Q, such\n*  that if ( UPPER ) then\n*\n*            U'*A*Q = U'*( A1 A2 )*Q = ( x  0  )\n*                        ( 0  A3 )     ( x  x  )\n*  and\n*            V'*B*Q = V'*( B1 B2 )*Q = ( x  0  )\n*                        ( 0  B3 )     ( x  x  )\n*\n*  or if ( .NOT.UPPER ) then\n*\n*            U'*A*Q = U'*( A1 0  )*Q = ( x  x  )\n*                        ( A2 A3 )     ( 0  x  )\n*  and\n*            V'*B*Q = V'*( B1 0  )*Q = ( x  x  )\n*                        ( B2 B3 )     ( 0  x  )\n*  where\n*\n*    U = (     CSU      SNU ), V = (     CSV     SNV ),\n*        ( -CONJG(SNU)  CSU )      ( -CONJG(SNV) CSV )\n*\n*    Q = (     CSQ      SNQ )\n*        ( -CONJG(SNQ)  CSQ )\n*\n*  Z' denotes the conjugate transpose of Z.\n*\n*  The rows of the transformed A and B are parallel. Moreover, if the\n*  input 2-by-2 matrix A is not zero, then the transformed (1,1) entry\n*  of A is not zero. If the input matrices A and B are both not zero,\n*  then the transformed (2,2) element of B is not zero, except when the\n*  first rows of input A and B are parallel and the second rows are\n*  zero.\n*\n\n*  Arguments\n*  =========\n*\n*  UPPER   (input) LOGICAL\n*          = .TRUE.: the input matrices A and B are upper triangular.\n*          = .FALSE.: the input matrices A and B are lower triangular.\n*\n*  A1      (input) DOUBLE PRECISION\n*  A2      (input) COMPLEX*16\n*  A3      (input) DOUBLE PRECISION\n*          On entry, A1, A2 and A3 are elements of the input 2-by-2\n*          upper (lower) triangular matrix A.\n*\n*  B1      (input) DOUBLE PRECISION\n*  B2      (input) COMPLEX*16\n*  B3      (input) DOUBLE PRECISION\n*          On entry, B1, B2 and B3 are elements of the input 2-by-2\n*          upper (lower) triangular matrix B.\n*\n*  CSU     (output) DOUBLE PRECISION\n*  SNU     (output) COMPLEX*16\n*          The desired unitary matrix U.\n*\n*  CSV     (output) DOUBLE PRECISION\n*  SNV     (output) COMPLEX*16\n*          The desired unitary matrix V.\n*\n*  CSQ     (output) DOUBLE PRECISION\n*  SNQ     (output) COMPLEX*16\n*          The desired unitary matrix Q.\n*\n\n*  =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  csu, snu, csv, snv, csq, snq = NumRu::Lapack.zlags2( upper, a1, a2, a3, b1, b2, b3, [:usage => usage, :help => help])\n");
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
  b2.r = NUM2DBL(rb_funcall(rblapack_b2, rb_intern("real"), 0));
  b2.i = NUM2DBL(rb_funcall(rblapack_b2, rb_intern("imag"), 0));
  a1 = NUM2DBL(rblapack_a1);
  b3 = NUM2DBL(rblapack_b3);
  a2.r = NUM2DBL(rb_funcall(rblapack_a2, rb_intern("real"), 0));
  a2.i = NUM2DBL(rb_funcall(rblapack_a2, rb_intern("imag"), 0));
  a3 = NUM2DBL(rblapack_a3);

  zlags2_(&upper, &a1, &a2, &a3, &b1, &b2, &b3, &csu, &snu, &csv, &snv, &csq, &snq);

  rblapack_csu = rb_float_new((double)csu);
  rblapack_snu = rb_funcall(rb_gv_get("Complex"), rb_intern("new"), 2, rb_float_new((double)(snu.r)), rb_float_new((double)(snu.i)));
  rblapack_csv = rb_float_new((double)csv);
  rblapack_snv = rb_funcall(rb_gv_get("Complex"), rb_intern("new"), 2, rb_float_new((double)(snv.r)), rb_float_new((double)(snv.i)));
  rblapack_csq = rb_float_new((double)csq);
  rblapack_snq = rb_funcall(rb_gv_get("Complex"), rb_intern("new"), 2, rb_float_new((double)(snq.r)), rb_float_new((double)(snq.i)));
  return rb_ary_new3(6, rblapack_csu, rblapack_snu, rblapack_csv, rblapack_snv, rblapack_csq, rblapack_snq);
}

void
init_lapack_zlags2(VALUE mLapack){
  rb_define_module_function(mLapack, "zlags2", rblapack_zlags2, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}

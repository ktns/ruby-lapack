#include "rb_lapack.h"

extern VOID slagv2_(real *a, integer *lda, real *b, integer *ldb, real *alphar, real *alphai, real *beta, real *csl, real *snl, real *csr, real *snr);

static VALUE sHelp, sUsage;

static VALUE
rblapack_slagv2(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_a;
  real *a; 
  VALUE rblapack_b;
  real *b; 
  VALUE rblapack_alphar;
  real *alphar; 
  VALUE rblapack_alphai;
  real *alphai; 
  VALUE rblapack_beta;
  real *beta; 
  VALUE rblapack_csl;
  real csl; 
  VALUE rblapack_snl;
  real snl; 
  VALUE rblapack_csr;
  real csr; 
  VALUE rblapack_snr;
  real snr; 
  VALUE rblapack_a_out__;
  real *a_out__;
  VALUE rblapack_b_out__;
  real *b_out__;

  integer lda;
  integer ldb;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  alphar, alphai, beta, csl, snl, csr, snr, a, b = NumRu::Lapack.slagv2( a, b, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE SLAGV2( A, LDA, B, LDB, ALPHAR, ALPHAI, BETA, CSL, SNL, CSR, SNR )\n\n*  Purpose\n*  =======\n*\n*  SLAGV2 computes the Generalized Schur factorization of a real 2-by-2\n*  matrix pencil (A,B) where B is upper triangular. This routine\n*  computes orthogonal (rotation) matrices given by CSL, SNL and CSR,\n*  SNR such that\n*\n*  1) if the pencil (A,B) has two real eigenvalues (include 0/0 or 1/0\n*     types), then\n*\n*     [ a11 a12 ] := [  CSL  SNL ] [ a11 a12 ] [  CSR -SNR ]\n*     [  0  a22 ]    [ -SNL  CSL ] [ a21 a22 ] [  SNR  CSR ]\n*\n*     [ b11 b12 ] := [  CSL  SNL ] [ b11 b12 ] [  CSR -SNR ]\n*     [  0  b22 ]    [ -SNL  CSL ] [  0  b22 ] [  SNR  CSR ],\n*\n*  2) if the pencil (A,B) has a pair of complex conjugate eigenvalues,\n*     then\n*\n*     [ a11 a12 ] := [  CSL  SNL ] [ a11 a12 ] [  CSR -SNR ]\n*     [ a21 a22 ]    [ -SNL  CSL ] [ a21 a22 ] [  SNR  CSR ]\n*\n*     [ b11  0  ] := [  CSL  SNL ] [ b11 b12 ] [  CSR -SNR ]\n*     [  0  b22 ]    [ -SNL  CSL ] [  0  b22 ] [  SNR  CSR ]\n*\n*     where b11 >= b22 > 0.\n*\n*\n\n*  Arguments\n*  =========\n*\n*  A       (input/output) REAL array, dimension (LDA, 2)\n*          On entry, the 2 x 2 matrix A.\n*          On exit, A is overwritten by the ``A-part'' of the\n*          generalized Schur form.\n*\n*  LDA     (input) INTEGER\n*          THe leading dimension of the array A.  LDA >= 2.\n*\n*  B       (input/output) REAL array, dimension (LDB, 2)\n*          On entry, the upper triangular 2 x 2 matrix B.\n*          On exit, B is overwritten by the ``B-part'' of the\n*          generalized Schur form.\n*\n*  LDB     (input) INTEGER\n*          THe leading dimension of the array B.  LDB >= 2.\n*\n*  ALPHAR  (output) REAL array, dimension (2)\n*  ALPHAI  (output) REAL array, dimension (2)\n*  BETA    (output) REAL array, dimension (2)\n*          (ALPHAR(k)+i*ALPHAI(k))/BETA(k) are the eigenvalues of the\n*          pencil (A,B), k=1,2, i = sqrt(-1).  Note that BETA(k) may\n*          be zero.\n*\n*  CSL     (output) REAL\n*          The cosine of the left rotation matrix.\n*\n*  SNL     (output) REAL\n*          The sine of the left rotation matrix.\n*\n*  CSR     (output) REAL\n*          The cosine of the right rotation matrix.\n*\n*  SNR     (output) REAL\n*          The sine of the right rotation matrix.\n*\n\n*  Further Details\n*  ===============\n*\n*  Based on contributions by\n*     Mark Fahey, Department of Mathematics, Univ. of Kentucky, USA\n*\n*  =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  alphar, alphai, beta, csl, snl, csr, snr, a, b = NumRu::Lapack.slagv2( a, b, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 2)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 2)", argc);
  rblapack_a = argv[0];
  rblapack_b = argv[1];
  if (rb_options != Qnil) {
  }

  if (!NA_IsNArray(rblapack_a))
    rb_raise(rb_eArgError, "a (1th argument) must be NArray");
  if (NA_RANK(rblapack_a) != 2)
    rb_raise(rb_eArgError, "rank of a (1th argument) must be %d", 2);
  if (NA_SHAPE1(rblapack_a) != (2))
    rb_raise(rb_eRuntimeError, "shape 1 of a must be %d", 2);
  lda = NA_SHAPE0(rblapack_a);
  if (NA_TYPE(rblapack_a) != NA_SFLOAT)
    rblapack_a = na_change_type(rblapack_a, NA_SFLOAT);
  a = NA_PTR_TYPE(rblapack_a, real*);
  if (!NA_IsNArray(rblapack_b))
    rb_raise(rb_eArgError, "b (2th argument) must be NArray");
  if (NA_RANK(rblapack_b) != 2)
    rb_raise(rb_eArgError, "rank of b (2th argument) must be %d", 2);
  if (NA_SHAPE1(rblapack_b) != (2))
    rb_raise(rb_eRuntimeError, "shape 1 of b must be %d", 2);
  ldb = NA_SHAPE0(rblapack_b);
  if (NA_TYPE(rblapack_b) != NA_SFLOAT)
    rblapack_b = na_change_type(rblapack_b, NA_SFLOAT);
  b = NA_PTR_TYPE(rblapack_b, real*);
  {
    int shape[1];
    shape[0] = 2;
    rblapack_alphar = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  alphar = NA_PTR_TYPE(rblapack_alphar, real*);
  {
    int shape[1];
    shape[0] = 2;
    rblapack_alphai = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  alphai = NA_PTR_TYPE(rblapack_alphai, real*);
  {
    int shape[1];
    shape[0] = 2;
    rblapack_beta = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  beta = NA_PTR_TYPE(rblapack_beta, real*);
  {
    int shape[2];
    shape[0] = lda;
    shape[1] = 2;
    rblapack_a_out__ = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  a_out__ = NA_PTR_TYPE(rblapack_a_out__, real*);
  MEMCPY(a_out__, a, real, NA_TOTAL(rblapack_a));
  rblapack_a = rblapack_a_out__;
  a = a_out__;
  {
    int shape[2];
    shape[0] = ldb;
    shape[1] = 2;
    rblapack_b_out__ = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  b_out__ = NA_PTR_TYPE(rblapack_b_out__, real*);
  MEMCPY(b_out__, b, real, NA_TOTAL(rblapack_b));
  rblapack_b = rblapack_b_out__;
  b = b_out__;

  slagv2_(a, &lda, b, &ldb, alphar, alphai, beta, &csl, &snl, &csr, &snr);

  rblapack_csl = rb_float_new((double)csl);
  rblapack_snl = rb_float_new((double)snl);
  rblapack_csr = rb_float_new((double)csr);
  rblapack_snr = rb_float_new((double)snr);
  return rb_ary_new3(9, rblapack_alphar, rblapack_alphai, rblapack_beta, rblapack_csl, rblapack_snl, rblapack_csr, rblapack_snr, rblapack_a, rblapack_b);
}

void
init_lapack_slagv2(VALUE mLapack){
  rb_define_module_function(mLapack, "slagv2", rblapack_slagv2, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}

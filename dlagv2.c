#include "rb_lapack.h"

static VALUE
rb_dlagv2(int argc, VALUE *argv, VALUE self){
  VALUE rb_a;
  doublereal *a; 
  VALUE rb_b;
  doublereal *b; 
  VALUE rb_alphar;
  doublereal *alphar; 
  VALUE rb_alphai;
  doublereal *alphai; 
  VALUE rb_beta;
  doublereal *beta; 
  VALUE rb_csl;
  doublereal csl; 
  VALUE rb_snl;
  doublereal snl; 
  VALUE rb_csr;
  doublereal csr; 
  VALUE rb_snr;
  doublereal snr; 
  VALUE rb_a_out__;
  doublereal *a_out__;
  VALUE rb_b_out__;
  doublereal *b_out__;

  integer lda;
  integer ldb;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  alphar, alphai, beta, csl, snl, csr, snr, a, b = NumRu::Lapack.dlagv2( a, b)\n    or\n  NumRu::Lapack.dlagv2  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE DLAGV2( A, LDA, B, LDB, ALPHAR, ALPHAI, BETA, CSL, SNL, CSR, SNR )\n\n*  Purpose\n*  =======\n*\n*  DLAGV2 computes the Generalized Schur factorization of a real 2-by-2\n*  matrix pencil (A,B) where B is upper triangular. This routine\n*  computes orthogonal (rotation) matrices given by CSL, SNL and CSR,\n*  SNR such that\n*\n*  1) if the pencil (A,B) has two real eigenvalues (include 0/0 or 1/0\n*     types), then\n*\n*     [ a11 a12 ] := [  CSL  SNL ] [ a11 a12 ] [  CSR -SNR ]\n*     [  0  a22 ]    [ -SNL  CSL ] [ a21 a22 ] [  SNR  CSR ]\n*\n*     [ b11 b12 ] := [  CSL  SNL ] [ b11 b12 ] [  CSR -SNR ]\n*     [  0  b22 ]    [ -SNL  CSL ] [  0  b22 ] [  SNR  CSR ],\n*\n*  2) if the pencil (A,B) has a pair of complex conjugate eigenvalues,\n*     then\n*\n*     [ a11 a12 ] := [  CSL  SNL ] [ a11 a12 ] [  CSR -SNR ]\n*     [ a21 a22 ]    [ -SNL  CSL ] [ a21 a22 ] [  SNR  CSR ]\n*\n*     [ b11  0  ] := [  CSL  SNL ] [ b11 b12 ] [  CSR -SNR ]\n*     [  0  b22 ]    [ -SNL  CSL ] [  0  b22 ] [  SNR  CSR ]\n*\n*     where b11 >= b22 > 0.\n*\n*\n\n*  Arguments\n*  =========\n*\n*  A       (input/output) DOUBLE PRECISION array, dimension (LDA, 2)\n*          On entry, the 2 x 2 matrix A.\n*          On exit, A is overwritten by the ``A-part'' of the\n*          generalized Schur form.\n*\n*  LDA     (input) INTEGER\n*          THe leading dimension of the array A.  LDA >= 2.\n*\n*  B       (input/output) DOUBLE PRECISION array, dimension (LDB, 2)\n*          On entry, the upper triangular 2 x 2 matrix B.\n*          On exit, B is overwritten by the ``B-part'' of the\n*          generalized Schur form.\n*\n*  LDB     (input) INTEGER\n*          THe leading dimension of the array B.  LDB >= 2.\n*\n*  ALPHAR  (output) DOUBLE PRECISION array, dimension (2)\n*  ALPHAI  (output) DOUBLE PRECISION array, dimension (2)\n*  BETA    (output) DOUBLE PRECISION array, dimension (2)\n*          (ALPHAR(k)+i*ALPHAI(k))/BETA(k) are the eigenvalues of the\n*          pencil (A,B), k=1,2, i = sqrt(-1).  Note that BETA(k) may\n*          be zero.\n*\n*  CSL     (output) DOUBLE PRECISION\n*          The cosine of the left rotation matrix.\n*\n*  SNL     (output) DOUBLE PRECISION\n*          The sine of the left rotation matrix.\n*\n*  CSR     (output) DOUBLE PRECISION\n*          The cosine of the right rotation matrix.\n*\n*  SNR     (output) DOUBLE PRECISION\n*          The sine of the right rotation matrix.\n*\n\n*  Further Details\n*  ===============\n*\n*  Based on contributions by\n*     Mark Fahey, Department of Mathematics, Univ. of Kentucky, USA\n*\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 2)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 2)", argc);
  rb_a = argv[0];
  rb_b = argv[1];

  if (!NA_IsNArray(rb_a))
    rb_raise(rb_eArgError, "a (1th argument) must be NArray");
  if (NA_RANK(rb_a) != 2)
    rb_raise(rb_eArgError, "rank of a (1th argument) must be %d", 2);
  lda = NA_SHAPE0(rb_a);
  if (NA_SHAPE1(rb_a) != (2))
    rb_raise(rb_eRuntimeError, "shape 1 of a must be %d", 2);
  if (NA_TYPE(rb_a) != NA_DFLOAT)
    rb_a = na_change_type(rb_a, NA_DFLOAT);
  a = NA_PTR_TYPE(rb_a, doublereal*);
  if (!NA_IsNArray(rb_b))
    rb_raise(rb_eArgError, "b (2th argument) must be NArray");
  if (NA_RANK(rb_b) != 2)
    rb_raise(rb_eArgError, "rank of b (2th argument) must be %d", 2);
  ldb = NA_SHAPE0(rb_b);
  if (NA_SHAPE1(rb_b) != (2))
    rb_raise(rb_eRuntimeError, "shape 1 of b must be %d", 2);
  if (NA_TYPE(rb_b) != NA_DFLOAT)
    rb_b = na_change_type(rb_b, NA_DFLOAT);
  b = NA_PTR_TYPE(rb_b, doublereal*);
  {
    int shape[1];
    shape[0] = 2;
    rb_alphar = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  alphar = NA_PTR_TYPE(rb_alphar, doublereal*);
  {
    int shape[1];
    shape[0] = 2;
    rb_alphai = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  alphai = NA_PTR_TYPE(rb_alphai, doublereal*);
  {
    int shape[1];
    shape[0] = 2;
    rb_beta = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  beta = NA_PTR_TYPE(rb_beta, doublereal*);
  {
    int shape[2];
    shape[0] = lda;
    shape[1] = 2;
    rb_a_out__ = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  a_out__ = NA_PTR_TYPE(rb_a_out__, doublereal*);
  MEMCPY(a_out__, a, doublereal, NA_TOTAL(rb_a));
  rb_a = rb_a_out__;
  a = a_out__;
  {
    int shape[2];
    shape[0] = ldb;
    shape[1] = 2;
    rb_b_out__ = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  b_out__ = NA_PTR_TYPE(rb_b_out__, doublereal*);
  MEMCPY(b_out__, b, doublereal, NA_TOTAL(rb_b));
  rb_b = rb_b_out__;
  b = b_out__;

  dlagv2_(a, &lda, b, &ldb, alphar, alphai, beta, &csl, &snl, &csr, &snr);

  rb_csl = rb_float_new((double)csl);
  rb_snl = rb_float_new((double)snl);
  rb_csr = rb_float_new((double)csr);
  rb_snr = rb_float_new((double)snr);
  return rb_ary_new3(9, rb_alphar, rb_alphai, rb_beta, rb_csl, rb_snl, rb_csr, rb_snr, rb_a, rb_b);
}

void
init_lapack_dlagv2(VALUE mLapack){
  rb_define_module_function(mLapack, "dlagv2", rb_dlagv2, -1);
}

#include "rb_lapack.h"

extern VOID slag2_(real *a, integer *lda, real *b, integer *ldb, real *safmin, real *scale1, real *scale2, real *wr1, real *wr2, real *wi);

static VALUE
rb_slag2(int argc, VALUE *argv, VALUE self){
  VALUE rb_a;
  real *a; 
  VALUE rb_b;
  real *b; 
  VALUE rb_safmin;
  real safmin; 
  VALUE rb_scale1;
  real scale1; 
  VALUE rb_scale2;
  real scale2; 
  VALUE rb_wr1;
  real wr1; 
  VALUE rb_wr2;
  real wr2; 
  VALUE rb_wi;
  real wi; 

  integer lda;
  integer ldb;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  scale1, scale2, wr1, wr2, wi = NumRu::Lapack.slag2( a, b, safmin)\n    or\n  NumRu::Lapack.slag2  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE SLAG2( A, LDA, B, LDB, SAFMIN, SCALE1, SCALE2, WR1, WR2, WI )\n\n*  Purpose\n*  =======\n*\n*  SLAG2 computes the eigenvalues of a 2 x 2 generalized eigenvalue\n*  problem  A - w B, with scaling as necessary to avoid over-/underflow.\n*\n*  The scaling factor \"s\" results in a modified eigenvalue equation\n*\n*      s A - w B\n*\n*  where  s  is a non-negative scaling factor chosen so that  w,  w B,\n*  and  s A  do not overflow and, if possible, do not underflow, either.\n*\n\n*  Arguments\n*  =========\n*\n*  A       (input) REAL array, dimension (LDA, 2)\n*          On entry, the 2 x 2 matrix A.  It is assumed that its 1-norm\n*          is less than 1/SAFMIN.  Entries less than\n*          sqrt(SAFMIN)*norm(A) are subject to being treated as zero.\n*\n*  LDA     (input) INTEGER\n*          The leading dimension of the array A.  LDA >= 2.\n*\n*  B       (input) REAL array, dimension (LDB, 2)\n*          On entry, the 2 x 2 upper triangular matrix B.  It is\n*          assumed that the one-norm of B is less than 1/SAFMIN.  The\n*          diagonals should be at least sqrt(SAFMIN) times the largest\n*          element of B (in absolute value); if a diagonal is smaller\n*          than that, then  +/- sqrt(SAFMIN) will be used instead of\n*          that diagonal.\n*\n*  LDB     (input) INTEGER\n*          The leading dimension of the array B.  LDB >= 2.\n*\n*  SAFMIN  (input) REAL\n*          The smallest positive number s.t. 1/SAFMIN does not\n*          overflow.  (This should always be SLAMCH('S') -- it is an\n*          argument in order to avoid having to call SLAMCH frequently.)\n*\n*  SCALE1  (output) REAL\n*          A scaling factor used to avoid over-/underflow in the\n*          eigenvalue equation which defines the first eigenvalue.  If\n*          the eigenvalues are complex, then the eigenvalues are\n*          ( WR1  +/-  WI i ) / SCALE1  (which may lie outside the\n*          exponent range of the machine), SCALE1=SCALE2, and SCALE1\n*          will always be positive.  If the eigenvalues are real, then\n*          the first (real) eigenvalue is  WR1 / SCALE1 , but this may\n*          overflow or underflow, and in fact, SCALE1 may be zero or\n*          less than the underflow threshhold if the exact eigenvalue\n*          is sufficiently large.\n*\n*  SCALE2  (output) REAL\n*          A scaling factor used to avoid over-/underflow in the\n*          eigenvalue equation which defines the second eigenvalue.  If\n*          the eigenvalues are complex, then SCALE2=SCALE1.  If the\n*          eigenvalues are real, then the second (real) eigenvalue is\n*          WR2 / SCALE2 , but this may overflow or underflow, and in\n*          fact, SCALE2 may be zero or less than the underflow\n*          threshhold if the exact eigenvalue is sufficiently large.\n*\n*  WR1     (output) REAL\n*          If the eigenvalue is real, then WR1 is SCALE1 times the\n*          eigenvalue closest to the (2,2) element of A B**(-1).  If the\n*          eigenvalue is complex, then WR1=WR2 is SCALE1 times the real\n*          part of the eigenvalues.\n*\n*  WR2     (output) REAL\n*          If the eigenvalue is real, then WR2 is SCALE2 times the\n*          other eigenvalue.  If the eigenvalue is complex, then\n*          WR1=WR2 is SCALE1 times the real part of the eigenvalues.\n*\n*  WI      (output) REAL\n*          If the eigenvalue is real, then WI is zero.  If the\n*          eigenvalue is complex, then WI is SCALE1 times the imaginary\n*          part of the eigenvalues.  WI will always be non-negative.\n*\n\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 3)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 3)", argc);
  rb_a = argv[0];
  rb_b = argv[1];
  rb_safmin = argv[2];

  if (!NA_IsNArray(rb_a))
    rb_raise(rb_eArgError, "a (1th argument) must be NArray");
  if (NA_RANK(rb_a) != 2)
    rb_raise(rb_eArgError, "rank of a (1th argument) must be %d", 2);
  if (NA_SHAPE1(rb_a) != (2))
    rb_raise(rb_eRuntimeError, "shape 1 of a must be %d", 2);
  lda = NA_SHAPE0(rb_a);
  if (NA_TYPE(rb_a) != NA_SFLOAT)
    rb_a = na_change_type(rb_a, NA_SFLOAT);
  a = NA_PTR_TYPE(rb_a, real*);
  if (!NA_IsNArray(rb_b))
    rb_raise(rb_eArgError, "b (2th argument) must be NArray");
  if (NA_RANK(rb_b) != 2)
    rb_raise(rb_eArgError, "rank of b (2th argument) must be %d", 2);
  if (NA_SHAPE1(rb_b) != (2))
    rb_raise(rb_eRuntimeError, "shape 1 of b must be %d", 2);
  ldb = NA_SHAPE0(rb_b);
  if (NA_TYPE(rb_b) != NA_SFLOAT)
    rb_b = na_change_type(rb_b, NA_SFLOAT);
  b = NA_PTR_TYPE(rb_b, real*);
  safmin = (real)NUM2DBL(rb_safmin);

  slag2_(a, &lda, b, &ldb, &safmin, &scale1, &scale2, &wr1, &wr2, &wi);

  rb_scale1 = rb_float_new((double)scale1);
  rb_scale2 = rb_float_new((double)scale2);
  rb_wr1 = rb_float_new((double)wr1);
  rb_wr2 = rb_float_new((double)wr2);
  rb_wi = rb_float_new((double)wi);
  return rb_ary_new3(5, rb_scale1, rb_scale2, rb_wr1, rb_wr2, rb_wi);
}

void
init_lapack_slag2(VALUE mLapack){
  rb_define_module_function(mLapack, "slag2", rb_slag2, -1);
}

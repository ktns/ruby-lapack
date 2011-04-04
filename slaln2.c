#include "rb_lapack.h"

extern VOID slaln2_(logical *ltrans, integer *na, integer *nw, real *smin, real *ca, real *a, integer *lda, real *d1, real *d2, real *b, integer *ldb, real *wr, real *wi, real *x, integer *ldx, real *scale, real *xnorm, integer *info);

static VALUE
rb_slaln2(int argc, VALUE *argv, VALUE self){
  VALUE rb_ltrans;
  logical ltrans; 
  VALUE rb_smin;
  real smin; 
  VALUE rb_ca;
  real ca; 
  VALUE rb_a;
  real *a; 
  VALUE rb_d1;
  real d1; 
  VALUE rb_d2;
  real d2; 
  VALUE rb_b;
  real *b; 
  VALUE rb_wr;
  real wr; 
  VALUE rb_wi;
  real wi; 
  VALUE rb_x;
  real *x; 
  VALUE rb_scale;
  real scale; 
  VALUE rb_xnorm;
  real xnorm; 
  VALUE rb_info;
  integer info; 

  integer lda;
  integer na;
  integer ldb;
  integer nw;
  integer ldx;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  x, scale, xnorm, info = NumRu::Lapack.slaln2( ltrans, smin, ca, a, d1, d2, b, wr, wi)\n    or\n  NumRu::Lapack.slaln2  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE SLALN2( LTRANS, NA, NW, SMIN, CA, A, LDA, D1, D2, B, LDB, WR, WI, X, LDX, SCALE, XNORM, INFO )\n\n*  Purpose\n*  =======\n*\n*  SLALN2 solves a system of the form  (ca A - w D ) X = s B\n*  or (ca A' - w D) X = s B   with possible scaling (\"s\") and\n*  perturbation of A.  (A' means A-transpose.)\n*\n*  A is an NA x NA real matrix, ca is a real scalar, D is an NA x NA\n*  real diagonal matrix, w is a real or complex value, and X and B are\n*  NA x 1 matrices -- real if w is real, complex if w is complex.  NA\n*  may be 1 or 2.\n*\n*  If w is complex, X and B are represented as NA x 2 matrices,\n*  the first column of each being the real part and the second\n*  being the imaginary part.\n*\n*  \"s\" is a scaling factor (.LE. 1), computed by SLALN2, which is\n*  so chosen that X can be computed without overflow.  X is further\n*  scaled if necessary to assure that norm(ca A - w D)*norm(X) is less\n*  than overflow.\n*\n*  If both singular values of (ca A - w D) are less than SMIN,\n*  SMIN*identity will be used instead of (ca A - w D).  If only one\n*  singular value is less than SMIN, one element of (ca A - w D) will be\n*  perturbed enough to make the smallest singular value roughly SMIN.\n*  If both singular values are at least SMIN, (ca A - w D) will not be\n*  perturbed.  In any case, the perturbation will be at most some small\n*  multiple of max( SMIN, ulp*norm(ca A - w D) ).  The singular values\n*  are computed by infinity-norm approximations, and thus will only be\n*  correct to a factor of 2 or so.\n*\n*  Note: all input quantities are assumed to be smaller than overflow\n*  by a reasonable factor.  (See BIGNUM.)\n*\n\n*  Arguments\n*  ==========\n*\n*  LTRANS  (input) LOGICAL\n*          =.TRUE.:  A-transpose will be used.\n*          =.FALSE.: A will be used (not transposed.)\n*\n*  NA      (input) INTEGER\n*          The size of the matrix A.  It may (only) be 1 or 2.\n*\n*  NW      (input) INTEGER\n*          1 if \"w\" is real, 2 if \"w\" is complex.  It may only be 1\n*          or 2.\n*\n*  SMIN    (input) REAL\n*          The desired lower bound on the singular values of A.  This\n*          should be a safe distance away from underflow or overflow,\n*          say, between (underflow/machine precision) and  (machine\n*          precision * overflow ).  (See BIGNUM and ULP.)\n*\n*  CA      (input) REAL\n*          The coefficient c, which A is multiplied by.\n*\n*  A       (input) REAL array, dimension (LDA,NA)\n*          The NA x NA matrix A.\n*\n*  LDA     (input) INTEGER\n*          The leading dimension of A.  It must be at least NA.\n*\n*  D1      (input) REAL\n*          The 1,1 element in the diagonal matrix D.\n*\n*  D2      (input) REAL\n*          The 2,2 element in the diagonal matrix D.  Not used if NW=1.\n*\n*  B       (input) REAL array, dimension (LDB,NW)\n*          The NA x NW matrix B (right-hand side).  If NW=2 (\"w\" is\n*          complex), column 1 contains the real part of B and column 2\n*          contains the imaginary part.\n*\n*  LDB     (input) INTEGER\n*          The leading dimension of B.  It must be at least NA.\n*\n*  WR      (input) REAL\n*          The real part of the scalar \"w\".\n*\n*  WI      (input) REAL\n*          The imaginary part of the scalar \"w\".  Not used if NW=1.\n*\n*  X       (output) REAL array, dimension (LDX,NW)\n*          The NA x NW matrix X (unknowns), as computed by SLALN2.\n*          If NW=2 (\"w\" is complex), on exit, column 1 will contain\n*          the real part of X and column 2 will contain the imaginary\n*          part.\n*\n*  LDX     (input) INTEGER\n*          The leading dimension of X.  It must be at least NA.\n*\n*  SCALE   (output) REAL\n*          The scale factor that B must be multiplied by to insure\n*          that overflow does not occur when computing X.  Thus,\n*          (ca A - w D) X  will be SCALE*B, not B (ignoring\n*          perturbations of A.)  It will be at most 1.\n*\n*  XNORM   (output) REAL\n*          The infinity-norm of X, when X is regarded as an NA x NW\n*          real matrix.\n*\n*  INFO    (output) INTEGER\n*          An error flag.  It will be set to zero if no error occurs,\n*          a negative number if an argument is in error, or a positive\n*          number if  ca A - w D  had to be perturbed.\n*          The possible values are:\n*          = 0: No error occurred, and (ca A - w D) did not have to be\n*                 perturbed.\n*          = 1: (ca A - w D) had to be perturbed to make its smallest\n*               (or only) singular value greater than SMIN.\n*          NOTE: In the interests of speed, this routine does not\n*                check the inputs for errors.\n*\n\n* =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 9)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 9)", argc);
  rb_ltrans = argv[0];
  rb_smin = argv[1];
  rb_ca = argv[2];
  rb_a = argv[3];
  rb_d1 = argv[4];
  rb_d2 = argv[5];
  rb_b = argv[6];
  rb_wr = argv[7];
  rb_wi = argv[8];

  smin = (real)NUM2DBL(rb_smin);
  if (!NA_IsNArray(rb_a))
    rb_raise(rb_eArgError, "a (4th argument) must be NArray");
  if (NA_RANK(rb_a) != 2)
    rb_raise(rb_eArgError, "rank of a (4th argument) must be %d", 2);
  na = NA_SHAPE1(rb_a);
  lda = NA_SHAPE0(rb_a);
  if (NA_TYPE(rb_a) != NA_SFLOAT)
    rb_a = na_change_type(rb_a, NA_SFLOAT);
  a = NA_PTR_TYPE(rb_a, real*);
  if (!NA_IsNArray(rb_b))
    rb_raise(rb_eArgError, "b (7th argument) must be NArray");
  if (NA_RANK(rb_b) != 2)
    rb_raise(rb_eArgError, "rank of b (7th argument) must be %d", 2);
  nw = NA_SHAPE1(rb_b);
  ldb = NA_SHAPE0(rb_b);
  if (NA_TYPE(rb_b) != NA_SFLOAT)
    rb_b = na_change_type(rb_b, NA_SFLOAT);
  b = NA_PTR_TYPE(rb_b, real*);
  d1 = (real)NUM2DBL(rb_d1);
  d2 = (real)NUM2DBL(rb_d2);
  ca = (real)NUM2DBL(rb_ca);
  ltrans = (rb_ltrans == Qtrue);
  wi = (real)NUM2DBL(rb_wi);
  wr = (real)NUM2DBL(rb_wr);
  ldx = na;
  {
    int shape[2];
    shape[0] = ldx;
    shape[1] = nw;
    rb_x = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  x = NA_PTR_TYPE(rb_x, real*);

  slaln2_(&ltrans, &na, &nw, &smin, &ca, a, &lda, &d1, &d2, b, &ldb, &wr, &wi, x, &ldx, &scale, &xnorm, &info);

  rb_scale = rb_float_new((double)scale);
  rb_xnorm = rb_float_new((double)xnorm);
  rb_info = INT2NUM(info);
  return rb_ary_new3(4, rb_x, rb_scale, rb_xnorm, rb_info);
}

void
init_lapack_slaln2(VALUE mLapack){
  rb_define_module_function(mLapack, "slaln2", rb_slaln2, -1);
}

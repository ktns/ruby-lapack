#include "rb_lapack.h"

extern VOID dlaln2_(logical *ltrans, integer *na, integer *nw, doublereal *smin, doublereal *ca, doublereal *a, integer *lda, doublereal *d1, doublereal *d2, doublereal *b, integer *ldb, doublereal *wr, doublereal *wi, doublereal *x, integer *ldx, doublereal *scale, doublereal *xnorm, integer *info);

static VALUE sHelp, sUsage;

static VALUE
rblapack_dlaln2(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_ltrans;
  logical ltrans; 
  VALUE rblapack_smin;
  doublereal smin; 
  VALUE rblapack_ca;
  doublereal ca; 
  VALUE rblapack_a;
  doublereal *a; 
  VALUE rblapack_d1;
  doublereal d1; 
  VALUE rblapack_d2;
  doublereal d2; 
  VALUE rblapack_b;
  doublereal *b; 
  VALUE rblapack_wr;
  doublereal wr; 
  VALUE rblapack_wi;
  doublereal wi; 
  VALUE rblapack_x;
  doublereal *x; 
  VALUE rblapack_scale;
  doublereal scale; 
  VALUE rblapack_xnorm;
  doublereal xnorm; 
  VALUE rblapack_info;
  integer info; 

  integer lda;
  integer na;
  integer ldb;
  integer nw;
  integer ldx;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  x, scale, xnorm, info = NumRu::Lapack.dlaln2( ltrans, smin, ca, a, d1, d2, b, wr, wi, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE DLALN2( LTRANS, NA, NW, SMIN, CA, A, LDA, D1, D2, B, LDB, WR, WI, X, LDX, SCALE, XNORM, INFO )\n\n*  Purpose\n*  =======\n*\n*  DLALN2 solves a system of the form  (ca A - w D ) X = s B\n*  or (ca A' - w D) X = s B   with possible scaling (\"s\") and\n*  perturbation of A.  (A' means A-transpose.)\n*\n*  A is an NA x NA real matrix, ca is a real scalar, D is an NA x NA\n*  real diagonal matrix, w is a real or complex value, and X and B are\n*  NA x 1 matrices -- real if w is real, complex if w is complex.  NA\n*  may be 1 or 2.\n*\n*  If w is complex, X and B are represented as NA x 2 matrices,\n*  the first column of each being the real part and the second\n*  being the imaginary part.\n*\n*  \"s\" is a scaling factor (.LE. 1), computed by DLALN2, which is\n*  so chosen that X can be computed without overflow.  X is further\n*  scaled if necessary to assure that norm(ca A - w D)*norm(X) is less\n*  than overflow.\n*\n*  If both singular values of (ca A - w D) are less than SMIN,\n*  SMIN*identity will be used instead of (ca A - w D).  If only one\n*  singular value is less than SMIN, one element of (ca A - w D) will be\n*  perturbed enough to make the smallest singular value roughly SMIN.\n*  If both singular values are at least SMIN, (ca A - w D) will not be\n*  perturbed.  In any case, the perturbation will be at most some small\n*  multiple of max( SMIN, ulp*norm(ca A - w D) ).  The singular values\n*  are computed by infinity-norm approximations, and thus will only be\n*  correct to a factor of 2 or so.\n*\n*  Note: all input quantities are assumed to be smaller than overflow\n*  by a reasonable factor.  (See BIGNUM.)\n*\n\n*  Arguments\n*  ==========\n*\n*  LTRANS  (input) LOGICAL\n*          =.TRUE.:  A-transpose will be used.\n*          =.FALSE.: A will be used (not transposed.)\n*\n*  NA      (input) INTEGER\n*          The size of the matrix A.  It may (only) be 1 or 2.\n*\n*  NW      (input) INTEGER\n*          1 if \"w\" is real, 2 if \"w\" is complex.  It may only be 1\n*          or 2.\n*\n*  SMIN    (input) DOUBLE PRECISION\n*          The desired lower bound on the singular values of A.  This\n*          should be a safe distance away from underflow or overflow,\n*          say, between (underflow/machine precision) and  (machine\n*          precision * overflow ).  (See BIGNUM and ULP.)\n*\n*  CA      (input) DOUBLE PRECISION\n*          The coefficient c, which A is multiplied by.\n*\n*  A       (input) DOUBLE PRECISION array, dimension (LDA,NA)\n*          The NA x NA matrix A.\n*\n*  LDA     (input) INTEGER\n*          The leading dimension of A.  It must be at least NA.\n*\n*  D1      (input) DOUBLE PRECISION\n*          The 1,1 element in the diagonal matrix D.\n*\n*  D2      (input) DOUBLE PRECISION\n*          The 2,2 element in the diagonal matrix D.  Not used if NW=1.\n*\n*  B       (input) DOUBLE PRECISION array, dimension (LDB,NW)\n*          The NA x NW matrix B (right-hand side).  If NW=2 (\"w\" is\n*          complex), column 1 contains the real part of B and column 2\n*          contains the imaginary part.\n*\n*  LDB     (input) INTEGER\n*          The leading dimension of B.  It must be at least NA.\n*\n*  WR      (input) DOUBLE PRECISION\n*          The real part of the scalar \"w\".\n*\n*  WI      (input) DOUBLE PRECISION\n*          The imaginary part of the scalar \"w\".  Not used if NW=1.\n*\n*  X       (output) DOUBLE PRECISION array, dimension (LDX,NW)\n*          The NA x NW matrix X (unknowns), as computed by DLALN2.\n*          If NW=2 (\"w\" is complex), on exit, column 1 will contain\n*          the real part of X and column 2 will contain the imaginary\n*          part.\n*\n*  LDX     (input) INTEGER\n*          The leading dimension of X.  It must be at least NA.\n*\n*  SCALE   (output) DOUBLE PRECISION\n*          The scale factor that B must be multiplied by to insure\n*          that overflow does not occur when computing X.  Thus,\n*          (ca A - w D) X  will be SCALE*B, not B (ignoring\n*          perturbations of A.)  It will be at most 1.\n*\n*  XNORM   (output) DOUBLE PRECISION\n*          The infinity-norm of X, when X is regarded as an NA x NW\n*          real matrix.\n*\n*  INFO    (output) INTEGER\n*          An error flag.  It will be set to zero if no error occurs,\n*          a negative number if an argument is in error, or a positive\n*          number if  ca A - w D  had to be perturbed.\n*          The possible values are:\n*          = 0: No error occurred, and (ca A - w D) did not have to be\n*                 perturbed.\n*          = 1: (ca A - w D) had to be perturbed to make its smallest\n*               (or only) singular value greater than SMIN.\n*          NOTE: In the interests of speed, this routine does not\n*                check the inputs for errors.\n*\n\n* =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  x, scale, xnorm, info = NumRu::Lapack.dlaln2( ltrans, smin, ca, a, d1, d2, b, wr, wi, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 9)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 9)", argc);
  rblapack_ltrans = argv[0];
  rblapack_smin = argv[1];
  rblapack_ca = argv[2];
  rblapack_a = argv[3];
  rblapack_d1 = argv[4];
  rblapack_d2 = argv[5];
  rblapack_b = argv[6];
  rblapack_wr = argv[7];
  rblapack_wi = argv[8];
  if (rb_options != Qnil) {
  }

  smin = NUM2DBL(rblapack_smin);
  if (!NA_IsNArray(rblapack_a))
    rb_raise(rb_eArgError, "a (4th argument) must be NArray");
  if (NA_RANK(rblapack_a) != 2)
    rb_raise(rb_eArgError, "rank of a (4th argument) must be %d", 2);
  na = NA_SHAPE1(rblapack_a);
  lda = NA_SHAPE0(rblapack_a);
  if (NA_TYPE(rblapack_a) != NA_DFLOAT)
    rblapack_a = na_change_type(rblapack_a, NA_DFLOAT);
  a = NA_PTR_TYPE(rblapack_a, doublereal*);
  if (!NA_IsNArray(rblapack_b))
    rb_raise(rb_eArgError, "b (7th argument) must be NArray");
  if (NA_RANK(rblapack_b) != 2)
    rb_raise(rb_eArgError, "rank of b (7th argument) must be %d", 2);
  nw = NA_SHAPE1(rblapack_b);
  ldb = NA_SHAPE0(rblapack_b);
  if (NA_TYPE(rblapack_b) != NA_DFLOAT)
    rblapack_b = na_change_type(rblapack_b, NA_DFLOAT);
  b = NA_PTR_TYPE(rblapack_b, doublereal*);
  d1 = NUM2DBL(rblapack_d1);
  d2 = NUM2DBL(rblapack_d2);
  ca = NUM2DBL(rblapack_ca);
  ltrans = (rblapack_ltrans == Qtrue);
  wi = NUM2DBL(rblapack_wi);
  wr = NUM2DBL(rblapack_wr);
  ldx = na;
  {
    int shape[2];
    shape[0] = ldx;
    shape[1] = nw;
    rblapack_x = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  x = NA_PTR_TYPE(rblapack_x, doublereal*);

  dlaln2_(&ltrans, &na, &nw, &smin, &ca, a, &lda, &d1, &d2, b, &ldb, &wr, &wi, x, &ldx, &scale, &xnorm, &info);

  rblapack_scale = rb_float_new((double)scale);
  rblapack_xnorm = rb_float_new((double)xnorm);
  rblapack_info = INT2NUM(info);
  return rb_ary_new3(4, rblapack_x, rblapack_scale, rblapack_xnorm, rblapack_info);
}

void
init_lapack_dlaln2(VALUE mLapack){
  rb_define_module_function(mLapack, "dlaln2", rblapack_dlaln2, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}

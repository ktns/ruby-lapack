#include "rb_lapack.h"

extern integer dlaneg_(integer *n, doublereal *d, doublereal *lld, doublereal *sigma, doublereal *pivmin, integer *r);

static VALUE
rb_dlaneg(int argc, VALUE *argv, VALUE self){
  VALUE rb_d;
  doublereal *d; 
  VALUE rb_lld;
  doublereal *lld; 
  VALUE rb_sigma;
  doublereal sigma; 
  VALUE rb_pivmin;
  doublereal pivmin; 
  VALUE rb_r;
  integer r; 
  VALUE rb___out__;
  integer __out__; 

  integer n;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  __out__ = NumRu::Lapack.dlaneg( d, lld, sigma, pivmin, r)\n    or\n  NumRu::Lapack.dlaneg  # print help\n\n\nFORTRAN MANUAL\n      INTEGER FUNCTION DLANEG( N, D, LLD, SIGMA, PIVMIN, R )\n\n*  Purpose\n*  =======\n*\n*  DLANEG computes the Sturm count, the number of negative pivots\n*  encountered while factoring tridiagonal T - sigma I = L D L^T.\n*  This implementation works directly on the factors without forming\n*  the tridiagonal matrix T.  The Sturm count is also the number of\n*  eigenvalues of T less than sigma.\n*\n*  This routine is called from DLARRB.\n*\n*  The current routine does not use the PIVMIN parameter but rather\n*  requires IEEE-754 propagation of Infinities and NaNs.  This\n*  routine also has no input range restrictions but does require\n*  default exception handling such that x/0 produces Inf when x is\n*  non-zero, and Inf/Inf produces NaN.  For more information, see:\n*\n*    Marques, Riedy, and Voemel, \"Benefits of IEEE-754 Features in\n*    Modern Symmetric Tridiagonal Eigensolvers,\" SIAM Journal on\n*    Scientific Computing, v28, n5, 2006.  DOI 10.1137/050641624\n*    (Tech report version in LAWN 172 with the same title.)\n*\n\n*  Arguments\n*  =========\n*\n*  N       (input) INTEGER\n*          The order of the matrix.\n*\n*  D       (input) DOUBLE PRECISION array, dimension (N)\n*          The N diagonal elements of the diagonal matrix D.\n*\n*  LLD     (input) DOUBLE PRECISION array, dimension (N-1)\n*          The (N-1) elements L(i)*L(i)*D(i).\n*\n*  SIGMA   (input) DOUBLE PRECISION\n*          Shift amount in T - sigma I = L D L^T.\n*\n*  PIVMIN  (input) DOUBLE PRECISION\n*          The minimum pivot in the Sturm sequence.  May be used\n*          when zero pivots are encountered on non-IEEE-754\n*          architectures.\n*\n*  R       (input) INTEGER\n*          The twist index for the twisted factorization that is used\n*          for the negcount.\n*\n\n*  Further Details\n*  ===============\n*\n*  Based on contributions by\n*     Osni Marques, LBNL/NERSC, USA\n*     Christof Voemel, University of California, Berkeley, USA\n*     Jason Riedy, University of California, Berkeley, USA\n*\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 5)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 5)", argc);
  rb_d = argv[0];
  rb_lld = argv[1];
  rb_sigma = argv[2];
  rb_pivmin = argv[3];
  rb_r = argv[4];

  pivmin = NUM2DBL(rb_pivmin);
  if (!NA_IsNArray(rb_d))
    rb_raise(rb_eArgError, "d (1th argument) must be NArray");
  if (NA_RANK(rb_d) != 1)
    rb_raise(rb_eArgError, "rank of d (1th argument) must be %d", 1);
  n = NA_SHAPE0(rb_d);
  if (NA_TYPE(rb_d) != NA_DFLOAT)
    rb_d = na_change_type(rb_d, NA_DFLOAT);
  d = NA_PTR_TYPE(rb_d, doublereal*);
  sigma = NUM2DBL(rb_sigma);
  r = NUM2INT(rb_r);
  if (!NA_IsNArray(rb_lld))
    rb_raise(rb_eArgError, "lld (2th argument) must be NArray");
  if (NA_RANK(rb_lld) != 1)
    rb_raise(rb_eArgError, "rank of lld (2th argument) must be %d", 1);
  if (NA_SHAPE0(rb_lld) != (n-1))
    rb_raise(rb_eRuntimeError, "shape 0 of lld must be %d", n-1);
  if (NA_TYPE(rb_lld) != NA_DFLOAT)
    rb_lld = na_change_type(rb_lld, NA_DFLOAT);
  lld = NA_PTR_TYPE(rb_lld, doublereal*);

  __out__ = dlaneg_(&n, d, lld, &sigma, &pivmin, &r);

  rb___out__ = INT2NUM(__out__);
  return rb___out__;
}

void
init_lapack_dlaneg(VALUE mLapack){
  rb_define_module_function(mLapack, "dlaneg", rb_dlaneg, -1);
}

#include "rb_lapack.h"

extern VOID dlar1v_(integer *n, integer *b1, integer *bn, doublereal *lambda, doublereal *d, doublereal *l, doublereal *ld, doublereal *lld, doublereal *pivmin, doublereal *gaptol, doublereal *z, logical *wantnc, integer *negcnt, doublereal *ztz, doublereal *mingma, integer *r, integer *isuppz, doublereal *nrminv, doublereal *resid, doublereal *rqcorr, doublereal *work);

static VALUE sHelp, sUsage;

static VALUE
rblapack_dlar1v(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_b1;
  integer b1; 
  VALUE rblapack_bn;
  integer bn; 
  VALUE rblapack_lambda;
  doublereal lambda; 
  VALUE rblapack_d;
  doublereal *d; 
  VALUE rblapack_l;
  doublereal *l; 
  VALUE rblapack_ld;
  doublereal *ld; 
  VALUE rblapack_lld;
  doublereal *lld; 
  VALUE rblapack_pivmin;
  doublereal pivmin; 
  VALUE rblapack_gaptol;
  doublereal gaptol; 
  VALUE rblapack_z;
  doublereal *z; 
  VALUE rblapack_wantnc;
  logical wantnc; 
  VALUE rblapack_r;
  integer r; 
  VALUE rblapack_negcnt;
  integer negcnt; 
  VALUE rblapack_ztz;
  doublereal ztz; 
  VALUE rblapack_mingma;
  doublereal mingma; 
  VALUE rblapack_isuppz;
  integer *isuppz; 
  VALUE rblapack_nrminv;
  doublereal nrminv; 
  VALUE rblapack_resid;
  doublereal resid; 
  VALUE rblapack_rqcorr;
  doublereal rqcorr; 
  VALUE rblapack_z_out__;
  doublereal *z_out__;
  doublereal *work;

  integer n;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  negcnt, ztz, mingma, isuppz, nrminv, resid, rqcorr, z, r = NumRu::Lapack.dlar1v( b1, bn, lambda, d, l, ld, lld, pivmin, gaptol, z, wantnc, r, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE DLAR1V( N, B1, BN, LAMBDA, D, L, LD, LLD, PIVMIN, GAPTOL, Z, WANTNC, NEGCNT, ZTZ, MINGMA, R, ISUPPZ, NRMINV, RESID, RQCORR, WORK )\n\n*  Purpose\n*  =======\n*\n*  DLAR1V computes the (scaled) r-th column of the inverse of\n*  the sumbmatrix in rows B1 through BN of the tridiagonal matrix\n*  L D L^T - sigma I. When sigma is close to an eigenvalue, the\n*  computed vector is an accurate eigenvector. Usually, r corresponds\n*  to the index where the eigenvector is largest in magnitude.\n*  The following steps accomplish this computation :\n*  (a) Stationary qd transform,  L D L^T - sigma I = L(+) D(+) L(+)^T,\n*  (b) Progressive qd transform, L D L^T - sigma I = U(-) D(-) U(-)^T,\n*  (c) Computation of the diagonal elements of the inverse of\n*      L D L^T - sigma I by combining the above transforms, and choosing\n*      r as the index where the diagonal of the inverse is (one of the)\n*      largest in magnitude.\n*  (d) Computation of the (scaled) r-th column of the inverse using the\n*      twisted factorization obtained by combining the top part of the\n*      the stationary and the bottom part of the progressive transform.\n*\n\n*  Arguments\n*  =========\n*\n*  N        (input) INTEGER\n*           The order of the matrix L D L^T.\n*\n*  B1       (input) INTEGER\n*           First index of the submatrix of L D L^T.\n*\n*  BN       (input) INTEGER\n*           Last index of the submatrix of L D L^T.\n*\n*  LAMBDA    (input) DOUBLE PRECISION\n*           The shift. In order to compute an accurate eigenvector,\n*           LAMBDA should be a good approximation to an eigenvalue\n*           of L D L^T.\n*\n*  L        (input) DOUBLE PRECISION array, dimension (N-1)\n*           The (n-1) subdiagonal elements of the unit bidiagonal matrix\n*           L, in elements 1 to N-1.\n*\n*  D        (input) DOUBLE PRECISION array, dimension (N)\n*           The n diagonal elements of the diagonal matrix D.\n*\n*  LD       (input) DOUBLE PRECISION array, dimension (N-1)\n*           The n-1 elements L(i)*D(i).\n*\n*  LLD      (input) DOUBLE PRECISION array, dimension (N-1)\n*           The n-1 elements L(i)*L(i)*D(i).\n*\n*  PIVMIN   (input) DOUBLE PRECISION\n*           The minimum pivot in the Sturm sequence.\n*\n*  GAPTOL   (input) DOUBLE PRECISION\n*           Tolerance that indicates when eigenvector entries are negligible\n*           w.r.t. their contribution to the residual.\n*\n*  Z        (input/output) DOUBLE PRECISION array, dimension (N)\n*           On input, all entries of Z must be set to 0.\n*           On output, Z contains the (scaled) r-th column of the\n*           inverse. The scaling is such that Z(R) equals 1.\n*\n*  WANTNC   (input) LOGICAL\n*           Specifies whether NEGCNT has to be computed.\n*\n*  NEGCNT   (output) INTEGER\n*           If WANTNC is .TRUE. then NEGCNT = the number of pivots < pivmin\n*           in the  matrix factorization L D L^T, and NEGCNT = -1 otherwise.\n*\n*  ZTZ      (output) DOUBLE PRECISION\n*           The square of the 2-norm of Z.\n*\n*  MINGMA   (output) DOUBLE PRECISION\n*           The reciprocal of the largest (in magnitude) diagonal\n*           element of the inverse of L D L^T - sigma I.\n*\n*  R        (input/output) INTEGER\n*           The twist index for the twisted factorization used to\n*           compute Z.\n*           On input, 0 <= R <= N. If R is input as 0, R is set to\n*           the index where (L D L^T - sigma I)^{-1} is largest\n*           in magnitude. If 1 <= R <= N, R is unchanged.\n*           On output, R contains the twist index used to compute Z.\n*           Ideally, R designates the position of the maximum entry in the\n*           eigenvector.\n*\n*  ISUPPZ   (output) INTEGER array, dimension (2)\n*           The support of the vector in Z, i.e., the vector Z is\n*           nonzero only in elements ISUPPZ(1) through ISUPPZ( 2 ).\n*\n*  NRMINV   (output) DOUBLE PRECISION\n*           NRMINV = 1/SQRT( ZTZ )\n*\n*  RESID    (output) DOUBLE PRECISION\n*           The residual of the FP vector.\n*           RESID = ABS( MINGMA )/SQRT( ZTZ )\n*\n*  RQCORR   (output) DOUBLE PRECISION\n*           The Rayleigh Quotient correction to LAMBDA.\n*           RQCORR = MINGMA*TMP\n*\n*  WORK     (workspace) DOUBLE PRECISION array, dimension (4*N)\n*\n\n*  Further Details\n*  ===============\n*\n*  Based on contributions by\n*     Beresford Parlett, University of California, Berkeley, USA\n*     Jim Demmel, University of California, Berkeley, USA\n*     Inderjit Dhillon, University of Texas, Austin, USA\n*     Osni Marques, LBNL/NERSC, USA\n*     Christof Voemel, University of California, Berkeley, USA\n*\n*  =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  negcnt, ztz, mingma, isuppz, nrminv, resid, rqcorr, z, r = NumRu::Lapack.dlar1v( b1, bn, lambda, d, l, ld, lld, pivmin, gaptol, z, wantnc, r, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 12)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 12)", argc);
  rblapack_b1 = argv[0];
  rblapack_bn = argv[1];
  rblapack_lambda = argv[2];
  rblapack_d = argv[3];
  rblapack_l = argv[4];
  rblapack_ld = argv[5];
  rblapack_lld = argv[6];
  rblapack_pivmin = argv[7];
  rblapack_gaptol = argv[8];
  rblapack_z = argv[9];
  rblapack_wantnc = argv[10];
  rblapack_r = argv[11];
  if (rb_options != Qnil) {
  }

  pivmin = NUM2DBL(rblapack_pivmin);
  bn = NUM2INT(rblapack_bn);
  lambda = NUM2DBL(rblapack_lambda);
  if (!NA_IsNArray(rblapack_z))
    rb_raise(rb_eArgError, "z (10th argument) must be NArray");
  if (NA_RANK(rblapack_z) != 1)
    rb_raise(rb_eArgError, "rank of z (10th argument) must be %d", 1);
  n = NA_SHAPE0(rblapack_z);
  if (NA_TYPE(rblapack_z) != NA_DFLOAT)
    rblapack_z = na_change_type(rblapack_z, NA_DFLOAT);
  z = NA_PTR_TYPE(rblapack_z, doublereal*);
  wantnc = (rblapack_wantnc == Qtrue);
  if (!NA_IsNArray(rblapack_d))
    rb_raise(rb_eArgError, "d (4th argument) must be NArray");
  if (NA_RANK(rblapack_d) != 1)
    rb_raise(rb_eArgError, "rank of d (4th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_d) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of d must be the same as shape 0 of z");
  if (NA_TYPE(rblapack_d) != NA_DFLOAT)
    rblapack_d = na_change_type(rblapack_d, NA_DFLOAT);
  d = NA_PTR_TYPE(rblapack_d, doublereal*);
  r = NUM2INT(rblapack_r);
  gaptol = NUM2DBL(rblapack_gaptol);
  b1 = NUM2INT(rblapack_b1);
  if (!NA_IsNArray(rblapack_lld))
    rb_raise(rb_eArgError, "lld (7th argument) must be NArray");
  if (NA_RANK(rblapack_lld) != 1)
    rb_raise(rb_eArgError, "rank of lld (7th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_lld) != (n-1))
    rb_raise(rb_eRuntimeError, "shape 0 of lld must be %d", n-1);
  if (NA_TYPE(rblapack_lld) != NA_DFLOAT)
    rblapack_lld = na_change_type(rblapack_lld, NA_DFLOAT);
  lld = NA_PTR_TYPE(rblapack_lld, doublereal*);
  if (!NA_IsNArray(rblapack_ld))
    rb_raise(rb_eArgError, "ld (6th argument) must be NArray");
  if (NA_RANK(rblapack_ld) != 1)
    rb_raise(rb_eArgError, "rank of ld (6th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_ld) != (n-1))
    rb_raise(rb_eRuntimeError, "shape 0 of ld must be %d", n-1);
  if (NA_TYPE(rblapack_ld) != NA_DFLOAT)
    rblapack_ld = na_change_type(rblapack_ld, NA_DFLOAT);
  ld = NA_PTR_TYPE(rblapack_ld, doublereal*);
  if (!NA_IsNArray(rblapack_l))
    rb_raise(rb_eArgError, "l (5th argument) must be NArray");
  if (NA_RANK(rblapack_l) != 1)
    rb_raise(rb_eArgError, "rank of l (5th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_l) != (n-1))
    rb_raise(rb_eRuntimeError, "shape 0 of l must be %d", n-1);
  if (NA_TYPE(rblapack_l) != NA_DFLOAT)
    rblapack_l = na_change_type(rblapack_l, NA_DFLOAT);
  l = NA_PTR_TYPE(rblapack_l, doublereal*);
  {
    int shape[1];
    shape[0] = 2;
    rblapack_isuppz = na_make_object(NA_LINT, 1, shape, cNArray);
  }
  isuppz = NA_PTR_TYPE(rblapack_isuppz, integer*);
  {
    int shape[1];
    shape[0] = n;
    rblapack_z_out__ = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  z_out__ = NA_PTR_TYPE(rblapack_z_out__, doublereal*);
  MEMCPY(z_out__, z, doublereal, NA_TOTAL(rblapack_z));
  rblapack_z = rblapack_z_out__;
  z = z_out__;
  work = ALLOC_N(doublereal, (4*n));

  dlar1v_(&n, &b1, &bn, &lambda, d, l, ld, lld, &pivmin, &gaptol, z, &wantnc, &negcnt, &ztz, &mingma, &r, isuppz, &nrminv, &resid, &rqcorr, work);

  free(work);
  rblapack_negcnt = INT2NUM(negcnt);
  rblapack_ztz = rb_float_new((double)ztz);
  rblapack_mingma = rb_float_new((double)mingma);
  rblapack_nrminv = rb_float_new((double)nrminv);
  rblapack_resid = rb_float_new((double)resid);
  rblapack_rqcorr = rb_float_new((double)rqcorr);
  rblapack_r = INT2NUM(r);
  return rb_ary_new3(9, rblapack_negcnt, rblapack_ztz, rblapack_mingma, rblapack_isuppz, rblapack_nrminv, rblapack_resid, rblapack_rqcorr, rblapack_z, rblapack_r);
}

void
init_lapack_dlar1v(VALUE mLapack){
  rb_define_module_function(mLapack, "dlar1v", rblapack_dlar1v, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}

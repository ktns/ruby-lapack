#include "rb_lapack.h"

static VALUE
rb_slar1v(int argc, VALUE *argv, VALUE self){
  VALUE rb_b1;
  integer b1; 
  VALUE rb_bn;
  integer bn; 
  VALUE rb_lambda;
  real lambda; 
  VALUE rb_d;
  real *d; 
  VALUE rb_l;
  real *l; 
  VALUE rb_ld;
  real *ld; 
  VALUE rb_lld;
  real *lld; 
  VALUE rb_pivmin;
  real pivmin; 
  VALUE rb_gaptol;
  real gaptol; 
  VALUE rb_z;
  real *z; 
  VALUE rb_wantnc;
  logical wantnc; 
  VALUE rb_r;
  integer r; 
  VALUE rb_negcnt;
  integer negcnt; 
  VALUE rb_ztz;
  real ztz; 
  VALUE rb_mingma;
  real mingma; 
  VALUE rb_isuppz;
  integer *isuppz; 
  VALUE rb_nrminv;
  real nrminv; 
  VALUE rb_resid;
  real resid; 
  VALUE rb_rqcorr;
  real rqcorr; 
  VALUE rb_z_out__;
  real *z_out__;
  real *work;

  integer n;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  negcnt, ztz, mingma, isuppz, nrminv, resid, rqcorr, z, r = NumRu::Lapack.slar1v( b1, bn, lambda, d, l, ld, lld, pivmin, gaptol, z, wantnc, r)\n    or\n  NumRu::Lapack.slar1v  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE SLAR1V( N, B1, BN, LAMBDA, D, L, LD, LLD, PIVMIN, GAPTOL, Z, WANTNC, NEGCNT, ZTZ, MINGMA, R, ISUPPZ, NRMINV, RESID, RQCORR, WORK )\n\n*  Purpose\n*  =======\n*\n*  SLAR1V computes the (scaled) r-th column of the inverse of\n*  the sumbmatrix in rows B1 through BN of the tridiagonal matrix\n*  L D L^T - sigma I. When sigma is close to an eigenvalue, the\n*  computed vector is an accurate eigenvector. Usually, r corresponds\n*  to the index where the eigenvector is largest in magnitude.\n*  The following steps accomplish this computation :\n*  (a) Stationary qd transform,  L D L^T - sigma I = L(+) D(+) L(+)^T,\n*  (b) Progressive qd transform, L D L^T - sigma I = U(-) D(-) U(-)^T,\n*  (c) Computation of the diagonal elements of the inverse of\n*      L D L^T - sigma I by combining the above transforms, and choosing\n*      r as the index where the diagonal of the inverse is (one of the)\n*      largest in magnitude.\n*  (d) Computation of the (scaled) r-th column of the inverse using the\n*      twisted factorization obtained by combining the top part of the\n*      the stationary and the bottom part of the progressive transform.\n*\n\n*  Arguments\n*  =========\n*\n*  N        (input) INTEGER\n*           The order of the matrix L D L^T.\n*\n*  B1       (input) INTEGER\n*           First index of the submatrix of L D L^T.\n*\n*  BN       (input) INTEGER\n*           Last index of the submatrix of L D L^T.\n*\n*  LAMBDA    (input) REAL            \n*           The shift. In order to compute an accurate eigenvector,\n*           LAMBDA should be a good approximation to an eigenvalue\n*           of L D L^T.\n*\n*  L        (input) REAL             array, dimension (N-1)\n*           The (n-1) subdiagonal elements of the unit bidiagonal matrix\n*           L, in elements 1 to N-1.\n*\n*  D        (input) REAL             array, dimension (N)\n*           The n diagonal elements of the diagonal matrix D.\n*\n*  LD       (input) REAL             array, dimension (N-1)\n*           The n-1 elements L(i)*D(i).\n*\n*  LLD      (input) REAL             array, dimension (N-1)\n*           The n-1 elements L(i)*L(i)*D(i).\n*\n*  PIVMIN   (input) REAL            \n*           The minimum pivot in the Sturm sequence.\n*\n*  GAPTOL   (input) REAL            \n*           Tolerance that indicates when eigenvector entries are negligible\n*           w.r.t. their contribution to the residual.\n*\n*  Z        (input/output) REAL             array, dimension (N)\n*           On input, all entries of Z must be set to 0.\n*           On output, Z contains the (scaled) r-th column of the\n*           inverse. The scaling is such that Z(R) equals 1.\n*\n*  WANTNC   (input) LOGICAL\n*           Specifies whether NEGCNT has to be computed.\n*\n*  NEGCNT   (output) INTEGER\n*           If WANTNC is .TRUE. then NEGCNT = the number of pivots < pivmin\n*           in the  matrix factorization L D L^T, and NEGCNT = -1 otherwise.\n*\n*  ZTZ      (output) REAL            \n*           The square of the 2-norm of Z.\n*\n*  MINGMA   (output) REAL            \n*           The reciprocal of the largest (in magnitude) diagonal\n*           element of the inverse of L D L^T - sigma I.\n*\n*  R        (input/output) INTEGER\n*           The twist index for the twisted factorization used to\n*           compute Z.\n*           On input, 0 <= R <= N. If R is input as 0, R is set to\n*           the index where (L D L^T - sigma I)^{-1} is largest\n*           in magnitude. If 1 <= R <= N, R is unchanged.\n*           On output, R contains the twist index used to compute Z.\n*           Ideally, R designates the position of the maximum entry in the\n*           eigenvector.\n*\n*  ISUPPZ   (output) INTEGER array, dimension (2)\n*           The support of the vector in Z, i.e., the vector Z is\n*           nonzero only in elements ISUPPZ(1) through ISUPPZ( 2 ).\n*\n*  NRMINV   (output) REAL            \n*           NRMINV = 1/SQRT( ZTZ )\n*\n*  RESID    (output) REAL            \n*           The residual of the FP vector.\n*           RESID = ABS( MINGMA )/SQRT( ZTZ )\n*\n*  RQCORR   (output) REAL            \n*           The Rayleigh Quotient correction to LAMBDA.\n*           RQCORR = MINGMA*TMP\n*\n*  WORK     (workspace) REAL             array, dimension (4*N)\n*\n\n*  Further Details\n*  ===============\n*\n*  Based on contributions by\n*     Beresford Parlett, University of California, Berkeley, USA\n*     Jim Demmel, University of California, Berkeley, USA\n*     Inderjit Dhillon, University of Texas, Austin, USA\n*     Osni Marques, LBNL/NERSC, USA\n*     Christof Voemel, University of California, Berkeley, USA\n*\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 12)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 12)", argc);
  rb_b1 = argv[0];
  rb_bn = argv[1];
  rb_lambda = argv[2];
  rb_d = argv[3];
  rb_l = argv[4];
  rb_ld = argv[5];
  rb_lld = argv[6];
  rb_pivmin = argv[7];
  rb_gaptol = argv[8];
  rb_z = argv[9];
  rb_wantnc = argv[10];
  rb_r = argv[11];

  b1 = NUM2INT(rb_b1);
  bn = NUM2INT(rb_bn);
  lambda = (real)NUM2DBL(rb_lambda);
  pivmin = (real)NUM2DBL(rb_pivmin);
  gaptol = (real)NUM2DBL(rb_gaptol);
  wantnc = (rb_wantnc == Qtrue);
  r = NUM2INT(rb_r);
  if (!NA_IsNArray(rb_d))
    rb_raise(rb_eArgError, "d (4th argument) must be NArray");
  if (NA_RANK(rb_d) != 1)
    rb_raise(rb_eArgError, "rank of d (4th argument) must be %d", 1);
  n = NA_SHAPE0(rb_d);
  if (NA_TYPE(rb_d) != NA_SFLOAT)
    rb_d = na_change_type(rb_d, NA_SFLOAT);
  d = NA_PTR_TYPE(rb_d, real*);
  if (!NA_IsNArray(rb_l))
    rb_raise(rb_eArgError, "l (5th argument) must be NArray");
  if (NA_RANK(rb_l) != 1)
    rb_raise(rb_eArgError, "rank of l (5th argument) must be %d", 1);
  if (NA_SHAPE0(rb_l) != (n-1))
    rb_raise(rb_eRuntimeError, "shape 0 of l must be %d", n-1);
  if (NA_TYPE(rb_l) != NA_SFLOAT)
    rb_l = na_change_type(rb_l, NA_SFLOAT);
  l = NA_PTR_TYPE(rb_l, real*);
  if (!NA_IsNArray(rb_ld))
    rb_raise(rb_eArgError, "ld (6th argument) must be NArray");
  if (NA_RANK(rb_ld) != 1)
    rb_raise(rb_eArgError, "rank of ld (6th argument) must be %d", 1);
  if (NA_SHAPE0(rb_ld) != (n-1))
    rb_raise(rb_eRuntimeError, "shape 0 of ld must be %d", n-1);
  if (NA_TYPE(rb_ld) != NA_SFLOAT)
    rb_ld = na_change_type(rb_ld, NA_SFLOAT);
  ld = NA_PTR_TYPE(rb_ld, real*);
  if (!NA_IsNArray(rb_lld))
    rb_raise(rb_eArgError, "lld (7th argument) must be NArray");
  if (NA_RANK(rb_lld) != 1)
    rb_raise(rb_eArgError, "rank of lld (7th argument) must be %d", 1);
  if (NA_SHAPE0(rb_lld) != (n-1))
    rb_raise(rb_eRuntimeError, "shape 0 of lld must be %d", n-1);
  if (NA_TYPE(rb_lld) != NA_SFLOAT)
    rb_lld = na_change_type(rb_lld, NA_SFLOAT);
  lld = NA_PTR_TYPE(rb_lld, real*);
  if (!NA_IsNArray(rb_z))
    rb_raise(rb_eArgError, "z (10th argument) must be NArray");
  if (NA_RANK(rb_z) != 1)
    rb_raise(rb_eArgError, "rank of z (10th argument) must be %d", 1);
  if (NA_SHAPE0(rb_z) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of z must be the same as shape 0 of d");
  if (NA_TYPE(rb_z) != NA_SFLOAT)
    rb_z = na_change_type(rb_z, NA_SFLOAT);
  z = NA_PTR_TYPE(rb_z, real*);
  {
    int shape[1];
    shape[0] = DIM_LEN(2);
    rb_isuppz = na_make_object(NA_LINT, 1, shape, cNArray);
  }
  isuppz = NA_PTR_TYPE(rb_isuppz, integer*);
  {
    int shape[1];
    shape[0] = n;
    rb_z_out__ = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  z_out__ = NA_PTR_TYPE(rb_z_out__, real*);
  MEMCPY(z_out__, z, real, NA_TOTAL(rb_z));
  rb_z = rb_z_out__;
  z = z_out__;
  work = ALLOC_N(real, (4*n));

  slar1v_(&n, &b1, &bn, &lambda, d, l, ld, lld, &pivmin, &gaptol, z, &wantnc, &negcnt, &ztz, &mingma, &r, isuppz, &nrminv, &resid, &rqcorr, work);

  free(work);
  rb_negcnt = INT2NUM(negcnt);
  rb_ztz = rb_float_new((double)ztz);
  rb_mingma = rb_float_new((double)mingma);
  rb_nrminv = rb_float_new((double)nrminv);
  rb_resid = rb_float_new((double)resid);
  rb_rqcorr = rb_float_new((double)rqcorr);
  rb_r = INT2NUM(r);
  return rb_ary_new3(9, rb_negcnt, rb_ztz, rb_mingma, rb_isuppz, rb_nrminv, rb_resid, rb_rqcorr, rb_z, rb_r);
}

void
init_lapack_slar1v(VALUE mLapack){
  rb_define_module_function(mLapack, "slar1v", rb_slar1v, -1);
}

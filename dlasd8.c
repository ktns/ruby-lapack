#include "rb_lapack.h"

static VALUE
rb_dlasd8(int argc, VALUE *argv, VALUE self){
  VALUE rb_icompq;
  integer icompq; 
  VALUE rb_z;
  doublereal *z; 
  VALUE rb_vf;
  doublereal *vf; 
  VALUE rb_vl;
  doublereal *vl; 
  VALUE rb_lddifr;
  integer lddifr; 
  VALUE rb_dsigma;
  doublereal *dsigma; 
  VALUE rb_d;
  doublereal *d; 
  VALUE rb_difl;
  doublereal *difl; 
  VALUE rb_difr;
  doublereal *difr; 
  VALUE rb_info;
  integer info; 
  VALUE rb_z_out__;
  doublereal *z_out__;
  VALUE rb_vf_out__;
  doublereal *vf_out__;
  VALUE rb_vl_out__;
  doublereal *vl_out__;
  VALUE rb_dsigma_out__;
  doublereal *dsigma_out__;
  doublereal *work;

  integer k;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  d, difl, difr, info, z, vf, vl, dsigma = NumRu::Lapack.dlasd8( icompq, z, vf, vl, lddifr, dsigma)\n    or\n  NumRu::Lapack.dlasd8  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE DLASD8( ICOMPQ, K, D, Z, VF, VL, DIFL, DIFR, LDDIFR, DSIGMA, WORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  DLASD8 finds the square roots of the roots of the secular equation,\n*  as defined by the values in DSIGMA and Z. It makes the appropriate\n*  calls to DLASD4, and stores, for each  element in D, the distance\n*  to its two nearest poles (elements in DSIGMA). It also updates\n*  the arrays VF and VL, the first and last components of all the\n*  right singular vectors of the original bidiagonal matrix.\n*\n*  DLASD8 is called from DLASD6.\n*\n\n*  Arguments\n*  =========\n*\n*  ICOMPQ  (input) INTEGER\n*          Specifies whether singular vectors are to be computed in\n*          factored form in the calling routine:\n*          = 0: Compute singular values only.\n*          = 1: Compute singular vectors in factored form as well.\n*\n*  K       (input) INTEGER\n*          The number of terms in the rational function to be solved\n*          by DLASD4.  K >= 1.\n*\n*  D       (output) DOUBLE PRECISION array, dimension ( K )\n*          On output, D contains the updated singular values.\n*\n*  Z       (input/output) DOUBLE PRECISION array, dimension ( K )\n*          On entry, the first K elements of this array contain the\n*          components of the deflation-adjusted updating row vector.\n*          On exit, Z is updated.\n*\n*  VF      (input/output) DOUBLE PRECISION array, dimension ( K )\n*          On entry, VF contains  information passed through DBEDE8.\n*          On exit, VF contains the first K components of the first\n*          components of all right singular vectors of the bidiagonal\n*          matrix.\n*\n*  VL      (input/output) DOUBLE PRECISION array, dimension ( K )\n*          On entry, VL contains  information passed through DBEDE8.\n*          On exit, VL contains the first K components of the last\n*          components of all right singular vectors of the bidiagonal\n*          matrix.\n*\n*  DIFL    (output) DOUBLE PRECISION array, dimension ( K )\n*          On exit, DIFL(I) = D(I) - DSIGMA(I).\n*\n*  DIFR    (output) DOUBLE PRECISION array,\n*                   dimension ( LDDIFR, 2 ) if ICOMPQ = 1 and\n*                   dimension ( K ) if ICOMPQ = 0.\n*          On exit, DIFR(I,1) = D(I) - DSIGMA(I+1), DIFR(K,1) is not\n*          defined and will not be referenced.\n*\n*          If ICOMPQ = 1, DIFR(1:K,2) is an array containing the\n*          normalizing factors for the right singular vector matrix.\n*\n*  LDDIFR  (input) INTEGER\n*          The leading dimension of DIFR, must be at least K.\n*\n*  DSIGMA  (input/output) DOUBLE PRECISION array, dimension ( K )\n*          On entry, the first K elements of this array contain the old\n*          roots of the deflated updating problem.  These are the poles\n*          of the secular equation.\n*          On exit, the elements of DSIGMA may be very slightly altered\n*          in value.\n*\n*  WORK    (workspace) DOUBLE PRECISION array, dimension at least 3 * K\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit.\n*          < 0:  if INFO = -i, the i-th argument had an illegal value.\n*          > 0:  if INFO = 1, an singular value did not converge\n*\n\n*  Further Details\n*  ===============\n*\n*  Based on contributions by\n*     Ming Gu and Huan Ren, Computer Science Division, University of\n*     California at Berkeley, USA\n*\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 6)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 6)", argc);
  rb_icompq = argv[0];
  rb_z = argv[1];
  rb_vf = argv[2];
  rb_vl = argv[3];
  rb_lddifr = argv[4];
  rb_dsigma = argv[5];

  icompq = NUM2INT(rb_icompq);
  lddifr = NUM2INT(rb_lddifr);
  if (!NA_IsNArray(rb_z))
    rb_raise(rb_eArgError, "z (2th argument) must be NArray");
  if (NA_RANK(rb_z) != 1)
    rb_raise(rb_eArgError, "rank of z (2th argument) must be %d", 1);
  k = NA_SHAPE0(rb_z);
  if (NA_TYPE(rb_z) != NA_DFLOAT)
    rb_z = na_change_type(rb_z, NA_DFLOAT);
  z = NA_PTR_TYPE(rb_z, doublereal*);
  if (!NA_IsNArray(rb_vf))
    rb_raise(rb_eArgError, "vf (3th argument) must be NArray");
  if (NA_RANK(rb_vf) != 1)
    rb_raise(rb_eArgError, "rank of vf (3th argument) must be %d", 1);
  if (NA_SHAPE0(rb_vf) != k)
    rb_raise(rb_eRuntimeError, "shape 0 of vf must be the same as shape 0 of z");
  if (NA_TYPE(rb_vf) != NA_DFLOAT)
    rb_vf = na_change_type(rb_vf, NA_DFLOAT);
  vf = NA_PTR_TYPE(rb_vf, doublereal*);
  if (!NA_IsNArray(rb_vl))
    rb_raise(rb_eArgError, "vl (4th argument) must be NArray");
  if (NA_RANK(rb_vl) != 1)
    rb_raise(rb_eArgError, "rank of vl (4th argument) must be %d", 1);
  if (NA_SHAPE0(rb_vl) != k)
    rb_raise(rb_eRuntimeError, "shape 0 of vl must be the same as shape 0 of z");
  if (NA_TYPE(rb_vl) != NA_DFLOAT)
    rb_vl = na_change_type(rb_vl, NA_DFLOAT);
  vl = NA_PTR_TYPE(rb_vl, doublereal*);
  if (!NA_IsNArray(rb_dsigma))
    rb_raise(rb_eArgError, "dsigma (6th argument) must be NArray");
  if (NA_RANK(rb_dsigma) != 1)
    rb_raise(rb_eArgError, "rank of dsigma (6th argument) must be %d", 1);
  if (NA_SHAPE0(rb_dsigma) != k)
    rb_raise(rb_eRuntimeError, "shape 0 of dsigma must be the same as shape 0 of z");
  if (NA_TYPE(rb_dsigma) != NA_DFLOAT)
    rb_dsigma = na_change_type(rb_dsigma, NA_DFLOAT);
  dsigma = NA_PTR_TYPE(rb_dsigma, doublereal*);
  {
    int shape[1];
    shape[0] = k;
    rb_d = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  d = NA_PTR_TYPE(rb_d, doublereal*);
  {
    int shape[1];
    shape[0] = k;
    rb_difl = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  difl = NA_PTR_TYPE(rb_difl, doublereal*);
  lddifr = k;
  {
    int shape[2];
    shape[0] = icompq == 1 ? lddifr : icompq == 0 ? k : 0;
    shape[1] = icompq == 1 ? 2 : 0;
    rb_difr = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  difr = NA_PTR_TYPE(rb_difr, doublereal*);
  {
    int shape[1];
    shape[0] = k;
    rb_z_out__ = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  z_out__ = NA_PTR_TYPE(rb_z_out__, doublereal*);
  MEMCPY(z_out__, z, doublereal, NA_TOTAL(rb_z));
  rb_z = rb_z_out__;
  z = z_out__;
  {
    int shape[1];
    shape[0] = k;
    rb_vf_out__ = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  vf_out__ = NA_PTR_TYPE(rb_vf_out__, doublereal*);
  MEMCPY(vf_out__, vf, doublereal, NA_TOTAL(rb_vf));
  rb_vf = rb_vf_out__;
  vf = vf_out__;
  {
    int shape[1];
    shape[0] = k;
    rb_vl_out__ = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  vl_out__ = NA_PTR_TYPE(rb_vl_out__, doublereal*);
  MEMCPY(vl_out__, vl, doublereal, NA_TOTAL(rb_vl));
  rb_vl = rb_vl_out__;
  vl = vl_out__;
  {
    int shape[1];
    shape[0] = k;
    rb_dsigma_out__ = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  dsigma_out__ = NA_PTR_TYPE(rb_dsigma_out__, doublereal*);
  MEMCPY(dsigma_out__, dsigma, doublereal, NA_TOTAL(rb_dsigma));
  rb_dsigma = rb_dsigma_out__;
  dsigma = dsigma_out__;
  work = ALLOC_N(doublereal, (3 * k));

  dlasd8_(&icompq, &k, d, z, vf, vl, difl, difr, &lddifr, dsigma, work, &info);

  free(work);
  rb_info = INT2NUM(info);
  return rb_ary_new3(8, rb_d, rb_difl, rb_difr, rb_info, rb_z, rb_vf, rb_vl, rb_dsigma);
}

void
init_lapack_dlasd8(VALUE mLapack){
  rb_define_module_function(mLapack, "dlasd8", rb_dlasd8, -1);
}

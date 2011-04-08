#include "rb_lapack.h"

extern VOID slasd8_(integer *icompq, integer *k, real *d, real *z, real *vf, real *vl, real *difl, real *difr, integer *lddifr, real *dsigma, real *work, integer *info);

static VALUE sHelp, sUsage;

static VALUE
rblapack_slasd8(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_icompq;
  integer icompq; 
  VALUE rblapack_z;
  real *z; 
  VALUE rblapack_vf;
  real *vf; 
  VALUE rblapack_vl;
  real *vl; 
  VALUE rblapack_lddifr;
  integer lddifr; 
  VALUE rblapack_dsigma;
  real *dsigma; 
  VALUE rblapack_d;
  real *d; 
  VALUE rblapack_difl;
  real *difl; 
  VALUE rblapack_difr;
  real *difr; 
  VALUE rblapack_info;
  integer info; 
  VALUE rblapack_z_out__;
  real *z_out__;
  VALUE rblapack_vf_out__;
  real *vf_out__;
  VALUE rblapack_vl_out__;
  real *vl_out__;
  VALUE rblapack_dsigma_out__;
  real *dsigma_out__;
  real *work;

  integer k;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  d, difl, difr, info, z, vf, vl, dsigma = NumRu::Lapack.slasd8( icompq, z, vf, vl, lddifr, dsigma, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE SLASD8( ICOMPQ, K, D, Z, VF, VL, DIFL, DIFR, LDDIFR, DSIGMA, WORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  SLASD8 finds the square roots of the roots of the secular equation,\n*  as defined by the values in DSIGMA and Z. It makes the appropriate\n*  calls to SLASD4, and stores, for each  element in D, the distance\n*  to its two nearest poles (elements in DSIGMA). It also updates\n*  the arrays VF and VL, the first and last components of all the\n*  right singular vectors of the original bidiagonal matrix.\n*\n*  SLASD8 is called from SLASD6.\n*\n\n*  Arguments\n*  =========\n*\n*  ICOMPQ  (input) INTEGER\n*          Specifies whether singular vectors are to be computed in\n*          factored form in the calling routine:\n*          = 0: Compute singular values only.\n*          = 1: Compute singular vectors in factored form as well.\n*\n*  K       (input) INTEGER\n*          The number of terms in the rational function to be solved\n*          by SLASD4.  K >= 1.\n*\n*  D       (output) REAL array, dimension ( K )\n*          On output, D contains the updated singular values.\n*\n*  Z       (input/output) REAL array, dimension ( K )\n*          On entry, the first K elements of this array contain the\n*          components of the deflation-adjusted updating row vector.\n*          On exit, Z is updated.\n*\n*  VF      (input/output) REAL array, dimension ( K )\n*          On entry, VF contains  information passed through DBEDE8.\n*          On exit, VF contains the first K components of the first\n*          components of all right singular vectors of the bidiagonal\n*          matrix.\n*\n*  VL      (input/output) REAL array, dimension ( K )\n*          On entry, VL contains  information passed through DBEDE8.\n*          On exit, VL contains the first K components of the last\n*          components of all right singular vectors of the bidiagonal\n*          matrix.\n*\n*  DIFL    (output) REAL array, dimension ( K )\n*          On exit, DIFL(I) = D(I) - DSIGMA(I).\n*\n*  DIFR    (output) REAL array,\n*                   dimension ( LDDIFR, 2 ) if ICOMPQ = 1 and\n*                   dimension ( K ) if ICOMPQ = 0.\n*          On exit, DIFR(I,1) = D(I) - DSIGMA(I+1), DIFR(K,1) is not\n*          defined and will not be referenced.\n*\n*          If ICOMPQ = 1, DIFR(1:K,2) is an array containing the\n*          normalizing factors for the right singular vector matrix.\n*\n*  LDDIFR  (input) INTEGER\n*          The leading dimension of DIFR, must be at least K.\n*\n*  DSIGMA  (input/output) REAL array, dimension ( K )\n*          On entry, the first K elements of this array contain the old\n*          roots of the deflated updating problem.  These are the poles\n*          of the secular equation.\n*          On exit, the elements of DSIGMA may be very slightly altered\n*          in value.\n*\n*  WORK    (workspace) REAL array, dimension at least 3 * K\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit.\n*          < 0:  if INFO = -i, the i-th argument had an illegal value.\n*          > 0:  if INFO = 1, a singular value did not converge\n*\n\n*  Further Details\n*  ===============\n*\n*  Based on contributions by\n*     Ming Gu and Huan Ren, Computer Science Division, University of\n*     California at Berkeley, USA\n*\n*  =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  d, difl, difr, info, z, vf, vl, dsigma = NumRu::Lapack.slasd8( icompq, z, vf, vl, lddifr, dsigma, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 6)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 6)", argc);
  rblapack_icompq = argv[0];
  rblapack_z = argv[1];
  rblapack_vf = argv[2];
  rblapack_vl = argv[3];
  rblapack_lddifr = argv[4];
  rblapack_dsigma = argv[5];
  if (rb_options != Qnil) {
  }

  icompq = NUM2INT(rblapack_icompq);
  if (!NA_IsNArray(rblapack_vl))
    rb_raise(rb_eArgError, "vl (4th argument) must be NArray");
  if (NA_RANK(rblapack_vl) != 1)
    rb_raise(rb_eArgError, "rank of vl (4th argument) must be %d", 1);
  k = NA_SHAPE0(rblapack_vl);
  if (NA_TYPE(rblapack_vl) != NA_SFLOAT)
    rblapack_vl = na_change_type(rblapack_vl, NA_SFLOAT);
  vl = NA_PTR_TYPE(rblapack_vl, real*);
  if (!NA_IsNArray(rblapack_dsigma))
    rb_raise(rb_eArgError, "dsigma (6th argument) must be NArray");
  if (NA_RANK(rblapack_dsigma) != 1)
    rb_raise(rb_eArgError, "rank of dsigma (6th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_dsigma) != k)
    rb_raise(rb_eRuntimeError, "shape 0 of dsigma must be the same as shape 0 of vl");
  if (NA_TYPE(rblapack_dsigma) != NA_SFLOAT)
    rblapack_dsigma = na_change_type(rblapack_dsigma, NA_SFLOAT);
  dsigma = NA_PTR_TYPE(rblapack_dsigma, real*);
  if (!NA_IsNArray(rblapack_vf))
    rb_raise(rb_eArgError, "vf (3th argument) must be NArray");
  if (NA_RANK(rblapack_vf) != 1)
    rb_raise(rb_eArgError, "rank of vf (3th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_vf) != k)
    rb_raise(rb_eRuntimeError, "shape 0 of vf must be the same as shape 0 of vl");
  if (NA_TYPE(rblapack_vf) != NA_SFLOAT)
    rblapack_vf = na_change_type(rblapack_vf, NA_SFLOAT);
  vf = NA_PTR_TYPE(rblapack_vf, real*);
  if (!NA_IsNArray(rblapack_z))
    rb_raise(rb_eArgError, "z (2th argument) must be NArray");
  if (NA_RANK(rblapack_z) != 1)
    rb_raise(rb_eArgError, "rank of z (2th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_z) != k)
    rb_raise(rb_eRuntimeError, "shape 0 of z must be the same as shape 0 of vl");
  if (NA_TYPE(rblapack_z) != NA_SFLOAT)
    rblapack_z = na_change_type(rblapack_z, NA_SFLOAT);
  z = NA_PTR_TYPE(rblapack_z, real*);
  lddifr = k;
  {
    int shape[1];
    shape[0] = k;
    rblapack_d = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  d = NA_PTR_TYPE(rblapack_d, real*);
  {
    int shape[1];
    shape[0] = k;
    rblapack_difl = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  difl = NA_PTR_TYPE(rblapack_difl, real*);
  {
    int shape[2];
    shape[0] = icompq == 1 ? lddifr : icompq == 0 ? k : 0;
    shape[1] = icompq == 1 ? 2 : 0;
    rblapack_difr = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  difr = NA_PTR_TYPE(rblapack_difr, real*);
  {
    int shape[1];
    shape[0] = k;
    rblapack_z_out__ = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  z_out__ = NA_PTR_TYPE(rblapack_z_out__, real*);
  MEMCPY(z_out__, z, real, NA_TOTAL(rblapack_z));
  rblapack_z = rblapack_z_out__;
  z = z_out__;
  {
    int shape[1];
    shape[0] = k;
    rblapack_vf_out__ = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  vf_out__ = NA_PTR_TYPE(rblapack_vf_out__, real*);
  MEMCPY(vf_out__, vf, real, NA_TOTAL(rblapack_vf));
  rblapack_vf = rblapack_vf_out__;
  vf = vf_out__;
  {
    int shape[1];
    shape[0] = k;
    rblapack_vl_out__ = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  vl_out__ = NA_PTR_TYPE(rblapack_vl_out__, real*);
  MEMCPY(vl_out__, vl, real, NA_TOTAL(rblapack_vl));
  rblapack_vl = rblapack_vl_out__;
  vl = vl_out__;
  {
    int shape[1];
    shape[0] = k;
    rblapack_dsigma_out__ = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  dsigma_out__ = NA_PTR_TYPE(rblapack_dsigma_out__, real*);
  MEMCPY(dsigma_out__, dsigma, real, NA_TOTAL(rblapack_dsigma));
  rblapack_dsigma = rblapack_dsigma_out__;
  dsigma = dsigma_out__;
  work = ALLOC_N(real, (3 * k));

  slasd8_(&icompq, &k, d, z, vf, vl, difl, difr, &lddifr, dsigma, work, &info);

  free(work);
  rblapack_info = INT2NUM(info);
  return rb_ary_new3(8, rblapack_d, rblapack_difl, rblapack_difr, rblapack_info, rblapack_z, rblapack_vf, rblapack_vl, rblapack_dsigma);
}

void
init_lapack_slasd8(VALUE mLapack){
  rb_define_module_function(mLapack, "slasd8", rblapack_slasd8, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}

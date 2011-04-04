#include "rb_lapack.h"

extern VOID slar1v_(integer *n, integer *b1, integer *bn, real *lambda, real *d, real *l, real *ld, real *lld, real *pivmin, real *gaptol, real *z, logical *wantnc, integer *negcnt, real *ztz, real *mingma, integer *r, integer *isuppz, real *nrminv, real *resid, real *rqcorr, real *work);

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
    printf("%s\n", "USAGE:\n  negcnt, ztz, mingma, isuppz, nrminv, resid, rqcorr, z, r = NumRu::Lapack.slar1v( b1, bn, lambda, d, l, ld, lld, pivmin, gaptol, z, wantnc, r)\n    or\n  NumRu::Lapack.slar1v  # print help\n\n\nFORTRAN MANUAL\n\n");
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

  pivmin = (real)NUM2DBL(rb_pivmin);
  bn = NUM2INT(rb_bn);
  lambda = (real)NUM2DBL(rb_lambda);
  if (!NA_IsNArray(rb_z))
    rb_raise(rb_eArgError, "z (10th argument) must be NArray");
  if (NA_RANK(rb_z) != 1)
    rb_raise(rb_eArgError, "rank of z (10th argument) must be %d", 1);
  n = NA_SHAPE0(rb_z);
  if (NA_TYPE(rb_z) != NA_SFLOAT)
    rb_z = na_change_type(rb_z, NA_SFLOAT);
  z = NA_PTR_TYPE(rb_z, real*);
  wantnc = (rb_wantnc == Qtrue);
  if (!NA_IsNArray(rb_d))
    rb_raise(rb_eArgError, "d (4th argument) must be NArray");
  if (NA_RANK(rb_d) != 1)
    rb_raise(rb_eArgError, "rank of d (4th argument) must be %d", 1);
  if (NA_SHAPE0(rb_d) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of d must be the same as shape 0 of z");
  if (NA_TYPE(rb_d) != NA_SFLOAT)
    rb_d = na_change_type(rb_d, NA_SFLOAT);
  d = NA_PTR_TYPE(rb_d, real*);
  r = NUM2INT(rb_r);
  gaptol = (real)NUM2DBL(rb_gaptol);
  b1 = NUM2INT(rb_b1);
  if (!NA_IsNArray(rb_lld))
    rb_raise(rb_eArgError, "lld (7th argument) must be NArray");
  if (NA_RANK(rb_lld) != 1)
    rb_raise(rb_eArgError, "rank of lld (7th argument) must be %d", 1);
  if (NA_SHAPE0(rb_lld) != (n-1))
    rb_raise(rb_eRuntimeError, "shape 0 of lld must be %d", n-1);
  if (NA_TYPE(rb_lld) != NA_SFLOAT)
    rb_lld = na_change_type(rb_lld, NA_SFLOAT);
  lld = NA_PTR_TYPE(rb_lld, real*);
  if (!NA_IsNArray(rb_ld))
    rb_raise(rb_eArgError, "ld (6th argument) must be NArray");
  if (NA_RANK(rb_ld) != 1)
    rb_raise(rb_eArgError, "rank of ld (6th argument) must be %d", 1);
  if (NA_SHAPE0(rb_ld) != (n-1))
    rb_raise(rb_eRuntimeError, "shape 0 of ld must be %d", n-1);
  if (NA_TYPE(rb_ld) != NA_SFLOAT)
    rb_ld = na_change_type(rb_ld, NA_SFLOAT);
  ld = NA_PTR_TYPE(rb_ld, real*);
  if (!NA_IsNArray(rb_l))
    rb_raise(rb_eArgError, "l (5th argument) must be NArray");
  if (NA_RANK(rb_l) != 1)
    rb_raise(rb_eArgError, "rank of l (5th argument) must be %d", 1);
  if (NA_SHAPE0(rb_l) != (n-1))
    rb_raise(rb_eRuntimeError, "shape 0 of l must be %d", n-1);
  if (NA_TYPE(rb_l) != NA_SFLOAT)
    rb_l = na_change_type(rb_l, NA_SFLOAT);
  l = NA_PTR_TYPE(rb_l, real*);
  {
    int shape[1];
    shape[0] = 2;
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

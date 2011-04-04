#include "rb_lapack.h"

extern VOID slasd8_(integer *icompq, integer *k, real *d, real *z, real *vf, real *vl, real *difl, real *difr, integer *lddifr, real *dsigma, real *work, integer *info);

static VALUE
rb_slasd8(int argc, VALUE *argv, VALUE self){
  VALUE rb_icompq;
  integer icompq; 
  VALUE rb_z;
  real *z; 
  VALUE rb_vf;
  real *vf; 
  VALUE rb_vl;
  real *vl; 
  VALUE rb_lddifr;
  integer lddifr; 
  VALUE rb_dsigma;
  real *dsigma; 
  VALUE rb_d;
  real *d; 
  VALUE rb_difl;
  real *difl; 
  VALUE rb_difr;
  real *difr; 
  VALUE rb_info;
  integer info; 
  VALUE rb_z_out__;
  real *z_out__;
  VALUE rb_vf_out__;
  real *vf_out__;
  VALUE rb_vl_out__;
  real *vl_out__;
  VALUE rb_dsigma_out__;
  real *dsigma_out__;
  real *work;

  integer k;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  d, difl, difr, info, z, vf, vl, dsigma = NumRu::Lapack.slasd8( icompq, z, vf, vl, lddifr, dsigma)\n    or\n  NumRu::Lapack.slasd8  # print help\n\n\nFORTRAN MANUAL\n\n");
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
  if (!NA_IsNArray(rb_vl))
    rb_raise(rb_eArgError, "vl (4th argument) must be NArray");
  if (NA_RANK(rb_vl) != 1)
    rb_raise(rb_eArgError, "rank of vl (4th argument) must be %d", 1);
  k = NA_SHAPE0(rb_vl);
  if (NA_TYPE(rb_vl) != NA_SFLOAT)
    rb_vl = na_change_type(rb_vl, NA_SFLOAT);
  vl = NA_PTR_TYPE(rb_vl, real*);
  if (!NA_IsNArray(rb_dsigma))
    rb_raise(rb_eArgError, "dsigma (6th argument) must be NArray");
  if (NA_RANK(rb_dsigma) != 1)
    rb_raise(rb_eArgError, "rank of dsigma (6th argument) must be %d", 1);
  if (NA_SHAPE0(rb_dsigma) != k)
    rb_raise(rb_eRuntimeError, "shape 0 of dsigma must be the same as shape 0 of vl");
  if (NA_TYPE(rb_dsigma) != NA_SFLOAT)
    rb_dsigma = na_change_type(rb_dsigma, NA_SFLOAT);
  dsigma = NA_PTR_TYPE(rb_dsigma, real*);
  if (!NA_IsNArray(rb_vf))
    rb_raise(rb_eArgError, "vf (3th argument) must be NArray");
  if (NA_RANK(rb_vf) != 1)
    rb_raise(rb_eArgError, "rank of vf (3th argument) must be %d", 1);
  if (NA_SHAPE0(rb_vf) != k)
    rb_raise(rb_eRuntimeError, "shape 0 of vf must be the same as shape 0 of vl");
  if (NA_TYPE(rb_vf) != NA_SFLOAT)
    rb_vf = na_change_type(rb_vf, NA_SFLOAT);
  vf = NA_PTR_TYPE(rb_vf, real*);
  if (!NA_IsNArray(rb_z))
    rb_raise(rb_eArgError, "z (2th argument) must be NArray");
  if (NA_RANK(rb_z) != 1)
    rb_raise(rb_eArgError, "rank of z (2th argument) must be %d", 1);
  if (NA_SHAPE0(rb_z) != k)
    rb_raise(rb_eRuntimeError, "shape 0 of z must be the same as shape 0 of vl");
  if (NA_TYPE(rb_z) != NA_SFLOAT)
    rb_z = na_change_type(rb_z, NA_SFLOAT);
  z = NA_PTR_TYPE(rb_z, real*);
  lddifr = k;
  {
    int shape[1];
    shape[0] = k;
    rb_d = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  d = NA_PTR_TYPE(rb_d, real*);
  {
    int shape[1];
    shape[0] = k;
    rb_difl = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  difl = NA_PTR_TYPE(rb_difl, real*);
  {
    int shape[2];
    shape[0] = icompq == 1 ? lddifr : icompq == 0 ? k : 0;
    shape[1] = icompq == 1 ? 2 : 0;
    rb_difr = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  difr = NA_PTR_TYPE(rb_difr, real*);
  {
    int shape[1];
    shape[0] = k;
    rb_z_out__ = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  z_out__ = NA_PTR_TYPE(rb_z_out__, real*);
  MEMCPY(z_out__, z, real, NA_TOTAL(rb_z));
  rb_z = rb_z_out__;
  z = z_out__;
  {
    int shape[1];
    shape[0] = k;
    rb_vf_out__ = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  vf_out__ = NA_PTR_TYPE(rb_vf_out__, real*);
  MEMCPY(vf_out__, vf, real, NA_TOTAL(rb_vf));
  rb_vf = rb_vf_out__;
  vf = vf_out__;
  {
    int shape[1];
    shape[0] = k;
    rb_vl_out__ = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  vl_out__ = NA_PTR_TYPE(rb_vl_out__, real*);
  MEMCPY(vl_out__, vl, real, NA_TOTAL(rb_vl));
  rb_vl = rb_vl_out__;
  vl = vl_out__;
  {
    int shape[1];
    shape[0] = k;
    rb_dsigma_out__ = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  dsigma_out__ = NA_PTR_TYPE(rb_dsigma_out__, real*);
  MEMCPY(dsigma_out__, dsigma, real, NA_TOTAL(rb_dsigma));
  rb_dsigma = rb_dsigma_out__;
  dsigma = dsigma_out__;
  work = ALLOC_N(real, (3 * k));

  slasd8_(&icompq, &k, d, z, vf, vl, difl, difr, &lddifr, dsigma, work, &info);

  free(work);
  rb_info = INT2NUM(info);
  return rb_ary_new3(8, rb_d, rb_difl, rb_difr, rb_info, rb_z, rb_vf, rb_vl, rb_dsigma);
}

void
init_lapack_slasd8(VALUE mLapack){
  rb_define_module_function(mLapack, "slasd8", rb_slasd8, -1);
}

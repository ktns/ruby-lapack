#include "rb_lapack.h"

extern VOID dhseqr_(char *job, char *compz, integer *n, integer *ilo, integer *ihi, doublereal *h, integer *ldh, doublereal *wr, doublereal *wi, doublereal *z, integer *ldz, doublereal *work, integer *lwork, integer *info);

static VALUE
rb_dhseqr(int argc, VALUE *argv, VALUE self){
  VALUE rb_job;
  char job; 
  VALUE rb_compz;
  char compz; 
  VALUE rb_ilo;
  integer ilo; 
  VALUE rb_ihi;
  integer ihi; 
  VALUE rb_h;
  doublereal *h; 
  VALUE rb_z;
  doublereal *z; 
  VALUE rb_ldz;
  integer ldz; 
  VALUE rb_lwork;
  integer lwork; 
  VALUE rb_wr;
  doublereal *wr; 
  VALUE rb_wi;
  doublereal *wi; 
  VALUE rb_work;
  doublereal *work; 
  VALUE rb_info;
  integer info; 
  VALUE rb_h_out__;
  doublereal *h_out__;
  VALUE rb_z_out__;
  doublereal *z_out__;

  integer ldh;
  integer n;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  wr, wi, work, info, h, z = NumRu::Lapack.dhseqr( job, compz, ilo, ihi, h, z, ldz, lwork)\n    or\n  NumRu::Lapack.dhseqr  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 8)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 8)", argc);
  rb_job = argv[0];
  rb_compz = argv[1];
  rb_ilo = argv[2];
  rb_ihi = argv[3];
  rb_h = argv[4];
  rb_z = argv[5];
  rb_ldz = argv[6];
  rb_lwork = argv[7];

  ilo = NUM2INT(rb_ilo);
  ldz = NUM2INT(rb_ldz);
  compz = StringValueCStr(rb_compz)[0];
  ihi = NUM2INT(rb_ihi);
  if (!NA_IsNArray(rb_h))
    rb_raise(rb_eArgError, "h (5th argument) must be NArray");
  if (NA_RANK(rb_h) != 2)
    rb_raise(rb_eArgError, "rank of h (5th argument) must be %d", 2);
  n = NA_SHAPE1(rb_h);
  ldh = NA_SHAPE0(rb_h);
  if (NA_TYPE(rb_h) != NA_DFLOAT)
    rb_h = na_change_type(rb_h, NA_DFLOAT);
  h = NA_PTR_TYPE(rb_h, doublereal*);
  job = StringValueCStr(rb_job)[0];
  lwork = NUM2INT(rb_lwork);
  if (!NA_IsNArray(rb_z))
    rb_raise(rb_eArgError, "z (6th argument) must be NArray");
  if (NA_RANK(rb_z) != 2)
    rb_raise(rb_eArgError, "rank of z (6th argument) must be %d", 2);
  if (NA_SHAPE1(rb_z) != (lsame_(&compz,"N") ? 0 : n))
    rb_raise(rb_eRuntimeError, "shape 1 of z must be %d", lsame_(&compz,"N") ? 0 : n);
  if (NA_SHAPE0(rb_z) != (lsame_(&compz,"N") ? 0 : ldz))
    rb_raise(rb_eRuntimeError, "shape 0 of z must be %d", lsame_(&compz,"N") ? 0 : ldz);
  if (NA_TYPE(rb_z) != NA_DFLOAT)
    rb_z = na_change_type(rb_z, NA_DFLOAT);
  z = NA_PTR_TYPE(rb_z, doublereal*);
  {
    int shape[1];
    shape[0] = n;
    rb_wr = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  wr = NA_PTR_TYPE(rb_wr, doublereal*);
  {
    int shape[1];
    shape[0] = n;
    rb_wi = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  wi = NA_PTR_TYPE(rb_wi, doublereal*);
  {
    int shape[1];
    shape[0] = MAX(1,lwork);
    rb_work = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  work = NA_PTR_TYPE(rb_work, doublereal*);
  {
    int shape[2];
    shape[0] = ldh;
    shape[1] = n;
    rb_h_out__ = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  h_out__ = NA_PTR_TYPE(rb_h_out__, doublereal*);
  MEMCPY(h_out__, h, doublereal, NA_TOTAL(rb_h));
  rb_h = rb_h_out__;
  h = h_out__;
  {
    int shape[2];
    shape[0] = lsame_(&compz,"N") ? 0 : ldz;
    shape[1] = lsame_(&compz,"N") ? 0 : n;
    rb_z_out__ = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  z_out__ = NA_PTR_TYPE(rb_z_out__, doublereal*);
  MEMCPY(z_out__, z, doublereal, NA_TOTAL(rb_z));
  rb_z = rb_z_out__;
  z = z_out__;

  dhseqr_(&job, &compz, &n, &ilo, &ihi, h, &ldh, wr, wi, z, &ldz, work, &lwork, &info);

  rb_info = INT2NUM(info);
  return rb_ary_new3(6, rb_wr, rb_wi, rb_work, rb_info, rb_h, rb_z);
}

void
init_lapack_dhseqr(VALUE mLapack){
  rb_define_module_function(mLapack, "dhseqr", rb_dhseqr, -1);
}

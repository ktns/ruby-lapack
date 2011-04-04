#include "rb_lapack.h"

extern VOID zhseqr_(char *job, char *compz, integer *n, integer *ilo, integer *ihi, doublecomplex *h, integer *ldh, doublecomplex *w, doublecomplex *z, integer *ldz, doublecomplex *work, integer *lwork, integer *info);

static VALUE
rb_zhseqr(int argc, VALUE *argv, VALUE self){
  VALUE rb_job;
  char job; 
  VALUE rb_compz;
  char compz; 
  VALUE rb_ilo;
  integer ilo; 
  VALUE rb_ihi;
  integer ihi; 
  VALUE rb_h;
  doublecomplex *h; 
  VALUE rb_z;
  doublecomplex *z; 
  VALUE rb_ldz;
  integer ldz; 
  VALUE rb_lwork;
  integer lwork; 
  VALUE rb_w;
  doublecomplex *w; 
  VALUE rb_work;
  doublecomplex *work; 
  VALUE rb_info;
  integer info; 
  VALUE rb_h_out__;
  doublecomplex *h_out__;
  VALUE rb_z_out__;
  doublecomplex *z_out__;

  integer ldh;
  integer n;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  w, work, info, h, z = NumRu::Lapack.zhseqr( job, compz, ilo, ihi, h, z, ldz, lwork)\n    or\n  NumRu::Lapack.zhseqr  # print help\n\n\nFORTRAN MANUAL\n\n");
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
  if (NA_TYPE(rb_h) != NA_DCOMPLEX)
    rb_h = na_change_type(rb_h, NA_DCOMPLEX);
  h = NA_PTR_TYPE(rb_h, doublecomplex*);
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
  if (NA_TYPE(rb_z) != NA_DCOMPLEX)
    rb_z = na_change_type(rb_z, NA_DCOMPLEX);
  z = NA_PTR_TYPE(rb_z, doublecomplex*);
  {
    int shape[1];
    shape[0] = n;
    rb_w = na_make_object(NA_DCOMPLEX, 1, shape, cNArray);
  }
  w = NA_PTR_TYPE(rb_w, doublecomplex*);
  {
    int shape[1];
    shape[0] = MAX(1,lwork);
    rb_work = na_make_object(NA_DCOMPLEX, 1, shape, cNArray);
  }
  work = NA_PTR_TYPE(rb_work, doublecomplex*);
  {
    int shape[2];
    shape[0] = ldh;
    shape[1] = n;
    rb_h_out__ = na_make_object(NA_DCOMPLEX, 2, shape, cNArray);
  }
  h_out__ = NA_PTR_TYPE(rb_h_out__, doublecomplex*);
  MEMCPY(h_out__, h, doublecomplex, NA_TOTAL(rb_h));
  rb_h = rb_h_out__;
  h = h_out__;
  {
    int shape[2];
    shape[0] = lsame_(&compz,"N") ? 0 : ldz;
    shape[1] = lsame_(&compz,"N") ? 0 : n;
    rb_z_out__ = na_make_object(NA_DCOMPLEX, 2, shape, cNArray);
  }
  z_out__ = NA_PTR_TYPE(rb_z_out__, doublecomplex*);
  MEMCPY(z_out__, z, doublecomplex, NA_TOTAL(rb_z));
  rb_z = rb_z_out__;
  z = z_out__;

  zhseqr_(&job, &compz, &n, &ilo, &ihi, h, &ldh, w, z, &ldz, work, &lwork, &info);

  rb_info = INT2NUM(info);
  return rb_ary_new3(5, rb_w, rb_work, rb_info, rb_h, rb_z);
}

void
init_lapack_zhseqr(VALUE mLapack){
  rb_define_module_function(mLapack, "zhseqr", rb_zhseqr, -1);
}

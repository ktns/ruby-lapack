#include "rb_lapack.h"

extern VOID slasd3_(integer *nl, integer *nr, integer *sqre, integer *k, real *d, real *q, integer *ldq, real *dsigma, real *u, integer *ldu, real *u2, integer *ldu2, real *vt, integer *ldvt, real *vt2, integer *ldvt2, integer *idxc, integer *ctot, real *z, integer *info);

static VALUE
rb_slasd3(int argc, VALUE *argv, VALUE self){
  VALUE rb_nl;
  integer nl; 
  VALUE rb_nr;
  integer nr; 
  VALUE rb_sqre;
  integer sqre; 
  VALUE rb_dsigma;
  real *dsigma; 
  VALUE rb_u2;
  real *u2; 
  VALUE rb_vt2;
  real *vt2; 
  VALUE rb_idxc;
  integer *idxc; 
  VALUE rb_ctot;
  integer *ctot; 
  VALUE rb_z;
  real *z; 
  VALUE rb_d;
  real *d; 
  VALUE rb_u;
  real *u; 
  VALUE rb_vt;
  real *vt; 
  VALUE rb_info;
  integer info; 
  VALUE rb_dsigma_out__;
  real *dsigma_out__;
  VALUE rb_vt2_out__;
  real *vt2_out__;
  VALUE rb_z_out__;
  real *z_out__;
  real *q;

  integer k;
  integer ldu2;
  integer n;
  integer ldvt2;
  integer ldu;
  integer ldvt;
  integer m;
  integer ldq;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  d, u, vt, info, dsigma, vt2, z = NumRu::Lapack.slasd3( nl, nr, sqre, dsigma, u2, vt2, idxc, ctot, z)\n    or\n  NumRu::Lapack.slasd3  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 9)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 9)", argc);
  rb_nl = argv[0];
  rb_nr = argv[1];
  rb_sqre = argv[2];
  rb_dsigma = argv[3];
  rb_u2 = argv[4];
  rb_vt2 = argv[5];
  rb_idxc = argv[6];
  rb_ctot = argv[7];
  rb_z = argv[8];

  if (!NA_IsNArray(rb_ctot))
    rb_raise(rb_eArgError, "ctot (8th argument) must be NArray");
  if (NA_RANK(rb_ctot) != 1)
    rb_raise(rb_eArgError, "rank of ctot (8th argument) must be %d", 1);
  if (NA_SHAPE0(rb_ctot) != (4))
    rb_raise(rb_eRuntimeError, "shape 0 of ctot must be %d", 4);
  if (NA_TYPE(rb_ctot) != NA_LINT)
    rb_ctot = na_change_type(rb_ctot, NA_LINT);
  ctot = NA_PTR_TYPE(rb_ctot, integer*);
  if (!NA_IsNArray(rb_vt2))
    rb_raise(rb_eArgError, "vt2 (6th argument) must be NArray");
  if (NA_RANK(rb_vt2) != 2)
    rb_raise(rb_eArgError, "rank of vt2 (6th argument) must be %d", 2);
  n = NA_SHAPE1(rb_vt2);
  if (n != (nl + nr + 1))
    rb_raise(rb_eRuntimeError, "shape 1 of vt2 must be %d", nl + nr + 1);
  ldvt2 = NA_SHAPE0(rb_vt2);
  if (ldvt2 != (n))
    rb_raise(rb_eRuntimeError, "shape 0 of vt2 must be %d", n);
  n = ldvt2;
  if (NA_TYPE(rb_vt2) != NA_SFLOAT)
    rb_vt2 = na_change_type(rb_vt2, NA_SFLOAT);
  vt2 = NA_PTR_TYPE(rb_vt2, real*);
  nl = NUM2INT(rb_nl);
  if (!NA_IsNArray(rb_idxc))
    rb_raise(rb_eArgError, "idxc (7th argument) must be NArray");
  if (NA_RANK(rb_idxc) != 1)
    rb_raise(rb_eArgError, "rank of idxc (7th argument) must be %d", 1);
  if (NA_SHAPE0(rb_idxc) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of idxc must be n");
  if (NA_TYPE(rb_idxc) != NA_LINT)
    rb_idxc = na_change_type(rb_idxc, NA_LINT);
  idxc = NA_PTR_TYPE(rb_idxc, integer*);
  if (!NA_IsNArray(rb_u2))
    rb_raise(rb_eArgError, "u2 (5th argument) must be NArray");
  if (NA_RANK(rb_u2) != 2)
    rb_raise(rb_eArgError, "rank of u2 (5th argument) must be %d", 2);
  if (NA_SHAPE1(rb_u2) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of u2 must be n");
  ldu2 = NA_SHAPE0(rb_u2);
  if (ldu2 != (n))
    rb_raise(rb_eRuntimeError, "shape 0 of u2 must be %d", n);
  n = ldu2;
  if (NA_TYPE(rb_u2) != NA_SFLOAT)
    rb_u2 = na_change_type(rb_u2, NA_SFLOAT);
  u2 = NA_PTR_TYPE(rb_u2, real*);
  if (!NA_IsNArray(rb_dsigma))
    rb_raise(rb_eArgError, "dsigma (4th argument) must be NArray");
  if (NA_RANK(rb_dsigma) != 1)
    rb_raise(rb_eArgError, "rank of dsigma (4th argument) must be %d", 1);
  k = NA_SHAPE0(rb_dsigma);
  if (NA_TYPE(rb_dsigma) != NA_SFLOAT)
    rb_dsigma = na_change_type(rb_dsigma, NA_SFLOAT);
  dsigma = NA_PTR_TYPE(rb_dsigma, real*);
  sqre = NUM2INT(rb_sqre);
  if (!NA_IsNArray(rb_z))
    rb_raise(rb_eArgError, "z (9th argument) must be NArray");
  if (NA_RANK(rb_z) != 1)
    rb_raise(rb_eArgError, "rank of z (9th argument) must be %d", 1);
  if (NA_SHAPE0(rb_z) != k)
    rb_raise(rb_eRuntimeError, "shape 0 of z must be the same as shape 0 of dsigma");
  if (NA_TYPE(rb_z) != NA_SFLOAT)
    rb_z = na_change_type(rb_z, NA_SFLOAT);
  z = NA_PTR_TYPE(rb_z, real*);
  nr = NUM2INT(rb_nr);
  ldu2 = n;
  ldu = n;
  ldvt = n;
  ldvt2 = n;
  ldq = k;
  m = n+sqre;
  {
    int shape[1];
    shape[0] = k;
    rb_d = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  d = NA_PTR_TYPE(rb_d, real*);
  {
    int shape[2];
    shape[0] = ldu;
    shape[1] = n;
    rb_u = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  u = NA_PTR_TYPE(rb_u, real*);
  {
    int shape[2];
    shape[0] = ldvt;
    shape[1] = m;
    rb_vt = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  vt = NA_PTR_TYPE(rb_vt, real*);
  {
    int shape[1];
    shape[0] = k;
    rb_dsigma_out__ = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  dsigma_out__ = NA_PTR_TYPE(rb_dsigma_out__, real*);
  MEMCPY(dsigma_out__, dsigma, real, NA_TOTAL(rb_dsigma));
  rb_dsigma = rb_dsigma_out__;
  dsigma = dsigma_out__;
  {
    int shape[2];
    shape[0] = ldvt2;
    shape[1] = n;
    rb_vt2_out__ = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  vt2_out__ = NA_PTR_TYPE(rb_vt2_out__, real*);
  MEMCPY(vt2_out__, vt2, real, NA_TOTAL(rb_vt2));
  rb_vt2 = rb_vt2_out__;
  vt2 = vt2_out__;
  {
    int shape[1];
    shape[0] = k;
    rb_z_out__ = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  z_out__ = NA_PTR_TYPE(rb_z_out__, real*);
  MEMCPY(z_out__, z, real, NA_TOTAL(rb_z));
  rb_z = rb_z_out__;
  z = z_out__;
  q = ALLOC_N(real, (ldq)*(k));

  slasd3_(&nl, &nr, &sqre, &k, d, q, &ldq, dsigma, u, &ldu, u2, &ldu2, vt, &ldvt, vt2, &ldvt2, idxc, ctot, z, &info);

  free(q);
  rb_info = INT2NUM(info);
  return rb_ary_new3(7, rb_d, rb_u, rb_vt, rb_info, rb_dsigma, rb_vt2, rb_z);
}

void
init_lapack_slasd3(VALUE mLapack){
  rb_define_module_function(mLapack, "slasd3", rb_slasd3, -1);
}

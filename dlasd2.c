#include "rb_lapack.h"

extern VOID dlasd2_(integer *nl, integer *nr, integer *sqre, integer *k, doublereal *d, doublereal *z, doublereal *alpha, doublereal *beta, doublereal *u, integer *ldu, doublereal *vt, integer *ldvt, doublereal *dsigma, doublereal *u2, integer *ldu2, doublereal *vt2, integer *ldvt2, integer *idxp, integer *idx, integer *idxc, integer *idxq, integer *coltyp, integer *info);

static VALUE
rb_dlasd2(int argc, VALUE *argv, VALUE self){
  VALUE rb_nl;
  integer nl; 
  VALUE rb_nr;
  integer nr; 
  VALUE rb_sqre;
  integer sqre; 
  VALUE rb_d;
  doublereal *d; 
  VALUE rb_alpha;
  doublereal alpha; 
  VALUE rb_beta;
  doublereal beta; 
  VALUE rb_u;
  doublereal *u; 
  VALUE rb_vt;
  doublereal *vt; 
  VALUE rb_idxq;
  integer *idxq; 
  VALUE rb_k;
  integer k; 
  VALUE rb_z;
  doublereal *z; 
  VALUE rb_dsigma;
  doublereal *dsigma; 
  VALUE rb_u2;
  doublereal *u2; 
  VALUE rb_vt2;
  doublereal *vt2; 
  VALUE rb_idxc;
  integer *idxc; 
  VALUE rb_coltyp;
  integer *coltyp; 
  VALUE rb_info;
  integer info; 
  VALUE rb_d_out__;
  doublereal *d_out__;
  VALUE rb_u_out__;
  doublereal *u_out__;
  VALUE rb_vt_out__;
  doublereal *vt_out__;
  VALUE rb_idxq_out__;
  integer *idxq_out__;
  integer *idxp;
  integer *idx;

  integer n;
  integer ldu;
  integer ldvt;
  integer m;
  integer ldu2;
  integer ldvt2;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  k, z, dsigma, u2, vt2, idxc, coltyp, info, d, u, vt, idxq = NumRu::Lapack.dlasd2( nl, nr, sqre, d, alpha, beta, u, vt, idxq)\n    or\n  NumRu::Lapack.dlasd2  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 9)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 9)", argc);
  rb_nl = argv[0];
  rb_nr = argv[1];
  rb_sqre = argv[2];
  rb_d = argv[3];
  rb_alpha = argv[4];
  rb_beta = argv[5];
  rb_u = argv[6];
  rb_vt = argv[7];
  rb_idxq = argv[8];

  if (!NA_IsNArray(rb_idxq))
    rb_raise(rb_eArgError, "idxq (9th argument) must be NArray");
  if (NA_RANK(rb_idxq) != 1)
    rb_raise(rb_eArgError, "rank of idxq (9th argument) must be %d", 1);
  n = NA_SHAPE0(rb_idxq);
  if (NA_TYPE(rb_idxq) != NA_LINT)
    rb_idxq = na_change_type(rb_idxq, NA_LINT);
  idxq = NA_PTR_TYPE(rb_idxq, integer*);
  if (!NA_IsNArray(rb_u))
    rb_raise(rb_eArgError, "u (7th argument) must be NArray");
  if (NA_RANK(rb_u) != 2)
    rb_raise(rb_eArgError, "rank of u (7th argument) must be %d", 2);
  if (NA_SHAPE1(rb_u) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of u must be the same as shape 0 of idxq");
  ldu = NA_SHAPE0(rb_u);
  if (NA_TYPE(rb_u) != NA_DFLOAT)
    rb_u = na_change_type(rb_u, NA_DFLOAT);
  u = NA_PTR_TYPE(rb_u, doublereal*);
  if (!NA_IsNArray(rb_d))
    rb_raise(rb_eArgError, "d (4th argument) must be NArray");
  if (NA_RANK(rb_d) != 1)
    rb_raise(rb_eArgError, "rank of d (4th argument) must be %d", 1);
  if (NA_SHAPE0(rb_d) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of d must be the same as shape 0 of idxq");
  if (NA_TYPE(rb_d) != NA_DFLOAT)
    rb_d = na_change_type(rb_d, NA_DFLOAT);
  d = NA_PTR_TYPE(rb_d, doublereal*);
  nr = NUM2INT(rb_nr);
  beta = NUM2DBL(rb_beta);
  nl = NUM2INT(rb_nl);
  sqre = NUM2INT(rb_sqre);
  if (!NA_IsNArray(rb_vt))
    rb_raise(rb_eArgError, "vt (8th argument) must be NArray");
  if (NA_RANK(rb_vt) != 2)
    rb_raise(rb_eArgError, "rank of vt (8th argument) must be %d", 2);
  m = NA_SHAPE1(rb_vt);
  ldvt = NA_SHAPE0(rb_vt);
  if (NA_TYPE(rb_vt) != NA_DFLOAT)
    rb_vt = na_change_type(rb_vt, NA_DFLOAT);
  vt = NA_PTR_TYPE(rb_vt, doublereal*);
  alpha = NUM2DBL(rb_alpha);
  ldu2 = n;
  ldvt2 = m;
  {
    int shape[1];
    shape[0] = n;
    rb_z = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  z = NA_PTR_TYPE(rb_z, doublereal*);
  {
    int shape[1];
    shape[0] = n;
    rb_dsigma = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  dsigma = NA_PTR_TYPE(rb_dsigma, doublereal*);
  {
    int shape[2];
    shape[0] = ldu2;
    shape[1] = n;
    rb_u2 = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  u2 = NA_PTR_TYPE(rb_u2, doublereal*);
  {
    int shape[2];
    shape[0] = ldvt2;
    shape[1] = n;
    rb_vt2 = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  vt2 = NA_PTR_TYPE(rb_vt2, doublereal*);
  {
    int shape[1];
    shape[0] = n;
    rb_idxc = na_make_object(NA_LINT, 1, shape, cNArray);
  }
  idxc = NA_PTR_TYPE(rb_idxc, integer*);
  {
    int shape[1];
    shape[0] = n;
    rb_coltyp = na_make_object(NA_LINT, 1, shape, cNArray);
  }
  coltyp = NA_PTR_TYPE(rb_coltyp, integer*);
  {
    int shape[1];
    shape[0] = n;
    rb_d_out__ = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  d_out__ = NA_PTR_TYPE(rb_d_out__, doublereal*);
  MEMCPY(d_out__, d, doublereal, NA_TOTAL(rb_d));
  rb_d = rb_d_out__;
  d = d_out__;
  {
    int shape[2];
    shape[0] = ldu;
    shape[1] = n;
    rb_u_out__ = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  u_out__ = NA_PTR_TYPE(rb_u_out__, doublereal*);
  MEMCPY(u_out__, u, doublereal, NA_TOTAL(rb_u));
  rb_u = rb_u_out__;
  u = u_out__;
  {
    int shape[2];
    shape[0] = ldvt;
    shape[1] = m;
    rb_vt_out__ = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  vt_out__ = NA_PTR_TYPE(rb_vt_out__, doublereal*);
  MEMCPY(vt_out__, vt, doublereal, NA_TOTAL(rb_vt));
  rb_vt = rb_vt_out__;
  vt = vt_out__;
  {
    int shape[1];
    shape[0] = n;
    rb_idxq_out__ = na_make_object(NA_LINT, 1, shape, cNArray);
  }
  idxq_out__ = NA_PTR_TYPE(rb_idxq_out__, integer*);
  MEMCPY(idxq_out__, idxq, integer, NA_TOTAL(rb_idxq));
  rb_idxq = rb_idxq_out__;
  idxq = idxq_out__;
  idxp = ALLOC_N(integer, (n));
  idx = ALLOC_N(integer, (n));

  dlasd2_(&nl, &nr, &sqre, &k, d, z, &alpha, &beta, u, &ldu, vt, &ldvt, dsigma, u2, &ldu2, vt2, &ldvt2, idxp, idx, idxc, idxq, coltyp, &info);

  free(idxp);
  free(idx);
  rb_k = INT2NUM(k);
  rb_info = INT2NUM(info);
  return rb_ary_new3(12, rb_k, rb_z, rb_dsigma, rb_u2, rb_vt2, rb_idxc, rb_coltyp, rb_info, rb_d, rb_u, rb_vt, rb_idxq);
}

void
init_lapack_dlasd2(VALUE mLapack){
  rb_define_module_function(mLapack, "dlasd2", rb_dlasd2, -1);
}

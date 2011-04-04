#include "rb_lapack.h"

extern VOID slasd6_(integer *icompq, integer *nl, integer *nr, integer *sqre, real *d, real *vf, real *vl, real *alpha, real *beta, integer *idxq, integer *perm, integer *givptr, integer *givcol, integer *ldgcol, real *givnum, integer *ldgnum, real *poles, real *difl, real *difr, real *z, integer *k, real *c, real *s, real *work, integer *iwork, integer *info);

static VALUE
rb_slasd6(int argc, VALUE *argv, VALUE self){
  VALUE rb_icompq;
  integer icompq; 
  VALUE rb_nl;
  integer nl; 
  VALUE rb_nr;
  integer nr; 
  VALUE rb_sqre;
  integer sqre; 
  VALUE rb_d;
  real *d; 
  VALUE rb_vf;
  real *vf; 
  VALUE rb_vl;
  real *vl; 
  VALUE rb_alpha;
  real alpha; 
  VALUE rb_beta;
  real beta; 
  VALUE rb_idxq;
  integer *idxq; 
  VALUE rb_perm;
  integer *perm; 
  VALUE rb_givptr;
  integer givptr; 
  VALUE rb_givcol;
  integer *givcol; 
  VALUE rb_givnum;
  real *givnum; 
  VALUE rb_poles;
  real *poles; 
  VALUE rb_difl;
  real *difl; 
  VALUE rb_difr;
  real *difr; 
  VALUE rb_z;
  real *z; 
  VALUE rb_k;
  integer k; 
  VALUE rb_c;
  real c; 
  VALUE rb_s;
  real s; 
  VALUE rb_info;
  integer info; 
  VALUE rb_d_out__;
  real *d_out__;
  VALUE rb_vf_out__;
  real *vf_out__;
  VALUE rb_vl_out__;
  real *vl_out__;
  real *work;
  integer *iwork;

  integer m;
  integer n;
  integer ldgcol;
  integer ldgnum;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  idxq, perm, givptr, givcol, givnum, poles, difl, difr, z, k, c, s, info, d, vf, vl, alpha, beta = NumRu::Lapack.slasd6( icompq, nl, nr, sqre, d, vf, vl, alpha, beta)\n    or\n  NumRu::Lapack.slasd6  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 9)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 9)", argc);
  rb_icompq = argv[0];
  rb_nl = argv[1];
  rb_nr = argv[2];
  rb_sqre = argv[3];
  rb_d = argv[4];
  rb_vf = argv[5];
  rb_vl = argv[6];
  rb_alpha = argv[7];
  rb_beta = argv[8];

  if (!NA_IsNArray(rb_vl))
    rb_raise(rb_eArgError, "vl (7th argument) must be NArray");
  if (NA_RANK(rb_vl) != 1)
    rb_raise(rb_eArgError, "rank of vl (7th argument) must be %d", 1);
  m = NA_SHAPE0(rb_vl);
  if (m != (n + sqre))
    rb_raise(rb_eRuntimeError, "shape 0 of vl must be %d", n + sqre);
  if (NA_TYPE(rb_vl) != NA_SFLOAT)
    rb_vl = na_change_type(rb_vl, NA_SFLOAT);
  vl = NA_PTR_TYPE(rb_vl, real*);
  alpha = (real)NUM2DBL(rb_alpha);
  nl = NUM2INT(rb_nl);
  sqre = NUM2INT(rb_sqre);
  nr = NUM2INT(rb_nr);
  icompq = NUM2INT(rb_icompq);
  if (!NA_IsNArray(rb_vf))
    rb_raise(rb_eArgError, "vf (6th argument) must be NArray");
  if (NA_RANK(rb_vf) != 1)
    rb_raise(rb_eArgError, "rank of vf (6th argument) must be %d", 1);
  if (NA_SHAPE0(rb_vf) != m)
    rb_raise(rb_eRuntimeError, "shape 0 of vf must be the same as shape 0 of vl");
  if (NA_TYPE(rb_vf) != NA_SFLOAT)
    rb_vf = na_change_type(rb_vf, NA_SFLOAT);
  vf = NA_PTR_TYPE(rb_vf, real*);
  beta = (real)NUM2DBL(rb_beta);
  n = nl + nr + 1;
  ldgnum = n;
  if (!NA_IsNArray(rb_d))
    rb_raise(rb_eArgError, "d (5th argument) must be NArray");
  if (NA_RANK(rb_d) != 1)
    rb_raise(rb_eArgError, "rank of d (5th argument) must be %d", 1);
  if (NA_SHAPE0(rb_d) != (nl+nr+1))
    rb_raise(rb_eRuntimeError, "shape 0 of d must be %d", nl+nr+1);
  if (NA_TYPE(rb_d) != NA_SFLOAT)
    rb_d = na_change_type(rb_d, NA_SFLOAT);
  d = NA_PTR_TYPE(rb_d, real*);
  ldgcol = n;
  m = n + sqre;
  {
    int shape[1];
    shape[0] = n;
    rb_idxq = na_make_object(NA_LINT, 1, shape, cNArray);
  }
  idxq = NA_PTR_TYPE(rb_idxq, integer*);
  {
    int shape[1];
    shape[0] = n;
    rb_perm = na_make_object(NA_LINT, 1, shape, cNArray);
  }
  perm = NA_PTR_TYPE(rb_perm, integer*);
  {
    int shape[2];
    shape[0] = ldgcol;
    shape[1] = 2;
    rb_givcol = na_make_object(NA_LINT, 2, shape, cNArray);
  }
  givcol = NA_PTR_TYPE(rb_givcol, integer*);
  {
    int shape[2];
    shape[0] = ldgnum;
    shape[1] = 2;
    rb_givnum = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  givnum = NA_PTR_TYPE(rb_givnum, real*);
  {
    int shape[2];
    shape[0] = ldgnum;
    shape[1] = 2;
    rb_poles = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  poles = NA_PTR_TYPE(rb_poles, real*);
  {
    int shape[1];
    shape[0] = n;
    rb_difl = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  difl = NA_PTR_TYPE(rb_difl, real*);
  {
    int shape[2];
    shape[0] = icompq == 1 ? ldgnum : icompq == 0 ? n : 0;
    shape[1] = icompq == 1 ? 2 : 0;
    rb_difr = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  difr = NA_PTR_TYPE(rb_difr, real*);
  {
    int shape[1];
    shape[0] = m;
    rb_z = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  z = NA_PTR_TYPE(rb_z, real*);
  {
    int shape[1];
    shape[0] = nl+nr+1;
    rb_d_out__ = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  d_out__ = NA_PTR_TYPE(rb_d_out__, real*);
  MEMCPY(d_out__, d, real, NA_TOTAL(rb_d));
  rb_d = rb_d_out__;
  d = d_out__;
  {
    int shape[1];
    shape[0] = m;
    rb_vf_out__ = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  vf_out__ = NA_PTR_TYPE(rb_vf_out__, real*);
  MEMCPY(vf_out__, vf, real, NA_TOTAL(rb_vf));
  rb_vf = rb_vf_out__;
  vf = vf_out__;
  {
    int shape[1];
    shape[0] = m;
    rb_vl_out__ = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  vl_out__ = NA_PTR_TYPE(rb_vl_out__, real*);
  MEMCPY(vl_out__, vl, real, NA_TOTAL(rb_vl));
  rb_vl = rb_vl_out__;
  vl = vl_out__;
  work = ALLOC_N(real, (4 * m));
  iwork = ALLOC_N(integer, (3 * n));

  slasd6_(&icompq, &nl, &nr, &sqre, d, vf, vl, &alpha, &beta, idxq, perm, &givptr, givcol, &ldgcol, givnum, &ldgnum, poles, difl, difr, z, &k, &c, &s, work, iwork, &info);

  free(work);
  free(iwork);
  rb_givptr = INT2NUM(givptr);
  rb_k = INT2NUM(k);
  rb_c = rb_float_new((double)c);
  rb_s = rb_float_new((double)s);
  rb_info = INT2NUM(info);
  rb_alpha = rb_float_new((double)alpha);
  rb_beta = rb_float_new((double)beta);
  return rb_ary_new3(18, rb_idxq, rb_perm, rb_givptr, rb_givcol, rb_givnum, rb_poles, rb_difl, rb_difr, rb_z, rb_k, rb_c, rb_s, rb_info, rb_d, rb_vf, rb_vl, rb_alpha, rb_beta);
}

void
init_lapack_slasd6(VALUE mLapack){
  rb_define_module_function(mLapack, "slasd6", rb_slasd6, -1);
}

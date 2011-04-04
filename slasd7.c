#include "rb_lapack.h"

extern VOID slasd7_(integer *icompq, integer *nl, integer *nr, integer *sqre, integer *k, real *d, real *z, real *zw, real *vf, real *vfw, real *vl, real *vlw, real *alpha, real *beta, real *dsigma, integer *idx, integer *idxp, integer *idxq, integer *perm, integer *givptr, integer *givcol, integer *ldgcol, real *givnum, integer *ldgnum, real *c, real *s, integer *info);

static VALUE
rb_slasd7(int argc, VALUE *argv, VALUE self){
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
  VALUE rb_k;
  integer k; 
  VALUE rb_z;
  real *z; 
  VALUE rb_dsigma;
  real *dsigma; 
  VALUE rb_perm;
  integer *perm; 
  VALUE rb_givptr;
  integer givptr; 
  VALUE rb_givcol;
  integer *givcol; 
  VALUE rb_givnum;
  real *givnum; 
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
  real *zw;
  real *vfw;
  real *vlw;
  integer *idx;
  integer *idxp;

  integer n;
  integer m;
  integer ldgcol;
  integer ldgnum;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  k, z, dsigma, perm, givptr, givcol, givnum, c, s, info, d, vf, vl = NumRu::Lapack.slasd7( icompq, nl, nr, sqre, d, vf, vl, alpha, beta, idxq)\n    or\n  NumRu::Lapack.slasd7  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 10)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 10)", argc);
  rb_icompq = argv[0];
  rb_nl = argv[1];
  rb_nr = argv[2];
  rb_sqre = argv[3];
  rb_d = argv[4];
  rb_vf = argv[5];
  rb_vl = argv[6];
  rb_alpha = argv[7];
  rb_beta = argv[8];
  rb_idxq = argv[9];

  if (!NA_IsNArray(rb_idxq))
    rb_raise(rb_eArgError, "idxq (10th argument) must be NArray");
  if (NA_RANK(rb_idxq) != 1)
    rb_raise(rb_eArgError, "rank of idxq (10th argument) must be %d", 1);
  n = NA_SHAPE0(rb_idxq);
  if (NA_TYPE(rb_idxq) != NA_LINT)
    rb_idxq = na_change_type(rb_idxq, NA_LINT);
  idxq = NA_PTR_TYPE(rb_idxq, integer*);
  if (!NA_IsNArray(rb_vl))
    rb_raise(rb_eArgError, "vl (7th argument) must be NArray");
  if (NA_RANK(rb_vl) != 1)
    rb_raise(rb_eArgError, "rank of vl (7th argument) must be %d", 1);
  m = NA_SHAPE0(rb_vl);
  if (NA_TYPE(rb_vl) != NA_SFLOAT)
    rb_vl = na_change_type(rb_vl, NA_SFLOAT);
  vl = NA_PTR_TYPE(rb_vl, real*);
  if (!NA_IsNArray(rb_d))
    rb_raise(rb_eArgError, "d (5th argument) must be NArray");
  if (NA_RANK(rb_d) != 1)
    rb_raise(rb_eArgError, "rank of d (5th argument) must be %d", 1);
  if (NA_SHAPE0(rb_d) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of d must be the same as shape 0 of idxq");
  if (NA_TYPE(rb_d) != NA_SFLOAT)
    rb_d = na_change_type(rb_d, NA_SFLOAT);
  d = NA_PTR_TYPE(rb_d, real*);
  nr = NUM2INT(rb_nr);
  alpha = (real)NUM2DBL(rb_alpha);
  beta = (real)NUM2DBL(rb_beta);
  nl = NUM2INT(rb_nl);
  icompq = NUM2INT(rb_icompq);
  sqre = NUM2INT(rb_sqre);
  if (!NA_IsNArray(rb_vf))
    rb_raise(rb_eArgError, "vf (6th argument) must be NArray");
  if (NA_RANK(rb_vf) != 1)
    rb_raise(rb_eArgError, "rank of vf (6th argument) must be %d", 1);
  if (NA_SHAPE0(rb_vf) != m)
    rb_raise(rb_eRuntimeError, "shape 0 of vf must be the same as shape 0 of vl");
  if (NA_TYPE(rb_vf) != NA_SFLOAT)
    rb_vf = na_change_type(rb_vf, NA_SFLOAT);
  vf = NA_PTR_TYPE(rb_vf, real*);
  ldgcol = n;
  ldgnum = n;
  {
    int shape[1];
    shape[0] = m;
    rb_z = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  z = NA_PTR_TYPE(rb_z, real*);
  {
    int shape[1];
    shape[0] = n;
    rb_dsigma = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  dsigma = NA_PTR_TYPE(rb_dsigma, real*);
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
    int shape[1];
    shape[0] = n;
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
  zw = ALLOC_N(real, (m));
  vfw = ALLOC_N(real, (m));
  vlw = ALLOC_N(real, (m));
  idx = ALLOC_N(integer, (n));
  idxp = ALLOC_N(integer, (n));

  slasd7_(&icompq, &nl, &nr, &sqre, &k, d, z, zw, vf, vfw, vl, vlw, &alpha, &beta, dsigma, idx, idxp, idxq, perm, &givptr, givcol, &ldgcol, givnum, &ldgnum, &c, &s, &info);

  free(zw);
  free(vfw);
  free(vlw);
  free(idx);
  free(idxp);
  rb_k = INT2NUM(k);
  rb_givptr = INT2NUM(givptr);
  rb_c = rb_float_new((double)c);
  rb_s = rb_float_new((double)s);
  rb_info = INT2NUM(info);
  return rb_ary_new3(13, rb_k, rb_z, rb_dsigma, rb_perm, rb_givptr, rb_givcol, rb_givnum, rb_c, rb_s, rb_info, rb_d, rb_vf, rb_vl);
}

void
init_lapack_slasd7(VALUE mLapack){
  rb_define_module_function(mLapack, "slasd7", rb_slasd7, -1);
}

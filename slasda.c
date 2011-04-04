#include "rb_lapack.h"

extern VOID slasda_(integer *icompq, integer *smlsiz, integer *n, integer *sqre, real *d, real *e, real *u, integer *ldu, real *vt, integer *k, real *difl, real *difr, real *z, real *poles, integer *givptr, integer *givcol, integer *ldgcol, integer *perm, real *givnum, real *c, real *s, real *work, integer *iwork, integer *info);

static VALUE
rb_slasda(int argc, VALUE *argv, VALUE self){
  VALUE rb_icompq;
  integer icompq; 
  VALUE rb_smlsiz;
  integer smlsiz; 
  VALUE rb_sqre;
  integer sqre; 
  VALUE rb_d;
  real *d; 
  VALUE rb_e;
  real *e; 
  VALUE rb_u;
  real *u; 
  VALUE rb_vt;
  real *vt; 
  VALUE rb_k;
  integer *k; 
  VALUE rb_difl;
  real *difl; 
  VALUE rb_difr;
  real *difr; 
  VALUE rb_z;
  real *z; 
  VALUE rb_poles;
  real *poles; 
  VALUE rb_givptr;
  integer *givptr; 
  VALUE rb_givcol;
  integer *givcol; 
  VALUE rb_perm;
  integer *perm; 
  VALUE rb_givnum;
  real *givnum; 
  VALUE rb_c;
  real *c; 
  VALUE rb_s;
  real *s; 
  VALUE rb_info;
  integer info; 
  VALUE rb_d_out__;
  real *d_out__;
  real *work;
  integer *iwork;

  integer n;
  integer ldu;
  integer nlvl;
  integer ldgcol;
  integer m;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  u, vt, k, difl, difr, z, poles, givptr, givcol, perm, givnum, c, s, info, d = NumRu::Lapack.slasda( icompq, smlsiz, sqre, d, e)\n    or\n  NumRu::Lapack.slasda  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 5)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 5)", argc);
  rb_icompq = argv[0];
  rb_smlsiz = argv[1];
  rb_sqre = argv[2];
  rb_d = argv[3];
  rb_e = argv[4];

  sqre = NUM2INT(rb_sqre);
  if (!NA_IsNArray(rb_d))
    rb_raise(rb_eArgError, "d (4th argument) must be NArray");
  if (NA_RANK(rb_d) != 1)
    rb_raise(rb_eArgError, "rank of d (4th argument) must be %d", 1);
  n = NA_SHAPE0(rb_d);
  if (NA_TYPE(rb_d) != NA_SFLOAT)
    rb_d = na_change_type(rb_d, NA_SFLOAT);
  d = NA_PTR_TYPE(rb_d, real*);
  smlsiz = NUM2INT(rb_smlsiz);
  icompq = NUM2INT(rb_icompq);
  m = sqre == 0 ? n : sqre == 1 ? n+1 : 0;
  ldu = n;
  nlvl = floor(1.0/log(2.0)*log((double)n/smlsiz));
  if (!NA_IsNArray(rb_e))
    rb_raise(rb_eArgError, "e (5th argument) must be NArray");
  if (NA_RANK(rb_e) != 1)
    rb_raise(rb_eArgError, "rank of e (5th argument) must be %d", 1);
  if (NA_SHAPE0(rb_e) != (m-1))
    rb_raise(rb_eRuntimeError, "shape 0 of e must be %d", m-1);
  if (NA_TYPE(rb_e) != NA_SFLOAT)
    rb_e = na_change_type(rb_e, NA_SFLOAT);
  e = NA_PTR_TYPE(rb_e, real*);
  ldgcol = n;
  {
    int shape[2];
    shape[0] = ldu;
    shape[1] = MAX(1,smlsiz);
    rb_u = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  u = NA_PTR_TYPE(rb_u, real*);
  {
    int shape[2];
    shape[0] = ldu;
    shape[1] = smlsiz+1;
    rb_vt = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  vt = NA_PTR_TYPE(rb_vt, real*);
  {
    int shape[1];
    shape[0] = icompq == 1 ? n : icompq == 0 ? 1 : 0;
    rb_k = na_make_object(NA_LINT, 1, shape, cNArray);
  }
  k = NA_PTR_TYPE(rb_k, integer*);
  {
    int shape[2];
    shape[0] = ldu;
    shape[1] = nlvl;
    rb_difl = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  difl = NA_PTR_TYPE(rb_difl, real*);
  {
    int shape[2];
    shape[0] = icompq == 1 ? ldu : icompq == 0 ? n : 0;
    shape[1] = icompq == 1 ? 2 * nlvl : 0;
    rb_difr = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  difr = NA_PTR_TYPE(rb_difr, real*);
  {
    int shape[2];
    shape[0] = icompq == 1 ? ldu : icompq == 0 ? n : 0;
    shape[1] = icompq == 1 ? nlvl : 0;
    rb_z = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  z = NA_PTR_TYPE(rb_z, real*);
  {
    int shape[2];
    shape[0] = ldu;
    shape[1] = 2 * nlvl;
    rb_poles = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  poles = NA_PTR_TYPE(rb_poles, real*);
  {
    int shape[1];
    shape[0] = n;
    rb_givptr = na_make_object(NA_LINT, 1, shape, cNArray);
  }
  givptr = NA_PTR_TYPE(rb_givptr, integer*);
  {
    int shape[2];
    shape[0] = ldgcol;
    shape[1] = 2 * nlvl;
    rb_givcol = na_make_object(NA_LINT, 2, shape, cNArray);
  }
  givcol = NA_PTR_TYPE(rb_givcol, integer*);
  {
    int shape[2];
    shape[0] = ldgcol;
    shape[1] = nlvl;
    rb_perm = na_make_object(NA_LINT, 2, shape, cNArray);
  }
  perm = NA_PTR_TYPE(rb_perm, integer*);
  {
    int shape[2];
    shape[0] = ldu;
    shape[1] = 2 * nlvl;
    rb_givnum = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  givnum = NA_PTR_TYPE(rb_givnum, real*);
  {
    int shape[1];
    shape[0] = icompq == 1 ? n : icompq == 0 ? 1 : 0;
    rb_c = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  c = NA_PTR_TYPE(rb_c, real*);
  {
    int shape[1];
    shape[0] = icompq==1 ? n : icompq==0 ? 1 : 0;
    rb_s = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  s = NA_PTR_TYPE(rb_s, real*);
  {
    int shape[1];
    shape[0] = n;
    rb_d_out__ = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  d_out__ = NA_PTR_TYPE(rb_d_out__, real*);
  MEMCPY(d_out__, d, real, NA_TOTAL(rb_d));
  rb_d = rb_d_out__;
  d = d_out__;
  work = ALLOC_N(real, (6 * n + (smlsiz + 1)*(smlsiz + 1)));
  iwork = ALLOC_N(integer, (7*n));

  slasda_(&icompq, &smlsiz, &n, &sqre, d, e, u, &ldu, vt, k, difl, difr, z, poles, givptr, givcol, &ldgcol, perm, givnum, c, s, work, iwork, &info);

  free(work);
  free(iwork);
  rb_info = INT2NUM(info);
  return rb_ary_new3(15, rb_u, rb_vt, rb_k, rb_difl, rb_difr, rb_z, rb_poles, rb_givptr, rb_givcol, rb_perm, rb_givnum, rb_c, rb_s, rb_info, rb_d);
}

void
init_lapack_slasda(VALUE mLapack){
  rb_define_module_function(mLapack, "slasda", rb_slasda, -1);
}

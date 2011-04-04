#include "rb_lapack.h"

extern VOID dlasda_(integer *icompq, integer *smlsiz, integer *n, integer *sqre, doublereal *d, doublereal *e, doublereal *u, integer *ldu, doublereal *vt, integer *k, doublereal *difl, doublereal *difr, doublereal *z, doublereal *poles, integer *givptr, integer *givcol, integer *ldgcol, integer *perm, doublereal *givnum, doublereal *c, doublereal *s, doublereal *work, integer *iwork, integer *info);

static VALUE
rb_dlasda(int argc, VALUE *argv, VALUE self){
  VALUE rb_icompq;
  integer icompq; 
  VALUE rb_smlsiz;
  integer smlsiz; 
  VALUE rb_sqre;
  integer sqre; 
  VALUE rb_d;
  doublereal *d; 
  VALUE rb_e;
  doublereal *e; 
  VALUE rb_u;
  doublereal *u; 
  VALUE rb_vt;
  doublereal *vt; 
  VALUE rb_k;
  integer *k; 
  VALUE rb_difl;
  doublereal *difl; 
  VALUE rb_difr;
  doublereal *difr; 
  VALUE rb_z;
  doublereal *z; 
  VALUE rb_poles;
  doublereal *poles; 
  VALUE rb_givptr;
  integer *givptr; 
  VALUE rb_givcol;
  integer *givcol; 
  VALUE rb_perm;
  integer *perm; 
  VALUE rb_givnum;
  doublereal *givnum; 
  VALUE rb_c;
  doublereal *c; 
  VALUE rb_s;
  doublereal *s; 
  VALUE rb_info;
  integer info; 
  VALUE rb_d_out__;
  doublereal *d_out__;
  doublereal *work;
  integer *iwork;

  integer n;
  integer ldu;
  integer nlvl;
  integer ldgcol;
  integer m;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  u, vt, k, difl, difr, z, poles, givptr, givcol, perm, givnum, c, s, info, d = NumRu::Lapack.dlasda( icompq, smlsiz, sqre, d, e)\n    or\n  NumRu::Lapack.dlasda  # print help\n\n\nFORTRAN MANUAL\n\n");
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
  if (NA_TYPE(rb_d) != NA_DFLOAT)
    rb_d = na_change_type(rb_d, NA_DFLOAT);
  d = NA_PTR_TYPE(rb_d, doublereal*);
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
  if (NA_TYPE(rb_e) != NA_DFLOAT)
    rb_e = na_change_type(rb_e, NA_DFLOAT);
  e = NA_PTR_TYPE(rb_e, doublereal*);
  ldgcol = n;
  {
    int shape[2];
    shape[0] = ldu;
    shape[1] = MAX(1,smlsiz);
    rb_u = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  u = NA_PTR_TYPE(rb_u, doublereal*);
  {
    int shape[2];
    shape[0] = ldu;
    shape[1] = smlsiz+1;
    rb_vt = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  vt = NA_PTR_TYPE(rb_vt, doublereal*);
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
    rb_difl = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  difl = NA_PTR_TYPE(rb_difl, doublereal*);
  {
    int shape[2];
    shape[0] = icompq == 1 ? ldu : icompq == 0 ? n : 0;
    shape[1] = icompq == 1 ? 2 * nlvl : 0;
    rb_difr = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  difr = NA_PTR_TYPE(rb_difr, doublereal*);
  {
    int shape[2];
    shape[0] = icompq == 1 ? ldu : icompq == 0 ? n : 0;
    shape[1] = icompq == 1 ? nlvl : 0;
    rb_z = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  z = NA_PTR_TYPE(rb_z, doublereal*);
  {
    int shape[2];
    shape[0] = ldu;
    shape[1] = 2 * nlvl;
    rb_poles = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  poles = NA_PTR_TYPE(rb_poles, doublereal*);
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
    rb_givnum = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  givnum = NA_PTR_TYPE(rb_givnum, doublereal*);
  {
    int shape[1];
    shape[0] = icompq == 1 ? n : icompq == 0 ? 1 : 0;
    rb_c = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  c = NA_PTR_TYPE(rb_c, doublereal*);
  {
    int shape[1];
    shape[0] = icompq==1 ? n : icompq==0 ? 1 : 0;
    rb_s = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  s = NA_PTR_TYPE(rb_s, doublereal*);
  {
    int shape[1];
    shape[0] = n;
    rb_d_out__ = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  d_out__ = NA_PTR_TYPE(rb_d_out__, doublereal*);
  MEMCPY(d_out__, d, doublereal, NA_TOTAL(rb_d));
  rb_d = rb_d_out__;
  d = d_out__;
  work = ALLOC_N(doublereal, (6 * n + (smlsiz + 1)*(smlsiz + 1)));
  iwork = ALLOC_N(integer, ((7 * n)));

  dlasda_(&icompq, &smlsiz, &n, &sqre, d, e, u, &ldu, vt, k, difl, difr, z, poles, givptr, givcol, &ldgcol, perm, givnum, c, s, work, iwork, &info);

  free(work);
  free(iwork);
  rb_info = INT2NUM(info);
  return rb_ary_new3(15, rb_u, rb_vt, rb_k, rb_difl, rb_difr, rb_z, rb_poles, rb_givptr, rb_givcol, rb_perm, rb_givnum, rb_c, rb_s, rb_info, rb_d);
}

void
init_lapack_dlasda(VALUE mLapack){
  rb_define_module_function(mLapack, "dlasda", rb_dlasda, -1);
}

#include "rb_lapack.h"

extern VOID zlalsa_(integer *icompq, integer *smlsiz, integer *n, integer *nrhs, doublecomplex *b, integer *ldb, doublecomplex *bx, integer *ldbx, doublereal *u, integer *ldu, doublereal *vt, integer *k, doublereal *difl, doublereal *difr, doublereal *z, doublereal *poles, integer *givptr, integer *givcol, integer *ldgcol, integer *perm, doublereal *givnum, doublereal *c, doublereal *s, doublereal *rwork, integer *iwork, integer *info);

static VALUE
rb_zlalsa(int argc, VALUE *argv, VALUE self){
  VALUE rb_icompq;
  integer icompq; 
  VALUE rb_b;
  doublecomplex *b; 
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
  VALUE rb_bx;
  doublecomplex *bx; 
  VALUE rb_info;
  integer info; 
  VALUE rb_b_out__;
  doublecomplex *b_out__;
  doublereal *rwork;
  integer *iwork;

  integer ldb;
  integer nrhs;
  integer ldu;
  integer smlsiz;
  integer n;
  integer nlvl;
  integer ldgcol;
  integer ldbx;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  bx, info, b = NumRu::Lapack.zlalsa( icompq, b, u, vt, k, difl, difr, z, poles, givptr, givcol, perm, givnum, c, s)\n    or\n  NumRu::Lapack.zlalsa  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 15)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 15)", argc);
  rb_icompq = argv[0];
  rb_b = argv[1];
  rb_u = argv[2];
  rb_vt = argv[3];
  rb_k = argv[4];
  rb_difl = argv[5];
  rb_difr = argv[6];
  rb_z = argv[7];
  rb_poles = argv[8];
  rb_givptr = argv[9];
  rb_givcol = argv[10];
  rb_perm = argv[11];
  rb_givnum = argv[12];
  rb_c = argv[13];
  rb_s = argv[14];

  if (!NA_IsNArray(rb_k))
    rb_raise(rb_eArgError, "k (5th argument) must be NArray");
  if (NA_RANK(rb_k) != 1)
    rb_raise(rb_eArgError, "rank of k (5th argument) must be %d", 1);
  n = NA_SHAPE0(rb_k);
  if (NA_TYPE(rb_k) != NA_LINT)
    rb_k = na_change_type(rb_k, NA_LINT);
  k = NA_PTR_TYPE(rb_k, integer*);
  if (!NA_IsNArray(rb_difl))
    rb_raise(rb_eArgError, "difl (6th argument) must be NArray");
  if (NA_RANK(rb_difl) != 2)
    rb_raise(rb_eArgError, "rank of difl (6th argument) must be %d", 2);
  nlvl = NA_SHAPE1(rb_difl);
  if (nlvl != ((int)(1.0/log(2.0)*log((double)n/(smlsiz+1))) + 1))
    rb_raise(rb_eRuntimeError, "shape 1 of difl must be %d", (int)(1.0/log(2.0)*log((double)n/(smlsiz+1))) + 1);
  ldu = NA_SHAPE0(rb_difl);
  if (NA_TYPE(rb_difl) != NA_DFLOAT)
    rb_difl = na_change_type(rb_difl, NA_DFLOAT);
  difl = NA_PTR_TYPE(rb_difl, doublereal*);
  if (!NA_IsNArray(rb_b))
    rb_raise(rb_eArgError, "b (2th argument) must be NArray");
  if (NA_RANK(rb_b) != 2)
    rb_raise(rb_eArgError, "rank of b (2th argument) must be %d", 2);
  nrhs = NA_SHAPE1(rb_b);
  ldb = NA_SHAPE0(rb_b);
  if (NA_TYPE(rb_b) != NA_DCOMPLEX)
    rb_b = na_change_type(rb_b, NA_DCOMPLEX);
  b = NA_PTR_TYPE(rb_b, doublecomplex*);
  if (!NA_IsNArray(rb_c))
    rb_raise(rb_eArgError, "c (14th argument) must be NArray");
  if (NA_RANK(rb_c) != 1)
    rb_raise(rb_eArgError, "rank of c (14th argument) must be %d", 1);
  if (NA_SHAPE0(rb_c) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of c must be the same as shape 0 of k");
  if (NA_TYPE(rb_c) != NA_DFLOAT)
    rb_c = na_change_type(rb_c, NA_DFLOAT);
  c = NA_PTR_TYPE(rb_c, doublereal*);
  if (!NA_IsNArray(rb_u))
    rb_raise(rb_eArgError, "u (3th argument) must be NArray");
  if (NA_RANK(rb_u) != 2)
    rb_raise(rb_eArgError, "rank of u (3th argument) must be %d", 2);
  smlsiz = NA_SHAPE1(rb_u);
  if (NA_SHAPE0(rb_u) != ldu)
    rb_raise(rb_eRuntimeError, "shape 0 of u must be the same as shape 0 of difl");
  if (NA_TYPE(rb_u) != NA_DFLOAT)
    rb_u = na_change_type(rb_u, NA_DFLOAT);
  u = NA_PTR_TYPE(rb_u, doublereal*);
  if (!NA_IsNArray(rb_z))
    rb_raise(rb_eArgError, "z (8th argument) must be NArray");
  if (NA_RANK(rb_z) != 2)
    rb_raise(rb_eArgError, "rank of z (8th argument) must be %d", 2);
  if (NA_SHAPE1(rb_z) != nlvl)
    rb_raise(rb_eRuntimeError, "shape 1 of z must be the same as shape 1 of difl");
  if (NA_SHAPE0(rb_z) != ldu)
    rb_raise(rb_eRuntimeError, "shape 0 of z must be the same as shape 0 of difl");
  if (NA_TYPE(rb_z) != NA_DFLOAT)
    rb_z = na_change_type(rb_z, NA_DFLOAT);
  z = NA_PTR_TYPE(rb_z, doublereal*);
  if (!NA_IsNArray(rb_s))
    rb_raise(rb_eArgError, "s (15th argument) must be NArray");
  if (NA_RANK(rb_s) != 1)
    rb_raise(rb_eArgError, "rank of s (15th argument) must be %d", 1);
  if (NA_SHAPE0(rb_s) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of s must be the same as shape 0 of k");
  if (NA_TYPE(rb_s) != NA_DFLOAT)
    rb_s = na_change_type(rb_s, NA_DFLOAT);
  s = NA_PTR_TYPE(rb_s, doublereal*);
  icompq = NUM2INT(rb_icompq);
  if (!NA_IsNArray(rb_perm))
    rb_raise(rb_eArgError, "perm (12th argument) must be NArray");
  if (NA_RANK(rb_perm) != 2)
    rb_raise(rb_eArgError, "rank of perm (12th argument) must be %d", 2);
  if (NA_SHAPE1(rb_perm) != nlvl)
    rb_raise(rb_eRuntimeError, "shape 1 of perm must be the same as shape 1 of difl");
  ldgcol = NA_SHAPE0(rb_perm);
  if (NA_TYPE(rb_perm) != NA_LINT)
    rb_perm = na_change_type(rb_perm, NA_LINT);
  perm = NA_PTR_TYPE(rb_perm, integer*);
  if (!NA_IsNArray(rb_givptr))
    rb_raise(rb_eArgError, "givptr (10th argument) must be NArray");
  if (NA_RANK(rb_givptr) != 1)
    rb_raise(rb_eArgError, "rank of givptr (10th argument) must be %d", 1);
  if (NA_SHAPE0(rb_givptr) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of givptr must be the same as shape 0 of k");
  if (NA_TYPE(rb_givptr) != NA_LINT)
    rb_givptr = na_change_type(rb_givptr, NA_LINT);
  givptr = NA_PTR_TYPE(rb_givptr, integer*);
  ldbx = n;
  if (!NA_IsNArray(rb_poles))
    rb_raise(rb_eArgError, "poles (9th argument) must be NArray");
  if (NA_RANK(rb_poles) != 2)
    rb_raise(rb_eArgError, "rank of poles (9th argument) must be %d", 2);
  if (NA_SHAPE1(rb_poles) != (2 * nlvl))
    rb_raise(rb_eRuntimeError, "shape 1 of poles must be %d", 2 * nlvl);
  if (NA_SHAPE0(rb_poles) != ldu)
    rb_raise(rb_eRuntimeError, "shape 0 of poles must be the same as shape 0 of difl");
  if (NA_TYPE(rb_poles) != NA_DFLOAT)
    rb_poles = na_change_type(rb_poles, NA_DFLOAT);
  poles = NA_PTR_TYPE(rb_poles, doublereal*);
  if (!NA_IsNArray(rb_difr))
    rb_raise(rb_eArgError, "difr (7th argument) must be NArray");
  if (NA_RANK(rb_difr) != 2)
    rb_raise(rb_eArgError, "rank of difr (7th argument) must be %d", 2);
  if (NA_SHAPE1(rb_difr) != (2 * nlvl))
    rb_raise(rb_eRuntimeError, "shape 1 of difr must be %d", 2 * nlvl);
  if (NA_SHAPE0(rb_difr) != ldu)
    rb_raise(rb_eRuntimeError, "shape 0 of difr must be the same as shape 0 of difl");
  if (NA_TYPE(rb_difr) != NA_DFLOAT)
    rb_difr = na_change_type(rb_difr, NA_DFLOAT);
  difr = NA_PTR_TYPE(rb_difr, doublereal*);
  if (!NA_IsNArray(rb_vt))
    rb_raise(rb_eArgError, "vt (4th argument) must be NArray");
  if (NA_RANK(rb_vt) != 2)
    rb_raise(rb_eArgError, "rank of vt (4th argument) must be %d", 2);
  if (NA_SHAPE1(rb_vt) != (smlsiz+1))
    rb_raise(rb_eRuntimeError, "shape 1 of vt must be %d", smlsiz+1);
  if (NA_SHAPE0(rb_vt) != ldu)
    rb_raise(rb_eRuntimeError, "shape 0 of vt must be the same as shape 0 of difl");
  if (NA_TYPE(rb_vt) != NA_DFLOAT)
    rb_vt = na_change_type(rb_vt, NA_DFLOAT);
  vt = NA_PTR_TYPE(rb_vt, doublereal*);
  nlvl = (int)(1.0/log(2.0)*log((double)n/(smlsiz+1))) + 1;
  if (!NA_IsNArray(rb_givnum))
    rb_raise(rb_eArgError, "givnum (13th argument) must be NArray");
  if (NA_RANK(rb_givnum) != 2)
    rb_raise(rb_eArgError, "rank of givnum (13th argument) must be %d", 2);
  if (NA_SHAPE1(rb_givnum) != (2 * nlvl))
    rb_raise(rb_eRuntimeError, "shape 1 of givnum must be %d", 2 * nlvl);
  if (NA_SHAPE0(rb_givnum) != ldu)
    rb_raise(rb_eRuntimeError, "shape 0 of givnum must be the same as shape 0 of difl");
  if (NA_TYPE(rb_givnum) != NA_DFLOAT)
    rb_givnum = na_change_type(rb_givnum, NA_DFLOAT);
  givnum = NA_PTR_TYPE(rb_givnum, doublereal*);
  if (!NA_IsNArray(rb_givcol))
    rb_raise(rb_eArgError, "givcol (11th argument) must be NArray");
  if (NA_RANK(rb_givcol) != 2)
    rb_raise(rb_eArgError, "rank of givcol (11th argument) must be %d", 2);
  if (NA_SHAPE1(rb_givcol) != (2 * nlvl))
    rb_raise(rb_eRuntimeError, "shape 1 of givcol must be %d", 2 * nlvl);
  if (NA_SHAPE0(rb_givcol) != ldgcol)
    rb_raise(rb_eRuntimeError, "shape 0 of givcol must be the same as shape 0 of perm");
  if (NA_TYPE(rb_givcol) != NA_LINT)
    rb_givcol = na_change_type(rb_givcol, NA_LINT);
  givcol = NA_PTR_TYPE(rb_givcol, integer*);
  {
    int shape[2];
    shape[0] = ldbx;
    shape[1] = nrhs;
    rb_bx = na_make_object(NA_DCOMPLEX, 2, shape, cNArray);
  }
  bx = NA_PTR_TYPE(rb_bx, doublecomplex*);
  {
    int shape[2];
    shape[0] = ldb;
    shape[1] = nrhs;
    rb_b_out__ = na_make_object(NA_DCOMPLEX, 2, shape, cNArray);
  }
  b_out__ = NA_PTR_TYPE(rb_b_out__, doublecomplex*);
  MEMCPY(b_out__, b, doublecomplex, NA_TOTAL(rb_b));
  rb_b = rb_b_out__;
  b = b_out__;
  rwork = ALLOC_N(doublereal, (MAX(n,(smlsiz+1)*nrhs*3)));
  iwork = ALLOC_N(integer, (3 * n));

  zlalsa_(&icompq, &smlsiz, &n, &nrhs, b, &ldb, bx, &ldbx, u, &ldu, vt, k, difl, difr, z, poles, givptr, givcol, &ldgcol, perm, givnum, c, s, rwork, iwork, &info);

  free(rwork);
  free(iwork);
  rb_info = INT2NUM(info);
  return rb_ary_new3(3, rb_bx, rb_info, rb_b);
}

void
init_lapack_zlalsa(VALUE mLapack){
  rb_define_module_function(mLapack, "zlalsa", rb_zlalsa, -1);
}

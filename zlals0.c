#include "rb_lapack.h"

extern VOID zlals0_(integer *icompq, integer *nl, integer *nr, integer *sqre, integer *nrhs, doublecomplex *b, integer *ldb, doublecomplex *bx, integer *ldbx, integer *perm, integer *givptr, integer *givcol, integer *ldgcol, doublereal *givnum, integer *ldgnum, doublereal *poles, doublereal *difl, doublereal *difr, doublereal *z, integer *k, doublereal *c, doublereal *s, doublereal *rwork, integer *info);

static VALUE
rb_zlals0(int argc, VALUE *argv, VALUE self){
  VALUE rb_icompq;
  integer icompq; 
  VALUE rb_nl;
  integer nl; 
  VALUE rb_nr;
  integer nr; 
  VALUE rb_sqre;
  integer sqre; 
  VALUE rb_b;
  doublecomplex *b; 
  VALUE rb_perm;
  integer *perm; 
  VALUE rb_givptr;
  integer givptr; 
  VALUE rb_givcol;
  integer *givcol; 
  VALUE rb_givnum;
  doublereal *givnum; 
  VALUE rb_poles;
  doublereal *poles; 
  VALUE rb_difl;
  doublereal *difl; 
  VALUE rb_difr;
  doublereal *difr; 
  VALUE rb_z;
  doublereal *z; 
  VALUE rb_c;
  doublereal c; 
  VALUE rb_s;
  doublereal s; 
  VALUE rb_info;
  integer info; 
  VALUE rb_b_out__;
  doublecomplex *b_out__;
  doublecomplex *bx;
  doublereal *rwork;

  integer ldb;
  integer nrhs;
  integer n;
  integer ldgcol;
  integer ldgnum;
  integer k;
  integer ldbx;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  info, b = NumRu::Lapack.zlals0( icompq, nl, nr, sqre, b, perm, givptr, givcol, givnum, poles, difl, difr, z, c, s)\n    or\n  NumRu::Lapack.zlals0  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 15)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 15)", argc);
  rb_icompq = argv[0];
  rb_nl = argv[1];
  rb_nr = argv[2];
  rb_sqre = argv[3];
  rb_b = argv[4];
  rb_perm = argv[5];
  rb_givptr = argv[6];
  rb_givcol = argv[7];
  rb_givnum = argv[8];
  rb_poles = argv[9];
  rb_difl = argv[10];
  rb_difr = argv[11];
  rb_z = argv[12];
  rb_c = argv[13];
  rb_s = argv[14];

  if (!NA_IsNArray(rb_difl))
    rb_raise(rb_eArgError, "difl (11th argument) must be NArray");
  if (NA_RANK(rb_difl) != 1)
    rb_raise(rb_eArgError, "rank of difl (11th argument) must be %d", 1);
  k = NA_SHAPE0(rb_difl);
  if (NA_TYPE(rb_difl) != NA_DFLOAT)
    rb_difl = na_change_type(rb_difl, NA_DFLOAT);
  difl = NA_PTR_TYPE(rb_difl, doublereal*);
  if (!NA_IsNArray(rb_b))
    rb_raise(rb_eArgError, "b (5th argument) must be NArray");
  if (NA_RANK(rb_b) != 2)
    rb_raise(rb_eArgError, "rank of b (5th argument) must be %d", 2);
  nrhs = NA_SHAPE1(rb_b);
  ldb = NA_SHAPE0(rb_b);
  if (NA_TYPE(rb_b) != NA_DCOMPLEX)
    rb_b = na_change_type(rb_b, NA_DCOMPLEX);
  b = NA_PTR_TYPE(rb_b, doublecomplex*);
  c = NUM2DBL(rb_c);
  if (!NA_IsNArray(rb_givcol))
    rb_raise(rb_eArgError, "givcol (8th argument) must be NArray");
  if (NA_RANK(rb_givcol) != 2)
    rb_raise(rb_eArgError, "rank of givcol (8th argument) must be %d", 2);
  if (NA_SHAPE1(rb_givcol) != (2))
    rb_raise(rb_eRuntimeError, "shape 1 of givcol must be %d", 2);
  ldgcol = NA_SHAPE0(rb_givcol);
  if (NA_TYPE(rb_givcol) != NA_LINT)
    rb_givcol = na_change_type(rb_givcol, NA_LINT);
  givcol = NA_PTR_TYPE(rb_givcol, integer*);
  if (!NA_IsNArray(rb_z))
    rb_raise(rb_eArgError, "z (13th argument) must be NArray");
  if (NA_RANK(rb_z) != 1)
    rb_raise(rb_eArgError, "rank of z (13th argument) must be %d", 1);
  if (NA_SHAPE0(rb_z) != k)
    rb_raise(rb_eRuntimeError, "shape 0 of z must be the same as shape 0 of difl");
  if (NA_TYPE(rb_z) != NA_DFLOAT)
    rb_z = na_change_type(rb_z, NA_DFLOAT);
  z = NA_PTR_TYPE(rb_z, doublereal*);
  nr = NUM2INT(rb_nr);
  if (!NA_IsNArray(rb_poles))
    rb_raise(rb_eArgError, "poles (10th argument) must be NArray");
  if (NA_RANK(rb_poles) != 2)
    rb_raise(rb_eArgError, "rank of poles (10th argument) must be %d", 2);
  if (NA_SHAPE1(rb_poles) != (2))
    rb_raise(rb_eRuntimeError, "shape 1 of poles must be %d", 2);
  ldgnum = NA_SHAPE0(rb_poles);
  if (NA_TYPE(rb_poles) != NA_DFLOAT)
    rb_poles = na_change_type(rb_poles, NA_DFLOAT);
  poles = NA_PTR_TYPE(rb_poles, doublereal*);
  icompq = NUM2INT(rb_icompq);
  nl = NUM2INT(rb_nl);
  sqre = NUM2INT(rb_sqre);
  if (!NA_IsNArray(rb_givnum))
    rb_raise(rb_eArgError, "givnum (9th argument) must be NArray");
  if (NA_RANK(rb_givnum) != 2)
    rb_raise(rb_eArgError, "rank of givnum (9th argument) must be %d", 2);
  if (NA_SHAPE1(rb_givnum) != (2))
    rb_raise(rb_eRuntimeError, "shape 1 of givnum must be %d", 2);
  if (NA_SHAPE0(rb_givnum) != ldgnum)
    rb_raise(rb_eRuntimeError, "shape 0 of givnum must be the same as shape 0 of poles");
  if (NA_TYPE(rb_givnum) != NA_DFLOAT)
    rb_givnum = na_change_type(rb_givnum, NA_DFLOAT);
  givnum = NA_PTR_TYPE(rb_givnum, doublereal*);
  if (!NA_IsNArray(rb_perm))
    rb_raise(rb_eArgError, "perm (6th argument) must be NArray");
  if (NA_RANK(rb_perm) != 1)
    rb_raise(rb_eArgError, "rank of perm (6th argument) must be %d", 1);
  n = NA_SHAPE0(rb_perm);
  if (NA_TYPE(rb_perm) != NA_LINT)
    rb_perm = na_change_type(rb_perm, NA_LINT);
  perm = NA_PTR_TYPE(rb_perm, integer*);
  s = NUM2DBL(rb_s);
  if (!NA_IsNArray(rb_difr))
    rb_raise(rb_eArgError, "difr (12th argument) must be NArray");
  if (NA_RANK(rb_difr) != 2)
    rb_raise(rb_eArgError, "rank of difr (12th argument) must be %d", 2);
  if (NA_SHAPE1(rb_difr) != (2))
    rb_raise(rb_eRuntimeError, "shape 1 of difr must be %d", 2);
  if (NA_SHAPE0(rb_difr) != ldgnum)
    rb_raise(rb_eRuntimeError, "shape 0 of difr must be the same as shape 0 of poles");
  if (NA_TYPE(rb_difr) != NA_DFLOAT)
    rb_difr = na_change_type(rb_difr, NA_DFLOAT);
  difr = NA_PTR_TYPE(rb_difr, doublereal*);
  givptr = NUM2INT(rb_givptr);
  ldbx = n;
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
  bx = ALLOC_N(doublecomplex, (ldbx)*(nrhs));
  rwork = ALLOC_N(doublereal, (k*(1+nrhs) + 2*nrhs));

  zlals0_(&icompq, &nl, &nr, &sqre, &nrhs, b, &ldb, bx, &ldbx, perm, &givptr, givcol, &ldgcol, givnum, &ldgnum, poles, difl, difr, z, &k, &c, &s, rwork, &info);

  free(bx);
  free(rwork);
  rb_info = INT2NUM(info);
  return rb_ary_new3(2, rb_info, rb_b);
}

void
init_lapack_zlals0(VALUE mLapack){
  rb_define_module_function(mLapack, "zlals0", rb_zlals0, -1);
}

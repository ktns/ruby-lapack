#include "rb_lapack.h"

extern VOID slals0_(integer *icompq, integer *nl, integer *nr, integer *sqre, integer *nrhs, real *b, integer *ldb, real *bx, integer *ldbx, integer *perm, integer *givptr, integer *givcol, integer *ldgcol, real *givnum, integer *ldgnum, real *poles, real *difl, real *difr, real *z, integer *k, real *c, real *s, real *work, integer *info);

static VALUE
rb_slals0(int argc, VALUE *argv, VALUE self){
  VALUE rb_icompq;
  integer icompq; 
  VALUE rb_nl;
  integer nl; 
  VALUE rb_nr;
  integer nr; 
  VALUE rb_sqre;
  integer sqre; 
  VALUE rb_b;
  real *b; 
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
  VALUE rb_c;
  real c; 
  VALUE rb_s;
  real s; 
  VALUE rb_info;
  integer info; 
  VALUE rb_b_out__;
  real *b_out__;
  real *bx;
  real *work;

  integer ldb;
  integer nrhs;
  integer n;
  integer ldgcol;
  integer ldgnum;
  integer k;
  integer ldbx;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  info, b = NumRu::Lapack.slals0( icompq, nl, nr, sqre, b, perm, givptr, givcol, givnum, poles, difl, difr, z, c, s)\n    or\n  NumRu::Lapack.slals0  # print help\n\n\nFORTRAN MANUAL\n\n");
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
  if (NA_TYPE(rb_difl) != NA_SFLOAT)
    rb_difl = na_change_type(rb_difl, NA_SFLOAT);
  difl = NA_PTR_TYPE(rb_difl, real*);
  if (!NA_IsNArray(rb_b))
    rb_raise(rb_eArgError, "b (5th argument) must be NArray");
  if (NA_RANK(rb_b) != 2)
    rb_raise(rb_eArgError, "rank of b (5th argument) must be %d", 2);
  nrhs = NA_SHAPE1(rb_b);
  ldb = NA_SHAPE0(rb_b);
  if (NA_TYPE(rb_b) != NA_SFLOAT)
    rb_b = na_change_type(rb_b, NA_SFLOAT);
  b = NA_PTR_TYPE(rb_b, real*);
  c = (real)NUM2DBL(rb_c);
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
  if (NA_TYPE(rb_z) != NA_SFLOAT)
    rb_z = na_change_type(rb_z, NA_SFLOAT);
  z = NA_PTR_TYPE(rb_z, real*);
  nr = NUM2INT(rb_nr);
  if (!NA_IsNArray(rb_poles))
    rb_raise(rb_eArgError, "poles (10th argument) must be NArray");
  if (NA_RANK(rb_poles) != 2)
    rb_raise(rb_eArgError, "rank of poles (10th argument) must be %d", 2);
  if (NA_SHAPE1(rb_poles) != (2))
    rb_raise(rb_eRuntimeError, "shape 1 of poles must be %d", 2);
  ldgnum = NA_SHAPE0(rb_poles);
  if (NA_TYPE(rb_poles) != NA_SFLOAT)
    rb_poles = na_change_type(rb_poles, NA_SFLOAT);
  poles = NA_PTR_TYPE(rb_poles, real*);
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
  if (NA_TYPE(rb_givnum) != NA_SFLOAT)
    rb_givnum = na_change_type(rb_givnum, NA_SFLOAT);
  givnum = NA_PTR_TYPE(rb_givnum, real*);
  if (!NA_IsNArray(rb_perm))
    rb_raise(rb_eArgError, "perm (6th argument) must be NArray");
  if (NA_RANK(rb_perm) != 1)
    rb_raise(rb_eArgError, "rank of perm (6th argument) must be %d", 1);
  n = NA_SHAPE0(rb_perm);
  if (NA_TYPE(rb_perm) != NA_LINT)
    rb_perm = na_change_type(rb_perm, NA_LINT);
  perm = NA_PTR_TYPE(rb_perm, integer*);
  s = (real)NUM2DBL(rb_s);
  if (!NA_IsNArray(rb_difr))
    rb_raise(rb_eArgError, "difr (12th argument) must be NArray");
  if (NA_RANK(rb_difr) != 2)
    rb_raise(rb_eArgError, "rank of difr (12th argument) must be %d", 2);
  if (NA_SHAPE1(rb_difr) != (2))
    rb_raise(rb_eRuntimeError, "shape 1 of difr must be %d", 2);
  if (NA_SHAPE0(rb_difr) != ldgnum)
    rb_raise(rb_eRuntimeError, "shape 0 of difr must be the same as shape 0 of poles");
  if (NA_TYPE(rb_difr) != NA_SFLOAT)
    rb_difr = na_change_type(rb_difr, NA_SFLOAT);
  difr = NA_PTR_TYPE(rb_difr, real*);
  givptr = NUM2INT(rb_givptr);
  ldbx = n;
  {
    int shape[2];
    shape[0] = ldb;
    shape[1] = nrhs;
    rb_b_out__ = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  b_out__ = NA_PTR_TYPE(rb_b_out__, real*);
  MEMCPY(b_out__, b, real, NA_TOTAL(rb_b));
  rb_b = rb_b_out__;
  b = b_out__;
  bx = ALLOC_N(real, (ldbx)*(nrhs));
  work = ALLOC_N(real, (k));

  slals0_(&icompq, &nl, &nr, &sqre, &nrhs, b, &ldb, bx, &ldbx, perm, &givptr, givcol, &ldgcol, givnum, &ldgnum, poles, difl, difr, z, &k, &c, &s, work, &info);

  free(bx);
  free(work);
  rb_info = INT2NUM(info);
  return rb_ary_new3(2, rb_info, rb_b);
}

void
init_lapack_slals0(VALUE mLapack){
  rb_define_module_function(mLapack, "slals0", rb_slals0, -1);
}

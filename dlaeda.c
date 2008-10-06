#include "rb_lapack.h"

static VALUE
rb_dlaeda(int argc, VALUE *argv, VALUE self){
  VALUE rb_tlvls;
  integer tlvls; 
  VALUE rb_curlvl;
  integer curlvl; 
  VALUE rb_curpbm;
  integer curpbm; 
  VALUE rb_prmptr;
  integer *prmptr; 
  VALUE rb_perm;
  integer *perm; 
  VALUE rb_givptr;
  integer *givptr; 
  VALUE rb_givcol;
  integer *givcol; 
  VALUE rb_givnum;
  doublereal *givnum; 
  VALUE rb_q;
  doublereal *q; 
  VALUE rb_qptr;
  integer *qptr; 
  VALUE rb_z;
  doublereal *z; 
  VALUE rb_info;
  integer info; 
  doublereal *ztemp;

  integer n;
  integer ldqptr;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  z, info = NumRu::Lapack.dlaeda( tlvls, curlvl, curpbm, prmptr, perm, givptr, givcol, givnum, q, qptr)\n    or\n  NumRu::Lapack.dlaeda  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE DLAEDA( N, TLVLS, CURLVL, CURPBM, PRMPTR, PERM, GIVPTR, GIVCOL, GIVNUM, Q, QPTR, Z, ZTEMP, INFO )\n\n*  Purpose\n*  =======\n*\n*  DLAEDA computes the Z vector corresponding to the merge step in the\n*  CURLVLth step of the merge process with TLVLS steps for the CURPBMth\n*  problem.\n*\n\n*  Arguments\n*  =========\n*\n*  N      (input) INTEGER\n*         The dimension of the symmetric tridiagonal matrix.  N >= 0.\n*\n*  TLVLS  (input) INTEGER\n*         The total number of merging levels in the overall divide and\n*         conquer tree.\n*\n*  CURLVL (input) INTEGER\n*         The current level in the overall merge routine,\n*         0 <= curlvl <= tlvls.\n*\n*  CURPBM (input) INTEGER\n*         The current problem in the current level in the overall\n*         merge routine (counting from upper left to lower right).\n*\n*  PRMPTR (input) INTEGER array, dimension (N lg N)\n*         Contains a list of pointers which indicate where in PERM a\n*         level's permutation is stored.  PRMPTR(i+1) - PRMPTR(i)\n*         indicates the size of the permutation and incidentally the\n*         size of the full, non-deflated problem.\n*\n*  PERM   (input) INTEGER array, dimension (N lg N)\n*         Contains the permutations (from deflation and sorting) to be\n*         applied to each eigenblock.\n*\n*  GIVPTR (input) INTEGER array, dimension (N lg N)\n*         Contains a list of pointers which indicate where in GIVCOL a\n*         level's Givens rotations are stored.  GIVPTR(i+1) - GIVPTR(i)\n*         indicates the number of Givens rotations.\n*\n*  GIVCOL (input) INTEGER array, dimension (2, N lg N)\n*         Each pair of numbers indicates a pair of columns to take place\n*         in a Givens rotation.\n*\n*  GIVNUM (input) DOUBLE PRECISION array, dimension (2, N lg N)\n*         Each number indicates the S value to be used in the\n*         corresponding Givens rotation.\n*\n*  Q      (input) DOUBLE PRECISION array, dimension (N**2)\n*         Contains the square eigenblocks from previous levels, the\n*         starting positions for blocks are given by QPTR.\n*\n*  QPTR   (input) INTEGER array, dimension (N+2)\n*         Contains a list of pointers which indicate where in Q an\n*         eigenblock is stored.  SQRT( QPTR(i+1) - QPTR(i) ) indicates\n*         the size of the block.\n*\n*  Z      (output) DOUBLE PRECISION array, dimension (N)\n*         On output this vector contains the updating vector (the last\n*         row of the first sub-eigenvector matrix and the first row of\n*         the second sub-eigenvector matrix).\n*\n*  ZTEMP  (workspace) DOUBLE PRECISION array, dimension (N)\n*\n*  INFO   (output) INTEGER\n*          = 0:  successful exit.\n*          < 0:  if INFO = -i, the i-th argument had an illegal value.\n*\n\n*  Further Details\n*  ===============\n*\n*  Based on contributions by\n*     Jeff Rutter, Computer Science Division, University of California\n*     at Berkeley, USA\n*\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 10)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 10)", argc);
  rb_tlvls = argv[0];
  rb_curlvl = argv[1];
  rb_curpbm = argv[2];
  rb_prmptr = argv[3];
  rb_perm = argv[4];
  rb_givptr = argv[5];
  rb_givcol = argv[6];
  rb_givnum = argv[7];
  rb_q = argv[8];
  rb_qptr = argv[9];

  tlvls = NUM2INT(rb_tlvls);
  curpbm = NUM2INT(rb_curpbm);
  curlvl = NUM2INT(rb_curlvl);
  if (!NA_IsNArray(rb_qptr))
    rb_raise(rb_eArgError, "qptr (1th argument) must be NArray");
  if (NA_RANK(rb_qptr) != 1)
    rb_raise(rb_eArgError, "rank of qptr (1th argument) must be %d", 1);
  ldqptr = NA_SHAPE0(rb_qptr);
  if (NA_TYPE(rb_qptr) != NA_LINT)
    rb_qptr = na_change_type(rb_qptr, NA_LINT);
  qptr = NA_PTR_TYPE(rb_qptr, integer*);
  if (!NA_IsNArray(rb_givptr))
    rb_raise(rb_eArgError, "givptr (2th argument) must be NArray");
  if (NA_RANK(rb_givptr) != 1)
    rb_raise(rb_eArgError, "rank of givptr (2th argument) must be %d", 1);
  n = ldqptr-2;
  if (NA_SHAPE0(rb_givptr) != (n*LG(n)))
    rb_raise(rb_eRuntimeError, "shape 0 of givptr must be %d", n*LG(n));
  if (NA_TYPE(rb_givptr) != NA_LINT)
    rb_givptr = na_change_type(rb_givptr, NA_LINT);
  givptr = NA_PTR_TYPE(rb_givptr, integer*);
  if (!NA_IsNArray(rb_givcol))
    rb_raise(rb_eArgError, "givcol (4th argument) must be NArray");
  if (NA_RANK(rb_givcol) != 2)
    rb_raise(rb_eArgError, "rank of givcol (4th argument) must be %d", 2);
  if (NA_SHAPE0(rb_givcol) != (2))
    rb_raise(rb_eRuntimeError, "shape 0 of givcol must be %d", 2);
  if (NA_SHAPE1(rb_givcol) != (n*LG(n)))
    rb_raise(rb_eRuntimeError, "shape 1 of givcol must be %d", n*LG(n));
  if (NA_TYPE(rb_givcol) != NA_LINT)
    rb_givcol = na_change_type(rb_givcol, NA_LINT);
  givcol = NA_PTR_TYPE(rb_givcol, integer*);
  if (!NA_IsNArray(rb_givnum))
    rb_raise(rb_eArgError, "givnum (6th argument) must be NArray");
  if (NA_RANK(rb_givnum) != 2)
    rb_raise(rb_eArgError, "rank of givnum (6th argument) must be %d", 2);
  if (NA_SHAPE0(rb_givnum) != (2))
    rb_raise(rb_eRuntimeError, "shape 0 of givnum must be %d", 2);
  if (NA_SHAPE1(rb_givnum) != (n*LG(n)))
    rb_raise(rb_eRuntimeError, "shape 1 of givnum must be %d", n*LG(n));
  if (NA_TYPE(rb_givnum) != NA_DFLOAT)
    rb_givnum = na_change_type(rb_givnum, NA_DFLOAT);
  givnum = NA_PTR_TYPE(rb_givnum, doublereal*);
  if (!NA_IsNArray(rb_perm))
    rb_raise(rb_eArgError, "perm (7th argument) must be NArray");
  if (NA_RANK(rb_perm) != 1)
    rb_raise(rb_eArgError, "rank of perm (7th argument) must be %d", 1);
  if (NA_SHAPE0(rb_perm) != (n*LG(n)))
    rb_raise(rb_eRuntimeError, "shape 0 of perm must be %d", n*LG(n));
  if (NA_TYPE(rb_perm) != NA_LINT)
    rb_perm = na_change_type(rb_perm, NA_LINT);
  perm = NA_PTR_TYPE(rb_perm, integer*);
  if (!NA_IsNArray(rb_q))
    rb_raise(rb_eArgError, "q (8th argument) must be NArray");
  if (NA_RANK(rb_q) != 1)
    rb_raise(rb_eArgError, "rank of q (8th argument) must be %d", 1);
  if (NA_SHAPE0(rb_q) != (pow(n,2)))
    rb_raise(rb_eRuntimeError, "shape 0 of q must be %d", pow(n,2));
  if (NA_TYPE(rb_q) != NA_DFLOAT)
    rb_q = na_change_type(rb_q, NA_DFLOAT);
  q = NA_PTR_TYPE(rb_q, doublereal*);
  if (!NA_IsNArray(rb_prmptr))
    rb_raise(rb_eArgError, "prmptr (9th argument) must be NArray");
  if (NA_RANK(rb_prmptr) != 1)
    rb_raise(rb_eArgError, "rank of prmptr (9th argument) must be %d", 1);
  if (NA_SHAPE0(rb_prmptr) != (n*LG(n)))
    rb_raise(rb_eRuntimeError, "shape 0 of prmptr must be %d", n*LG(n));
  if (NA_TYPE(rb_prmptr) != NA_LINT)
    rb_prmptr = na_change_type(rb_prmptr, NA_LINT);
  prmptr = NA_PTR_TYPE(rb_prmptr, integer*);
  {
    int shape[1];
    shape[0] = n;
    rb_z = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  z = NA_PTR_TYPE(rb_z, doublereal*);
  ztemp = ALLOC_N(doublereal, (n));

  dlaeda_(&n, &tlvls, &curlvl, &curpbm, prmptr, perm, givptr, givcol, givnum, q, qptr, z, ztemp, &info);

  free(ztemp);
  rb_info = INT2NUM(info);
  return rb_ary_new3(2, rb_z, rb_info);
}

void
init_lapack_dlaeda(VALUE mLapack){
  rb_define_module_function(mLapack, "dlaeda", rb_dlaeda, -1);
}

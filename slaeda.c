#include "rb_lapack.h"

extern VOID slaeda_(integer *n, integer *tlvls, integer *curlvl, integer *curpbm, integer *prmptr, integer *perm, integer *givptr, integer *givcol, real *givnum, real *q, integer *qptr, real *z, real *ztemp, integer *info);

static VALUE sHelp, sUsage;

static VALUE
rblapack_slaeda(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_tlvls;
  integer tlvls; 
  VALUE rblapack_curlvl;
  integer curlvl; 
  VALUE rblapack_curpbm;
  integer curpbm; 
  VALUE rblapack_prmptr;
  integer *prmptr; 
  VALUE rblapack_perm;
  integer *perm; 
  VALUE rblapack_givptr;
  integer *givptr; 
  VALUE rblapack_givcol;
  integer *givcol; 
  VALUE rblapack_givnum;
  real *givnum; 
  VALUE rblapack_q;
  real *q; 
  VALUE rblapack_qptr;
  integer *qptr; 
  VALUE rblapack_z;
  real *z; 
  VALUE rblapack_info;
  integer info; 
  real *ztemp;

  integer ldqptr;
  integer n;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  z, info = NumRu::Lapack.slaeda( tlvls, curlvl, curpbm, prmptr, perm, givptr, givcol, givnum, q, qptr, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE SLAEDA( N, TLVLS, CURLVL, CURPBM, PRMPTR, PERM, GIVPTR, GIVCOL, GIVNUM, Q, QPTR, Z, ZTEMP, INFO )\n\n*  Purpose\n*  =======\n*\n*  SLAEDA computes the Z vector corresponding to the merge step in the\n*  CURLVLth step of the merge process with TLVLS steps for the CURPBMth\n*  problem.\n*\n\n*  Arguments\n*  =========\n*\n*  N      (input) INTEGER\n*         The dimension of the symmetric tridiagonal matrix.  N >= 0.\n*\n*  TLVLS  (input) INTEGER\n*         The total number of merging levels in the overall divide and\n*         conquer tree.\n*\n*  CURLVL (input) INTEGER\n*         The current level in the overall merge routine,\n*         0 <= curlvl <= tlvls.\n*\n*  CURPBM (input) INTEGER\n*         The current problem in the current level in the overall\n*         merge routine (counting from upper left to lower right).\n*\n*  PRMPTR (input) INTEGER array, dimension (N lg N)\n*         Contains a list of pointers which indicate where in PERM a\n*         level's permutation is stored.  PRMPTR(i+1) - PRMPTR(i)\n*         indicates the size of the permutation and incidentally the\n*         size of the full, non-deflated problem.\n*\n*  PERM   (input) INTEGER array, dimension (N lg N)\n*         Contains the permutations (from deflation and sorting) to be\n*         applied to each eigenblock.\n*\n*  GIVPTR (input) INTEGER array, dimension (N lg N)\n*         Contains a list of pointers which indicate where in GIVCOL a\n*         level's Givens rotations are stored.  GIVPTR(i+1) - GIVPTR(i)\n*         indicates the number of Givens rotations.\n*\n*  GIVCOL (input) INTEGER array, dimension (2, N lg N)\n*         Each pair of numbers indicates a pair of columns to take place\n*         in a Givens rotation.\n*\n*  GIVNUM (input) REAL array, dimension (2, N lg N)\n*         Each number indicates the S value to be used in the\n*         corresponding Givens rotation.\n*\n*  Q      (input) REAL array, dimension (N**2)\n*         Contains the square eigenblocks from previous levels, the\n*         starting positions for blocks are given by QPTR.\n*\n*  QPTR   (input) INTEGER array, dimension (N+2)\n*         Contains a list of pointers which indicate where in Q an\n*         eigenblock is stored.  SQRT( QPTR(i+1) - QPTR(i) ) indicates\n*         the size of the block.\n*\n*  Z      (output) REAL array, dimension (N)\n*         On output this vector contains the updating vector (the last\n*         row of the first sub-eigenvector matrix and the first row of\n*         the second sub-eigenvector matrix).\n*\n*  ZTEMP  (workspace) REAL array, dimension (N)\n*\n*  INFO   (output) INTEGER\n*          = 0:  successful exit.\n*          < 0:  if INFO = -i, the i-th argument had an illegal value.\n*\n\n*  Further Details\n*  ===============\n*\n*  Based on contributions by\n*     Jeff Rutter, Computer Science Division, University of California\n*     at Berkeley, USA\n*\n*  =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  z, info = NumRu::Lapack.slaeda( tlvls, curlvl, curpbm, prmptr, perm, givptr, givcol, givnum, q, qptr, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 10)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 10)", argc);
  rblapack_tlvls = argv[0];
  rblapack_curlvl = argv[1];
  rblapack_curpbm = argv[2];
  rblapack_prmptr = argv[3];
  rblapack_perm = argv[4];
  rblapack_givptr = argv[5];
  rblapack_givcol = argv[6];
  rblapack_givnum = argv[7];
  rblapack_q = argv[8];
  rblapack_qptr = argv[9];
  if (rb_options != Qnil) {
  }

  curpbm = NUM2INT(rblapack_curpbm);
  if (!NA_IsNArray(rblapack_qptr))
    rb_raise(rb_eArgError, "qptr (10th argument) must be NArray");
  if (NA_RANK(rblapack_qptr) != 1)
    rb_raise(rb_eArgError, "rank of qptr (10th argument) must be %d", 1);
  ldqptr = NA_SHAPE0(rblapack_qptr);
  if (NA_TYPE(rblapack_qptr) != NA_LINT)
    rblapack_qptr = na_change_type(rblapack_qptr, NA_LINT);
  qptr = NA_PTR_TYPE(rblapack_qptr, integer*);
  tlvls = NUM2INT(rblapack_tlvls);
  curlvl = NUM2INT(rblapack_curlvl);
  n = ldqptr-2;
  if (!NA_IsNArray(rblapack_perm))
    rb_raise(rb_eArgError, "perm (5th argument) must be NArray");
  if (NA_RANK(rblapack_perm) != 1)
    rb_raise(rb_eArgError, "rank of perm (5th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_perm) != (n*LG(n)))
    rb_raise(rb_eRuntimeError, "shape 0 of perm must be %d", n*LG(n));
  if (NA_TYPE(rblapack_perm) != NA_LINT)
    rblapack_perm = na_change_type(rblapack_perm, NA_LINT);
  perm = NA_PTR_TYPE(rblapack_perm, integer*);
  if (!NA_IsNArray(rblapack_prmptr))
    rb_raise(rb_eArgError, "prmptr (4th argument) must be NArray");
  if (NA_RANK(rblapack_prmptr) != 1)
    rb_raise(rb_eArgError, "rank of prmptr (4th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_prmptr) != (n*LG(n)))
    rb_raise(rb_eRuntimeError, "shape 0 of prmptr must be %d", n*LG(n));
  if (NA_TYPE(rblapack_prmptr) != NA_LINT)
    rblapack_prmptr = na_change_type(rblapack_prmptr, NA_LINT);
  prmptr = NA_PTR_TYPE(rblapack_prmptr, integer*);
  if (!NA_IsNArray(rblapack_givptr))
    rb_raise(rb_eArgError, "givptr (6th argument) must be NArray");
  if (NA_RANK(rblapack_givptr) != 1)
    rb_raise(rb_eArgError, "rank of givptr (6th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_givptr) != (n*LG(n)))
    rb_raise(rb_eRuntimeError, "shape 0 of givptr must be %d", n*LG(n));
  if (NA_TYPE(rblapack_givptr) != NA_LINT)
    rblapack_givptr = na_change_type(rblapack_givptr, NA_LINT);
  givptr = NA_PTR_TYPE(rblapack_givptr, integer*);
  if (!NA_IsNArray(rblapack_q))
    rb_raise(rb_eArgError, "q (9th argument) must be NArray");
  if (NA_RANK(rblapack_q) != 1)
    rb_raise(rb_eArgError, "rank of q (9th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_q) != (pow(n,2)))
    rb_raise(rb_eRuntimeError, "shape 0 of q must be %d", pow(n,2));
  if (NA_TYPE(rblapack_q) != NA_SFLOAT)
    rblapack_q = na_change_type(rblapack_q, NA_SFLOAT);
  q = NA_PTR_TYPE(rblapack_q, real*);
  if (!NA_IsNArray(rblapack_givcol))
    rb_raise(rb_eArgError, "givcol (7th argument) must be NArray");
  if (NA_RANK(rblapack_givcol) != 2)
    rb_raise(rb_eArgError, "rank of givcol (7th argument) must be %d", 2);
  if (NA_SHAPE1(rblapack_givcol) != (n*LG(n)))
    rb_raise(rb_eRuntimeError, "shape 1 of givcol must be %d", n*LG(n));
  if (NA_SHAPE0(rblapack_givcol) != (2))
    rb_raise(rb_eRuntimeError, "shape 0 of givcol must be %d", 2);
  if (NA_TYPE(rblapack_givcol) != NA_LINT)
    rblapack_givcol = na_change_type(rblapack_givcol, NA_LINT);
  givcol = NA_PTR_TYPE(rblapack_givcol, integer*);
  if (!NA_IsNArray(rblapack_givnum))
    rb_raise(rb_eArgError, "givnum (8th argument) must be NArray");
  if (NA_RANK(rblapack_givnum) != 2)
    rb_raise(rb_eArgError, "rank of givnum (8th argument) must be %d", 2);
  if (NA_SHAPE1(rblapack_givnum) != (n*LG(n)))
    rb_raise(rb_eRuntimeError, "shape 1 of givnum must be %d", n*LG(n));
  if (NA_SHAPE0(rblapack_givnum) != (2))
    rb_raise(rb_eRuntimeError, "shape 0 of givnum must be %d", 2);
  if (NA_TYPE(rblapack_givnum) != NA_SFLOAT)
    rblapack_givnum = na_change_type(rblapack_givnum, NA_SFLOAT);
  givnum = NA_PTR_TYPE(rblapack_givnum, real*);
  {
    int shape[1];
    shape[0] = n;
    rblapack_z = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  z = NA_PTR_TYPE(rblapack_z, real*);
  ztemp = ALLOC_N(real, (n));

  slaeda_(&n, &tlvls, &curlvl, &curpbm, prmptr, perm, givptr, givcol, givnum, q, qptr, z, ztemp, &info);

  free(ztemp);
  rblapack_info = INT2NUM(info);
  return rb_ary_new3(2, rblapack_z, rblapack_info);
}

void
init_lapack_slaeda(VALUE mLapack){
  rb_define_module_function(mLapack, "slaeda", rblapack_slaeda, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}

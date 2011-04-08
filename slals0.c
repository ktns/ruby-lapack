#include "rb_lapack.h"

extern VOID slals0_(integer *icompq, integer *nl, integer *nr, integer *sqre, integer *nrhs, real *b, integer *ldb, real *bx, integer *ldbx, integer *perm, integer *givptr, integer *givcol, integer *ldgcol, real *givnum, integer *ldgnum, real *poles, real *difl, real *difr, real *z, integer *k, real *c, real *s, real *work, integer *info);

static VALUE sHelp, sUsage;

static VALUE
rblapack_slals0(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_icompq;
  integer icompq; 
  VALUE rblapack_nl;
  integer nl; 
  VALUE rblapack_nr;
  integer nr; 
  VALUE rblapack_sqre;
  integer sqre; 
  VALUE rblapack_b;
  real *b; 
  VALUE rblapack_perm;
  integer *perm; 
  VALUE rblapack_givptr;
  integer givptr; 
  VALUE rblapack_givcol;
  integer *givcol; 
  VALUE rblapack_givnum;
  real *givnum; 
  VALUE rblapack_poles;
  real *poles; 
  VALUE rblapack_difl;
  real *difl; 
  VALUE rblapack_difr;
  real *difr; 
  VALUE rblapack_z;
  real *z; 
  VALUE rblapack_c;
  real c; 
  VALUE rblapack_s;
  real s; 
  VALUE rblapack_info;
  integer info; 
  VALUE rblapack_b_out__;
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

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  info, b = NumRu::Lapack.slals0( icompq, nl, nr, sqre, b, perm, givptr, givcol, givnum, poles, difl, difr, z, c, s, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE SLALS0( ICOMPQ, NL, NR, SQRE, NRHS, B, LDB, BX, LDBX, PERM, GIVPTR, GIVCOL, LDGCOL, GIVNUM, LDGNUM, POLES, DIFL, DIFR, Z, K, C, S, WORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  SLALS0 applies back the multiplying factors of either the left or the\n*  right singular vector matrix of a diagonal matrix appended by a row\n*  to the right hand side matrix B in solving the least squares problem\n*  using the divide-and-conquer SVD approach.\n*\n*  For the left singular vector matrix, three types of orthogonal\n*  matrices are involved:\n*\n*  (1L) Givens rotations: the number of such rotations is GIVPTR; the\n*       pairs of columns/rows they were applied to are stored in GIVCOL;\n*       and the C- and S-values of these rotations are stored in GIVNUM.\n*\n*  (2L) Permutation. The (NL+1)-st row of B is to be moved to the first\n*       row, and for J=2:N, PERM(J)-th row of B is to be moved to the\n*       J-th row.\n*\n*  (3L) The left singular vector matrix of the remaining matrix.\n*\n*  For the right singular vector matrix, four types of orthogonal\n*  matrices are involved:\n*\n*  (1R) The right singular vector matrix of the remaining matrix.\n*\n*  (2R) If SQRE = 1, one extra Givens rotation to generate the right\n*       null space.\n*\n*  (3R) The inverse transformation of (2L).\n*\n*  (4R) The inverse transformation of (1L).\n*\n\n*  Arguments\n*  =========\n*\n*  ICOMPQ (input) INTEGER\n*         Specifies whether singular vectors are to be computed in\n*         factored form:\n*         = 0: Left singular vector matrix.\n*         = 1: Right singular vector matrix.\n*\n*  NL     (input) INTEGER\n*         The row dimension of the upper block. NL >= 1.\n*\n*  NR     (input) INTEGER\n*         The row dimension of the lower block. NR >= 1.\n*\n*  SQRE   (input) INTEGER\n*         = 0: the lower block is an NR-by-NR square matrix.\n*         = 1: the lower block is an NR-by-(NR+1) rectangular matrix.\n*\n*         The bidiagonal matrix has row dimension N = NL + NR + 1,\n*         and column dimension M = N + SQRE.\n*\n*  NRHS   (input) INTEGER\n*         The number of columns of B and BX. NRHS must be at least 1.\n*\n*  B      (input/output) REAL array, dimension ( LDB, NRHS )\n*         On input, B contains the right hand sides of the least\n*         squares problem in rows 1 through M. On output, B contains\n*         the solution X in rows 1 through N.\n*\n*  LDB    (input) INTEGER\n*         The leading dimension of B. LDB must be at least\n*         max(1,MAX( M, N ) ).\n*\n*  BX     (workspace) REAL array, dimension ( LDBX, NRHS )\n*\n*  LDBX   (input) INTEGER\n*         The leading dimension of BX.\n*\n*  PERM   (input) INTEGER array, dimension ( N )\n*         The permutations (from deflation and sorting) applied\n*         to the two blocks.\n*\n*  GIVPTR (input) INTEGER\n*         The number of Givens rotations which took place in this\n*         subproblem.\n*\n*  GIVCOL (input) INTEGER array, dimension ( LDGCOL, 2 )\n*         Each pair of numbers indicates a pair of rows/columns\n*         involved in a Givens rotation.\n*\n*  LDGCOL (input) INTEGER\n*         The leading dimension of GIVCOL, must be at least N.\n*\n*  GIVNUM (input) REAL array, dimension ( LDGNUM, 2 )\n*         Each number indicates the C or S value used in the\n*         corresponding Givens rotation.\n*\n*  LDGNUM (input) INTEGER\n*         The leading dimension of arrays DIFR, POLES and\n*         GIVNUM, must be at least K.\n*\n*  POLES  (input) REAL array, dimension ( LDGNUM, 2 )\n*         On entry, POLES(1:K, 1) contains the new singular\n*         values obtained from solving the secular equation, and\n*         POLES(1:K, 2) is an array containing the poles in the secular\n*         equation.\n*\n*  DIFL   (input) REAL array, dimension ( K ).\n*         On entry, DIFL(I) is the distance between I-th updated\n*         (undeflated) singular value and the I-th (undeflated) old\n*         singular value.\n*\n*  DIFR   (input) REAL array, dimension ( LDGNUM, 2 ).\n*         On entry, DIFR(I, 1) contains the distances between I-th\n*         updated (undeflated) singular value and the I+1-th\n*         (undeflated) old singular value. And DIFR(I, 2) is the\n*         normalizing factor for the I-th right singular vector.\n*\n*  Z      (input) REAL array, dimension ( K )\n*         Contain the components of the deflation-adjusted updating row\n*         vector.\n*\n*  K      (input) INTEGER\n*         Contains the dimension of the non-deflated matrix,\n*         This is the order of the related secular equation. 1 <= K <=N.\n*\n*  C      (input) REAL\n*         C contains garbage if SQRE =0 and the C-value of a Givens\n*         rotation related to the right null space if SQRE = 1.\n*\n*  S      (input) REAL\n*         S contains garbage if SQRE =0 and the S-value of a Givens\n*         rotation related to the right null space if SQRE = 1.\n*\n*  WORK   (workspace) REAL array, dimension ( K )\n*\n*  INFO   (output) INTEGER\n*          = 0:  successful exit.\n*          < 0:  if INFO = -i, the i-th argument had an illegal value.\n*\n\n*  Further Details\n*  ===============\n*\n*  Based on contributions by\n*     Ming Gu and Ren-Cang Li, Computer Science Division, University of\n*       California at Berkeley, USA\n*     Osni Marques, LBNL/NERSC, USA\n*\n*  =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  info, b = NumRu::Lapack.slals0( icompq, nl, nr, sqre, b, perm, givptr, givcol, givnum, poles, difl, difr, z, c, s, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 15)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 15)", argc);
  rblapack_icompq = argv[0];
  rblapack_nl = argv[1];
  rblapack_nr = argv[2];
  rblapack_sqre = argv[3];
  rblapack_b = argv[4];
  rblapack_perm = argv[5];
  rblapack_givptr = argv[6];
  rblapack_givcol = argv[7];
  rblapack_givnum = argv[8];
  rblapack_poles = argv[9];
  rblapack_difl = argv[10];
  rblapack_difr = argv[11];
  rblapack_z = argv[12];
  rblapack_c = argv[13];
  rblapack_s = argv[14];
  if (rb_options != Qnil) {
  }

  if (!NA_IsNArray(rblapack_difl))
    rb_raise(rb_eArgError, "difl (11th argument) must be NArray");
  if (NA_RANK(rblapack_difl) != 1)
    rb_raise(rb_eArgError, "rank of difl (11th argument) must be %d", 1);
  k = NA_SHAPE0(rblapack_difl);
  if (NA_TYPE(rblapack_difl) != NA_SFLOAT)
    rblapack_difl = na_change_type(rblapack_difl, NA_SFLOAT);
  difl = NA_PTR_TYPE(rblapack_difl, real*);
  if (!NA_IsNArray(rblapack_b))
    rb_raise(rb_eArgError, "b (5th argument) must be NArray");
  if (NA_RANK(rblapack_b) != 2)
    rb_raise(rb_eArgError, "rank of b (5th argument) must be %d", 2);
  nrhs = NA_SHAPE1(rblapack_b);
  ldb = NA_SHAPE0(rblapack_b);
  if (NA_TYPE(rblapack_b) != NA_SFLOAT)
    rblapack_b = na_change_type(rblapack_b, NA_SFLOAT);
  b = NA_PTR_TYPE(rblapack_b, real*);
  c = (real)NUM2DBL(rblapack_c);
  if (!NA_IsNArray(rblapack_givcol))
    rb_raise(rb_eArgError, "givcol (8th argument) must be NArray");
  if (NA_RANK(rblapack_givcol) != 2)
    rb_raise(rb_eArgError, "rank of givcol (8th argument) must be %d", 2);
  if (NA_SHAPE1(rblapack_givcol) != (2))
    rb_raise(rb_eRuntimeError, "shape 1 of givcol must be %d", 2);
  ldgcol = NA_SHAPE0(rblapack_givcol);
  if (NA_TYPE(rblapack_givcol) != NA_LINT)
    rblapack_givcol = na_change_type(rblapack_givcol, NA_LINT);
  givcol = NA_PTR_TYPE(rblapack_givcol, integer*);
  if (!NA_IsNArray(rblapack_z))
    rb_raise(rb_eArgError, "z (13th argument) must be NArray");
  if (NA_RANK(rblapack_z) != 1)
    rb_raise(rb_eArgError, "rank of z (13th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_z) != k)
    rb_raise(rb_eRuntimeError, "shape 0 of z must be the same as shape 0 of difl");
  if (NA_TYPE(rblapack_z) != NA_SFLOAT)
    rblapack_z = na_change_type(rblapack_z, NA_SFLOAT);
  z = NA_PTR_TYPE(rblapack_z, real*);
  nr = NUM2INT(rblapack_nr);
  if (!NA_IsNArray(rblapack_poles))
    rb_raise(rb_eArgError, "poles (10th argument) must be NArray");
  if (NA_RANK(rblapack_poles) != 2)
    rb_raise(rb_eArgError, "rank of poles (10th argument) must be %d", 2);
  if (NA_SHAPE1(rblapack_poles) != (2))
    rb_raise(rb_eRuntimeError, "shape 1 of poles must be %d", 2);
  ldgnum = NA_SHAPE0(rblapack_poles);
  if (NA_TYPE(rblapack_poles) != NA_SFLOAT)
    rblapack_poles = na_change_type(rblapack_poles, NA_SFLOAT);
  poles = NA_PTR_TYPE(rblapack_poles, real*);
  icompq = NUM2INT(rblapack_icompq);
  nl = NUM2INT(rblapack_nl);
  sqre = NUM2INT(rblapack_sqre);
  if (!NA_IsNArray(rblapack_givnum))
    rb_raise(rb_eArgError, "givnum (9th argument) must be NArray");
  if (NA_RANK(rblapack_givnum) != 2)
    rb_raise(rb_eArgError, "rank of givnum (9th argument) must be %d", 2);
  if (NA_SHAPE1(rblapack_givnum) != (2))
    rb_raise(rb_eRuntimeError, "shape 1 of givnum must be %d", 2);
  if (NA_SHAPE0(rblapack_givnum) != ldgnum)
    rb_raise(rb_eRuntimeError, "shape 0 of givnum must be the same as shape 0 of poles");
  if (NA_TYPE(rblapack_givnum) != NA_SFLOAT)
    rblapack_givnum = na_change_type(rblapack_givnum, NA_SFLOAT);
  givnum = NA_PTR_TYPE(rblapack_givnum, real*);
  if (!NA_IsNArray(rblapack_perm))
    rb_raise(rb_eArgError, "perm (6th argument) must be NArray");
  if (NA_RANK(rblapack_perm) != 1)
    rb_raise(rb_eArgError, "rank of perm (6th argument) must be %d", 1);
  n = NA_SHAPE0(rblapack_perm);
  if (NA_TYPE(rblapack_perm) != NA_LINT)
    rblapack_perm = na_change_type(rblapack_perm, NA_LINT);
  perm = NA_PTR_TYPE(rblapack_perm, integer*);
  s = (real)NUM2DBL(rblapack_s);
  if (!NA_IsNArray(rblapack_difr))
    rb_raise(rb_eArgError, "difr (12th argument) must be NArray");
  if (NA_RANK(rblapack_difr) != 2)
    rb_raise(rb_eArgError, "rank of difr (12th argument) must be %d", 2);
  if (NA_SHAPE1(rblapack_difr) != (2))
    rb_raise(rb_eRuntimeError, "shape 1 of difr must be %d", 2);
  if (NA_SHAPE0(rblapack_difr) != ldgnum)
    rb_raise(rb_eRuntimeError, "shape 0 of difr must be the same as shape 0 of poles");
  if (NA_TYPE(rblapack_difr) != NA_SFLOAT)
    rblapack_difr = na_change_type(rblapack_difr, NA_SFLOAT);
  difr = NA_PTR_TYPE(rblapack_difr, real*);
  givptr = NUM2INT(rblapack_givptr);
  ldbx = n;
  {
    int shape[2];
    shape[0] = ldb;
    shape[1] = nrhs;
    rblapack_b_out__ = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  b_out__ = NA_PTR_TYPE(rblapack_b_out__, real*);
  MEMCPY(b_out__, b, real, NA_TOTAL(rblapack_b));
  rblapack_b = rblapack_b_out__;
  b = b_out__;
  bx = ALLOC_N(real, (ldbx)*(nrhs));
  work = ALLOC_N(real, (k));

  slals0_(&icompq, &nl, &nr, &sqre, &nrhs, b, &ldb, bx, &ldbx, perm, &givptr, givcol, &ldgcol, givnum, &ldgnum, poles, difl, difr, z, &k, &c, &s, work, &info);

  free(bx);
  free(work);
  rblapack_info = INT2NUM(info);
  return rb_ary_new3(2, rblapack_info, rblapack_b);
}

void
init_lapack_slals0(VALUE mLapack){
  rb_define_module_function(mLapack, "slals0", rblapack_slals0, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}

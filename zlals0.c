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
    printf("%s\n", "USAGE:\n  info, b = NumRu::Lapack.zlals0( icompq, nl, nr, sqre, b, perm, givptr, givcol, givnum, poles, difl, difr, z, c, s)\n    or\n  NumRu::Lapack.zlals0  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE ZLALS0( ICOMPQ, NL, NR, SQRE, NRHS, B, LDB, BX, LDBX, PERM, GIVPTR, GIVCOL, LDGCOL, GIVNUM, LDGNUM, POLES, DIFL, DIFR, Z, K, C, S, RWORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  ZLALS0 applies back the multiplying factors of either the left or the\n*  right singular vector matrix of a diagonal matrix appended by a row\n*  to the right hand side matrix B in solving the least squares problem\n*  using the divide-and-conquer SVD approach.\n*\n*  For the left singular vector matrix, three types of orthogonal\n*  matrices are involved:\n*\n*  (1L) Givens rotations: the number of such rotations is GIVPTR; the\n*       pairs of columns/rows they were applied to are stored in GIVCOL;\n*       and the C- and S-values of these rotations are stored in GIVNUM.\n*\n*  (2L) Permutation. The (NL+1)-st row of B is to be moved to the first\n*       row, and for J=2:N, PERM(J)-th row of B is to be moved to the\n*       J-th row.\n*\n*  (3L) The left singular vector matrix of the remaining matrix.\n*\n*  For the right singular vector matrix, four types of orthogonal\n*  matrices are involved:\n*\n*  (1R) The right singular vector matrix of the remaining matrix.\n*\n*  (2R) If SQRE = 1, one extra Givens rotation to generate the right\n*       null space.\n*\n*  (3R) The inverse transformation of (2L).\n*\n*  (4R) The inverse transformation of (1L).\n*\n\n*  Arguments\n*  =========\n*\n*  ICOMPQ (input) INTEGER\n*         Specifies whether singular vectors are to be computed in\n*         factored form:\n*         = 0: Left singular vector matrix.\n*         = 1: Right singular vector matrix.\n*\n*  NL     (input) INTEGER\n*         The row dimension of the upper block. NL >= 1.\n*\n*  NR     (input) INTEGER\n*         The row dimension of the lower block. NR >= 1.\n*\n*  SQRE   (input) INTEGER\n*         = 0: the lower block is an NR-by-NR square matrix.\n*         = 1: the lower block is an NR-by-(NR+1) rectangular matrix.\n*\n*         The bidiagonal matrix has row dimension N = NL + NR + 1,\n*         and column dimension M = N + SQRE.\n*\n*  NRHS   (input) INTEGER\n*         The number of columns of B and BX. NRHS must be at least 1.\n*\n*  B      (input/output) COMPLEX*16 array, dimension ( LDB, NRHS )\n*         On input, B contains the right hand sides of the least\n*         squares problem in rows 1 through M. On output, B contains\n*         the solution X in rows 1 through N.\n*\n*  LDB    (input) INTEGER\n*         The leading dimension of B. LDB must be at least\n*         max(1,MAX( M, N ) ).\n*\n*  BX     (workspace) COMPLEX*16 array, dimension ( LDBX, NRHS )\n*\n*  LDBX   (input) INTEGER\n*         The leading dimension of BX.\n*\n*  PERM   (input) INTEGER array, dimension ( N )\n*         The permutations (from deflation and sorting) applied\n*         to the two blocks.\n*\n*  GIVPTR (input) INTEGER\n*         The number of Givens rotations which took place in this\n*         subproblem.\n*\n*  GIVCOL (input) INTEGER array, dimension ( LDGCOL, 2 )\n*         Each pair of numbers indicates a pair of rows/columns\n*         involved in a Givens rotation.\n*\n*  LDGCOL (input) INTEGER\n*         The leading dimension of GIVCOL, must be at least N.\n*\n*  GIVNUM (input) DOUBLE PRECISION array, dimension ( LDGNUM, 2 )\n*         Each number indicates the C or S value used in the\n*         corresponding Givens rotation.\n*\n*  LDGNUM (input) INTEGER\n*         The leading dimension of arrays DIFR, POLES and\n*         GIVNUM, must be at least K.\n*\n*  POLES  (input) DOUBLE PRECISION array, dimension ( LDGNUM, 2 )\n*         On entry, POLES(1:K, 1) contains the new singular\n*         values obtained from solving the secular equation, and\n*         POLES(1:K, 2) is an array containing the poles in the secular\n*         equation.\n*\n*  DIFL   (input) DOUBLE PRECISION array, dimension ( K ).\n*         On entry, DIFL(I) is the distance between I-th updated\n*         (undeflated) singular value and the I-th (undeflated) old\n*         singular value.\n*\n*  DIFR   (input) DOUBLE PRECISION array, dimension ( LDGNUM, 2 ).\n*         On entry, DIFR(I, 1) contains the distances between I-th\n*         updated (undeflated) singular value and the I+1-th\n*         (undeflated) old singular value. And DIFR(I, 2) is the\n*         normalizing factor for the I-th right singular vector.\n*\n*  Z      (input) DOUBLE PRECISION array, dimension ( K )\n*         Contain the components of the deflation-adjusted updating row\n*         vector.\n*\n*  K      (input) INTEGER\n*         Contains the dimension of the non-deflated matrix,\n*         This is the order of the related secular equation. 1 <= K <=N.\n*\n*  C      (input) DOUBLE PRECISION\n*         C contains garbage if SQRE =0 and the C-value of a Givens\n*         rotation related to the right null space if SQRE = 1.\n*\n*  S      (input) DOUBLE PRECISION\n*         S contains garbage if SQRE =0 and the S-value of a Givens\n*         rotation related to the right null space if SQRE = 1.\n*\n*  RWORK  (workspace) DOUBLE PRECISION array, dimension\n*         ( K*(1+NRHS) + 2*NRHS )\n*\n*  INFO   (output) INTEGER\n*          = 0:  successful exit.\n*          < 0:  if INFO = -i, the i-th argument had an illegal value.\n*\n\n*  Further Details\n*  ===============\n*\n*  Based on contributions by\n*     Ming Gu and Ren-Cang Li, Computer Science Division, University of\n*       California at Berkeley, USA\n*     Osni Marques, LBNL/NERSC, USA\n*\n*  =====================================================================\n*\n\n");
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

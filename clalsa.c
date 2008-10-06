#include "rb_lapack.h"

static VALUE
rb_clalsa(int argc, VALUE *argv, VALUE self){
  VALUE rb_icompq;
  integer icompq; 
  VALUE rb_b;
  complex *b; 
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
  VALUE rb_bx;
  complex *bx; 
  VALUE rb_info;
  integer info; 
  VALUE rb_b_out__;
  complex *b_out__;
  real *rwork;
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
    printf("%s\n", "USAGE:\n  bx, info, b = NumRu::Lapack.clalsa( icompq, b, u, vt, k, difl, difr, z, poles, givptr, givcol, perm, givnum, c, s)\n    or\n  NumRu::Lapack.clalsa  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE CLALSA( ICOMPQ, SMLSIZ, N, NRHS, B, LDB, BX, LDBX, U, LDU, VT, K, DIFL, DIFR, Z, POLES, GIVPTR, GIVCOL, LDGCOL, PERM, GIVNUM, C, S, RWORK, IWORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  CLALSA is an itermediate step in solving the least squares problem\n*  by computing the SVD of the coefficient matrix in compact form (The\n*  singular vectors are computed as products of simple orthorgonal\n*  matrices.).\n*\n*  If ICOMPQ = 0, CLALSA applies the inverse of the left singular vector\n*  matrix of an upper bidiagonal matrix to the right hand side; and if\n*  ICOMPQ = 1, CLALSA applies the right singular vector matrix to the\n*  right hand side. The singular vector matrices were generated in\n*  compact form by CLALSA.\n*\n\n*  Arguments\n*  =========\n*\n*  ICOMPQ (input) INTEGER\n*         Specifies whether the left or the right singular vector\n*         matrix is involved.\n*         = 0: Left singular vector matrix\n*         = 1: Right singular vector matrix\n*\n*  SMLSIZ (input) INTEGER\n*         The maximum size of the subproblems at the bottom of the\n*         computation tree.\n*\n*  N      (input) INTEGER\n*         The row and column dimensions of the upper bidiagonal matrix.\n*\n*  NRHS   (input) INTEGER\n*         The number of columns of B and BX. NRHS must be at least 1.\n*\n*  B      (input/output) COMPLEX array, dimension ( LDB, NRHS )\n*         On input, B contains the right hand sides of the least\n*         squares problem in rows 1 through M.\n*         On output, B contains the solution X in rows 1 through N.\n*\n*  LDB    (input) INTEGER\n*         The leading dimension of B in the calling subprogram.\n*         LDB must be at least max(1,MAX( M, N ) ).\n*\n*  BX     (output) COMPLEX array, dimension ( LDBX, NRHS )\n*         On exit, the result of applying the left or right singular\n*         vector matrix to B.\n*\n*  LDBX   (input) INTEGER\n*         The leading dimension of BX.\n*\n*  U      (input) REAL array, dimension ( LDU, SMLSIZ ).\n*         On entry, U contains the left singular vector matrices of all\n*         subproblems at the bottom level.\n*\n*  LDU    (input) INTEGER, LDU = > N.\n*         The leading dimension of arrays U, VT, DIFL, DIFR,\n*         POLES, GIVNUM, and Z.\n*\n*  VT     (input) REAL array, dimension ( LDU, SMLSIZ+1 ).\n*         On entry, VT' contains the right singular vector matrices of\n*         all subproblems at the bottom level.\n*\n*  K      (input) INTEGER array, dimension ( N ).\n*\n*  DIFL   (input) REAL array, dimension ( LDU, NLVL ).\n*         where NLVL = INT(log_2 (N/(SMLSIZ+1))) + 1.\n*\n*  DIFR   (input) REAL array, dimension ( LDU, 2 * NLVL ).\n*         On entry, DIFL(*, I) and DIFR(*, 2 * I -1) record\n*         distances between singular values on the I-th level and\n*         singular values on the (I -1)-th level, and DIFR(*, 2 * I)\n*         record the normalizing factors of the right singular vectors\n*         matrices of subproblems on I-th level.\n*\n*  Z      (input) REAL array, dimension ( LDU, NLVL ).\n*         On entry, Z(1, I) contains the components of the deflation-\n*         adjusted updating row vector for subproblems on the I-th\n*         level.\n*\n*  POLES  (input) REAL array, dimension ( LDU, 2 * NLVL ).\n*         On entry, POLES(*, 2 * I -1: 2 * I) contains the new and old\n*         singular values involved in the secular equations on the I-th\n*         level.\n*\n*  GIVPTR (input) INTEGER array, dimension ( N ).\n*         On entry, GIVPTR( I ) records the number of Givens\n*         rotations performed on the I-th problem on the computation\n*         tree.\n*\n*  GIVCOL (input) INTEGER array, dimension ( LDGCOL, 2 * NLVL ).\n*         On entry, for each I, GIVCOL(*, 2 * I - 1: 2 * I) records the\n*         locations of Givens rotations performed on the I-th level on\n*         the computation tree.\n*\n*  LDGCOL (input) INTEGER, LDGCOL = > N.\n*         The leading dimension of arrays GIVCOL and PERM.\n*\n*  PERM   (input) INTEGER array, dimension ( LDGCOL, NLVL ).\n*         On entry, PERM(*, I) records permutations done on the I-th\n*         level of the computation tree.\n*\n*  GIVNUM (input) REAL array, dimension ( LDU, 2 * NLVL ).\n*         On entry, GIVNUM(*, 2 *I -1 : 2 * I) records the C- and S-\n*         values of Givens rotations performed on the I-th level on the\n*         computation tree.\n*\n*  C      (input) REAL array, dimension ( N ).\n*         On entry, if the I-th subproblem is not square,\n*         C( I ) contains the C-value of a Givens rotation related to\n*         the right null space of the I-th subproblem.\n*\n*  S      (input) REAL array, dimension ( N ).\n*         On entry, if the I-th subproblem is not square,\n*         S( I ) contains the S-value of a Givens rotation related to\n*         the right null space of the I-th subproblem.\n*\n*  RWORK  (workspace) REAL array, dimension at least\n*         max ( N, (SMLSZ+1)*NRHS*3 ).\n*\n*  IWORK  (workspace) INTEGER array.\n*         The dimension must be at least 3 * N\n*\n*  INFO   (output) INTEGER\n*          = 0:  successful exit.\n*          < 0:  if INFO = -i, the i-th argument had an illegal value.\n*\n\n*  Further Details\n*  ===============\n*\n*  Based on contributions by\n*     Ming Gu and Ren-Cang Li, Computer Science Division, University of\n*       California at Berkeley, USA\n*     Osni Marques, LBNL/NERSC, USA\n*\n*  =====================================================================\n*\n\n");
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

  icompq = NUM2INT(rb_icompq);
  if (!NA_IsNArray(rb_b))
    rb_raise(rb_eArgError, "b (2th argument) must be NArray");
  if (NA_RANK(rb_b) != 2)
    rb_raise(rb_eArgError, "rank of b (2th argument) must be %d", 2);
  ldb = NA_SHAPE0(rb_b);
  nrhs = NA_SHAPE1(rb_b);
  if (NA_TYPE(rb_b) != NA_SCOMPLEX)
    rb_b = na_change_type(rb_b, NA_SCOMPLEX);
  b = NA_PTR_TYPE(rb_b, complex*);
  if (!NA_IsNArray(rb_u))
    rb_raise(rb_eArgError, "u (3th argument) must be NArray");
  if (NA_RANK(rb_u) != 2)
    rb_raise(rb_eArgError, "rank of u (3th argument) must be %d", 2);
  ldu = NA_SHAPE0(rb_u);
  smlsiz = NA_SHAPE1(rb_u);
  if (NA_TYPE(rb_u) != NA_SFLOAT)
    rb_u = na_change_type(rb_u, NA_SFLOAT);
  u = NA_PTR_TYPE(rb_u, real*);
  if (!NA_IsNArray(rb_vt))
    rb_raise(rb_eArgError, "vt (4th argument) must be NArray");
  if (NA_RANK(rb_vt) != 2)
    rb_raise(rb_eArgError, "rank of vt (4th argument) must be %d", 2);
  if (NA_SHAPE0(rb_vt) != ldu)
    rb_raise(rb_eRuntimeError, "shape 0 of vt must be the same as shape 0 of u");
  if (NA_SHAPE1(rb_vt) != (smlsiz+1))
    rb_raise(rb_eRuntimeError, "shape 1 of vt must be %d", smlsiz+1);
  if (NA_TYPE(rb_vt) != NA_SFLOAT)
    rb_vt = na_change_type(rb_vt, NA_SFLOAT);
  vt = NA_PTR_TYPE(rb_vt, real*);
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
  nlvl = (int)(1.0/log(2.0)*log((double)n/(smlsiz+1))) + 1;
  if (NA_SHAPE0(rb_difl) != ldu)
    rb_raise(rb_eRuntimeError, "shape 0 of difl must be the same as shape 0 of u");
  if (NA_SHAPE1(rb_difl) != nlvl)
    rb_raise(rb_eRuntimeError, "shape 1 of difl must be nlvl");
  if (NA_TYPE(rb_difl) != NA_SFLOAT)
    rb_difl = na_change_type(rb_difl, NA_SFLOAT);
  difl = NA_PTR_TYPE(rb_difl, real*);
  if (!NA_IsNArray(rb_difr))
    rb_raise(rb_eArgError, "difr (7th argument) must be NArray");
  if (NA_RANK(rb_difr) != 2)
    rb_raise(rb_eArgError, "rank of difr (7th argument) must be %d", 2);
  if (NA_SHAPE0(rb_difr) != ldu)
    rb_raise(rb_eRuntimeError, "shape 0 of difr must be the same as shape 0 of u");
  if (NA_SHAPE1(rb_difr) != (2 * nlvl))
    rb_raise(rb_eRuntimeError, "shape 1 of difr must be %d", 2 * nlvl);
  if (NA_TYPE(rb_difr) != NA_SFLOAT)
    rb_difr = na_change_type(rb_difr, NA_SFLOAT);
  difr = NA_PTR_TYPE(rb_difr, real*);
  if (!NA_IsNArray(rb_z))
    rb_raise(rb_eArgError, "z (8th argument) must be NArray");
  if (NA_RANK(rb_z) != 2)
    rb_raise(rb_eArgError, "rank of z (8th argument) must be %d", 2);
  if (NA_SHAPE0(rb_z) != ldu)
    rb_raise(rb_eRuntimeError, "shape 0 of z must be the same as shape 0 of u");
  if (NA_SHAPE1(rb_z) != nlvl)
    rb_raise(rb_eRuntimeError, "shape 1 of z must be nlvl");
  if (NA_TYPE(rb_z) != NA_SFLOAT)
    rb_z = na_change_type(rb_z, NA_SFLOAT);
  z = NA_PTR_TYPE(rb_z, real*);
  if (!NA_IsNArray(rb_poles))
    rb_raise(rb_eArgError, "poles (9th argument) must be NArray");
  if (NA_RANK(rb_poles) != 2)
    rb_raise(rb_eArgError, "rank of poles (9th argument) must be %d", 2);
  if (NA_SHAPE0(rb_poles) != ldu)
    rb_raise(rb_eRuntimeError, "shape 0 of poles must be the same as shape 0 of u");
  if (NA_SHAPE1(rb_poles) != (2 * nlvl))
    rb_raise(rb_eRuntimeError, "shape 1 of poles must be %d", 2 * nlvl);
  if (NA_TYPE(rb_poles) != NA_SFLOAT)
    rb_poles = na_change_type(rb_poles, NA_SFLOAT);
  poles = NA_PTR_TYPE(rb_poles, real*);
  if (!NA_IsNArray(rb_givptr))
    rb_raise(rb_eArgError, "givptr (10th argument) must be NArray");
  if (NA_RANK(rb_givptr) != 1)
    rb_raise(rb_eArgError, "rank of givptr (10th argument) must be %d", 1);
  if (NA_SHAPE0(rb_givptr) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of givptr must be the same as shape 0 of k");
  if (NA_TYPE(rb_givptr) != NA_LINT)
    rb_givptr = na_change_type(rb_givptr, NA_LINT);
  givptr = NA_PTR_TYPE(rb_givptr, integer*);
  if (!NA_IsNArray(rb_givcol))
    rb_raise(rb_eArgError, "givcol (11th argument) must be NArray");
  if (NA_RANK(rb_givcol) != 2)
    rb_raise(rb_eArgError, "rank of givcol (11th argument) must be %d", 2);
  ldgcol = NA_SHAPE0(rb_givcol);
  if (NA_SHAPE1(rb_givcol) != (2 * nlvl))
    rb_raise(rb_eRuntimeError, "shape 1 of givcol must be %d", 2 * nlvl);
  if (NA_TYPE(rb_givcol) != NA_LINT)
    rb_givcol = na_change_type(rb_givcol, NA_LINT);
  givcol = NA_PTR_TYPE(rb_givcol, integer*);
  if (!NA_IsNArray(rb_perm))
    rb_raise(rb_eArgError, "perm (12th argument) must be NArray");
  if (NA_RANK(rb_perm) != 2)
    rb_raise(rb_eArgError, "rank of perm (12th argument) must be %d", 2);
  if (NA_SHAPE0(rb_perm) != ldgcol)
    rb_raise(rb_eRuntimeError, "shape 0 of perm must be the same as shape 0 of givcol");
  if (NA_SHAPE1(rb_perm) != nlvl)
    rb_raise(rb_eRuntimeError, "shape 1 of perm must be nlvl");
  if (NA_TYPE(rb_perm) != NA_LINT)
    rb_perm = na_change_type(rb_perm, NA_LINT);
  perm = NA_PTR_TYPE(rb_perm, integer*);
  if (!NA_IsNArray(rb_givnum))
    rb_raise(rb_eArgError, "givnum (13th argument) must be NArray");
  if (NA_RANK(rb_givnum) != 2)
    rb_raise(rb_eArgError, "rank of givnum (13th argument) must be %d", 2);
  if (NA_SHAPE0(rb_givnum) != ldu)
    rb_raise(rb_eRuntimeError, "shape 0 of givnum must be the same as shape 0 of u");
  if (NA_SHAPE1(rb_givnum) != (2 * nlvl))
    rb_raise(rb_eRuntimeError, "shape 1 of givnum must be %d", 2 * nlvl);
  if (NA_TYPE(rb_givnum) != NA_SFLOAT)
    rb_givnum = na_change_type(rb_givnum, NA_SFLOAT);
  givnum = NA_PTR_TYPE(rb_givnum, real*);
  if (!NA_IsNArray(rb_c))
    rb_raise(rb_eArgError, "c (14th argument) must be NArray");
  if (NA_RANK(rb_c) != 1)
    rb_raise(rb_eArgError, "rank of c (14th argument) must be %d", 1);
  if (NA_SHAPE0(rb_c) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of c must be the same as shape 0 of k");
  if (NA_TYPE(rb_c) != NA_SFLOAT)
    rb_c = na_change_type(rb_c, NA_SFLOAT);
  c = NA_PTR_TYPE(rb_c, real*);
  if (!NA_IsNArray(rb_s))
    rb_raise(rb_eArgError, "s (15th argument) must be NArray");
  if (NA_RANK(rb_s) != 1)
    rb_raise(rb_eArgError, "rank of s (15th argument) must be %d", 1);
  if (NA_SHAPE0(rb_s) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of s must be the same as shape 0 of k");
  if (NA_TYPE(rb_s) != NA_SFLOAT)
    rb_s = na_change_type(rb_s, NA_SFLOAT);
  s = NA_PTR_TYPE(rb_s, real*);
  ldbx = n;
  {
    int shape[2];
    shape[0] = ldbx;
    shape[1] = nrhs;
    rb_bx = na_make_object(NA_SCOMPLEX, 2, shape, cNArray);
  }
  bx = NA_PTR_TYPE(rb_bx, complex*);
  {
    int shape[2];
    shape[0] = ldb;
    shape[1] = nrhs;
    rb_b_out__ = na_make_object(NA_SCOMPLEX, 2, shape, cNArray);
  }
  b_out__ = NA_PTR_TYPE(rb_b_out__, complex*);
  MEMCPY(b_out__, b, complex, NA_TOTAL(rb_b));
  rb_b = rb_b_out__;
  b = b_out__;
  rwork = ALLOC_N(real, (MAX(n,(smlsiz+1)*nrhs*3)));
  iwork = ALLOC_N(integer, (3 * n));

  clalsa_(&icompq, &smlsiz, &n, &nrhs, b, &ldb, bx, &ldbx, u, &ldu, vt, k, difl, difr, z, poles, givptr, givcol, &ldgcol, perm, givnum, c, s, rwork, iwork, &info);

  free(rwork);
  free(iwork);
  rb_info = INT2NUM(info);
  return rb_ary_new3(3, rb_bx, rb_info, rb_b);
}

void
init_lapack_clalsa(VALUE mLapack){
  rb_define_module_function(mLapack, "clalsa", rb_clalsa, -1);
}

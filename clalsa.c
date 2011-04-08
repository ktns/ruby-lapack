#include "rb_lapack.h"

extern VOID clalsa_(integer *icompq, integer *smlsiz, integer *n, integer *nrhs, complex *b, integer *ldb, complex *bx, integer *ldbx, real *u, integer *ldu, real *vt, integer *k, real *difl, real *difr, real *z, real *poles, integer *givptr, integer *givcol, integer *ldgcol, integer *perm, real *givnum, real *c, real *s, real *rwork, integer *iwork, integer *info);

static VALUE sHelp, sUsage;

static VALUE
rblapack_clalsa(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_icompq;
  integer icompq; 
  VALUE rblapack_b;
  complex *b; 
  VALUE rblapack_u;
  real *u; 
  VALUE rblapack_vt;
  real *vt; 
  VALUE rblapack_k;
  integer *k; 
  VALUE rblapack_difl;
  real *difl; 
  VALUE rblapack_difr;
  real *difr; 
  VALUE rblapack_z;
  real *z; 
  VALUE rblapack_poles;
  real *poles; 
  VALUE rblapack_givptr;
  integer *givptr; 
  VALUE rblapack_givcol;
  integer *givcol; 
  VALUE rblapack_perm;
  integer *perm; 
  VALUE rblapack_givnum;
  real *givnum; 
  VALUE rblapack_c;
  real *c; 
  VALUE rblapack_s;
  real *s; 
  VALUE rblapack_bx;
  complex *bx; 
  VALUE rblapack_info;
  integer info; 
  VALUE rblapack_b_out__;
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

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  bx, info, b = NumRu::Lapack.clalsa( icompq, b, u, vt, k, difl, difr, z, poles, givptr, givcol, perm, givnum, c, s, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE CLALSA( ICOMPQ, SMLSIZ, N, NRHS, B, LDB, BX, LDBX, U, LDU, VT, K, DIFL, DIFR, Z, POLES, GIVPTR, GIVCOL, LDGCOL, PERM, GIVNUM, C, S, RWORK, IWORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  CLALSA is an itermediate step in solving the least squares problem\n*  by computing the SVD of the coefficient matrix in compact form (The\n*  singular vectors are computed as products of simple orthorgonal\n*  matrices.).\n*\n*  If ICOMPQ = 0, CLALSA applies the inverse of the left singular vector\n*  matrix of an upper bidiagonal matrix to the right hand side; and if\n*  ICOMPQ = 1, CLALSA applies the right singular vector matrix to the\n*  right hand side. The singular vector matrices were generated in\n*  compact form by CLALSA.\n*\n\n*  Arguments\n*  =========\n*\n*  ICOMPQ (input) INTEGER\n*         Specifies whether the left or the right singular vector\n*         matrix is involved.\n*         = 0: Left singular vector matrix\n*         = 1: Right singular vector matrix\n*\n*  SMLSIZ (input) INTEGER\n*         The maximum size of the subproblems at the bottom of the\n*         computation tree.\n*\n*  N      (input) INTEGER\n*         The row and column dimensions of the upper bidiagonal matrix.\n*\n*  NRHS   (input) INTEGER\n*         The number of columns of B and BX. NRHS must be at least 1.\n*\n*  B      (input/output) COMPLEX array, dimension ( LDB, NRHS )\n*         On input, B contains the right hand sides of the least\n*         squares problem in rows 1 through M.\n*         On output, B contains the solution X in rows 1 through N.\n*\n*  LDB    (input) INTEGER\n*         The leading dimension of B in the calling subprogram.\n*         LDB must be at least max(1,MAX( M, N ) ).\n*\n*  BX     (output) COMPLEX array, dimension ( LDBX, NRHS )\n*         On exit, the result of applying the left or right singular\n*         vector matrix to B.\n*\n*  LDBX   (input) INTEGER\n*         The leading dimension of BX.\n*\n*  U      (input) REAL array, dimension ( LDU, SMLSIZ ).\n*         On entry, U contains the left singular vector matrices of all\n*         subproblems at the bottom level.\n*\n*  LDU    (input) INTEGER, LDU = > N.\n*         The leading dimension of arrays U, VT, DIFL, DIFR,\n*         POLES, GIVNUM, and Z.\n*\n*  VT     (input) REAL array, dimension ( LDU, SMLSIZ+1 ).\n*         On entry, VT' contains the right singular vector matrices of\n*         all subproblems at the bottom level.\n*\n*  K      (input) INTEGER array, dimension ( N ).\n*\n*  DIFL   (input) REAL array, dimension ( LDU, NLVL ).\n*         where NLVL = INT(log_2 (N/(SMLSIZ+1))) + 1.\n*\n*  DIFR   (input) REAL array, dimension ( LDU, 2 * NLVL ).\n*         On entry, DIFL(*, I) and DIFR(*, 2 * I -1) record\n*         distances between singular values on the I-th level and\n*         singular values on the (I -1)-th level, and DIFR(*, 2 * I)\n*         record the normalizing factors of the right singular vectors\n*         matrices of subproblems on I-th level.\n*\n*  Z      (input) REAL array, dimension ( LDU, NLVL ).\n*         On entry, Z(1, I) contains the components of the deflation-\n*         adjusted updating row vector for subproblems on the I-th\n*         level.\n*\n*  POLES  (input) REAL array, dimension ( LDU, 2 * NLVL ).\n*         On entry, POLES(*, 2 * I -1: 2 * I) contains the new and old\n*         singular values involved in the secular equations on the I-th\n*         level.\n*\n*  GIVPTR (input) INTEGER array, dimension ( N ).\n*         On entry, GIVPTR( I ) records the number of Givens\n*         rotations performed on the I-th problem on the computation\n*         tree.\n*\n*  GIVCOL (input) INTEGER array, dimension ( LDGCOL, 2 * NLVL ).\n*         On entry, for each I, GIVCOL(*, 2 * I - 1: 2 * I) records the\n*         locations of Givens rotations performed on the I-th level on\n*         the computation tree.\n*\n*  LDGCOL (input) INTEGER, LDGCOL = > N.\n*         The leading dimension of arrays GIVCOL and PERM.\n*\n*  PERM   (input) INTEGER array, dimension ( LDGCOL, NLVL ).\n*         On entry, PERM(*, I) records permutations done on the I-th\n*         level of the computation tree.\n*\n*  GIVNUM (input) REAL array, dimension ( LDU, 2 * NLVL ).\n*         On entry, GIVNUM(*, 2 *I -1 : 2 * I) records the C- and S-\n*         values of Givens rotations performed on the I-th level on the\n*         computation tree.\n*\n*  C      (input) REAL array, dimension ( N ).\n*         On entry, if the I-th subproblem is not square,\n*         C( I ) contains the C-value of a Givens rotation related to\n*         the right null space of the I-th subproblem.\n*\n*  S      (input) REAL array, dimension ( N ).\n*         On entry, if the I-th subproblem is not square,\n*         S( I ) contains the S-value of a Givens rotation related to\n*         the right null space of the I-th subproblem.\n*\n*  RWORK  (workspace) REAL array, dimension at least\n*         MAX( (SMLSZ+1)*NRHS*3, N*(1+NRHS) + 2*NRHS ).\n*\n*  IWORK  (workspace) INTEGER array.\n*         The dimension must be at least 3 * N\n*\n*  INFO   (output) INTEGER\n*          = 0:  successful exit.\n*          < 0:  if INFO = -i, the i-th argument had an illegal value.\n*\n\n*  Further Details\n*  ===============\n*\n*  Based on contributions by\n*     Ming Gu and Ren-Cang Li, Computer Science Division, University of\n*       California at Berkeley, USA\n*     Osni Marques, LBNL/NERSC, USA\n*\n*  =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  bx, info, b = NumRu::Lapack.clalsa( icompq, b, u, vt, k, difl, difr, z, poles, givptr, givcol, perm, givnum, c, s, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 15)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 15)", argc);
  rblapack_icompq = argv[0];
  rblapack_b = argv[1];
  rblapack_u = argv[2];
  rblapack_vt = argv[3];
  rblapack_k = argv[4];
  rblapack_difl = argv[5];
  rblapack_difr = argv[6];
  rblapack_z = argv[7];
  rblapack_poles = argv[8];
  rblapack_givptr = argv[9];
  rblapack_givcol = argv[10];
  rblapack_perm = argv[11];
  rblapack_givnum = argv[12];
  rblapack_c = argv[13];
  rblapack_s = argv[14];
  if (rb_options != Qnil) {
  }

  if (!NA_IsNArray(rblapack_k))
    rb_raise(rb_eArgError, "k (5th argument) must be NArray");
  if (NA_RANK(rblapack_k) != 1)
    rb_raise(rb_eArgError, "rank of k (5th argument) must be %d", 1);
  n = NA_SHAPE0(rblapack_k);
  if (NA_TYPE(rblapack_k) != NA_LINT)
    rblapack_k = na_change_type(rblapack_k, NA_LINT);
  k = NA_PTR_TYPE(rblapack_k, integer*);
  if (!NA_IsNArray(rblapack_difl))
    rb_raise(rb_eArgError, "difl (6th argument) must be NArray");
  if (NA_RANK(rblapack_difl) != 2)
    rb_raise(rb_eArgError, "rank of difl (6th argument) must be %d", 2);
  nlvl = NA_SHAPE1(rblapack_difl);
  if (nlvl != ((int)(1.0/log(2.0)*log((double)n/(smlsiz+1))) + 1))
    rb_raise(rb_eRuntimeError, "shape 1 of difl must be %d", (int)(1.0/log(2.0)*log((double)n/(smlsiz+1))) + 1);
  ldu = NA_SHAPE0(rblapack_difl);
  if (NA_TYPE(rblapack_difl) != NA_SFLOAT)
    rblapack_difl = na_change_type(rblapack_difl, NA_SFLOAT);
  difl = NA_PTR_TYPE(rblapack_difl, real*);
  if (!NA_IsNArray(rblapack_b))
    rb_raise(rb_eArgError, "b (2th argument) must be NArray");
  if (NA_RANK(rblapack_b) != 2)
    rb_raise(rb_eArgError, "rank of b (2th argument) must be %d", 2);
  nrhs = NA_SHAPE1(rblapack_b);
  ldb = NA_SHAPE0(rblapack_b);
  if (NA_TYPE(rblapack_b) != NA_SCOMPLEX)
    rblapack_b = na_change_type(rblapack_b, NA_SCOMPLEX);
  b = NA_PTR_TYPE(rblapack_b, complex*);
  if (!NA_IsNArray(rblapack_c))
    rb_raise(rb_eArgError, "c (14th argument) must be NArray");
  if (NA_RANK(rblapack_c) != 1)
    rb_raise(rb_eArgError, "rank of c (14th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_c) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of c must be the same as shape 0 of k");
  if (NA_TYPE(rblapack_c) != NA_SFLOAT)
    rblapack_c = na_change_type(rblapack_c, NA_SFLOAT);
  c = NA_PTR_TYPE(rblapack_c, real*);
  if (!NA_IsNArray(rblapack_u))
    rb_raise(rb_eArgError, "u (3th argument) must be NArray");
  if (NA_RANK(rblapack_u) != 2)
    rb_raise(rb_eArgError, "rank of u (3th argument) must be %d", 2);
  smlsiz = NA_SHAPE1(rblapack_u);
  if (NA_SHAPE0(rblapack_u) != ldu)
    rb_raise(rb_eRuntimeError, "shape 0 of u must be the same as shape 0 of difl");
  if (NA_TYPE(rblapack_u) != NA_SFLOAT)
    rblapack_u = na_change_type(rblapack_u, NA_SFLOAT);
  u = NA_PTR_TYPE(rblapack_u, real*);
  if (!NA_IsNArray(rblapack_z))
    rb_raise(rb_eArgError, "z (8th argument) must be NArray");
  if (NA_RANK(rblapack_z) != 2)
    rb_raise(rb_eArgError, "rank of z (8th argument) must be %d", 2);
  if (NA_SHAPE1(rblapack_z) != nlvl)
    rb_raise(rb_eRuntimeError, "shape 1 of z must be the same as shape 1 of difl");
  if (NA_SHAPE0(rblapack_z) != ldu)
    rb_raise(rb_eRuntimeError, "shape 0 of z must be the same as shape 0 of difl");
  if (NA_TYPE(rblapack_z) != NA_SFLOAT)
    rblapack_z = na_change_type(rblapack_z, NA_SFLOAT);
  z = NA_PTR_TYPE(rblapack_z, real*);
  if (!NA_IsNArray(rblapack_s))
    rb_raise(rb_eArgError, "s (15th argument) must be NArray");
  if (NA_RANK(rblapack_s) != 1)
    rb_raise(rb_eArgError, "rank of s (15th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_s) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of s must be the same as shape 0 of k");
  if (NA_TYPE(rblapack_s) != NA_SFLOAT)
    rblapack_s = na_change_type(rblapack_s, NA_SFLOAT);
  s = NA_PTR_TYPE(rblapack_s, real*);
  icompq = NUM2INT(rblapack_icompq);
  if (!NA_IsNArray(rblapack_perm))
    rb_raise(rb_eArgError, "perm (12th argument) must be NArray");
  if (NA_RANK(rblapack_perm) != 2)
    rb_raise(rb_eArgError, "rank of perm (12th argument) must be %d", 2);
  if (NA_SHAPE1(rblapack_perm) != nlvl)
    rb_raise(rb_eRuntimeError, "shape 1 of perm must be the same as shape 1 of difl");
  ldgcol = NA_SHAPE0(rblapack_perm);
  if (NA_TYPE(rblapack_perm) != NA_LINT)
    rblapack_perm = na_change_type(rblapack_perm, NA_LINT);
  perm = NA_PTR_TYPE(rblapack_perm, integer*);
  if (!NA_IsNArray(rblapack_givptr))
    rb_raise(rb_eArgError, "givptr (10th argument) must be NArray");
  if (NA_RANK(rblapack_givptr) != 1)
    rb_raise(rb_eArgError, "rank of givptr (10th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_givptr) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of givptr must be the same as shape 0 of k");
  if (NA_TYPE(rblapack_givptr) != NA_LINT)
    rblapack_givptr = na_change_type(rblapack_givptr, NA_LINT);
  givptr = NA_PTR_TYPE(rblapack_givptr, integer*);
  ldbx = n;
  if (!NA_IsNArray(rblapack_poles))
    rb_raise(rb_eArgError, "poles (9th argument) must be NArray");
  if (NA_RANK(rblapack_poles) != 2)
    rb_raise(rb_eArgError, "rank of poles (9th argument) must be %d", 2);
  if (NA_SHAPE1(rblapack_poles) != (2 * nlvl))
    rb_raise(rb_eRuntimeError, "shape 1 of poles must be %d", 2 * nlvl);
  if (NA_SHAPE0(rblapack_poles) != ldu)
    rb_raise(rb_eRuntimeError, "shape 0 of poles must be the same as shape 0 of difl");
  if (NA_TYPE(rblapack_poles) != NA_SFLOAT)
    rblapack_poles = na_change_type(rblapack_poles, NA_SFLOAT);
  poles = NA_PTR_TYPE(rblapack_poles, real*);
  if (!NA_IsNArray(rblapack_difr))
    rb_raise(rb_eArgError, "difr (7th argument) must be NArray");
  if (NA_RANK(rblapack_difr) != 2)
    rb_raise(rb_eArgError, "rank of difr (7th argument) must be %d", 2);
  if (NA_SHAPE1(rblapack_difr) != (2 * nlvl))
    rb_raise(rb_eRuntimeError, "shape 1 of difr must be %d", 2 * nlvl);
  if (NA_SHAPE0(rblapack_difr) != ldu)
    rb_raise(rb_eRuntimeError, "shape 0 of difr must be the same as shape 0 of difl");
  if (NA_TYPE(rblapack_difr) != NA_SFLOAT)
    rblapack_difr = na_change_type(rblapack_difr, NA_SFLOAT);
  difr = NA_PTR_TYPE(rblapack_difr, real*);
  if (!NA_IsNArray(rblapack_vt))
    rb_raise(rb_eArgError, "vt (4th argument) must be NArray");
  if (NA_RANK(rblapack_vt) != 2)
    rb_raise(rb_eArgError, "rank of vt (4th argument) must be %d", 2);
  if (NA_SHAPE1(rblapack_vt) != (smlsiz+1))
    rb_raise(rb_eRuntimeError, "shape 1 of vt must be %d", smlsiz+1);
  if (NA_SHAPE0(rblapack_vt) != ldu)
    rb_raise(rb_eRuntimeError, "shape 0 of vt must be the same as shape 0 of difl");
  if (NA_TYPE(rblapack_vt) != NA_SFLOAT)
    rblapack_vt = na_change_type(rblapack_vt, NA_SFLOAT);
  vt = NA_PTR_TYPE(rblapack_vt, real*);
  nlvl = (int)(1.0/log(2.0)*log((double)n/(smlsiz+1))) + 1;
  if (!NA_IsNArray(rblapack_givnum))
    rb_raise(rb_eArgError, "givnum (13th argument) must be NArray");
  if (NA_RANK(rblapack_givnum) != 2)
    rb_raise(rb_eArgError, "rank of givnum (13th argument) must be %d", 2);
  if (NA_SHAPE1(rblapack_givnum) != (2 * nlvl))
    rb_raise(rb_eRuntimeError, "shape 1 of givnum must be %d", 2 * nlvl);
  if (NA_SHAPE0(rblapack_givnum) != ldu)
    rb_raise(rb_eRuntimeError, "shape 0 of givnum must be the same as shape 0 of difl");
  if (NA_TYPE(rblapack_givnum) != NA_SFLOAT)
    rblapack_givnum = na_change_type(rblapack_givnum, NA_SFLOAT);
  givnum = NA_PTR_TYPE(rblapack_givnum, real*);
  if (!NA_IsNArray(rblapack_givcol))
    rb_raise(rb_eArgError, "givcol (11th argument) must be NArray");
  if (NA_RANK(rblapack_givcol) != 2)
    rb_raise(rb_eArgError, "rank of givcol (11th argument) must be %d", 2);
  if (NA_SHAPE1(rblapack_givcol) != (2 * nlvl))
    rb_raise(rb_eRuntimeError, "shape 1 of givcol must be %d", 2 * nlvl);
  if (NA_SHAPE0(rblapack_givcol) != ldgcol)
    rb_raise(rb_eRuntimeError, "shape 0 of givcol must be the same as shape 0 of perm");
  if (NA_TYPE(rblapack_givcol) != NA_LINT)
    rblapack_givcol = na_change_type(rblapack_givcol, NA_LINT);
  givcol = NA_PTR_TYPE(rblapack_givcol, integer*);
  {
    int shape[2];
    shape[0] = ldbx;
    shape[1] = nrhs;
    rblapack_bx = na_make_object(NA_SCOMPLEX, 2, shape, cNArray);
  }
  bx = NA_PTR_TYPE(rblapack_bx, complex*);
  {
    int shape[2];
    shape[0] = ldb;
    shape[1] = nrhs;
    rblapack_b_out__ = na_make_object(NA_SCOMPLEX, 2, shape, cNArray);
  }
  b_out__ = NA_PTR_TYPE(rblapack_b_out__, complex*);
  MEMCPY(b_out__, b, complex, NA_TOTAL(rblapack_b));
  rblapack_b = rblapack_b_out__;
  b = b_out__;
  rwork = ALLOC_N(real, (MAX(n,(smlsiz+1)*nrhs*3)));
  iwork = ALLOC_N(integer, (3 * n));

  clalsa_(&icompq, &smlsiz, &n, &nrhs, b, &ldb, bx, &ldbx, u, &ldu, vt, k, difl, difr, z, poles, givptr, givcol, &ldgcol, perm, givnum, c, s, rwork, iwork, &info);

  free(rwork);
  free(iwork);
  rblapack_info = INT2NUM(info);
  return rb_ary_new3(3, rblapack_bx, rblapack_info, rblapack_b);
}

void
init_lapack_clalsa(VALUE mLapack){
  rb_define_module_function(mLapack, "clalsa", rblapack_clalsa, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}

#include "rb_lapack.h"

extern VOID dlasda_(integer *icompq, integer *smlsiz, integer *n, integer *sqre, doublereal *d, doublereal *e, doublereal *u, integer *ldu, doublereal *vt, integer *k, doublereal *difl, doublereal *difr, doublereal *z, doublereal *poles, integer *givptr, integer *givcol, integer *ldgcol, integer *perm, doublereal *givnum, doublereal *c, doublereal *s, doublereal *work, integer *iwork, integer *info);

static VALUE sHelp, sUsage;

static VALUE
rblapack_dlasda(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_icompq;
  integer icompq; 
  VALUE rblapack_smlsiz;
  integer smlsiz; 
  VALUE rblapack_sqre;
  integer sqre; 
  VALUE rblapack_d;
  doublereal *d; 
  VALUE rblapack_e;
  doublereal *e; 
  VALUE rblapack_u;
  doublereal *u; 
  VALUE rblapack_vt;
  doublereal *vt; 
  VALUE rblapack_k;
  integer *k; 
  VALUE rblapack_difl;
  doublereal *difl; 
  VALUE rblapack_difr;
  doublereal *difr; 
  VALUE rblapack_z;
  doublereal *z; 
  VALUE rblapack_poles;
  doublereal *poles; 
  VALUE rblapack_givptr;
  integer *givptr; 
  VALUE rblapack_givcol;
  integer *givcol; 
  VALUE rblapack_perm;
  integer *perm; 
  VALUE rblapack_givnum;
  doublereal *givnum; 
  VALUE rblapack_c;
  doublereal *c; 
  VALUE rblapack_s;
  doublereal *s; 
  VALUE rblapack_info;
  integer info; 
  VALUE rblapack_d_out__;
  doublereal *d_out__;
  doublereal *work;
  integer *iwork;

  integer n;
  integer ldu;
  integer nlvl;
  integer ldgcol;
  integer m;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  u, vt, k, difl, difr, z, poles, givptr, givcol, perm, givnum, c, s, info, d = NumRu::Lapack.dlasda( icompq, smlsiz, sqre, d, e, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE DLASDA( ICOMPQ, SMLSIZ, N, SQRE, D, E, U, LDU, VT, K, DIFL, DIFR, Z, POLES, GIVPTR, GIVCOL, LDGCOL, PERM, GIVNUM, C, S, WORK, IWORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  Using a divide and conquer approach, DLASDA computes the singular\n*  value decomposition (SVD) of a real upper bidiagonal N-by-M matrix\n*  B with diagonal D and offdiagonal E, where M = N + SQRE. The\n*  algorithm computes the singular values in the SVD B = U * S * VT.\n*  The orthogonal matrices U and VT are optionally computed in\n*  compact form.\n*\n*  A related subroutine, DLASD0, computes the singular values and\n*  the singular vectors in explicit form.\n*\n\n*  Arguments\n*  =========\n*\n*  ICOMPQ (input) INTEGER\n*         Specifies whether singular vectors are to be computed\n*         in compact form, as follows\n*         = 0: Compute singular values only.\n*         = 1: Compute singular vectors of upper bidiagonal\n*              matrix in compact form.\n*\n*  SMLSIZ (input) INTEGER\n*         The maximum size of the subproblems at the bottom of the\n*         computation tree.\n*\n*  N      (input) INTEGER\n*         The row dimension of the upper bidiagonal matrix. This is\n*         also the dimension of the main diagonal array D.\n*\n*  SQRE   (input) INTEGER\n*         Specifies the column dimension of the bidiagonal matrix.\n*         = 0: The bidiagonal matrix has column dimension M = N;\n*         = 1: The bidiagonal matrix has column dimension M = N + 1.\n*\n*  D      (input/output) DOUBLE PRECISION array, dimension ( N )\n*         On entry D contains the main diagonal of the bidiagonal\n*         matrix. On exit D, if INFO = 0, contains its singular values.\n*\n*  E      (input) DOUBLE PRECISION array, dimension ( M-1 )\n*         Contains the subdiagonal entries of the bidiagonal matrix.\n*         On exit, E has been destroyed.\n*\n*  U      (output) DOUBLE PRECISION array,\n*         dimension ( LDU, SMLSIZ ) if ICOMPQ = 1, and not referenced\n*         if ICOMPQ = 0. If ICOMPQ = 1, on exit, U contains the left\n*         singular vector matrices of all subproblems at the bottom\n*         level.\n*\n*  LDU    (input) INTEGER, LDU = > N.\n*         The leading dimension of arrays U, VT, DIFL, DIFR, POLES,\n*         GIVNUM, and Z.\n*\n*  VT     (output) DOUBLE PRECISION array,\n*         dimension ( LDU, SMLSIZ+1 ) if ICOMPQ = 1, and not referenced\n*         if ICOMPQ = 0. If ICOMPQ = 1, on exit, VT' contains the right\n*         singular vector matrices of all subproblems at the bottom\n*         level.\n*\n*  K      (output) INTEGER array,\n*         dimension ( N ) if ICOMPQ = 1 and dimension 1 if ICOMPQ = 0.\n*         If ICOMPQ = 1, on exit, K(I) is the dimension of the I-th\n*         secular equation on the computation tree.\n*\n*  DIFL   (output) DOUBLE PRECISION array, dimension ( LDU, NLVL ),\n*         where NLVL = floor(log_2 (N/SMLSIZ))).\n*\n*  DIFR   (output) DOUBLE PRECISION array,\n*                  dimension ( LDU, 2 * NLVL ) if ICOMPQ = 1 and\n*                  dimension ( N ) if ICOMPQ = 0.\n*         If ICOMPQ = 1, on exit, DIFL(1:N, I) and DIFR(1:N, 2 * I - 1)\n*         record distances between singular values on the I-th\n*         level and singular values on the (I -1)-th level, and\n*         DIFR(1:N, 2 * I ) contains the normalizing factors for\n*         the right singular vector matrix. See DLASD8 for details.\n*\n*  Z      (output) DOUBLE PRECISION array,\n*                  dimension ( LDU, NLVL ) if ICOMPQ = 1 and\n*                  dimension ( N ) if ICOMPQ = 0.\n*         The first K elements of Z(1, I) contain the components of\n*         the deflation-adjusted updating row vector for subproblems\n*         on the I-th level.\n*\n*  POLES  (output) DOUBLE PRECISION array,\n*         dimension ( LDU, 2 * NLVL ) if ICOMPQ = 1, and not referenced\n*         if ICOMPQ = 0. If ICOMPQ = 1, on exit, POLES(1, 2*I - 1) and\n*         POLES(1, 2*I) contain  the new and old singular values\n*         involved in the secular equations on the I-th level.\n*\n*  GIVPTR (output) INTEGER array,\n*         dimension ( N ) if ICOMPQ = 1, and not referenced if\n*         ICOMPQ = 0. If ICOMPQ = 1, on exit, GIVPTR( I ) records\n*         the number of Givens rotations performed on the I-th\n*         problem on the computation tree.\n*\n*  GIVCOL (output) INTEGER array,\n*         dimension ( LDGCOL, 2 * NLVL ) if ICOMPQ = 1, and not\n*         referenced if ICOMPQ = 0. If ICOMPQ = 1, on exit, for each I,\n*         GIVCOL(1, 2 *I - 1) and GIVCOL(1, 2 *I) record the locations\n*         of Givens rotations performed on the I-th level on the\n*         computation tree.\n*\n*  LDGCOL (input) INTEGER, LDGCOL = > N.\n*         The leading dimension of arrays GIVCOL and PERM.\n*\n*  PERM   (output) INTEGER array,\n*         dimension ( LDGCOL, NLVL ) if ICOMPQ = 1, and not referenced\n*         if ICOMPQ = 0. If ICOMPQ = 1, on exit, PERM(1, I) records\n*         permutations done on the I-th level of the computation tree.\n*\n*  GIVNUM (output) DOUBLE PRECISION array,\n*         dimension ( LDU,  2 * NLVL ) if ICOMPQ = 1, and not\n*         referenced if ICOMPQ = 0. If ICOMPQ = 1, on exit, for each I,\n*         GIVNUM(1, 2 *I - 1) and GIVNUM(1, 2 *I) record the C- and S-\n*         values of Givens rotations performed on the I-th level on\n*         the computation tree.\n*\n*  C      (output) DOUBLE PRECISION array,\n*         dimension ( N ) if ICOMPQ = 1, and dimension 1 if ICOMPQ = 0.\n*         If ICOMPQ = 1 and the I-th subproblem is not square, on exit,\n*         C( I ) contains the C-value of a Givens rotation related to\n*         the right null space of the I-th subproblem.\n*\n*  S      (output) DOUBLE PRECISION array, dimension ( N ) if\n*         ICOMPQ = 1, and dimension 1 if ICOMPQ = 0. If ICOMPQ = 1\n*         and the I-th subproblem is not square, on exit, S( I )\n*         contains the S-value of a Givens rotation related to\n*         the right null space of the I-th subproblem.\n*\n*  WORK   (workspace) DOUBLE PRECISION array, dimension\n*         (6 * N + (SMLSIZ + 1)*(SMLSIZ + 1)).\n*\n*  IWORK  (workspace) INTEGER array.\n*         Dimension must be at least (7 * N).\n*\n*  INFO   (output) INTEGER\n*          = 0:  successful exit.\n*          < 0:  if INFO = -i, the i-th argument had an illegal value.\n*          > 0:  if INFO = 1, a singular value did not converge\n*\n\n*  Further Details\n*  ===============\n*\n*  Based on contributions by\n*     Ming Gu and Huan Ren, Computer Science Division, University of\n*     California at Berkeley, USA\n*\n*  =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  u, vt, k, difl, difr, z, poles, givptr, givcol, perm, givnum, c, s, info, d = NumRu::Lapack.dlasda( icompq, smlsiz, sqre, d, e, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 5)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 5)", argc);
  rblapack_icompq = argv[0];
  rblapack_smlsiz = argv[1];
  rblapack_sqre = argv[2];
  rblapack_d = argv[3];
  rblapack_e = argv[4];
  if (rb_options != Qnil) {
  }

  sqre = NUM2INT(rblapack_sqre);
  if (!NA_IsNArray(rblapack_d))
    rb_raise(rb_eArgError, "d (4th argument) must be NArray");
  if (NA_RANK(rblapack_d) != 1)
    rb_raise(rb_eArgError, "rank of d (4th argument) must be %d", 1);
  n = NA_SHAPE0(rblapack_d);
  if (NA_TYPE(rblapack_d) != NA_DFLOAT)
    rblapack_d = na_change_type(rblapack_d, NA_DFLOAT);
  d = NA_PTR_TYPE(rblapack_d, doublereal*);
  smlsiz = NUM2INT(rblapack_smlsiz);
  icompq = NUM2INT(rblapack_icompq);
  m = sqre == 0 ? n : sqre == 1 ? n+1 : 0;
  ldu = n;
  nlvl = floor(1.0/log(2.0)*log((double)n/smlsiz));
  if (!NA_IsNArray(rblapack_e))
    rb_raise(rb_eArgError, "e (5th argument) must be NArray");
  if (NA_RANK(rblapack_e) != 1)
    rb_raise(rb_eArgError, "rank of e (5th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_e) != (m-1))
    rb_raise(rb_eRuntimeError, "shape 0 of e must be %d", m-1);
  if (NA_TYPE(rblapack_e) != NA_DFLOAT)
    rblapack_e = na_change_type(rblapack_e, NA_DFLOAT);
  e = NA_PTR_TYPE(rblapack_e, doublereal*);
  ldgcol = n;
  {
    int shape[2];
    shape[0] = ldu;
    shape[1] = MAX(1,smlsiz);
    rblapack_u = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  u = NA_PTR_TYPE(rblapack_u, doublereal*);
  {
    int shape[2];
    shape[0] = ldu;
    shape[1] = smlsiz+1;
    rblapack_vt = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  vt = NA_PTR_TYPE(rblapack_vt, doublereal*);
  {
    int shape[1];
    shape[0] = icompq == 1 ? n : icompq == 0 ? 1 : 0;
    rblapack_k = na_make_object(NA_LINT, 1, shape, cNArray);
  }
  k = NA_PTR_TYPE(rblapack_k, integer*);
  {
    int shape[2];
    shape[0] = ldu;
    shape[1] = nlvl;
    rblapack_difl = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  difl = NA_PTR_TYPE(rblapack_difl, doublereal*);
  {
    int shape[2];
    shape[0] = icompq == 1 ? ldu : icompq == 0 ? n : 0;
    shape[1] = icompq == 1 ? 2 * nlvl : 0;
    rblapack_difr = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  difr = NA_PTR_TYPE(rblapack_difr, doublereal*);
  {
    int shape[2];
    shape[0] = icompq == 1 ? ldu : icompq == 0 ? n : 0;
    shape[1] = icompq == 1 ? nlvl : 0;
    rblapack_z = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  z = NA_PTR_TYPE(rblapack_z, doublereal*);
  {
    int shape[2];
    shape[0] = ldu;
    shape[1] = 2 * nlvl;
    rblapack_poles = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  poles = NA_PTR_TYPE(rblapack_poles, doublereal*);
  {
    int shape[1];
    shape[0] = n;
    rblapack_givptr = na_make_object(NA_LINT, 1, shape, cNArray);
  }
  givptr = NA_PTR_TYPE(rblapack_givptr, integer*);
  {
    int shape[2];
    shape[0] = ldgcol;
    shape[1] = 2 * nlvl;
    rblapack_givcol = na_make_object(NA_LINT, 2, shape, cNArray);
  }
  givcol = NA_PTR_TYPE(rblapack_givcol, integer*);
  {
    int shape[2];
    shape[0] = ldgcol;
    shape[1] = nlvl;
    rblapack_perm = na_make_object(NA_LINT, 2, shape, cNArray);
  }
  perm = NA_PTR_TYPE(rblapack_perm, integer*);
  {
    int shape[2];
    shape[0] = ldu;
    shape[1] = 2 * nlvl;
    rblapack_givnum = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  givnum = NA_PTR_TYPE(rblapack_givnum, doublereal*);
  {
    int shape[1];
    shape[0] = icompq == 1 ? n : icompq == 0 ? 1 : 0;
    rblapack_c = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  c = NA_PTR_TYPE(rblapack_c, doublereal*);
  {
    int shape[1];
    shape[0] = icompq==1 ? n : icompq==0 ? 1 : 0;
    rblapack_s = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  s = NA_PTR_TYPE(rblapack_s, doublereal*);
  {
    int shape[1];
    shape[0] = n;
    rblapack_d_out__ = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  d_out__ = NA_PTR_TYPE(rblapack_d_out__, doublereal*);
  MEMCPY(d_out__, d, doublereal, NA_TOTAL(rblapack_d));
  rblapack_d = rblapack_d_out__;
  d = d_out__;
  work = ALLOC_N(doublereal, (6 * n + (smlsiz + 1)*(smlsiz + 1)));
  iwork = ALLOC_N(integer, ((7 * n)));

  dlasda_(&icompq, &smlsiz, &n, &sqre, d, e, u, &ldu, vt, k, difl, difr, z, poles, givptr, givcol, &ldgcol, perm, givnum, c, s, work, iwork, &info);

  free(work);
  free(iwork);
  rblapack_info = INT2NUM(info);
  return rb_ary_new3(15, rblapack_u, rblapack_vt, rblapack_k, rblapack_difl, rblapack_difr, rblapack_z, rblapack_poles, rblapack_givptr, rblapack_givcol, rblapack_perm, rblapack_givnum, rblapack_c, rblapack_s, rblapack_info, rblapack_d);
}

void
init_lapack_dlasda(VALUE mLapack){
  rb_define_module_function(mLapack, "dlasda", rblapack_dlasda, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}

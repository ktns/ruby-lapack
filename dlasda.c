#include "rb_lapack.h"

static VALUE
rb_dlasda(int argc, VALUE *argv, VALUE self){
  VALUE rb_icompq;
  integer icompq; 
  VALUE rb_smlsiz;
  integer smlsiz; 
  VALUE rb_sqre;
  integer sqre; 
  VALUE rb_d;
  doublereal *d; 
  VALUE rb_e;
  doublereal *e; 
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
  VALUE rb_info;
  integer info; 
  VALUE rb_d_out__;
  doublereal *d_out__;
  doublereal *work;
  integer *iwork;

  integer n;
  integer m;
  integer ldu;
  integer nlvl;
  integer ldgcol;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  u, vt, k, difl, difr, z, poles, givptr, givcol, perm, givnum, c, s, info, d = NumRu::Lapack.dlasda( icompq, smlsiz, sqre, d, e)\n    or\n  NumRu::Lapack.dlasda  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE DLASDA( ICOMPQ, SMLSIZ, N, SQRE, D, E, U, LDU, VT, K, DIFL, DIFR, Z, POLES, GIVPTR, GIVCOL, LDGCOL, PERM, GIVNUM, C, S, WORK, IWORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  Using a divide and conquer approach, DLASDA computes the singular\n*  value decomposition (SVD) of a real upper bidiagonal N-by-M matrix\n*  B with diagonal D and offdiagonal E, where M = N + SQRE. The\n*  algorithm computes the singular values in the SVD B = U * S * VT.\n*  The orthogonal matrices U and VT are optionally computed in\n*  compact form.\n*\n*  A related subroutine, DLASD0, computes the singular values and\n*  the singular vectors in explicit form.\n*\n\n*  Arguments\n*  =========\n*\n*  ICOMPQ (input) INTEGER\n*         Specifies whether singular vectors are to be computed\n*         in compact form, as follows\n*         = 0: Compute singular values only.\n*         = 1: Compute singular vectors of upper bidiagonal\n*              matrix in compact form.\n*\n*  SMLSIZ (input) INTEGER\n*         The maximum size of the subproblems at the bottom of the\n*         computation tree.\n*\n*  N      (input) INTEGER\n*         The row dimension of the upper bidiagonal matrix. This is\n*         also the dimension of the main diagonal array D.\n*\n*  SQRE   (input) INTEGER\n*         Specifies the column dimension of the bidiagonal matrix.\n*         = 0: The bidiagonal matrix has column dimension M = N;\n*         = 1: The bidiagonal matrix has column dimension M = N + 1.\n*\n*  D      (input/output) DOUBLE PRECISION array, dimension ( N )\n*         On entry D contains the main diagonal of the bidiagonal\n*         matrix. On exit D, if INFO = 0, contains its singular values.\n*\n*  E      (input) DOUBLE PRECISION array, dimension ( M-1 )\n*         Contains the subdiagonal entries of the bidiagonal matrix.\n*         On exit, E has been destroyed.\n*\n*  U      (output) DOUBLE PRECISION array,\n*         dimension ( LDU, SMLSIZ ) if ICOMPQ = 1, and not referenced\n*         if ICOMPQ = 0. If ICOMPQ = 1, on exit, U contains the left\n*         singular vector matrices of all subproblems at the bottom\n*         level.\n*\n*  LDU    (input) INTEGER, LDU = > N.\n*         The leading dimension of arrays U, VT, DIFL, DIFR, POLES,\n*         GIVNUM, and Z.\n*\n*  VT     (output) DOUBLE PRECISION array,\n*         dimension ( LDU, SMLSIZ+1 ) if ICOMPQ = 1, and not referenced\n*         if ICOMPQ = 0. If ICOMPQ = 1, on exit, VT' contains the right\n*         singular vector matrices of all subproblems at the bottom\n*         level.\n*\n*  K      (output) INTEGER array,\n*         dimension ( N ) if ICOMPQ = 1 and dimension 1 if ICOMPQ = 0.\n*         If ICOMPQ = 1, on exit, K(I) is the dimension of the I-th\n*         secular equation on the computation tree.\n*\n*  DIFL   (output) DOUBLE PRECISION array, dimension ( LDU, NLVL ),\n*         where NLVL = floor(log_2 (N/SMLSIZ))).\n*\n*  DIFR   (output) DOUBLE PRECISION array,\n*                  dimension ( LDU, 2 * NLVL ) if ICOMPQ = 1 and\n*                  dimension ( N ) if ICOMPQ = 0.\n*         If ICOMPQ = 1, on exit, DIFL(1:N, I) and DIFR(1:N, 2 * I - 1)\n*         record distances between singular values on the I-th\n*         level and singular values on the (I -1)-th level, and\n*         DIFR(1:N, 2 * I ) contains the normalizing factors for\n*         the right singular vector matrix. See DLASD8 for details.\n*\n*  Z      (output) DOUBLE PRECISION array,\n*                  dimension ( LDU, NLVL ) if ICOMPQ = 1 and\n*                  dimension ( N ) if ICOMPQ = 0.\n*         The first K elements of Z(1, I) contain the components of\n*         the deflation-adjusted updating row vector for subproblems\n*         on the I-th level.\n*\n*  POLES  (output) DOUBLE PRECISION array,\n*         dimension ( LDU, 2 * NLVL ) if ICOMPQ = 1, and not referenced\n*         if ICOMPQ = 0. If ICOMPQ = 1, on exit, POLES(1, 2*I - 1) and\n*         POLES(1, 2*I) contain  the new and old singular values\n*         involved in the secular equations on the I-th level.\n*\n*  GIVPTR (output) INTEGER array,\n*         dimension ( N ) if ICOMPQ = 1, and not referenced if\n*         ICOMPQ = 0. If ICOMPQ = 1, on exit, GIVPTR( I ) records\n*         the number of Givens rotations performed on the I-th\n*         problem on the computation tree.\n*\n*  GIVCOL (output) INTEGER array,\n*         dimension ( LDGCOL, 2 * NLVL ) if ICOMPQ = 1, and not\n*         referenced if ICOMPQ = 0. If ICOMPQ = 1, on exit, for each I,\n*         GIVCOL(1, 2 *I - 1) and GIVCOL(1, 2 *I) record the locations\n*         of Givens rotations performed on the I-th level on the\n*         computation tree.\n*\n*  LDGCOL (input) INTEGER, LDGCOL = > N.\n*         The leading dimension of arrays GIVCOL and PERM.\n*\n*  PERM   (output) INTEGER array,\n*         dimension ( LDGCOL, NLVL ) if ICOMPQ = 1, and not referenced\n*         if ICOMPQ = 0. If ICOMPQ = 1, on exit, PERM(1, I) records\n*         permutations done on the I-th level of the computation tree.\n*\n*  GIVNUM (output) DOUBLE PRECISION array,\n*         dimension ( LDU,  2 * NLVL ) if ICOMPQ = 1, and not\n*         referenced if ICOMPQ = 0. If ICOMPQ = 1, on exit, for each I,\n*         GIVNUM(1, 2 *I - 1) and GIVNUM(1, 2 *I) record the C- and S-\n*         values of Givens rotations performed on the I-th level on\n*         the computation tree.\n*\n*  C      (output) DOUBLE PRECISION array,\n*         dimension ( N ) if ICOMPQ = 1, and dimension 1 if ICOMPQ = 0.\n*         If ICOMPQ = 1 and the I-th subproblem is not square, on exit,\n*         C( I ) contains the C-value of a Givens rotation related to\n*         the right null space of the I-th subproblem.\n*\n*  S      (output) DOUBLE PRECISION array, dimension ( N ) if\n*         ICOMPQ = 1, and dimension 1 if ICOMPQ = 0. If ICOMPQ = 1\n*         and the I-th subproblem is not square, on exit, S( I )\n*         contains the S-value of a Givens rotation related to\n*         the right null space of the I-th subproblem.\n*\n*  WORK   (workspace) DOUBLE PRECISION array, dimension\n*         (6 * N + (SMLSIZ + 1)*(SMLSIZ + 1)).\n*\n*  IWORK  (workspace) INTEGER array.\n*         Dimension must be at least (7 * N).\n*\n*  INFO   (output) INTEGER\n*          = 0:  successful exit.\n*          < 0:  if INFO = -i, the i-th argument had an illegal value.\n*          > 0:  if INFO = 1, an singular value did not converge\n*\n\n*  Further Details\n*  ===============\n*\n*  Based on contributions by\n*     Ming Gu and Huan Ren, Computer Science Division, University of\n*     California at Berkeley, USA\n*\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 5)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 5)", argc);
  rb_icompq = argv[0];
  rb_smlsiz = argv[1];
  rb_sqre = argv[2];
  rb_d = argv[3];
  rb_e = argv[4];

  icompq = NUM2INT(rb_icompq);
  smlsiz = NUM2INT(rb_smlsiz);
  sqre = NUM2INT(rb_sqre);
  if (!NA_IsNArray(rb_d))
    rb_raise(rb_eArgError, "d (4th argument) must be NArray");
  if (NA_RANK(rb_d) != 1)
    rb_raise(rb_eArgError, "rank of d (4th argument) must be %d", 1);
  n = NA_SHAPE0(rb_d);
  if (NA_TYPE(rb_d) != NA_DFLOAT)
    rb_d = na_change_type(rb_d, NA_DFLOAT);
  d = NA_PTR_TYPE(rb_d, doublereal*);
  if (!NA_IsNArray(rb_e))
    rb_raise(rb_eArgError, "e (5th argument) must be NArray");
  if (NA_RANK(rb_e) != 1)
    rb_raise(rb_eArgError, "rank of e (5th argument) must be %d", 1);
  m = sqre == 0 ? n : sqre == 1 ? n+1 : 0;
  ldu = n;
  if (NA_SHAPE0(rb_e) != (m-1))
    rb_raise(rb_eRuntimeError, "shape 0 of e must be %d", m-1);
  if (NA_TYPE(rb_e) != NA_DFLOAT)
    rb_e = na_change_type(rb_e, NA_DFLOAT);
  e = NA_PTR_TYPE(rb_e, doublereal*);
  ldu = n;
  {
    int shape[2];
    shape[0] = ldu;
    shape[1] = MAX(1,smlsiz);
    rb_u = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  u = NA_PTR_TYPE(rb_u, doublereal*);
  {
    int shape[2];
    shape[0] = ldu;
    shape[1] = smlsiz+1;
    rb_vt = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  vt = NA_PTR_TYPE(rb_vt, doublereal*);
  {
    int shape[1];
    shape[0] = icompq == 1 ? n : icompq == 0 ? 1 : 0;
    rb_k = na_make_object(NA_LINT, 1, shape, cNArray);
  }
  k = NA_PTR_TYPE(rb_k, integer*);
  nlvl = floor(1.0/log(2.0)*log((double)n/smlsiz));
  {
    int shape[2];
    shape[0] = ldu;
    shape[1] = nlvl;
    rb_difl = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  difl = NA_PTR_TYPE(rb_difl, doublereal*);
  {
    int shape[2];
    shape[0] = icompq == 1 ? ldu : icompq == 0 ? n : 0;
    shape[1] = icompq == 1 ? 2 * nlvl : 0;
    rb_difr = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  difr = NA_PTR_TYPE(rb_difr, doublereal*);
  {
    int shape[2];
    shape[0] = icompq == 1 ? ldu : icompq == 0 ? n : 0;
    shape[1] = icompq == 1 ? nlvl : 0;
    rb_z = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  z = NA_PTR_TYPE(rb_z, doublereal*);
  {
    int shape[2];
    shape[0] = ldu;
    shape[1] = 2 * nlvl;
    rb_poles = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  poles = NA_PTR_TYPE(rb_poles, doublereal*);
  {
    int shape[1];
    shape[0] = n;
    rb_givptr = na_make_object(NA_LINT, 1, shape, cNArray);
  }
  givptr = NA_PTR_TYPE(rb_givptr, integer*);
  ldgcol = n;
  {
    int shape[2];
    shape[0] = ldgcol;
    shape[1] = 2 * nlvl;
    rb_givcol = na_make_object(NA_LINT, 2, shape, cNArray);
  }
  givcol = NA_PTR_TYPE(rb_givcol, integer*);
  {
    int shape[2];
    shape[0] = ldgcol;
    shape[1] = nlvl;
    rb_perm = na_make_object(NA_LINT, 2, shape, cNArray);
  }
  perm = NA_PTR_TYPE(rb_perm, integer*);
  {
    int shape[2];
    shape[0] = ldu;
    shape[1] = 2 * nlvl;
    rb_givnum = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  givnum = NA_PTR_TYPE(rb_givnum, doublereal*);
  {
    int shape[1];
    shape[0] = icompq == 1 ? n : icompq == 0 ? 1 : 0;
    rb_c = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  c = NA_PTR_TYPE(rb_c, doublereal*);
  {
    int shape[1];
    shape[0] = icompq==1 ? n : icompq==0 ? 1 : 0;
    rb_s = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  s = NA_PTR_TYPE(rb_s, doublereal*);
  {
    int shape[1];
    shape[0] = n;
    rb_d_out__ = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  d_out__ = NA_PTR_TYPE(rb_d_out__, doublereal*);
  MEMCPY(d_out__, d, doublereal, NA_TOTAL(rb_d));
  rb_d = rb_d_out__;
  d = d_out__;
  work = ALLOC_N(doublereal, (6 * n + (smlsiz + 1)*(smlsiz + 1)));
  iwork = ALLOC_N(integer, ((7 * n)));

  dlasda_(&icompq, &smlsiz, &n, &sqre, d, e, u, &ldu, vt, k, difl, difr, z, poles, givptr, givcol, &ldgcol, perm, givnum, c, s, work, iwork, &info);

  free(work);
  free(iwork);
  rb_info = INT2NUM(info);
  return rb_ary_new3(15, rb_u, rb_vt, rb_k, rb_difl, rb_difr, rb_z, rb_poles, rb_givptr, rb_givcol, rb_perm, rb_givnum, rb_c, rb_s, rb_info, rb_d);
}

void
init_lapack_dlasda(VALUE mLapack){
  rb_define_module_function(mLapack, "dlasda", rb_dlasda, -1);
}

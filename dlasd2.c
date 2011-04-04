#include "rb_lapack.h"

extern VOID dlasd2_(integer *nl, integer *nr, integer *sqre, integer *k, doublereal *d, doublereal *z, doublereal *alpha, doublereal *beta, doublereal *u, integer *ldu, doublereal *vt, integer *ldvt, doublereal *dsigma, doublereal *u2, integer *ldu2, doublereal *vt2, integer *ldvt2, integer *idxp, integer *idx, integer *idxc, integer *idxq, integer *coltyp, integer *info);

static VALUE
rb_dlasd2(int argc, VALUE *argv, VALUE self){
  VALUE rb_nl;
  integer nl; 
  VALUE rb_nr;
  integer nr; 
  VALUE rb_sqre;
  integer sqre; 
  VALUE rb_d;
  doublereal *d; 
  VALUE rb_alpha;
  doublereal alpha; 
  VALUE rb_beta;
  doublereal beta; 
  VALUE rb_u;
  doublereal *u; 
  VALUE rb_vt;
  doublereal *vt; 
  VALUE rb_idxq;
  integer *idxq; 
  VALUE rb_k;
  integer k; 
  VALUE rb_z;
  doublereal *z; 
  VALUE rb_dsigma;
  doublereal *dsigma; 
  VALUE rb_u2;
  doublereal *u2; 
  VALUE rb_vt2;
  doublereal *vt2; 
  VALUE rb_idxc;
  integer *idxc; 
  VALUE rb_coltyp;
  integer *coltyp; 
  VALUE rb_info;
  integer info; 
  VALUE rb_d_out__;
  doublereal *d_out__;
  VALUE rb_u_out__;
  doublereal *u_out__;
  VALUE rb_vt_out__;
  doublereal *vt_out__;
  VALUE rb_idxq_out__;
  integer *idxq_out__;
  integer *idxp;
  integer *idx;

  integer n;
  integer ldu;
  integer ldvt;
  integer m;
  integer ldu2;
  integer ldvt2;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  k, z, dsigma, u2, vt2, idxc, coltyp, info, d, u, vt, idxq = NumRu::Lapack.dlasd2( nl, nr, sqre, d, alpha, beta, u, vt, idxq)\n    or\n  NumRu::Lapack.dlasd2  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE DLASD2( NL, NR, SQRE, K, D, Z, ALPHA, BETA, U, LDU, VT, LDVT, DSIGMA, U2, LDU2, VT2, LDVT2, IDXP, IDX, IDXC, IDXQ, COLTYP, INFO )\n\n*  Purpose\n*  =======\n*\n*  DLASD2 merges the two sets of singular values together into a single\n*  sorted set.  Then it tries to deflate the size of the problem.\n*  There are two ways in which deflation can occur:  when two or more\n*  singular values are close together or if there is a tiny entry in the\n*  Z vector.  For each such occurrence the order of the related secular\n*  equation problem is reduced by one.\n*\n*  DLASD2 is called from DLASD1.\n*\n\n*  Arguments\n*  =========\n*\n*  NL     (input) INTEGER\n*         The row dimension of the upper block.  NL >= 1.\n*\n*  NR     (input) INTEGER\n*         The row dimension of the lower block.  NR >= 1.\n*\n*  SQRE   (input) INTEGER\n*         = 0: the lower block is an NR-by-NR square matrix.\n*         = 1: the lower block is an NR-by-(NR+1) rectangular matrix.\n*\n*         The bidiagonal matrix has N = NL + NR + 1 rows and\n*         M = N + SQRE >= N columns.\n*\n*  K      (output) INTEGER\n*         Contains the dimension of the non-deflated matrix,\n*         This is the order of the related secular equation. 1 <= K <=N.\n*\n*  D      (input/output) DOUBLE PRECISION array, dimension(N)\n*         On entry D contains the singular values of the two submatrices\n*         to be combined.  On exit D contains the trailing (N-K) updated\n*         singular values (those which were deflated) sorted into\n*         increasing order.\n*\n*  Z      (output) DOUBLE PRECISION array, dimension(N)\n*         On exit Z contains the updating row vector in the secular\n*         equation.\n*\n*  ALPHA  (input) DOUBLE PRECISION\n*         Contains the diagonal element associated with the added row.\n*\n*  BETA   (input) DOUBLE PRECISION\n*         Contains the off-diagonal element associated with the added\n*         row.\n*\n*  U      (input/output) DOUBLE PRECISION array, dimension(LDU,N)\n*         On entry U contains the left singular vectors of two\n*         submatrices in the two square blocks with corners at (1,1),\n*         (NL, NL), and (NL+2, NL+2), (N,N).\n*         On exit U contains the trailing (N-K) updated left singular\n*         vectors (those which were deflated) in its last N-K columns.\n*\n*  LDU    (input) INTEGER\n*         The leading dimension of the array U.  LDU >= N.\n*\n*  VT     (input/output) DOUBLE PRECISION array, dimension(LDVT,M)\n*         On entry VT' contains the right singular vectors of two\n*         submatrices in the two square blocks with corners at (1,1),\n*         (NL+1, NL+1), and (NL+2, NL+2), (M,M).\n*         On exit VT' contains the trailing (N-K) updated right singular\n*         vectors (those which were deflated) in its last N-K columns.\n*         In case SQRE =1, the last row of VT spans the right null\n*         space.\n*\n*  LDVT   (input) INTEGER\n*         The leading dimension of the array VT.  LDVT >= M.\n*\n*  DSIGMA (output) DOUBLE PRECISION array, dimension (N)\n*         Contains a copy of the diagonal elements (K-1 singular values\n*         and one zero) in the secular equation.\n*\n*  U2     (output) DOUBLE PRECISION array, dimension(LDU2,N)\n*         Contains a copy of the first K-1 left singular vectors which\n*         will be used by DLASD3 in a matrix multiply (DGEMM) to solve\n*         for the new left singular vectors. U2 is arranged into four\n*         blocks. The first block contains a column with 1 at NL+1 and\n*         zero everywhere else; the second block contains non-zero\n*         entries only at and above NL; the third contains non-zero\n*         entries only below NL+1; and the fourth is dense.\n*\n*  LDU2   (input) INTEGER\n*         The leading dimension of the array U2.  LDU2 >= N.\n*\n*  VT2    (output) DOUBLE PRECISION array, dimension(LDVT2,N)\n*         VT2' contains a copy of the first K right singular vectors\n*         which will be used by DLASD3 in a matrix multiply (DGEMM) to\n*         solve for the new right singular vectors. VT2 is arranged into\n*         three blocks. The first block contains a row that corresponds\n*         to the special 0 diagonal element in SIGMA; the second block\n*         contains non-zeros only at and before NL +1; the third block\n*         contains non-zeros only at and after  NL +2.\n*\n*  LDVT2  (input) INTEGER\n*         The leading dimension of the array VT2.  LDVT2 >= M.\n*\n*  IDXP   (workspace) INTEGER array dimension(N)\n*         This will contain the permutation used to place deflated\n*         values of D at the end of the array. On output IDXP(2:K)\n*         points to the nondeflated D-values and IDXP(K+1:N)\n*         points to the deflated singular values.\n*\n*  IDX    (workspace) INTEGER array dimension(N)\n*         This will contain the permutation used to sort the contents of\n*         D into ascending order.\n*\n*  IDXC   (output) INTEGER array dimension(N)\n*         This will contain the permutation used to arrange the columns\n*         of the deflated U matrix into three groups:  the first group\n*         contains non-zero entries only at and above NL, the second\n*         contains non-zero entries only below NL+2, and the third is\n*         dense.\n*\n*  IDXQ   (input/output) INTEGER array dimension(N)\n*         This contains the permutation which separately sorts the two\n*         sub-problems in D into ascending order.  Note that entries in\n*         the first hlaf of this permutation must first be moved one\n*         position backward; and entries in the second half\n*         must first have NL+1 added to their values.\n*\n*  COLTYP (workspace/output) INTEGER array dimension(N)\n*         As workspace, this will contain a label which will indicate\n*         which of the following types a column in the U2 matrix or a\n*         row in the VT2 matrix is:\n*         1 : non-zero in the upper half only\n*         2 : non-zero in the lower half only\n*         3 : dense\n*         4 : deflated\n*\n*         On exit, it is an array of dimension 4, with COLTYP(I) being\n*         the dimension of the I-th type columns.\n*\n*  INFO   (output) INTEGER\n*          = 0:  successful exit.\n*          < 0:  if INFO = -i, the i-th argument had an illegal value.\n*\n\n*  Further Details\n*  ===============\n*\n*  Based on contributions by\n*     Ming Gu and Huan Ren, Computer Science Division, University of\n*     California at Berkeley, USA\n*\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 9)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 9)", argc);
  rb_nl = argv[0];
  rb_nr = argv[1];
  rb_sqre = argv[2];
  rb_d = argv[3];
  rb_alpha = argv[4];
  rb_beta = argv[5];
  rb_u = argv[6];
  rb_vt = argv[7];
  rb_idxq = argv[8];

  if (!NA_IsNArray(rb_idxq))
    rb_raise(rb_eArgError, "idxq (9th argument) must be NArray");
  if (NA_RANK(rb_idxq) != 1)
    rb_raise(rb_eArgError, "rank of idxq (9th argument) must be %d", 1);
  n = NA_SHAPE0(rb_idxq);
  if (NA_TYPE(rb_idxq) != NA_LINT)
    rb_idxq = na_change_type(rb_idxq, NA_LINT);
  idxq = NA_PTR_TYPE(rb_idxq, integer*);
  if (!NA_IsNArray(rb_u))
    rb_raise(rb_eArgError, "u (7th argument) must be NArray");
  if (NA_RANK(rb_u) != 2)
    rb_raise(rb_eArgError, "rank of u (7th argument) must be %d", 2);
  if (NA_SHAPE1(rb_u) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of u must be the same as shape 0 of idxq");
  ldu = NA_SHAPE0(rb_u);
  if (NA_TYPE(rb_u) != NA_DFLOAT)
    rb_u = na_change_type(rb_u, NA_DFLOAT);
  u = NA_PTR_TYPE(rb_u, doublereal*);
  if (!NA_IsNArray(rb_d))
    rb_raise(rb_eArgError, "d (4th argument) must be NArray");
  if (NA_RANK(rb_d) != 1)
    rb_raise(rb_eArgError, "rank of d (4th argument) must be %d", 1);
  if (NA_SHAPE0(rb_d) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of d must be the same as shape 0 of idxq");
  if (NA_TYPE(rb_d) != NA_DFLOAT)
    rb_d = na_change_type(rb_d, NA_DFLOAT);
  d = NA_PTR_TYPE(rb_d, doublereal*);
  nr = NUM2INT(rb_nr);
  beta = NUM2DBL(rb_beta);
  nl = NUM2INT(rb_nl);
  sqre = NUM2INT(rb_sqre);
  if (!NA_IsNArray(rb_vt))
    rb_raise(rb_eArgError, "vt (8th argument) must be NArray");
  if (NA_RANK(rb_vt) != 2)
    rb_raise(rb_eArgError, "rank of vt (8th argument) must be %d", 2);
  m = NA_SHAPE1(rb_vt);
  ldvt = NA_SHAPE0(rb_vt);
  if (NA_TYPE(rb_vt) != NA_DFLOAT)
    rb_vt = na_change_type(rb_vt, NA_DFLOAT);
  vt = NA_PTR_TYPE(rb_vt, doublereal*);
  alpha = NUM2DBL(rb_alpha);
  ldu2 = n;
  ldvt2 = m;
  {
    int shape[1];
    shape[0] = n;
    rb_z = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  z = NA_PTR_TYPE(rb_z, doublereal*);
  {
    int shape[1];
    shape[0] = n;
    rb_dsigma = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  dsigma = NA_PTR_TYPE(rb_dsigma, doublereal*);
  {
    int shape[2];
    shape[0] = ldu2;
    shape[1] = n;
    rb_u2 = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  u2 = NA_PTR_TYPE(rb_u2, doublereal*);
  {
    int shape[2];
    shape[0] = ldvt2;
    shape[1] = n;
    rb_vt2 = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  vt2 = NA_PTR_TYPE(rb_vt2, doublereal*);
  {
    int shape[1];
    shape[0] = n;
    rb_idxc = na_make_object(NA_LINT, 1, shape, cNArray);
  }
  idxc = NA_PTR_TYPE(rb_idxc, integer*);
  {
    int shape[1];
    shape[0] = n;
    rb_coltyp = na_make_object(NA_LINT, 1, shape, cNArray);
  }
  coltyp = NA_PTR_TYPE(rb_coltyp, integer*);
  {
    int shape[1];
    shape[0] = n;
    rb_d_out__ = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  d_out__ = NA_PTR_TYPE(rb_d_out__, doublereal*);
  MEMCPY(d_out__, d, doublereal, NA_TOTAL(rb_d));
  rb_d = rb_d_out__;
  d = d_out__;
  {
    int shape[2];
    shape[0] = ldu;
    shape[1] = n;
    rb_u_out__ = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  u_out__ = NA_PTR_TYPE(rb_u_out__, doublereal*);
  MEMCPY(u_out__, u, doublereal, NA_TOTAL(rb_u));
  rb_u = rb_u_out__;
  u = u_out__;
  {
    int shape[2];
    shape[0] = ldvt;
    shape[1] = m;
    rb_vt_out__ = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  vt_out__ = NA_PTR_TYPE(rb_vt_out__, doublereal*);
  MEMCPY(vt_out__, vt, doublereal, NA_TOTAL(rb_vt));
  rb_vt = rb_vt_out__;
  vt = vt_out__;
  {
    int shape[1];
    shape[0] = n;
    rb_idxq_out__ = na_make_object(NA_LINT, 1, shape, cNArray);
  }
  idxq_out__ = NA_PTR_TYPE(rb_idxq_out__, integer*);
  MEMCPY(idxq_out__, idxq, integer, NA_TOTAL(rb_idxq));
  rb_idxq = rb_idxq_out__;
  idxq = idxq_out__;
  idxp = ALLOC_N(integer, (n));
  idx = ALLOC_N(integer, (n));

  dlasd2_(&nl, &nr, &sqre, &k, d, z, &alpha, &beta, u, &ldu, vt, &ldvt, dsigma, u2, &ldu2, vt2, &ldvt2, idxp, idx, idxc, idxq, coltyp, &info);

  free(idxp);
  free(idx);
  rb_k = INT2NUM(k);
  rb_info = INT2NUM(info);
  return rb_ary_new3(12, rb_k, rb_z, rb_dsigma, rb_u2, rb_vt2, rb_idxc, rb_coltyp, rb_info, rb_d, rb_u, rb_vt, rb_idxq);
}

void
init_lapack_dlasd2(VALUE mLapack){
  rb_define_module_function(mLapack, "dlasd2", rb_dlasd2, -1);
}

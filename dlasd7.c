#include "rb_lapack.h"

static VALUE
rb_dlasd7(int argc, VALUE *argv, VALUE self){
  VALUE rb_icompq;
  integer icompq; 
  VALUE rb_nl;
  integer nl; 
  VALUE rb_nr;
  integer nr; 
  VALUE rb_sqre;
  integer sqre; 
  VALUE rb_d;
  doublereal *d; 
  VALUE rb_vf;
  doublereal *vf; 
  VALUE rb_vl;
  doublereal *vl; 
  VALUE rb_alpha;
  doublereal alpha; 
  VALUE rb_beta;
  doublereal beta; 
  VALUE rb_idxq;
  integer *idxq; 
  VALUE rb_k;
  integer k; 
  VALUE rb_z;
  doublereal *z; 
  VALUE rb_dsigma;
  doublereal *dsigma; 
  VALUE rb_perm;
  integer *perm; 
  VALUE rb_givptr;
  integer givptr; 
  VALUE rb_givcol;
  integer *givcol; 
  VALUE rb_givnum;
  doublereal *givnum; 
  VALUE rb_c;
  doublereal c; 
  VALUE rb_s;
  doublereal s; 
  VALUE rb_info;
  integer info; 
  VALUE rb_d_out__;
  doublereal *d_out__;
  VALUE rb_vf_out__;
  doublereal *vf_out__;
  VALUE rb_vl_out__;
  doublereal *vl_out__;
  doublereal *zw;
  doublereal *vfw;
  doublereal *vlw;
  integer *idx;
  integer *idxp;

  integer n;
  integer m;
  integer ldgcol;
  integer ldgnum;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  k, z, dsigma, perm, givptr, givcol, givnum, c, s, info, d, vf, vl = NumRu::Lapack.dlasd7( icompq, nl, nr, sqre, d, vf, vl, alpha, beta, idxq)\n    or\n  NumRu::Lapack.dlasd7  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE DLASD7( ICOMPQ, NL, NR, SQRE, K, D, Z, ZW, VF, VFW, VL, VLW, ALPHA, BETA, DSIGMA, IDX, IDXP, IDXQ, PERM, GIVPTR, GIVCOL, LDGCOL, GIVNUM, LDGNUM, C, S, INFO )\n\n*  Purpose\n*  =======\n*\n*  DLASD7 merges the two sets of singular values together into a single\n*  sorted set. Then it tries to deflate the size of the problem. There\n*  are two ways in which deflation can occur:  when two or more singular\n*  values are close together or if there is a tiny entry in the Z\n*  vector. For each such occurrence the order of the related\n*  secular equation problem is reduced by one.\n*\n*  DLASD7 is called from DLASD6.\n*\n\n*  Arguments\n*  =========\n*\n*  ICOMPQ  (input) INTEGER\n*          Specifies whether singular vectors are to be computed\n*          in compact form, as follows:\n*          = 0: Compute singular values only.\n*          = 1: Compute singular vectors of upper\n*               bidiagonal matrix in compact form.\n*\n*  NL     (input) INTEGER\n*         The row dimension of the upper block. NL >= 1.\n*\n*  NR     (input) INTEGER\n*         The row dimension of the lower block. NR >= 1.\n*\n*  SQRE   (input) INTEGER\n*         = 0: the lower block is an NR-by-NR square matrix.\n*         = 1: the lower block is an NR-by-(NR+1) rectangular matrix.\n*\n*         The bidiagonal matrix has\n*         N = NL + NR + 1 rows and\n*         M = N + SQRE >= N columns.\n*\n*  K      (output) INTEGER\n*         Contains the dimension of the non-deflated matrix, this is\n*         the order of the related secular equation. 1 <= K <=N.\n*\n*  D      (input/output) DOUBLE PRECISION array, dimension ( N )\n*         On entry D contains the singular values of the two submatrices\n*         to be combined. On exit D contains the trailing (N-K) updated\n*         singular values (those which were deflated) sorted into\n*         increasing order.\n*\n*  Z      (output) DOUBLE PRECISION array, dimension ( M )\n*         On exit Z contains the updating row vector in the secular\n*         equation.\n*\n*  ZW     (workspace) DOUBLE PRECISION array, dimension ( M )\n*         Workspace for Z.\n*\n*  VF     (input/output) DOUBLE PRECISION array, dimension ( M )\n*         On entry, VF(1:NL+1) contains the first components of all\n*         right singular vectors of the upper block; and VF(NL+2:M)\n*         contains the first components of all right singular vectors\n*         of the lower block. On exit, VF contains the first components\n*         of all right singular vectors of the bidiagonal matrix.\n*\n*  VFW    (workspace) DOUBLE PRECISION array, dimension ( M )\n*         Workspace for VF.\n*\n*  VL     (input/output) DOUBLE PRECISION array, dimension ( M )\n*         On entry, VL(1:NL+1) contains the  last components of all\n*         right singular vectors of the upper block; and VL(NL+2:M)\n*         contains the last components of all right singular vectors\n*         of the lower block. On exit, VL contains the last components\n*         of all right singular vectors of the bidiagonal matrix.\n*\n*  VLW    (workspace) DOUBLE PRECISION array, dimension ( M )\n*         Workspace for VL.\n*\n*  ALPHA  (input) DOUBLE PRECISION\n*         Contains the diagonal element associated with the added row.\n*\n*  BETA   (input) DOUBLE PRECISION\n*         Contains the off-diagonal element associated with the added\n*         row.\n*\n*  DSIGMA (output) DOUBLE PRECISION array, dimension ( N )\n*         Contains a copy of the diagonal elements (K-1 singular values\n*         and one zero) in the secular equation.\n*\n*  IDX    (workspace) INTEGER array, dimension ( N )\n*         This will contain the permutation used to sort the contents of\n*         D into ascending order.\n*\n*  IDXP   (workspace) INTEGER array, dimension ( N )\n*         This will contain the permutation used to place deflated\n*         values of D at the end of the array. On output IDXP(2:K)\n*         points to the nondeflated D-values and IDXP(K+1:N)\n*         points to the deflated singular values.\n*\n*  IDXQ   (input) INTEGER array, dimension ( N )\n*         This contains the permutation which separately sorts the two\n*         sub-problems in D into ascending order.  Note that entries in\n*         the first half of this permutation must first be moved one\n*         position backward; and entries in the second half\n*         must first have NL+1 added to their values.\n*\n*  PERM   (output) INTEGER array, dimension ( N )\n*         The permutations (from deflation and sorting) to be applied\n*         to each singular block. Not referenced if ICOMPQ = 0.\n*\n*  GIVPTR (output) INTEGER\n*         The number of Givens rotations which took place in this\n*         subproblem. Not referenced if ICOMPQ = 0.\n*\n*  GIVCOL (output) INTEGER array, dimension ( LDGCOL, 2 )\n*         Each pair of numbers indicates a pair of columns to take place\n*         in a Givens rotation. Not referenced if ICOMPQ = 0.\n*\n*  LDGCOL (input) INTEGER\n*         The leading dimension of GIVCOL, must be at least N.\n*\n*  GIVNUM (output) DOUBLE PRECISION array, dimension ( LDGNUM, 2 )\n*         Each number indicates the C or S value to be used in the\n*         corresponding Givens rotation. Not referenced if ICOMPQ = 0.\n*\n*  LDGNUM (input) INTEGER\n*         The leading dimension of GIVNUM, must be at least N.\n*\n*  C      (output) DOUBLE PRECISION\n*         C contains garbage if SQRE =0 and the C-value of a Givens\n*         rotation related to the right null space if SQRE = 1.\n*\n*  S      (output) DOUBLE PRECISION\n*         S contains garbage if SQRE =0 and the S-value of a Givens\n*         rotation related to the right null space if SQRE = 1.\n*\n*  INFO   (output) INTEGER\n*         = 0:  successful exit.\n*         < 0:  if INFO = -i, the i-th argument had an illegal value.\n*\n\n*  Further Details\n*  ===============\n*\n*  Based on contributions by\n*     Ming Gu and Huan Ren, Computer Science Division, University of\n*     California at Berkeley, USA\n*\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 10)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 10)", argc);
  rb_icompq = argv[0];
  rb_nl = argv[1];
  rb_nr = argv[2];
  rb_sqre = argv[3];
  rb_d = argv[4];
  rb_vf = argv[5];
  rb_vl = argv[6];
  rb_alpha = argv[7];
  rb_beta = argv[8];
  rb_idxq = argv[9];

  icompq = NUM2INT(rb_icompq);
  nl = NUM2INT(rb_nl);
  nr = NUM2INT(rb_nr);
  sqre = NUM2INT(rb_sqre);
  alpha = NUM2DBL(rb_alpha);
  beta = NUM2DBL(rb_beta);
  if (!NA_IsNArray(rb_d))
    rb_raise(rb_eArgError, "d (5th argument) must be NArray");
  if (NA_RANK(rb_d) != 1)
    rb_raise(rb_eArgError, "rank of d (5th argument) must be %d", 1);
  n = NA_SHAPE0(rb_d);
  if (NA_TYPE(rb_d) != NA_DFLOAT)
    rb_d = na_change_type(rb_d, NA_DFLOAT);
  d = NA_PTR_TYPE(rb_d, doublereal*);
  if (!NA_IsNArray(rb_vf))
    rb_raise(rb_eArgError, "vf (6th argument) must be NArray");
  if (NA_RANK(rb_vf) != 1)
    rb_raise(rb_eArgError, "rank of vf (6th argument) must be %d", 1);
  m = NA_SHAPE0(rb_vf);
  if (NA_TYPE(rb_vf) != NA_DFLOAT)
    rb_vf = na_change_type(rb_vf, NA_DFLOAT);
  vf = NA_PTR_TYPE(rb_vf, doublereal*);
  if (!NA_IsNArray(rb_vl))
    rb_raise(rb_eArgError, "vl (7th argument) must be NArray");
  if (NA_RANK(rb_vl) != 1)
    rb_raise(rb_eArgError, "rank of vl (7th argument) must be %d", 1);
  if (NA_SHAPE0(rb_vl) != m)
    rb_raise(rb_eRuntimeError, "shape 0 of vl must be the same as shape 0 of vf");
  if (NA_TYPE(rb_vl) != NA_DFLOAT)
    rb_vl = na_change_type(rb_vl, NA_DFLOAT);
  vl = NA_PTR_TYPE(rb_vl, doublereal*);
  if (!NA_IsNArray(rb_idxq))
    rb_raise(rb_eArgError, "idxq (10th argument) must be NArray");
  if (NA_RANK(rb_idxq) != 1)
    rb_raise(rb_eArgError, "rank of idxq (10th argument) must be %d", 1);
  if (NA_SHAPE0(rb_idxq) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of idxq must be the same as shape 0 of d");
  if (NA_TYPE(rb_idxq) != NA_LINT)
    rb_idxq = na_change_type(rb_idxq, NA_LINT);
  idxq = NA_PTR_TYPE(rb_idxq, integer*);
  {
    int shape[1];
    shape[0] = m;
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
    int shape[1];
    shape[0] = n;
    rb_perm = na_make_object(NA_LINT, 1, shape, cNArray);
  }
  perm = NA_PTR_TYPE(rb_perm, integer*);
  ldgcol = n;
  {
    int shape[2];
    shape[0] = ldgcol;
    shape[1] = 2;
    rb_givcol = na_make_object(NA_LINT, 2, shape, cNArray);
  }
  givcol = NA_PTR_TYPE(rb_givcol, integer*);
  ldgnum = n;
  {
    int shape[2];
    shape[0] = ldgnum;
    shape[1] = 2;
    rb_givnum = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  givnum = NA_PTR_TYPE(rb_givnum, doublereal*);
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
    int shape[1];
    shape[0] = m;
    rb_vf_out__ = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  vf_out__ = NA_PTR_TYPE(rb_vf_out__, doublereal*);
  MEMCPY(vf_out__, vf, doublereal, NA_TOTAL(rb_vf));
  rb_vf = rb_vf_out__;
  vf = vf_out__;
  {
    int shape[1];
    shape[0] = m;
    rb_vl_out__ = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  vl_out__ = NA_PTR_TYPE(rb_vl_out__, doublereal*);
  MEMCPY(vl_out__, vl, doublereal, NA_TOTAL(rb_vl));
  rb_vl = rb_vl_out__;
  vl = vl_out__;
  zw = ALLOC_N(doublereal, (m));
  vfw = ALLOC_N(doublereal, (m));
  vlw = ALLOC_N(doublereal, (m));
  idx = ALLOC_N(integer, (n));
  idxp = ALLOC_N(integer, (n));

  dlasd7_(&icompq, &nl, &nr, &sqre, &k, d, z, zw, vf, vfw, vl, vlw, &alpha, &beta, dsigma, idx, idxp, idxq, perm, &givptr, givcol, &ldgcol, givnum, &ldgnum, &c, &s, &info);

  free(zw);
  free(vfw);
  free(vlw);
  free(idx);
  free(idxp);
  rb_k = INT2NUM(k);
  rb_givptr = INT2NUM(givptr);
  rb_c = rb_float_new((double)c);
  rb_s = rb_float_new((double)s);
  rb_info = INT2NUM(info);
  return rb_ary_new3(13, rb_k, rb_z, rb_dsigma, rb_perm, rb_givptr, rb_givcol, rb_givnum, rb_c, rb_s, rb_info, rb_d, rb_vf, rb_vl);
}

void
init_lapack_dlasd7(VALUE mLapack){
  rb_define_module_function(mLapack, "dlasd7", rb_dlasd7, -1);
}

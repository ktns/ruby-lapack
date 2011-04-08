#include "rb_lapack.h"

extern VOID slasd2_(integer *nl, integer *nr, integer *sqre, integer *k, real *d, real *z, real *alpha, real *beta, real *u, integer *ldu, real *vt, integer *ldvt, real *dsigma, real *u2, integer *ldu2, real *vt2, integer *ldvt2, integer *idxp, integer *idx, integer *idxc, integer *idxq, integer *coltyp, integer *info);

static VALUE sHelp, sUsage;

static VALUE
rblapack_slasd2(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_nl;
  integer nl; 
  VALUE rblapack_nr;
  integer nr; 
  VALUE rblapack_sqre;
  integer sqre; 
  VALUE rblapack_d;
  real *d; 
  VALUE rblapack_alpha;
  real alpha; 
  VALUE rblapack_beta;
  real beta; 
  VALUE rblapack_u;
  real *u; 
  VALUE rblapack_vt;
  real *vt; 
  VALUE rblapack_idxq;
  integer *idxq; 
  VALUE rblapack_k;
  integer k; 
  VALUE rblapack_z;
  real *z; 
  VALUE rblapack_dsigma;
  real *dsigma; 
  VALUE rblapack_u2;
  real *u2; 
  VALUE rblapack_vt2;
  real *vt2; 
  VALUE rblapack_idxc;
  integer *idxc; 
  VALUE rblapack_coltyp;
  integer *coltyp; 
  VALUE rblapack_info;
  integer info; 
  VALUE rblapack_d_out__;
  real *d_out__;
  VALUE rblapack_u_out__;
  real *u_out__;
  VALUE rblapack_vt_out__;
  real *vt_out__;
  VALUE rblapack_idxq_out__;
  integer *idxq_out__;
  integer *idxp;
  integer *idx;

  integer n;
  integer ldu;
  integer ldvt;
  integer m;
  integer ldu2;
  integer ldvt2;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  k, z, dsigma, u2, vt2, idxc, coltyp, info, d, u, vt, idxq = NumRu::Lapack.slasd2( nl, nr, sqre, d, alpha, beta, u, vt, idxq, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE SLASD2( NL, NR, SQRE, K, D, Z, ALPHA, BETA, U, LDU, VT, LDVT, DSIGMA, U2, LDU2, VT2, LDVT2, IDXP, IDX, IDXC, IDXQ, COLTYP, INFO )\n\n*  Purpose\n*  =======\n*\n*  SLASD2 merges the two sets of singular values together into a single\n*  sorted set.  Then it tries to deflate the size of the problem.\n*  There are two ways in which deflation can occur:  when two or more\n*  singular values are close together or if there is a tiny entry in the\n*  Z vector.  For each such occurrence the order of the related secular\n*  equation problem is reduced by one.\n*\n*  SLASD2 is called from SLASD1.\n*\n\n*  Arguments\n*  =========\n*\n*  NL     (input) INTEGER\n*         The row dimension of the upper block.  NL >= 1.\n*\n*  NR     (input) INTEGER\n*         The row dimension of the lower block.  NR >= 1.\n*\n*  SQRE   (input) INTEGER\n*         = 0: the lower block is an NR-by-NR square matrix.\n*         = 1: the lower block is an NR-by-(NR+1) rectangular matrix.\n*\n*         The bidiagonal matrix has N = NL + NR + 1 rows and\n*         M = N + SQRE >= N columns.\n*\n*  K      (output) INTEGER\n*         Contains the dimension of the non-deflated matrix,\n*         This is the order of the related secular equation. 1 <= K <=N.\n*\n*  D      (input/output) REAL array, dimension (N)\n*         On entry D contains the singular values of the two submatrices\n*         to be combined.  On exit D contains the trailing (N-K) updated\n*         singular values (those which were deflated) sorted into\n*         increasing order.\n*\n*  Z      (output) REAL array, dimension (N)\n*         On exit Z contains the updating row vector in the secular\n*         equation.\n*\n*  ALPHA  (input) REAL\n*         Contains the diagonal element associated with the added row.\n*\n*  BETA   (input) REAL\n*         Contains the off-diagonal element associated with the added\n*         row.\n*\n*  U      (input/output) REAL array, dimension (LDU,N)\n*         On entry U contains the left singular vectors of two\n*         submatrices in the two square blocks with corners at (1,1),\n*         (NL, NL), and (NL+2, NL+2), (N,N).\n*         On exit U contains the trailing (N-K) updated left singular\n*         vectors (those which were deflated) in its last N-K columns.\n*\n*  LDU    (input) INTEGER\n*         The leading dimension of the array U.  LDU >= N.\n*\n*  VT     (input/output) REAL array, dimension (LDVT,M)\n*         On entry VT' contains the right singular vectors of two\n*         submatrices in the two square blocks with corners at (1,1),\n*         (NL+1, NL+1), and (NL+2, NL+2), (M,M).\n*         On exit VT' contains the trailing (N-K) updated right singular\n*         vectors (those which were deflated) in its last N-K columns.\n*         In case SQRE =1, the last row of VT spans the right null\n*         space.\n*\n*  LDVT   (input) INTEGER\n*         The leading dimension of the array VT.  LDVT >= M.\n*\n*  DSIGMA (output) REAL array, dimension (N)\n*         Contains a copy of the diagonal elements (K-1 singular values\n*         and one zero) in the secular equation.\n*\n*  U2     (output) REAL array, dimension (LDU2,N)\n*         Contains a copy of the first K-1 left singular vectors which\n*         will be used by SLASD3 in a matrix multiply (SGEMM) to solve\n*         for the new left singular vectors. U2 is arranged into four\n*         blocks. The first block contains a column with 1 at NL+1 and\n*         zero everywhere else; the second block contains non-zero\n*         entries only at and above NL; the third contains non-zero\n*         entries only below NL+1; and the fourth is dense.\n*\n*  LDU2   (input) INTEGER\n*         The leading dimension of the array U2.  LDU2 >= N.\n*\n*  VT2    (output) REAL array, dimension (LDVT2,N)\n*         VT2' contains a copy of the first K right singular vectors\n*         which will be used by SLASD3 in a matrix multiply (SGEMM) to\n*         solve for the new right singular vectors. VT2 is arranged into\n*         three blocks. The first block contains a row that corresponds\n*         to the special 0 diagonal element in SIGMA; the second block\n*         contains non-zeros only at and before NL +1; the third block\n*         contains non-zeros only at and after  NL +2.\n*\n*  LDVT2  (input) INTEGER\n*         The leading dimension of the array VT2.  LDVT2 >= M.\n*\n*  IDXP   (workspace) INTEGER array, dimension (N)\n*         This will contain the permutation used to place deflated\n*         values of D at the end of the array. On output IDXP(2:K)\n*         points to the nondeflated D-values and IDXP(K+1:N)\n*         points to the deflated singular values.\n*\n*  IDX    (workspace) INTEGER array, dimension (N)\n*         This will contain the permutation used to sort the contents of\n*         D into ascending order.\n*\n*  IDXC   (output) INTEGER array, dimension (N)\n*         This will contain the permutation used to arrange the columns\n*         of the deflated U matrix into three groups:  the first group\n*         contains non-zero entries only at and above NL, the second\n*         contains non-zero entries only below NL+2, and the third is\n*         dense.\n*\n*  IDXQ   (input/output) INTEGER array, dimension (N)\n*         This contains the permutation which separately sorts the two\n*         sub-problems in D into ascending order.  Note that entries in\n*         the first hlaf of this permutation must first be moved one\n*         position backward; and entries in the second half\n*         must first have NL+1 added to their values.\n*\n*  COLTYP (workspace/output) INTEGER array, dimension (N)\n*         As workspace, this will contain a label which will indicate\n*         which of the following types a column in the U2 matrix or a\n*         row in the VT2 matrix is:\n*         1 : non-zero in the upper half only\n*         2 : non-zero in the lower half only\n*         3 : dense\n*         4 : deflated\n*\n*         On exit, it is an array of dimension 4, with COLTYP(I) being\n*         the dimension of the I-th type columns.\n*\n*  INFO   (output) INTEGER\n*          = 0:  successful exit.\n*          < 0:  if INFO = -i, the i-th argument had an illegal value.\n*\n\n*  Further Details\n*  ===============\n*\n*  Based on contributions by\n*     Ming Gu and Huan Ren, Computer Science Division, University of\n*     California at Berkeley, USA\n*\n*  =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  k, z, dsigma, u2, vt2, idxc, coltyp, info, d, u, vt, idxq = NumRu::Lapack.slasd2( nl, nr, sqre, d, alpha, beta, u, vt, idxq, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 9)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 9)", argc);
  rblapack_nl = argv[0];
  rblapack_nr = argv[1];
  rblapack_sqre = argv[2];
  rblapack_d = argv[3];
  rblapack_alpha = argv[4];
  rblapack_beta = argv[5];
  rblapack_u = argv[6];
  rblapack_vt = argv[7];
  rblapack_idxq = argv[8];
  if (rb_options != Qnil) {
  }

  if (!NA_IsNArray(rblapack_idxq))
    rb_raise(rb_eArgError, "idxq (9th argument) must be NArray");
  if (NA_RANK(rblapack_idxq) != 1)
    rb_raise(rb_eArgError, "rank of idxq (9th argument) must be %d", 1);
  n = NA_SHAPE0(rblapack_idxq);
  if (NA_TYPE(rblapack_idxq) != NA_LINT)
    rblapack_idxq = na_change_type(rblapack_idxq, NA_LINT);
  idxq = NA_PTR_TYPE(rblapack_idxq, integer*);
  if (!NA_IsNArray(rblapack_u))
    rb_raise(rb_eArgError, "u (7th argument) must be NArray");
  if (NA_RANK(rblapack_u) != 2)
    rb_raise(rb_eArgError, "rank of u (7th argument) must be %d", 2);
  if (NA_SHAPE1(rblapack_u) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of u must be the same as shape 0 of idxq");
  ldu = NA_SHAPE0(rblapack_u);
  if (NA_TYPE(rblapack_u) != NA_SFLOAT)
    rblapack_u = na_change_type(rblapack_u, NA_SFLOAT);
  u = NA_PTR_TYPE(rblapack_u, real*);
  if (!NA_IsNArray(rblapack_d))
    rb_raise(rb_eArgError, "d (4th argument) must be NArray");
  if (NA_RANK(rblapack_d) != 1)
    rb_raise(rb_eArgError, "rank of d (4th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_d) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of d must be the same as shape 0 of idxq");
  if (NA_TYPE(rblapack_d) != NA_SFLOAT)
    rblapack_d = na_change_type(rblapack_d, NA_SFLOAT);
  d = NA_PTR_TYPE(rblapack_d, real*);
  nr = NUM2INT(rblapack_nr);
  beta = (real)NUM2DBL(rblapack_beta);
  nl = NUM2INT(rblapack_nl);
  sqre = NUM2INT(rblapack_sqre);
  if (!NA_IsNArray(rblapack_vt))
    rb_raise(rb_eArgError, "vt (8th argument) must be NArray");
  if (NA_RANK(rblapack_vt) != 2)
    rb_raise(rb_eArgError, "rank of vt (8th argument) must be %d", 2);
  m = NA_SHAPE1(rblapack_vt);
  ldvt = NA_SHAPE0(rblapack_vt);
  if (NA_TYPE(rblapack_vt) != NA_SFLOAT)
    rblapack_vt = na_change_type(rblapack_vt, NA_SFLOAT);
  vt = NA_PTR_TYPE(rblapack_vt, real*);
  alpha = (real)NUM2DBL(rblapack_alpha);
  ldu2 = n;
  ldvt2 = m;
  {
    int shape[1];
    shape[0] = n;
    rblapack_z = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  z = NA_PTR_TYPE(rblapack_z, real*);
  {
    int shape[1];
    shape[0] = n;
    rblapack_dsigma = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  dsigma = NA_PTR_TYPE(rblapack_dsigma, real*);
  {
    int shape[2];
    shape[0] = ldu2;
    shape[1] = n;
    rblapack_u2 = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  u2 = NA_PTR_TYPE(rblapack_u2, real*);
  {
    int shape[2];
    shape[0] = ldvt2;
    shape[1] = n;
    rblapack_vt2 = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  vt2 = NA_PTR_TYPE(rblapack_vt2, real*);
  {
    int shape[1];
    shape[0] = n;
    rblapack_idxc = na_make_object(NA_LINT, 1, shape, cNArray);
  }
  idxc = NA_PTR_TYPE(rblapack_idxc, integer*);
  {
    int shape[1];
    shape[0] = n;
    rblapack_coltyp = na_make_object(NA_LINT, 1, shape, cNArray);
  }
  coltyp = NA_PTR_TYPE(rblapack_coltyp, integer*);
  {
    int shape[1];
    shape[0] = n;
    rblapack_d_out__ = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  d_out__ = NA_PTR_TYPE(rblapack_d_out__, real*);
  MEMCPY(d_out__, d, real, NA_TOTAL(rblapack_d));
  rblapack_d = rblapack_d_out__;
  d = d_out__;
  {
    int shape[2];
    shape[0] = ldu;
    shape[1] = n;
    rblapack_u_out__ = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  u_out__ = NA_PTR_TYPE(rblapack_u_out__, real*);
  MEMCPY(u_out__, u, real, NA_TOTAL(rblapack_u));
  rblapack_u = rblapack_u_out__;
  u = u_out__;
  {
    int shape[2];
    shape[0] = ldvt;
    shape[1] = m;
    rblapack_vt_out__ = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  vt_out__ = NA_PTR_TYPE(rblapack_vt_out__, real*);
  MEMCPY(vt_out__, vt, real, NA_TOTAL(rblapack_vt));
  rblapack_vt = rblapack_vt_out__;
  vt = vt_out__;
  {
    int shape[1];
    shape[0] = n;
    rblapack_idxq_out__ = na_make_object(NA_LINT, 1, shape, cNArray);
  }
  idxq_out__ = NA_PTR_TYPE(rblapack_idxq_out__, integer*);
  MEMCPY(idxq_out__, idxq, integer, NA_TOTAL(rblapack_idxq));
  rblapack_idxq = rblapack_idxq_out__;
  idxq = idxq_out__;
  idxp = ALLOC_N(integer, (n));
  idx = ALLOC_N(integer, (n));

  slasd2_(&nl, &nr, &sqre, &k, d, z, &alpha, &beta, u, &ldu, vt, &ldvt, dsigma, u2, &ldu2, vt2, &ldvt2, idxp, idx, idxc, idxq, coltyp, &info);

  free(idxp);
  free(idx);
  rblapack_k = INT2NUM(k);
  rblapack_info = INT2NUM(info);
  return rb_ary_new3(12, rblapack_k, rblapack_z, rblapack_dsigma, rblapack_u2, rblapack_vt2, rblapack_idxc, rblapack_coltyp, rblapack_info, rblapack_d, rblapack_u, rblapack_vt, rblapack_idxq);
}

void
init_lapack_slasd2(VALUE mLapack){
  rb_define_module_function(mLapack, "slasd2", rblapack_slasd2, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}

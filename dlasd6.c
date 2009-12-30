#include "rb_lapack.h"

static VALUE
rb_dlasd6(int argc, VALUE *argv, VALUE self){
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
  VALUE rb_k;
  integer k; 
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
  doublereal *work;
  integer *iwork;

  integer m;
  integer n;
  integer ldgcol;
  integer ldgnum;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  idxq, perm, givptr, givcol, givnum, poles, difl, difr, z, k, c, s, info, d, vf, vl, alpha, beta = NumRu::Lapack.dlasd6( icompq, nl, nr, sqre, d, vf, vl, alpha, beta)\n    or\n  NumRu::Lapack.dlasd6  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE DLASD6( ICOMPQ, NL, NR, SQRE, D, VF, VL, ALPHA, BETA, IDXQ, PERM, GIVPTR, GIVCOL, LDGCOL, GIVNUM, LDGNUM, POLES, DIFL, DIFR, Z, K, C, S, WORK, IWORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  DLASD6 computes the SVD of an updated upper bidiagonal matrix B\n*  obtained by merging two smaller ones by appending a row. This\n*  routine is used only for the problem which requires all singular\n*  values and optionally singular vector matrices in factored form.\n*  B is an N-by-M matrix with N = NL + NR + 1 and M = N + SQRE.\n*  A related subroutine, DLASD1, handles the case in which all singular\n*  values and singular vectors of the bidiagonal matrix are desired.\n*\n*  DLASD6 computes the SVD as follows:\n*\n*                ( D1(in)  0    0     0 )\n*    B = U(in) * (   Z1'   a   Z2'    b ) * VT(in)\n*                (   0     0   D2(in) 0 )\n*\n*      = U(out) * ( D(out) 0) * VT(out)\n*\n*  where Z' = (Z1' a Z2' b) = u' VT', and u is a vector of dimension M\n*  with ALPHA and BETA in the NL+1 and NL+2 th entries and zeros\n*  elsewhere; and the entry b is empty if SQRE = 0.\n*\n*  The singular values of B can be computed using D1, D2, the first\n*  components of all the right singular vectors of the lower block, and\n*  the last components of all the right singular vectors of the upper\n*  block. These components are stored and updated in VF and VL,\n*  respectively, in DLASD6. Hence U and VT are not explicitly\n*  referenced.\n*\n*  The singular values are stored in D. The algorithm consists of two\n*  stages:\n*\n*        The first stage consists of deflating the size of the problem\n*        when there are multiple singular values or if there is a zero\n*        in the Z vector. For each such occurence the dimension of the\n*        secular equation problem is reduced by one. This stage is\n*        performed by the routine DLASD7.\n*\n*        The second stage consists of calculating the updated\n*        singular values. This is done by finding the roots of the\n*        secular equation via the routine DLASD4 (as called by DLASD8).\n*        This routine also updates VF and VL and computes the distances\n*        between the updated singular values and the old singular\n*        values.\n*\n*  DLASD6 is called from DLASDA.\n*\n\n*  Arguments\n*  =========\n*\n*  ICOMPQ (input) INTEGER\n*         Specifies whether singular vectors are to be computed in\n*         factored form:\n*         = 0: Compute singular values only.\n*         = 1: Compute singular vectors in factored form as well.\n*\n*  NL     (input) INTEGER\n*         The row dimension of the upper block.  NL >= 1.\n*\n*  NR     (input) INTEGER\n*         The row dimension of the lower block.  NR >= 1.\n*\n*  SQRE   (input) INTEGER\n*         = 0: the lower block is an NR-by-NR square matrix.\n*         = 1: the lower block is an NR-by-(NR+1) rectangular matrix.\n*\n*         The bidiagonal matrix has row dimension N = NL + NR + 1,\n*         and column dimension M = N + SQRE.\n*\n*  D      (input/output) DOUBLE PRECISION array, dimension ( NL+NR+1 ).\n*         On entry D(1:NL,1:NL) contains the singular values of the\n*         upper block, and D(NL+2:N) contains the singular values\n*         of the lower block. On exit D(1:N) contains the singular\n*         values of the modified matrix.\n*\n*  VF     (input/output) DOUBLE PRECISION array, dimension ( M )\n*         On entry, VF(1:NL+1) contains the first components of all\n*         right singular vectors of the upper block; and VF(NL+2:M)\n*         contains the first components of all right singular vectors\n*         of the lower block. On exit, VF contains the first components\n*         of all right singular vectors of the bidiagonal matrix.\n*\n*  VL     (input/output) DOUBLE PRECISION array, dimension ( M )\n*         On entry, VL(1:NL+1) contains the  last components of all\n*         right singular vectors of the upper block; and VL(NL+2:M)\n*         contains the last components of all right singular vectors of\n*         the lower block. On exit, VL contains the last components of\n*         all right singular vectors of the bidiagonal matrix.\n*\n*  ALPHA  (input/output) DOUBLE PRECISION\n*         Contains the diagonal element associated with the added row.\n*\n*  BETA   (input/output) DOUBLE PRECISION\n*         Contains the off-diagonal element associated with the added\n*         row.\n*\n*  IDXQ   (output) INTEGER array, dimension ( N )\n*         This contains the permutation which will reintegrate the\n*         subproblem just solved back into sorted order, i.e.\n*         D( IDXQ( I = 1, N ) ) will be in ascending order.\n*\n*  PERM   (output) INTEGER array, dimension ( N )\n*         The permutations (from deflation and sorting) to be applied\n*         to each block. Not referenced if ICOMPQ = 0.\n*\n*  GIVPTR (output) INTEGER\n*         The number of Givens rotations which took place in this\n*         subproblem. Not referenced if ICOMPQ = 0.\n*\n*  GIVCOL (output) INTEGER array, dimension ( LDGCOL, 2 )\n*         Each pair of numbers indicates a pair of columns to take place\n*         in a Givens rotation. Not referenced if ICOMPQ = 0.\n*\n*  LDGCOL (input) INTEGER\n*         leading dimension of GIVCOL, must be at least N.\n*\n*  GIVNUM (output) DOUBLE PRECISION array, dimension ( LDGNUM, 2 )\n*         Each number indicates the C or S value to be used in the\n*         corresponding Givens rotation. Not referenced if ICOMPQ = 0.\n*\n*  LDGNUM (input) INTEGER\n*         The leading dimension of GIVNUM and POLES, must be at least N.\n*\n*  POLES  (output) DOUBLE PRECISION array, dimension ( LDGNUM, 2 )\n*         On exit, POLES(1,*) is an array containing the new singular\n*         values obtained from solving the secular equation, and\n*         POLES(2,*) is an array containing the poles in the secular\n*         equation. Not referenced if ICOMPQ = 0.\n*\n*  DIFL   (output) DOUBLE PRECISION array, dimension ( N )\n*         On exit, DIFL(I) is the distance between I-th updated\n*         (undeflated) singular value and the I-th (undeflated) old\n*         singular value.\n*\n*  DIFR   (output) DOUBLE PRECISION array,\n*                  dimension ( LDGNUM, 2 ) if ICOMPQ = 1 and\n*                  dimension ( N ) if ICOMPQ = 0.\n*         On exit, DIFR(I, 1) is the distance between I-th updated\n*         (undeflated) singular value and the I+1-th (undeflated) old\n*         singular value.\n*\n*         If ICOMPQ = 1, DIFR(1:K,2) is an array containing the\n*         normalizing factors for the right singular vector matrix.\n*\n*         See DLASD8 for details on DIFL and DIFR.\n*\n*  Z      (output) DOUBLE PRECISION array, dimension ( M )\n*         The first elements of this array contain the components\n*         of the deflation-adjusted updating row vector.\n*\n*  K      (output) INTEGER\n*         Contains the dimension of the non-deflated matrix,\n*         This is the order of the related secular equation. 1 <= K <=N.\n*\n*  C      (output) DOUBLE PRECISION\n*         C contains garbage if SQRE =0 and the C-value of a Givens\n*         rotation related to the right null space if SQRE = 1.\n*\n*  S      (output) DOUBLE PRECISION\n*         S contains garbage if SQRE =0 and the S-value of a Givens\n*         rotation related to the right null space if SQRE = 1.\n*\n*  WORK   (workspace) DOUBLE PRECISION array, dimension ( 4 * M )\n*\n*  IWORK  (workspace) INTEGER array, dimension ( 3 * N )\n*\n*  INFO   (output) INTEGER\n*          = 0:  successful exit.\n*          < 0:  if INFO = -i, the i-th argument had an illegal value.\n*          > 0:  if INFO = 1, an singular value did not converge\n*\n\n*  Further Details\n*  ===============\n*\n*  Based on contributions by\n*     Ming Gu and Huan Ren, Computer Science Division, University of\n*     California at Berkeley, USA\n*\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 9)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 9)", argc);
  rb_icompq = argv[0];
  rb_nl = argv[1];
  rb_nr = argv[2];
  rb_sqre = argv[3];
  rb_d = argv[4];
  rb_vf = argv[5];
  rb_vl = argv[6];
  rb_alpha = argv[7];
  rb_beta = argv[8];

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
  if (NA_SHAPE0(rb_d) != (nl+nr+1))
    rb_raise(rb_eRuntimeError, "shape 0 of d must be %d", nl+nr+1);
  if (NA_TYPE(rb_d) != NA_DFLOAT)
    rb_d = na_change_type(rb_d, NA_DFLOAT);
  d = NA_PTR_TYPE(rb_d, doublereal*);
  if (!NA_IsNArray(rb_vf))
    rb_raise(rb_eArgError, "vf (6th argument) must be NArray");
  if (NA_RANK(rb_vf) != 1)
    rb_raise(rb_eArgError, "rank of vf (6th argument) must be %d", 1);
  n = nl + nr + 1;
  m = n + sqre;
  if (NA_SHAPE0(rb_vf) != m)
    rb_raise(rb_eRuntimeError, "shape 0 of vf must be m");
  if (NA_TYPE(rb_vf) != NA_DFLOAT)
    rb_vf = na_change_type(rb_vf, NA_DFLOAT);
  vf = NA_PTR_TYPE(rb_vf, doublereal*);
  if (!NA_IsNArray(rb_vl))
    rb_raise(rb_eArgError, "vl (7th argument) must be NArray");
  if (NA_RANK(rb_vl) != 1)
    rb_raise(rb_eArgError, "rank of vl (7th argument) must be %d", 1);
  if (NA_SHAPE0(rb_vl) != m)
    rb_raise(rb_eRuntimeError, "shape 0 of vl must be m");
  if (NA_TYPE(rb_vl) != NA_DFLOAT)
    rb_vl = na_change_type(rb_vl, NA_DFLOAT);
  vl = NA_PTR_TYPE(rb_vl, doublereal*);
  {
    int shape[1];
    shape[0] = DIM_LEN(n);
    rb_idxq = na_make_object(NA_LINT, 1, shape, cNArray);
  }
  idxq = NA_PTR_TYPE(rb_idxq, integer*);
  {
    int shape[1];
    shape[0] = DIM_LEN(n);
    rb_perm = na_make_object(NA_LINT, 1, shape, cNArray);
  }
  perm = NA_PTR_TYPE(rb_perm, integer*);
  ldgcol = n;
  {
    int shape[2];
    shape[0] = DIM_LEN(ldgcol);
    shape[1] = DIM_LEN(2);
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
    int shape[2];
    shape[0] = ldgnum;
    shape[1] = 2;
    rb_poles = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  poles = NA_PTR_TYPE(rb_poles, doublereal*);
  {
    int shape[1];
    shape[0] = n;
    rb_difl = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  difl = NA_PTR_TYPE(rb_difl, doublereal*);
  {
    int shape[2];
    shape[0] = icompq == 1 ? ldgnum : icompq == 0 ? n : 0;
    shape[1] = icompq == 1 ? 2 : 0;
    rb_difr = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  difr = NA_PTR_TYPE(rb_difr, doublereal*);
  {
    int shape[1];
    shape[0] = m;
    rb_z = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  z = NA_PTR_TYPE(rb_z, doublereal*);
  {
    int shape[1];
    shape[0] = nl+nr+1;
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
  work = ALLOC_N(doublereal, (4 * m));
  iwork = ALLOC_N(integer, (3 * n));

  dlasd6_(&icompq, &nl, &nr, &sqre, d, vf, vl, &alpha, &beta, idxq, perm, &givptr, givcol, &ldgcol, givnum, &ldgnum, poles, difl, difr, z, &k, &c, &s, work, iwork, &info);

  free(work);
  free(iwork);
  rb_givptr = INT2NUM(givptr);
  rb_k = INT2NUM(k);
  rb_c = rb_float_new((double)c);
  rb_s = rb_float_new((double)s);
  rb_info = INT2NUM(info);
  rb_alpha = rb_float_new((double)alpha);
  rb_beta = rb_float_new((double)beta);
  return rb_ary_new3(18, rb_idxq, rb_perm, rb_givptr, rb_givcol, rb_givnum, rb_poles, rb_difl, rb_difr, rb_z, rb_k, rb_c, rb_s, rb_info, rb_d, rb_vf, rb_vl, rb_alpha, rb_beta);
}

void
init_lapack_dlasd6(VALUE mLapack){
  rb_define_module_function(mLapack, "dlasd6", rb_dlasd6, -1);
}

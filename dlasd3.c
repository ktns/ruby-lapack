#include "rb_lapack.h"

extern VOID dlasd3_(integer *nl, integer *nr, integer *sqre, integer *k, doublereal *d, doublereal *q, integer *ldq, doublereal *dsigma, doublereal *u, integer *ldu, doublereal *u2, integer *ldu2, doublereal *vt, integer *ldvt, doublereal *vt2, integer *ldvt2, integer *idxc, integer *ctot, doublereal *z, integer *info);

static VALUE
rb_dlasd3(int argc, VALUE *argv, VALUE self){
  VALUE rb_nl;
  integer nl; 
  VALUE rb_nr;
  integer nr; 
  VALUE rb_sqre;
  integer sqre; 
  VALUE rb_dsigma;
  doublereal *dsigma; 
  VALUE rb_u2;
  doublereal *u2; 
  VALUE rb_vt2;
  doublereal *vt2; 
  VALUE rb_idxc;
  integer *idxc; 
  VALUE rb_ctot;
  integer *ctot; 
  VALUE rb_z;
  doublereal *z; 
  VALUE rb_d;
  doublereal *d; 
  VALUE rb_u;
  doublereal *u; 
  VALUE rb_vt;
  doublereal *vt; 
  VALUE rb_info;
  integer info; 
  VALUE rb_u2_out__;
  doublereal *u2_out__;
  VALUE rb_vt2_out__;
  doublereal *vt2_out__;
  doublereal *q;

  integer k;
  integer ldu2;
  integer n;
  integer ldvt2;
  integer ldu;
  integer ldvt;
  integer m;
  integer ldq;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  d, u, vt, info, u2, vt2 = NumRu::Lapack.dlasd3( nl, nr, sqre, dsigma, u2, vt2, idxc, ctot, z)\n    or\n  NumRu::Lapack.dlasd3  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE DLASD3( NL, NR, SQRE, K, D, Q, LDQ, DSIGMA, U, LDU, U2, LDU2, VT, LDVT, VT2, LDVT2, IDXC, CTOT, Z, INFO )\n\n*  Purpose\n*  =======\n*\n*  DLASD3 finds all the square roots of the roots of the secular\n*  equation, as defined by the values in D and Z.  It makes the\n*  appropriate calls to DLASD4 and then updates the singular\n*  vectors by matrix multiplication.\n*\n*  This code makes very mild assumptions about floating point\n*  arithmetic. It will work on machines with a guard digit in\n*  add/subtract, or on those binary machines without guard digits\n*  which subtract like the Cray XMP, Cray YMP, Cray C 90, or Cray 2.\n*  It could conceivably fail on hexadecimal or decimal machines\n*  without guard digits, but we know of none.\n*\n*  DLASD3 is called from DLASD1.\n*\n\n*  Arguments\n*  =========\n*\n*  NL     (input) INTEGER\n*         The row dimension of the upper block.  NL >= 1.\n*\n*  NR     (input) INTEGER\n*         The row dimension of the lower block.  NR >= 1.\n*\n*  SQRE   (input) INTEGER\n*         = 0: the lower block is an NR-by-NR square matrix.\n*         = 1: the lower block is an NR-by-(NR+1) rectangular matrix.\n*\n*         The bidiagonal matrix has N = NL + NR + 1 rows and\n*         M = N + SQRE >= N columns.\n*\n*  K      (input) INTEGER\n*         The size of the secular equation, 1 =< K = < N.\n*\n*  D      (output) DOUBLE PRECISION array, dimension(K)\n*         On exit the square roots of the roots of the secular equation,\n*         in ascending order.\n*\n*  Q      (workspace) DOUBLE PRECISION array,\n*                     dimension at least (LDQ,K).\n*\n*  LDQ    (input) INTEGER\n*         The leading dimension of the array Q.  LDQ >= K.\n*\n*  DSIGMA (input) DOUBLE PRECISION array, dimension(K)\n*         The first K elements of this array contain the old roots\n*         of the deflated updating problem.  These are the poles\n*         of the secular equation.\n*\n*  U      (output) DOUBLE PRECISION array, dimension (LDU, N)\n*         The last N - K columns of this matrix contain the deflated\n*         left singular vectors.\n*\n*  LDU    (input) INTEGER\n*         The leading dimension of the array U.  LDU >= N.\n*\n*  U2     (input/output) DOUBLE PRECISION array, dimension (LDU2, N)\n*         The first K columns of this matrix contain the non-deflated\n*         left singular vectors for the split problem.\n*\n*  LDU2   (input) INTEGER\n*         The leading dimension of the array U2.  LDU2 >= N.\n*\n*  VT     (output) DOUBLE PRECISION array, dimension (LDVT, M)\n*         The last M - K columns of VT' contain the deflated\n*         right singular vectors.\n*\n*  LDVT   (input) INTEGER\n*         The leading dimension of the array VT.  LDVT >= N.\n*\n*  VT2    (input/output) DOUBLE PRECISION array, dimension (LDVT2, N)\n*         The first K columns of VT2' contain the non-deflated\n*         right singular vectors for the split problem.\n*\n*  LDVT2  (input) INTEGER\n*         The leading dimension of the array VT2.  LDVT2 >= N.\n*\n*  IDXC   (input) INTEGER array, dimension ( N )\n*         The permutation used to arrange the columns of U (and rows of\n*         VT) into three groups:  the first group contains non-zero\n*         entries only at and above (or before) NL +1; the second\n*         contains non-zero entries only at and below (or after) NL+2;\n*         and the third is dense. The first column of U and the row of\n*         VT are treated separately, however.\n*\n*         The rows of the singular vectors found by DLASD4\n*         must be likewise permuted before the matrix multiplies can\n*         take place.\n*\n*  CTOT   (input) INTEGER array, dimension ( 4 )\n*         A count of the total number of the various types of columns\n*         in U (or rows in VT), as described in IDXC. The fourth column\n*         type is any column which has been deflated.\n*\n*  Z      (input) DOUBLE PRECISION array, dimension (K)\n*         The first K elements of this array contain the components\n*         of the deflation-adjusted updating row vector.\n*\n*  INFO   (output) INTEGER\n*         = 0:  successful exit.\n*         < 0:  if INFO = -i, the i-th argument had an illegal value.\n*         > 0:  if INFO = 1, a singular value did not converge\n*\n\n*  Further Details\n*  ===============\n*\n*  Based on contributions by\n*     Ming Gu and Huan Ren, Computer Science Division, University of\n*     California at Berkeley, USA\n*\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 9)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 9)", argc);
  rb_nl = argv[0];
  rb_nr = argv[1];
  rb_sqre = argv[2];
  rb_dsigma = argv[3];
  rb_u2 = argv[4];
  rb_vt2 = argv[5];
  rb_idxc = argv[6];
  rb_ctot = argv[7];
  rb_z = argv[8];

  if (!NA_IsNArray(rb_ctot))
    rb_raise(rb_eArgError, "ctot (8th argument) must be NArray");
  if (NA_RANK(rb_ctot) != 1)
    rb_raise(rb_eArgError, "rank of ctot (8th argument) must be %d", 1);
  if (NA_SHAPE0(rb_ctot) != (4))
    rb_raise(rb_eRuntimeError, "shape 0 of ctot must be %d", 4);
  if (NA_TYPE(rb_ctot) != NA_LINT)
    rb_ctot = na_change_type(rb_ctot, NA_LINT);
  ctot = NA_PTR_TYPE(rb_ctot, integer*);
  if (!NA_IsNArray(rb_vt2))
    rb_raise(rb_eArgError, "vt2 (6th argument) must be NArray");
  if (NA_RANK(rb_vt2) != 2)
    rb_raise(rb_eArgError, "rank of vt2 (6th argument) must be %d", 2);
  n = NA_SHAPE1(rb_vt2);
  if (n != (nl + nr + 1))
    rb_raise(rb_eRuntimeError, "shape 1 of vt2 must be %d", nl + nr + 1);
  ldvt2 = NA_SHAPE0(rb_vt2);
  if (ldvt2 != (n))
    rb_raise(rb_eRuntimeError, "shape 0 of vt2 must be %d", n);
  n = ldvt2;
  if (NA_TYPE(rb_vt2) != NA_DFLOAT)
    rb_vt2 = na_change_type(rb_vt2, NA_DFLOAT);
  vt2 = NA_PTR_TYPE(rb_vt2, doublereal*);
  nl = NUM2INT(rb_nl);
  if (!NA_IsNArray(rb_idxc))
    rb_raise(rb_eArgError, "idxc (7th argument) must be NArray");
  if (NA_RANK(rb_idxc) != 1)
    rb_raise(rb_eArgError, "rank of idxc (7th argument) must be %d", 1);
  if (NA_SHAPE0(rb_idxc) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of idxc must be n");
  if (NA_TYPE(rb_idxc) != NA_LINT)
    rb_idxc = na_change_type(rb_idxc, NA_LINT);
  idxc = NA_PTR_TYPE(rb_idxc, integer*);
  if (!NA_IsNArray(rb_u2))
    rb_raise(rb_eArgError, "u2 (5th argument) must be NArray");
  if (NA_RANK(rb_u2) != 2)
    rb_raise(rb_eArgError, "rank of u2 (5th argument) must be %d", 2);
  if (NA_SHAPE1(rb_u2) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of u2 must be n");
  ldu2 = NA_SHAPE0(rb_u2);
  if (ldu2 != (n))
    rb_raise(rb_eRuntimeError, "shape 0 of u2 must be %d", n);
  n = ldu2;
  if (NA_TYPE(rb_u2) != NA_DFLOAT)
    rb_u2 = na_change_type(rb_u2, NA_DFLOAT);
  u2 = NA_PTR_TYPE(rb_u2, doublereal*);
  if (!NA_IsNArray(rb_dsigma))
    rb_raise(rb_eArgError, "dsigma (4th argument) must be NArray");
  if (NA_RANK(rb_dsigma) != 1)
    rb_raise(rb_eArgError, "rank of dsigma (4th argument) must be %d", 1);
  k = NA_SHAPE0(rb_dsigma);
  if (NA_TYPE(rb_dsigma) != NA_DFLOAT)
    rb_dsigma = na_change_type(rb_dsigma, NA_DFLOAT);
  dsigma = NA_PTR_TYPE(rb_dsigma, doublereal*);
  sqre = NUM2INT(rb_sqre);
  if (!NA_IsNArray(rb_z))
    rb_raise(rb_eArgError, "z (9th argument) must be NArray");
  if (NA_RANK(rb_z) != 1)
    rb_raise(rb_eArgError, "rank of z (9th argument) must be %d", 1);
  if (NA_SHAPE0(rb_z) != k)
    rb_raise(rb_eRuntimeError, "shape 0 of z must be the same as shape 0 of dsigma");
  if (NA_TYPE(rb_z) != NA_DFLOAT)
    rb_z = na_change_type(rb_z, NA_DFLOAT);
  z = NA_PTR_TYPE(rb_z, doublereal*);
  nr = NUM2INT(rb_nr);
  ldu2 = n;
  ldu = n;
  ldvt = n;
  ldvt2 = n;
  ldq = k;
  m = n + sqre;
  {
    int shape[1];
    shape[0] = k;
    rb_d = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  d = NA_PTR_TYPE(rb_d, doublereal*);
  {
    int shape[2];
    shape[0] = ldu;
    shape[1] = n;
    rb_u = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  u = NA_PTR_TYPE(rb_u, doublereal*);
  {
    int shape[2];
    shape[0] = ldvt;
    shape[1] = m;
    rb_vt = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  vt = NA_PTR_TYPE(rb_vt, doublereal*);
  {
    int shape[2];
    shape[0] = ldu2;
    shape[1] = n;
    rb_u2_out__ = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  u2_out__ = NA_PTR_TYPE(rb_u2_out__, doublereal*);
  MEMCPY(u2_out__, u2, doublereal, NA_TOTAL(rb_u2));
  rb_u2 = rb_u2_out__;
  u2 = u2_out__;
  {
    int shape[2];
    shape[0] = ldvt2;
    shape[1] = n;
    rb_vt2_out__ = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  vt2_out__ = NA_PTR_TYPE(rb_vt2_out__, doublereal*);
  MEMCPY(vt2_out__, vt2, doublereal, NA_TOTAL(rb_vt2));
  rb_vt2 = rb_vt2_out__;
  vt2 = vt2_out__;
  q = ALLOC_N(doublereal, (ldq)*(k));

  dlasd3_(&nl, &nr, &sqre, &k, d, q, &ldq, dsigma, u, &ldu, u2, &ldu2, vt, &ldvt, vt2, &ldvt2, idxc, ctot, z, &info);

  free(q);
  rb_info = INT2NUM(info);
  return rb_ary_new3(6, rb_d, rb_u, rb_vt, rb_info, rb_u2, rb_vt2);
}

void
init_lapack_dlasd3(VALUE mLapack){
  rb_define_module_function(mLapack, "dlasd3", rb_dlasd3, -1);
}

#include "rb_lapack.h"

static VALUE
rb_dbdsdc(int argc, VALUE *argv, VALUE self){
  VALUE rb_uplo;
  char uplo; 
  VALUE rb_compq;
  char compq; 
  VALUE rb_d;
  doublereal *d; 
  VALUE rb_e;
  doublereal *e; 
  VALUE rb_ldu;
  integer ldu; 
  VALUE rb_ldvt;
  integer ldvt; 
  VALUE rb_u;
  doublereal *u; 
  VALUE rb_vt;
  doublereal *vt; 
  VALUE rb_q;
  doublereal *q; 
  VALUE rb_iq;
  integer *iq; 
  VALUE rb_info;
  integer info; 
  VALUE rb_d_out__;
  doublereal *d_out__;
  VALUE rb_e_out__;
  doublereal *e_out__;
  doublereal *work;
  integer *iwork;

  integer n;
  integer smlsiz;
  integer c__9;
  integer c__0;
  integer ldq;
  integer ldiq;
  integer lwork;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  u, vt, q, iq, info, d, e = NumRu::Lapack.dbdsdc( uplo, compq, d, e, ldu, ldvt)\n    or\n  NumRu::Lapack.dbdsdc  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE DBDSDC( UPLO, COMPQ, N, D, E, U, LDU, VT, LDVT, Q, IQ, WORK, IWORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  DBDSDC computes the singular value decomposition (SVD) of a real\n*  N-by-N (upper or lower) bidiagonal matrix B:  B = U * S * VT,\n*  using a divide and conquer method, where S is a diagonal matrix\n*  with non-negative diagonal elements (the singular values of B), and\n*  U and VT are orthogonal matrices of left and right singular vectors,\n*  respectively. DBDSDC can be used to compute all singular values,\n*  and optionally, singular vectors or singular vectors in compact form.\n*\n*  This code makes very mild assumptions about floating point\n*  arithmetic. It will work on machines with a guard digit in\n*  add/subtract, or on those binary machines without guard digits\n*  which subtract like the Cray X-MP, Cray Y-MP, Cray C-90, or Cray-2.\n*  It could conceivably fail on hexadecimal or decimal machines\n*  without guard digits, but we know of none.  See DLASD3 for details.\n*\n*  The code currently calls DLASDQ if singular values only are desired.\n*  However, it can be slightly modified to compute singular values\n*  using the divide and conquer method.\n*\n\n*  Arguments\n*  =========\n*\n*  UPLO    (input) CHARACTER*1\n*          = 'U':  B is upper bidiagonal.\n*          = 'L':  B is lower bidiagonal.\n*\n*  COMPQ   (input) CHARACTER*1\n*          Specifies whether singular vectors are to be computed\n*          as follows:\n*          = 'N':  Compute singular values only;\n*          = 'P':  Compute singular values and compute singular\n*                  vectors in compact form;\n*          = 'I':  Compute singular values and singular vectors.\n*\n*  N       (input) INTEGER\n*          The order of the matrix B.  N >= 0.\n*\n*  D       (input/output) DOUBLE PRECISION array, dimension (N)\n*          On entry, the n diagonal elements of the bidiagonal matrix B.\n*          On exit, if INFO=0, the singular values of B.\n*\n*  E       (input/output) DOUBLE PRECISION array, dimension (N-1)\n*          On entry, the elements of E contain the offdiagonal\n*          elements of the bidiagonal matrix whose SVD is desired.\n*          On exit, E has been destroyed.\n*\n*  U       (output) DOUBLE PRECISION array, dimension (LDU,N)\n*          If  COMPQ = 'I', then:\n*             On exit, if INFO = 0, U contains the left singular vectors\n*             of the bidiagonal matrix.\n*          For other values of COMPQ, U is not referenced.\n*\n*  LDU     (input) INTEGER\n*          The leading dimension of the array U.  LDU >= 1.\n*          If singular vectors are desired, then LDU >= max( 1, N ).\n*\n*  VT      (output) DOUBLE PRECISION array, dimension (LDVT,N)\n*          If  COMPQ = 'I', then:\n*             On exit, if INFO = 0, VT' contains the right singular\n*             vectors of the bidiagonal matrix.\n*          For other values of COMPQ, VT is not referenced.\n*\n*  LDVT    (input) INTEGER\n*          The leading dimension of the array VT.  LDVT >= 1.\n*          If singular vectors are desired, then LDVT >= max( 1, N ).\n*\n*  Q       (output) DOUBLE PRECISION array, dimension (LDQ)\n*          If  COMPQ = 'P', then:\n*             On exit, if INFO = 0, Q and IQ contain the left\n*             and right singular vectors in a compact form,\n*             requiring O(N log N) space instead of 2*N**2.\n*             In particular, Q contains all the DOUBLE PRECISION data in\n*             LDQ >= N*(11 + 2*SMLSIZ + 8*INT(LOG_2(N/(SMLSIZ+1))))\n*             words of memory, where SMLSIZ is returned by ILAENV and\n*             is equal to the maximum size of the subproblems at the\n*             bottom of the computation tree (usually about 25).\n*          For other values of COMPQ, Q is not referenced.\n*\n*  IQ      (output) INTEGER array, dimension (LDIQ)\n*          If  COMPQ = 'P', then:\n*             On exit, if INFO = 0, Q and IQ contain the left\n*             and right singular vectors in a compact form,\n*             requiring O(N log N) space instead of 2*N**2.\n*             In particular, IQ contains all INTEGER data in\n*             LDIQ >= N*(3 + 3*INT(LOG_2(N/(SMLSIZ+1))))\n*             words of memory, where SMLSIZ is returned by ILAENV and\n*             is equal to the maximum size of the subproblems at the\n*             bottom of the computation tree (usually about 25).\n*          For other values of COMPQ, IQ is not referenced.\n*\n*  WORK    (workspace) DOUBLE PRECISION array, dimension (MAX(1,LWORK))\n*          If COMPQ = 'N' then LWORK >= (4 * N).\n*          If COMPQ = 'P' then LWORK >= (6 * N).\n*          If COMPQ = 'I' then LWORK >= (3 * N**2 + 4 * N).\n*\n*  IWORK   (workspace) INTEGER array, dimension (8*N)\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit.\n*          < 0:  if INFO = -i, the i-th argument had an illegal value.\n*          > 0:  The algorithm failed to compute an singular value.\n*                The update process of divide and conquer failed.\n*\n\n*  Further Details\n*  ===============\n*\n*  Based on contributions by\n*     Ming Gu and Huan Ren, Computer Science Division, University of\n*     California at Berkeley, USA\n*\n*  =====================================================================\n*  Changed dimension statement in comment describing E from (N) to\n*  (N-1).  Sven, 17 Feb 05.\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 6)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 6)", argc);
  rb_uplo = argv[0];
  rb_compq = argv[1];
  rb_d = argv[2];
  rb_e = argv[3];
  rb_ldu = argv[4];
  rb_ldvt = argv[5];

  uplo = StringValueCStr(rb_uplo)[0];
  compq = StringValueCStr(rb_compq)[0];
  ldu = NUM2INT(rb_ldu);
  ldvt = NUM2INT(rb_ldvt);
  if (!NA_IsNArray(rb_d))
    rb_raise(rb_eArgError, "d (3th argument) must be NArray");
  if (NA_RANK(rb_d) != 1)
    rb_raise(rb_eArgError, "rank of d (3th argument) must be %d", 1);
  n = NA_SHAPE0(rb_d);
  if (NA_TYPE(rb_d) != NA_DFLOAT)
    rb_d = na_change_type(rb_d, NA_DFLOAT);
  d = NA_PTR_TYPE(rb_d, doublereal*);
  if (!NA_IsNArray(rb_e))
    rb_raise(rb_eArgError, "e (4th argument) must be NArray");
  if (NA_RANK(rb_e) != 1)
    rb_raise(rb_eArgError, "rank of e (4th argument) must be %d", 1);
  if (NA_SHAPE0(rb_e) != (n-1))
    rb_raise(rb_eRuntimeError, "shape 0 of e must be %d", n-1);
  if (NA_TYPE(rb_e) != NA_DFLOAT)
    rb_e = na_change_type(rb_e, NA_DFLOAT);
  e = NA_PTR_TYPE(rb_e, doublereal*);
  ldu = lsame_(&compq,"I") ? MAX(1,n) : 0;
  {
    int shape[2];
    shape[0] = lsame_(&compq,"I") ? ldu : 0;
    shape[1] = lsame_(&compq,"I") ? n : 0;
    rb_u = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  u = NA_PTR_TYPE(rb_u, doublereal*);
  ldvt = lsame_(&compq,"I") ? MAX(1,n) : 0;
  {
    int shape[2];
    shape[0] = lsame_(&compq,"I") ? ldvt : 0;
    shape[1] = lsame_(&compq,"I") ? n : 0;
    rb_vt = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  vt = NA_PTR_TYPE(rb_vt, doublereal*);
  c__0 = 0;
  c__9 = 9;
  smlsiz = ilaenv_(&c__9, "DBDSDC", " ", &c__0, &c__0, &c__0, &c__0, (ftnlen)6, (ftnlen)1);
  ldq = lsame_(&compq,"P") ? n*(11+2*smlsiz+8*(int)(log(((double)n)/(smlsiz+1))/log(2.0))) : 0;
  {
    int shape[1];
    shape[0] = lsame_(&compq,"I") ? ldq : 0;
    rb_q = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  q = NA_PTR_TYPE(rb_q, doublereal*);
  ldiq = lsame_(&compq,"P") ? n*(3+3*(int)(log(((double)n)/(smlsiz+1))/log(2.0))) : 0;
  {
    int shape[1];
    shape[0] = lsame_(&compq,"I") ? ldiq : 0;
    rb_iq = na_make_object(NA_LINT, 1, shape, cNArray);
  }
  iq = NA_PTR_TYPE(rb_iq, integer*);
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
    shape[0] = n-1;
    rb_e_out__ = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  e_out__ = NA_PTR_TYPE(rb_e_out__, doublereal*);
  MEMCPY(e_out__, e, doublereal, NA_TOTAL(rb_e));
  rb_e = rb_e_out__;
  e = e_out__;
  lwork = lsame_(&compq,"N") ? 4*n : lsame_(&compq,"P") ? 6*n : lsame_(&compq,"I") ? 3*n*n+4*n : 0;
  work = ALLOC_N(doublereal, (MAX(1,lwork)));
  iwork = ALLOC_N(integer, (8*n));

  dbdsdc_(&uplo, &compq, &n, d, e, u, &ldu, vt, &ldvt, q, iq, work, iwork, &info);

  free(work);
  free(iwork);
  rb_info = INT2NUM(info);
  return rb_ary_new3(7, rb_u, rb_vt, rb_q, rb_iq, rb_info, rb_d, rb_e);
}

void
init_lapack_dbdsdc(VALUE mLapack){
  rb_define_module_function(mLapack, "dbdsdc", rb_dbdsdc, -1);
}

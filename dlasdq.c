#include "rb_lapack.h"

extern VOID dlasdq_(char *uplo, integer *sqre, integer *n, integer *ncvt, integer *nru, integer *ncc, doublereal *d, doublereal *e, doublereal *vt, integer *ldvt, doublereal *u, integer *ldu, doublereal *c, integer *ldc, doublereal *work, integer *info);

static VALUE
rb_dlasdq(int argc, VALUE *argv, VALUE self){
  VALUE rb_uplo;
  char uplo; 
  VALUE rb_sqre;
  integer sqre; 
  VALUE rb_nru;
  integer nru; 
  VALUE rb_d;
  doublereal *d; 
  VALUE rb_e;
  doublereal *e; 
  VALUE rb_vt;
  doublereal *vt; 
  VALUE rb_u;
  doublereal *u; 
  VALUE rb_c;
  doublereal *c; 
  VALUE rb_info;
  integer info; 
  VALUE rb_d_out__;
  doublereal *d_out__;
  VALUE rb_e_out__;
  doublereal *e_out__;
  VALUE rb_vt_out__;
  doublereal *vt_out__;
  VALUE rb_u_out__;
  doublereal *u_out__;
  VALUE rb_c_out__;
  doublereal *c_out__;
  doublereal *work;

  integer n;
  integer ldvt;
  integer ncvt;
  integer ldu;
  integer ldc;
  integer ncc;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  info, d, e, vt, u, c = NumRu::Lapack.dlasdq( uplo, sqre, nru, d, e, vt, u, c)\n    or\n  NumRu::Lapack.dlasdq  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE DLASDQ( UPLO, SQRE, N, NCVT, NRU, NCC, D, E, VT, LDVT, U, LDU, C, LDC, WORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  DLASDQ computes the singular value decomposition (SVD) of a real\n*  (upper or lower) bidiagonal matrix with diagonal D and offdiagonal\n*  E, accumulating the transformations if desired. Letting B denote\n*  the input bidiagonal matrix, the algorithm computes orthogonal\n*  matrices Q and P such that B = Q * S * P' (P' denotes the transpose\n*  of P). The singular values S are overwritten on D.\n*\n*  The input matrix U  is changed to U  * Q  if desired.\n*  The input matrix VT is changed to P' * VT if desired.\n*  The input matrix C  is changed to Q' * C  if desired.\n*\n*  See \"Computing  Small Singular Values of Bidiagonal Matrices With\n*  Guaranteed High Relative Accuracy,\" by J. Demmel and W. Kahan,\n*  LAPACK Working Note #3, for a detailed description of the algorithm.\n*\n\n*  Arguments\n*  =========\n*\n*  UPLO  (input) CHARACTER*1\n*        On entry, UPLO specifies whether the input bidiagonal matrix\n*        is upper or lower bidiagonal, and wether it is square are\n*        not.\n*           UPLO = 'U' or 'u'   B is upper bidiagonal.\n*           UPLO = 'L' or 'l'   B is lower bidiagonal.\n*\n*  SQRE  (input) INTEGER\n*        = 0: then the input matrix is N-by-N.\n*        = 1: then the input matrix is N-by-(N+1) if UPLU = 'U' and\n*             (N+1)-by-N if UPLU = 'L'.\n*\n*        The bidiagonal matrix has\n*        N = NL + NR + 1 rows and\n*        M = N + SQRE >= N columns.\n*\n*  N     (input) INTEGER\n*        On entry, N specifies the number of rows and columns\n*        in the matrix. N must be at least 0.\n*\n*  NCVT  (input) INTEGER\n*        On entry, NCVT specifies the number of columns of\n*        the matrix VT. NCVT must be at least 0.\n*\n*  NRU   (input) INTEGER\n*        On entry, NRU specifies the number of rows of\n*        the matrix U. NRU must be at least 0.\n*\n*  NCC   (input) INTEGER\n*        On entry, NCC specifies the number of columns of\n*        the matrix C. NCC must be at least 0.\n*\n*  D     (input/output) DOUBLE PRECISION array, dimension (N)\n*        On entry, D contains the diagonal entries of the\n*        bidiagonal matrix whose SVD is desired. On normal exit,\n*        D contains the singular values in ascending order.\n*\n*  E     (input/output) DOUBLE PRECISION array.\n*        dimension is (N-1) if SQRE = 0 and N if SQRE = 1.\n*        On entry, the entries of E contain the offdiagonal entries\n*        of the bidiagonal matrix whose SVD is desired. On normal\n*        exit, E will contain 0. If the algorithm does not converge,\n*        D and E will contain the diagonal and superdiagonal entries\n*        of a bidiagonal matrix orthogonally equivalent to the one\n*        given as input.\n*\n*  VT    (input/output) DOUBLE PRECISION array, dimension (LDVT, NCVT)\n*        On entry, contains a matrix which on exit has been\n*        premultiplied by P', dimension N-by-NCVT if SQRE = 0\n*        and (N+1)-by-NCVT if SQRE = 1 (not referenced if NCVT=0).\n*\n*  LDVT  (input) INTEGER\n*        On entry, LDVT specifies the leading dimension of VT as\n*        declared in the calling (sub) program. LDVT must be at\n*        least 1. If NCVT is nonzero LDVT must also be at least N.\n*\n*  U     (input/output) DOUBLE PRECISION array, dimension (LDU, N)\n*        On entry, contains a  matrix which on exit has been\n*        postmultiplied by Q, dimension NRU-by-N if SQRE = 0\n*        and NRU-by-(N+1) if SQRE = 1 (not referenced if NRU=0).\n*\n*  LDU   (input) INTEGER\n*        On entry, LDU  specifies the leading dimension of U as\n*        declared in the calling (sub) program. LDU must be at\n*        least max( 1, NRU ) .\n*\n*  C     (input/output) DOUBLE PRECISION array, dimension (LDC, NCC)\n*        On entry, contains an N-by-NCC matrix which on exit\n*        has been premultiplied by Q'  dimension N-by-NCC if SQRE = 0\n*        and (N+1)-by-NCC if SQRE = 1 (not referenced if NCC=0).\n*\n*  LDC   (input) INTEGER\n*        On entry, LDC  specifies the leading dimension of C as\n*        declared in the calling (sub) program. LDC must be at\n*        least 1. If NCC is nonzero, LDC must also be at least N.\n*\n*  WORK  (workspace) DOUBLE PRECISION array, dimension (4*N)\n*        Workspace. Only referenced if one of NCVT, NRU, or NCC is\n*        nonzero, and if N is at least 2.\n*\n*  INFO  (output) INTEGER\n*        On exit, a value of 0 indicates a successful exit.\n*        If INFO < 0, argument number -INFO is illegal.\n*        If INFO > 0, the algorithm did not converge, and INFO\n*        specifies how many superdiagonals did not converge.\n*\n\n*  Further Details\n*  ===============\n*\n*  Based on contributions by\n*     Ming Gu and Huan Ren, Computer Science Division, University of\n*     California at Berkeley, USA\n*\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 8)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 8)", argc);
  rb_uplo = argv[0];
  rb_sqre = argv[1];
  rb_nru = argv[2];
  rb_d = argv[3];
  rb_e = argv[4];
  rb_vt = argv[5];
  rb_u = argv[6];
  rb_c = argv[7];

  if (!NA_IsNArray(rb_c))
    rb_raise(rb_eArgError, "c (8th argument) must be NArray");
  if (NA_RANK(rb_c) != 2)
    rb_raise(rb_eArgError, "rank of c (8th argument) must be %d", 2);
  ncc = NA_SHAPE1(rb_c);
  ldc = NA_SHAPE0(rb_c);
  if (NA_TYPE(rb_c) != NA_DFLOAT)
    rb_c = na_change_type(rb_c, NA_DFLOAT);
  c = NA_PTR_TYPE(rb_c, doublereal*);
  if (!NA_IsNArray(rb_d))
    rb_raise(rb_eArgError, "d (4th argument) must be NArray");
  if (NA_RANK(rb_d) != 1)
    rb_raise(rb_eArgError, "rank of d (4th argument) must be %d", 1);
  n = NA_SHAPE0(rb_d);
  if (NA_TYPE(rb_d) != NA_DFLOAT)
    rb_d = na_change_type(rb_d, NA_DFLOAT);
  d = NA_PTR_TYPE(rb_d, doublereal*);
  uplo = StringValueCStr(rb_uplo)[0];
  nru = NUM2INT(rb_nru);
  if (!NA_IsNArray(rb_u))
    rb_raise(rb_eArgError, "u (7th argument) must be NArray");
  if (NA_RANK(rb_u) != 2)
    rb_raise(rb_eArgError, "rank of u (7th argument) must be %d", 2);
  if (NA_SHAPE1(rb_u) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of u must be the same as shape 0 of d");
  ldu = NA_SHAPE0(rb_u);
  if (NA_TYPE(rb_u) != NA_DFLOAT)
    rb_u = na_change_type(rb_u, NA_DFLOAT);
  u = NA_PTR_TYPE(rb_u, doublereal*);
  if (!NA_IsNArray(rb_vt))
    rb_raise(rb_eArgError, "vt (6th argument) must be NArray");
  if (NA_RANK(rb_vt) != 2)
    rb_raise(rb_eArgError, "rank of vt (6th argument) must be %d", 2);
  ncvt = NA_SHAPE1(rb_vt);
  ldvt = NA_SHAPE0(rb_vt);
  if (NA_TYPE(rb_vt) != NA_DFLOAT)
    rb_vt = na_change_type(rb_vt, NA_DFLOAT);
  vt = NA_PTR_TYPE(rb_vt, doublereal*);
  sqre = NUM2INT(rb_sqre);
  if (!NA_IsNArray(rb_e))
    rb_raise(rb_eArgError, "e (5th argument) must be NArray");
  if (NA_RANK(rb_e) != 1)
    rb_raise(rb_eArgError, "rank of e (5th argument) must be %d", 1);
  if (NA_SHAPE0(rb_e) != (sqre==0 ? n-1 : sqre==1 ? n : 0))
    rb_raise(rb_eRuntimeError, "shape 0 of e must be %d", sqre==0 ? n-1 : sqre==1 ? n : 0);
  if (NA_TYPE(rb_e) != NA_DFLOAT)
    rb_e = na_change_type(rb_e, NA_DFLOAT);
  e = NA_PTR_TYPE(rb_e, doublereal*);
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
    shape[0] = sqre==0 ? n-1 : sqre==1 ? n : 0;
    rb_e_out__ = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  e_out__ = NA_PTR_TYPE(rb_e_out__, doublereal*);
  MEMCPY(e_out__, e, doublereal, NA_TOTAL(rb_e));
  rb_e = rb_e_out__;
  e = e_out__;
  {
    int shape[2];
    shape[0] = ldvt;
    shape[1] = ncvt;
    rb_vt_out__ = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  vt_out__ = NA_PTR_TYPE(rb_vt_out__, doublereal*);
  MEMCPY(vt_out__, vt, doublereal, NA_TOTAL(rb_vt));
  rb_vt = rb_vt_out__;
  vt = vt_out__;
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
    shape[0] = ldc;
    shape[1] = ncc;
    rb_c_out__ = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  c_out__ = NA_PTR_TYPE(rb_c_out__, doublereal*);
  MEMCPY(c_out__, c, doublereal, NA_TOTAL(rb_c));
  rb_c = rb_c_out__;
  c = c_out__;
  work = ALLOC_N(doublereal, (4*n));

  dlasdq_(&uplo, &sqre, &n, &ncvt, &nru, &ncc, d, e, vt, &ldvt, u, &ldu, c, &ldc, work, &info);

  free(work);
  rb_info = INT2NUM(info);
  return rb_ary_new3(6, rb_info, rb_d, rb_e, rb_vt, rb_u, rb_c);
}

void
init_lapack_dlasdq(VALUE mLapack){
  rb_define_module_function(mLapack, "dlasdq", rb_dlasdq, -1);
}

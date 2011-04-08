#include "rb_lapack.h"

extern VOID slasdq_(char *uplo, integer *sqre, integer *n, integer *ncvt, integer *nru, integer *ncc, real *d, real *e, real *vt, integer *ldvt, real *u, integer *ldu, real *c, integer *ldc, real *work, integer *info);

static VALUE sHelp, sUsage;

static VALUE
rblapack_slasdq(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_uplo;
  char uplo; 
  VALUE rblapack_sqre;
  integer sqre; 
  VALUE rblapack_nru;
  integer nru; 
  VALUE rblapack_d;
  real *d; 
  VALUE rblapack_e;
  real *e; 
  VALUE rblapack_vt;
  real *vt; 
  VALUE rblapack_u;
  real *u; 
  VALUE rblapack_c;
  real *c; 
  VALUE rblapack_info;
  integer info; 
  VALUE rblapack_d_out__;
  real *d_out__;
  VALUE rblapack_e_out__;
  real *e_out__;
  VALUE rblapack_vt_out__;
  real *vt_out__;
  VALUE rblapack_u_out__;
  real *u_out__;
  VALUE rblapack_c_out__;
  real *c_out__;
  real *work;

  integer n;
  integer ldvt;
  integer ncvt;
  integer ldu;
  integer ldc;
  integer ncc;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  info, d, e, vt, u, c = NumRu::Lapack.slasdq( uplo, sqre, nru, d, e, vt, u, c, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE SLASDQ( UPLO, SQRE, N, NCVT, NRU, NCC, D, E, VT, LDVT, U, LDU, C, LDC, WORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  SLASDQ computes the singular value decomposition (SVD) of a real\n*  (upper or lower) bidiagonal matrix with diagonal D and offdiagonal\n*  E, accumulating the transformations if desired. Letting B denote\n*  the input bidiagonal matrix, the algorithm computes orthogonal\n*  matrices Q and P such that B = Q * S * P' (P' denotes the transpose\n*  of P). The singular values S are overwritten on D.\n*\n*  The input matrix U  is changed to U  * Q  if desired.\n*  The input matrix VT is changed to P' * VT if desired.\n*  The input matrix C  is changed to Q' * C  if desired.\n*\n*  See \"Computing  Small Singular Values of Bidiagonal Matrices With\n*  Guaranteed High Relative Accuracy,\" by J. Demmel and W. Kahan,\n*  LAPACK Working Note #3, for a detailed description of the algorithm.\n*\n\n*  Arguments\n*  =========\n*\n*  UPLO  (input) CHARACTER*1\n*        On entry, UPLO specifies whether the input bidiagonal matrix\n*        is upper or lower bidiagonal, and wether it is square are\n*        not.\n*           UPLO = 'U' or 'u'   B is upper bidiagonal.\n*           UPLO = 'L' or 'l'   B is lower bidiagonal.\n*\n*  SQRE  (input) INTEGER\n*        = 0: then the input matrix is N-by-N.\n*        = 1: then the input matrix is N-by-(N+1) if UPLU = 'U' and\n*             (N+1)-by-N if UPLU = 'L'.\n*\n*        The bidiagonal matrix has\n*        N = NL + NR + 1 rows and\n*        M = N + SQRE >= N columns.\n*\n*  N     (input) INTEGER\n*        On entry, N specifies the number of rows and columns\n*        in the matrix. N must be at least 0.\n*\n*  NCVT  (input) INTEGER\n*        On entry, NCVT specifies the number of columns of\n*        the matrix VT. NCVT must be at least 0.\n*\n*  NRU   (input) INTEGER\n*        On entry, NRU specifies the number of rows of\n*        the matrix U. NRU must be at least 0.\n*\n*  NCC   (input) INTEGER\n*        On entry, NCC specifies the number of columns of\n*        the matrix C. NCC must be at least 0.\n*\n*  D     (input/output) REAL array, dimension (N)\n*        On entry, D contains the diagonal entries of the\n*        bidiagonal matrix whose SVD is desired. On normal exit,\n*        D contains the singular values in ascending order.\n*\n*  E     (input/output) REAL array.\n*        dimension is (N-1) if SQRE = 0 and N if SQRE = 1.\n*        On entry, the entries of E contain the offdiagonal entries\n*        of the bidiagonal matrix whose SVD is desired. On normal\n*        exit, E will contain 0. If the algorithm does not converge,\n*        D and E will contain the diagonal and superdiagonal entries\n*        of a bidiagonal matrix orthogonally equivalent to the one\n*        given as input.\n*\n*  VT    (input/output) REAL array, dimension (LDVT, NCVT)\n*        On entry, contains a matrix which on exit has been\n*        premultiplied by P', dimension N-by-NCVT if SQRE = 0\n*        and (N+1)-by-NCVT if SQRE = 1 (not referenced if NCVT=0).\n*\n*  LDVT  (input) INTEGER\n*        On entry, LDVT specifies the leading dimension of VT as\n*        declared in the calling (sub) program. LDVT must be at\n*        least 1. If NCVT is nonzero LDVT must also be at least N.\n*\n*  U     (input/output) REAL array, dimension (LDU, N)\n*        On entry, contains a  matrix which on exit has been\n*        postmultiplied by Q, dimension NRU-by-N if SQRE = 0\n*        and NRU-by-(N+1) if SQRE = 1 (not referenced if NRU=0).\n*\n*  LDU   (input) INTEGER\n*        On entry, LDU  specifies the leading dimension of U as\n*        declared in the calling (sub) program. LDU must be at\n*        least max( 1, NRU ) .\n*\n*  C     (input/output) REAL array, dimension (LDC, NCC)\n*        On entry, contains an N-by-NCC matrix which on exit\n*        has been premultiplied by Q'  dimension N-by-NCC if SQRE = 0\n*        and (N+1)-by-NCC if SQRE = 1 (not referenced if NCC=0).\n*\n*  LDC   (input) INTEGER\n*        On entry, LDC  specifies the leading dimension of C as\n*        declared in the calling (sub) program. LDC must be at\n*        least 1. If NCC is nonzero, LDC must also be at least N.\n*\n*  WORK  (workspace) REAL array, dimension (4*N)\n*        Workspace. Only referenced if one of NCVT, NRU, or NCC is\n*        nonzero, and if N is at least 2.\n*\n*  INFO  (output) INTEGER\n*        On exit, a value of 0 indicates a successful exit.\n*        If INFO < 0, argument number -INFO is illegal.\n*        If INFO > 0, the algorithm did not converge, and INFO\n*        specifies how many superdiagonals did not converge.\n*\n\n*  Further Details\n*  ===============\n*\n*  Based on contributions by\n*     Ming Gu and Huan Ren, Computer Science Division, University of\n*     California at Berkeley, USA\n*\n*  =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  info, d, e, vt, u, c = NumRu::Lapack.slasdq( uplo, sqre, nru, d, e, vt, u, c, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 8)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 8)", argc);
  rblapack_uplo = argv[0];
  rblapack_sqre = argv[1];
  rblapack_nru = argv[2];
  rblapack_d = argv[3];
  rblapack_e = argv[4];
  rblapack_vt = argv[5];
  rblapack_u = argv[6];
  rblapack_c = argv[7];
  if (rb_options != Qnil) {
  }

  if (!NA_IsNArray(rblapack_c))
    rb_raise(rb_eArgError, "c (8th argument) must be NArray");
  if (NA_RANK(rblapack_c) != 2)
    rb_raise(rb_eArgError, "rank of c (8th argument) must be %d", 2);
  ncc = NA_SHAPE1(rblapack_c);
  ldc = NA_SHAPE0(rblapack_c);
  if (NA_TYPE(rblapack_c) != NA_SFLOAT)
    rblapack_c = na_change_type(rblapack_c, NA_SFLOAT);
  c = NA_PTR_TYPE(rblapack_c, real*);
  if (!NA_IsNArray(rblapack_d))
    rb_raise(rb_eArgError, "d (4th argument) must be NArray");
  if (NA_RANK(rblapack_d) != 1)
    rb_raise(rb_eArgError, "rank of d (4th argument) must be %d", 1);
  n = NA_SHAPE0(rblapack_d);
  if (NA_TYPE(rblapack_d) != NA_SFLOAT)
    rblapack_d = na_change_type(rblapack_d, NA_SFLOAT);
  d = NA_PTR_TYPE(rblapack_d, real*);
  uplo = StringValueCStr(rblapack_uplo)[0];
  nru = NUM2INT(rblapack_nru);
  if (!NA_IsNArray(rblapack_u))
    rb_raise(rb_eArgError, "u (7th argument) must be NArray");
  if (NA_RANK(rblapack_u) != 2)
    rb_raise(rb_eArgError, "rank of u (7th argument) must be %d", 2);
  if (NA_SHAPE1(rblapack_u) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of u must be the same as shape 0 of d");
  ldu = NA_SHAPE0(rblapack_u);
  if (NA_TYPE(rblapack_u) != NA_SFLOAT)
    rblapack_u = na_change_type(rblapack_u, NA_SFLOAT);
  u = NA_PTR_TYPE(rblapack_u, real*);
  if (!NA_IsNArray(rblapack_vt))
    rb_raise(rb_eArgError, "vt (6th argument) must be NArray");
  if (NA_RANK(rblapack_vt) != 2)
    rb_raise(rb_eArgError, "rank of vt (6th argument) must be %d", 2);
  ncvt = NA_SHAPE1(rblapack_vt);
  ldvt = NA_SHAPE0(rblapack_vt);
  if (NA_TYPE(rblapack_vt) != NA_SFLOAT)
    rblapack_vt = na_change_type(rblapack_vt, NA_SFLOAT);
  vt = NA_PTR_TYPE(rblapack_vt, real*);
  sqre = NUM2INT(rblapack_sqre);
  if (!NA_IsNArray(rblapack_e))
    rb_raise(rb_eArgError, "e (5th argument) must be NArray");
  if (NA_RANK(rblapack_e) != 1)
    rb_raise(rb_eArgError, "rank of e (5th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_e) != (sqre==0 ? n-1 : sqre==1 ? n : 0))
    rb_raise(rb_eRuntimeError, "shape 0 of e must be %d", sqre==0 ? n-1 : sqre==1 ? n : 0);
  if (NA_TYPE(rblapack_e) != NA_SFLOAT)
    rblapack_e = na_change_type(rblapack_e, NA_SFLOAT);
  e = NA_PTR_TYPE(rblapack_e, real*);
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
    int shape[1];
    shape[0] = sqre==0 ? n-1 : sqre==1 ? n : 0;
    rblapack_e_out__ = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  e_out__ = NA_PTR_TYPE(rblapack_e_out__, real*);
  MEMCPY(e_out__, e, real, NA_TOTAL(rblapack_e));
  rblapack_e = rblapack_e_out__;
  e = e_out__;
  {
    int shape[2];
    shape[0] = ldvt;
    shape[1] = ncvt;
    rblapack_vt_out__ = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  vt_out__ = NA_PTR_TYPE(rblapack_vt_out__, real*);
  MEMCPY(vt_out__, vt, real, NA_TOTAL(rblapack_vt));
  rblapack_vt = rblapack_vt_out__;
  vt = vt_out__;
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
    shape[0] = ldc;
    shape[1] = ncc;
    rblapack_c_out__ = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  c_out__ = NA_PTR_TYPE(rblapack_c_out__, real*);
  MEMCPY(c_out__, c, real, NA_TOTAL(rblapack_c));
  rblapack_c = rblapack_c_out__;
  c = c_out__;
  work = ALLOC_N(real, (4*n));

  slasdq_(&uplo, &sqre, &n, &ncvt, &nru, &ncc, d, e, vt, &ldvt, u, &ldu, c, &ldc, work, &info);

  free(work);
  rblapack_info = INT2NUM(info);
  return rb_ary_new3(6, rblapack_info, rblapack_d, rblapack_e, rblapack_vt, rblapack_u, rblapack_c);
}

void
init_lapack_slasdq(VALUE mLapack){
  rb_define_module_function(mLapack, "slasdq", rblapack_slasdq, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}

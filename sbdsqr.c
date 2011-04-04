#include "rb_lapack.h"

extern VOID sbdsqr_(char *uplo, integer *n, integer *ncvt, integer *nru, integer *ncc, real *d, real *e, real *vt, integer *ldvt, real *u, integer *ldu, real *c, integer *ldc, real *work, integer *info);

static VALUE
rb_sbdsqr(int argc, VALUE *argv, VALUE self){
  VALUE rb_uplo;
  char uplo; 
  VALUE rb_nru;
  integer nru; 
  VALUE rb_d;
  real *d; 
  VALUE rb_e;
  real *e; 
  VALUE rb_vt;
  real *vt; 
  VALUE rb_u;
  real *u; 
  VALUE rb_c;
  real *c; 
  VALUE rb_info;
  integer info; 
  VALUE rb_d_out__;
  real *d_out__;
  VALUE rb_e_out__;
  real *e_out__;
  VALUE rb_vt_out__;
  real *vt_out__;
  VALUE rb_u_out__;
  real *u_out__;
  VALUE rb_c_out__;
  real *c_out__;
  real *work;

  integer n;
  integer ldvt;
  integer ncvt;
  integer ldu;
  integer ldc;
  integer ncc;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  info, d, e, vt, u, c = NumRu::Lapack.sbdsqr( uplo, nru, d, e, vt, u, c)\n    or\n  NumRu::Lapack.sbdsqr  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE SBDSQR( UPLO, N, NCVT, NRU, NCC, D, E, VT, LDVT, U, LDU, C, LDC, WORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  SBDSQR computes the singular values and, optionally, the right and/or\n*  left singular vectors from the singular value decomposition (SVD) of\n*  a real N-by-N (upper or lower) bidiagonal matrix B using the implicit\n*  zero-shift QR algorithm.  The SVD of B has the form\n*  \n*     B = Q * S * P**T\n*  \n*  where S is the diagonal matrix of singular values, Q is an orthogonal\n*  matrix of left singular vectors, and P is an orthogonal matrix of\n*  right singular vectors.  If left singular vectors are requested, this\n*  subroutine actually returns U*Q instead of Q, and, if right singular\n*  vectors are requested, this subroutine returns P**T*VT instead of\n*  P**T, for given real input matrices U and VT.  When U and VT are the\n*  orthogonal matrices that reduce a general matrix A to bidiagonal\n*  form:  A = U*B*VT, as computed by SGEBRD, then\n* \n*     A = (U*Q) * S * (P**T*VT)\n* \n*  is the SVD of A.  Optionally, the subroutine may also compute Q**T*C\n*  for a given real input matrix C.\n*\n*  See \"Computing  Small Singular Values of Bidiagonal Matrices With\n*  Guaranteed High Relative Accuracy,\" by J. Demmel and W. Kahan,\n*  LAPACK Working Note #3 (or SIAM J. Sci. Statist. Comput. vol. 11,\n*  no. 5, pp. 873-912, Sept 1990) and\n*  \"Accurate singular values and differential qd algorithms,\" by\n*  B. Parlett and V. Fernando, Technical Report CPAM-554, Mathematics\n*  Department, University of California at Berkeley, July 1992\n*  for a detailed description of the algorithm.\n*\n\n*  Arguments\n*  =========\n*\n*  UPLO    (input) CHARACTER*1\n*          = 'U':  B is upper bidiagonal;\n*          = 'L':  B is lower bidiagonal.\n*\n*  N       (input) INTEGER\n*          The order of the matrix B.  N >= 0.\n*\n*  NCVT    (input) INTEGER\n*          The number of columns of the matrix VT. NCVT >= 0.\n*\n*  NRU     (input) INTEGER\n*          The number of rows of the matrix U. NRU >= 0.\n*\n*  NCC     (input) INTEGER\n*          The number of columns of the matrix C. NCC >= 0.\n*\n*  D       (input/output) REAL array, dimension (N)\n*          On entry, the n diagonal elements of the bidiagonal matrix B.\n*          On exit, if INFO=0, the singular values of B in decreasing\n*          order.\n*\n*  E       (input/output) REAL array, dimension (N-1)\n*          On entry, the N-1 offdiagonal elements of the bidiagonal\n*          matrix B.\n*          On exit, if INFO = 0, E is destroyed; if INFO > 0, D and E\n*          will contain the diagonal and superdiagonal elements of a\n*          bidiagonal matrix orthogonally equivalent to the one given\n*          as input.\n*\n*  VT      (input/output) REAL array, dimension (LDVT, NCVT)\n*          On entry, an N-by-NCVT matrix VT.\n*          On exit, VT is overwritten by P**T * VT.\n*          Not referenced if NCVT = 0.\n*\n*  LDVT    (input) INTEGER\n*          The leading dimension of the array VT.\n*          LDVT >= max(1,N) if NCVT > 0; LDVT >= 1 if NCVT = 0.\n*\n*  U       (input/output) REAL array, dimension (LDU, N)\n*          On entry, an NRU-by-N matrix U.\n*          On exit, U is overwritten by U * Q.\n*          Not referenced if NRU = 0.\n*\n*  LDU     (input) INTEGER\n*          The leading dimension of the array U.  LDU >= max(1,NRU).\n*\n*  C       (input/output) REAL array, dimension (LDC, NCC)\n*          On entry, an N-by-NCC matrix C.\n*          On exit, C is overwritten by Q**T * C.\n*          Not referenced if NCC = 0.\n*\n*  LDC     (input) INTEGER\n*          The leading dimension of the array C.\n*          LDC >= max(1,N) if NCC > 0; LDC >=1 if NCC = 0.\n*\n*  WORK    (workspace) REAL array, dimension (4*N)\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit\n*          < 0:  If INFO = -i, the i-th argument had an illegal value\n*          > 0:\n*             if NCVT = NRU = NCC = 0,\n*                = 1, a split was marked by a positive value in E\n*                = 2, current block of Z not diagonalized after 30*N\n*                     iterations (in inner while loop)\n*                = 3, termination criterion of outer while loop not met \n*                     (program created more than N unreduced blocks)\n*             else NCVT = NRU = NCC = 0,\n*                   the algorithm did not converge; D and E contain the\n*                   elements of a bidiagonal matrix which is orthogonally\n*                   similar to the input matrix B;  if INFO = i, i\n*                   elements of E have not converged to zero.\n*\n*  Internal Parameters\n*  ===================\n*\n*  TOLMUL  REAL, default = max(10,min(100,EPS**(-1/8)))\n*          TOLMUL controls the convergence criterion of the QR loop.\n*          If it is positive, TOLMUL*EPS is the desired relative\n*             precision in the computed singular values.\n*          If it is negative, abs(TOLMUL*EPS*sigma_max) is the\n*             desired absolute accuracy in the computed singular\n*             values (corresponds to relative accuracy\n*             abs(TOLMUL*EPS) in the largest singular value.\n*          abs(TOLMUL) should be between 1 and 1/EPS, and preferably\n*             between 10 (for fast convergence) and .1/EPS\n*             (for there to be some accuracy in the results).\n*          Default is to lose at either one eighth or 2 of the\n*             available decimal digits in each computed singular value\n*             (whichever is smaller).\n*\n*  MAXITR  INTEGER, default = 6\n*          MAXITR controls the maximum number of passes of the\n*          algorithm through its inner loop. The algorithms stops\n*          (and so fails to converge) if the number of passes\n*          through the inner loop exceeds MAXITR*N**2.\n*\n\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 7)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 7)", argc);
  rb_uplo = argv[0];
  rb_nru = argv[1];
  rb_d = argv[2];
  rb_e = argv[3];
  rb_vt = argv[4];
  rb_u = argv[5];
  rb_c = argv[6];

  if (!NA_IsNArray(rb_c))
    rb_raise(rb_eArgError, "c (7th argument) must be NArray");
  if (NA_RANK(rb_c) != 2)
    rb_raise(rb_eArgError, "rank of c (7th argument) must be %d", 2);
  ncc = NA_SHAPE1(rb_c);
  ldc = NA_SHAPE0(rb_c);
  if (NA_TYPE(rb_c) != NA_SFLOAT)
    rb_c = na_change_type(rb_c, NA_SFLOAT);
  c = NA_PTR_TYPE(rb_c, real*);
  if (!NA_IsNArray(rb_d))
    rb_raise(rb_eArgError, "d (3th argument) must be NArray");
  if (NA_RANK(rb_d) != 1)
    rb_raise(rb_eArgError, "rank of d (3th argument) must be %d", 1);
  n = NA_SHAPE0(rb_d);
  if (NA_TYPE(rb_d) != NA_SFLOAT)
    rb_d = na_change_type(rb_d, NA_SFLOAT);
  d = NA_PTR_TYPE(rb_d, real*);
  uplo = StringValueCStr(rb_uplo)[0];
  if (!NA_IsNArray(rb_u))
    rb_raise(rb_eArgError, "u (6th argument) must be NArray");
  if (NA_RANK(rb_u) != 2)
    rb_raise(rb_eArgError, "rank of u (6th argument) must be %d", 2);
  if (NA_SHAPE1(rb_u) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of u must be the same as shape 0 of d");
  ldu = NA_SHAPE0(rb_u);
  if (NA_TYPE(rb_u) != NA_SFLOAT)
    rb_u = na_change_type(rb_u, NA_SFLOAT);
  u = NA_PTR_TYPE(rb_u, real*);
  if (!NA_IsNArray(rb_vt))
    rb_raise(rb_eArgError, "vt (5th argument) must be NArray");
  if (NA_RANK(rb_vt) != 2)
    rb_raise(rb_eArgError, "rank of vt (5th argument) must be %d", 2);
  ncvt = NA_SHAPE1(rb_vt);
  ldvt = NA_SHAPE0(rb_vt);
  if (NA_TYPE(rb_vt) != NA_SFLOAT)
    rb_vt = na_change_type(rb_vt, NA_SFLOAT);
  vt = NA_PTR_TYPE(rb_vt, real*);
  nru = NUM2INT(rb_nru);
  if (!NA_IsNArray(rb_e))
    rb_raise(rb_eArgError, "e (4th argument) must be NArray");
  if (NA_RANK(rb_e) != 1)
    rb_raise(rb_eArgError, "rank of e (4th argument) must be %d", 1);
  if (NA_SHAPE0(rb_e) != (n-1))
    rb_raise(rb_eRuntimeError, "shape 0 of e must be %d", n-1);
  if (NA_TYPE(rb_e) != NA_SFLOAT)
    rb_e = na_change_type(rb_e, NA_SFLOAT);
  e = NA_PTR_TYPE(rb_e, real*);
  {
    int shape[1];
    shape[0] = n;
    rb_d_out__ = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  d_out__ = NA_PTR_TYPE(rb_d_out__, real*);
  MEMCPY(d_out__, d, real, NA_TOTAL(rb_d));
  rb_d = rb_d_out__;
  d = d_out__;
  {
    int shape[1];
    shape[0] = n-1;
    rb_e_out__ = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  e_out__ = NA_PTR_TYPE(rb_e_out__, real*);
  MEMCPY(e_out__, e, real, NA_TOTAL(rb_e));
  rb_e = rb_e_out__;
  e = e_out__;
  {
    int shape[2];
    shape[0] = ldvt;
    shape[1] = ncvt;
    rb_vt_out__ = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  vt_out__ = NA_PTR_TYPE(rb_vt_out__, real*);
  MEMCPY(vt_out__, vt, real, NA_TOTAL(rb_vt));
  rb_vt = rb_vt_out__;
  vt = vt_out__;
  {
    int shape[2];
    shape[0] = ldu;
    shape[1] = n;
    rb_u_out__ = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  u_out__ = NA_PTR_TYPE(rb_u_out__, real*);
  MEMCPY(u_out__, u, real, NA_TOTAL(rb_u));
  rb_u = rb_u_out__;
  u = u_out__;
  {
    int shape[2];
    shape[0] = ldc;
    shape[1] = ncc;
    rb_c_out__ = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  c_out__ = NA_PTR_TYPE(rb_c_out__, real*);
  MEMCPY(c_out__, c, real, NA_TOTAL(rb_c));
  rb_c = rb_c_out__;
  c = c_out__;
  work = ALLOC_N(real, (4*n));

  sbdsqr_(&uplo, &n, &ncvt, &nru, &ncc, d, e, vt, &ldvt, u, &ldu, c, &ldc, work, &info);

  free(work);
  rb_info = INT2NUM(info);
  return rb_ary_new3(6, rb_info, rb_d, rb_e, rb_vt, rb_u, rb_c);
}

void
init_lapack_sbdsqr(VALUE mLapack){
  rb_define_module_function(mLapack, "sbdsqr", rb_sbdsqr, -1);
}

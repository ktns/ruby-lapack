#include "rb_lapack.h"

extern VOID cbdsqr_(char *uplo, integer *n, integer *ncvt, integer *nru, integer *ncc, real *d, real *e, complex *vt, integer *ldvt, complex *u, integer *ldu, complex *c, integer *ldc, real *rwork, integer *info);

static VALUE sHelp, sUsage;

static VALUE
rblapack_cbdsqr(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_uplo;
  char uplo; 
  VALUE rblapack_nru;
  integer nru; 
  VALUE rblapack_d;
  real *d; 
  VALUE rblapack_e;
  real *e; 
  VALUE rblapack_vt;
  complex *vt; 
  VALUE rblapack_u;
  complex *u; 
  VALUE rblapack_c;
  complex *c; 
  VALUE rblapack_info;
  integer info; 
  VALUE rblapack_d_out__;
  real *d_out__;
  VALUE rblapack_e_out__;
  real *e_out__;
  VALUE rblapack_vt_out__;
  complex *vt_out__;
  VALUE rblapack_u_out__;
  complex *u_out__;
  VALUE rblapack_c_out__;
  complex *c_out__;
  real *rwork;

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
      printf("%s\n", "USAGE:\n  info, d, e, vt, u, c = NumRu::Lapack.cbdsqr( uplo, nru, d, e, vt, u, c, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE CBDSQR( UPLO, N, NCVT, NRU, NCC, D, E, VT, LDVT, U, LDU, C, LDC, RWORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  CBDSQR computes the singular values and, optionally, the right and/or\n*  left singular vectors from the singular value decomposition (SVD) of\n*  a real N-by-N (upper or lower) bidiagonal matrix B using the implicit\n*  zero-shift QR algorithm.  The SVD of B has the form\n*  \n*     B = Q * S * P**H\n*  \n*  where S is the diagonal matrix of singular values, Q is an orthogonal\n*  matrix of left singular vectors, and P is an orthogonal matrix of\n*  right singular vectors.  If left singular vectors are requested, this\n*  subroutine actually returns U*Q instead of Q, and, if right singular\n*  vectors are requested, this subroutine returns P**H*VT instead of\n*  P**H, for given complex input matrices U and VT.  When U and VT are\n*  the unitary matrices that reduce a general matrix A to bidiagonal\n*  form: A = U*B*VT, as computed by CGEBRD, then\n*  \n*     A = (U*Q) * S * (P**H*VT)\n*  \n*  is the SVD of A.  Optionally, the subroutine may also compute Q**H*C\n*  for a given complex input matrix C.\n*\n*  See \"Computing  Small Singular Values of Bidiagonal Matrices With\n*  Guaranteed High Relative Accuracy,\" by J. Demmel and W. Kahan,\n*  LAPACK Working Note #3 (or SIAM J. Sci. Statist. Comput. vol. 11,\n*  no. 5, pp. 873-912, Sept 1990) and\n*  \"Accurate singular values and differential qd algorithms,\" by\n*  B. Parlett and V. Fernando, Technical Report CPAM-554, Mathematics\n*  Department, University of California at Berkeley, July 1992\n*  for a detailed description of the algorithm.\n*\n\n*  Arguments\n*  =========\n*\n*  UPLO    (input) CHARACTER*1\n*          = 'U':  B is upper bidiagonal;\n*          = 'L':  B is lower bidiagonal.\n*\n*  N       (input) INTEGER\n*          The order of the matrix B.  N >= 0.\n*\n*  NCVT    (input) INTEGER\n*          The number of columns of the matrix VT. NCVT >= 0.\n*\n*  NRU     (input) INTEGER\n*          The number of rows of the matrix U. NRU >= 0.\n*\n*  NCC     (input) INTEGER\n*          The number of columns of the matrix C. NCC >= 0.\n*\n*  D       (input/output) REAL array, dimension (N)\n*          On entry, the n diagonal elements of the bidiagonal matrix B.\n*          On exit, if INFO=0, the singular values of B in decreasing\n*          order.\n*\n*  E       (input/output) REAL array, dimension (N-1)\n*          On entry, the N-1 offdiagonal elements of the bidiagonal\n*          matrix B.\n*          On exit, if INFO = 0, E is destroyed; if INFO > 0, D and E\n*          will contain the diagonal and superdiagonal elements of a\n*          bidiagonal matrix orthogonally equivalent to the one given\n*          as input.\n*\n*  VT      (input/output) COMPLEX array, dimension (LDVT, NCVT)\n*          On entry, an N-by-NCVT matrix VT.\n*          On exit, VT is overwritten by P**H * VT.\n*          Not referenced if NCVT = 0.\n*\n*  LDVT    (input) INTEGER\n*          The leading dimension of the array VT.\n*          LDVT >= max(1,N) if NCVT > 0; LDVT >= 1 if NCVT = 0.\n*\n*  U       (input/output) COMPLEX array, dimension (LDU, N)\n*          On entry, an NRU-by-N matrix U.\n*          On exit, U is overwritten by U * Q.\n*          Not referenced if NRU = 0.\n*\n*  LDU     (input) INTEGER\n*          The leading dimension of the array U.  LDU >= max(1,NRU).\n*\n*  C       (input/output) COMPLEX array, dimension (LDC, NCC)\n*          On entry, an N-by-NCC matrix C.\n*          On exit, C is overwritten by Q**H * C.\n*          Not referenced if NCC = 0.\n*\n*  LDC     (input) INTEGER\n*          The leading dimension of the array C.\n*          LDC >= max(1,N) if NCC > 0; LDC >=1 if NCC = 0.\n*\n*  RWORK   (workspace) REAL array, dimension (2*N) \n*          if NCVT = NRU = NCC = 0, (max(1, 4*N-4)) otherwise\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit\n*          < 0:  If INFO = -i, the i-th argument had an illegal value\n*          > 0:  the algorithm did not converge; D and E contain the\n*                elements of a bidiagonal matrix which is orthogonally\n*                similar to the input matrix B;  if INFO = i, i\n*                elements of E have not converged to zero.\n*\n*  Internal Parameters\n*  ===================\n*\n*  TOLMUL  REAL, default = max(10,min(100,EPS**(-1/8)))\n*          TOLMUL controls the convergence criterion of the QR loop.\n*          If it is positive, TOLMUL*EPS is the desired relative\n*             precision in the computed singular values.\n*          If it is negative, abs(TOLMUL*EPS*sigma_max) is the\n*             desired absolute accuracy in the computed singular\n*             values (corresponds to relative accuracy\n*             abs(TOLMUL*EPS) in the largest singular value.\n*          abs(TOLMUL) should be between 1 and 1/EPS, and preferably\n*             between 10 (for fast convergence) and .1/EPS\n*             (for there to be some accuracy in the results).\n*          Default is to lose at either one eighth or 2 of the\n*             available decimal digits in each computed singular value\n*             (whichever is smaller).\n*\n*  MAXITR  INTEGER, default = 6\n*          MAXITR controls the maximum number of passes of the\n*          algorithm through its inner loop. The algorithms stops\n*          (and so fails to converge) if the number of passes\n*          through the inner loop exceeds MAXITR*N**2.\n*\n\n*  =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  info, d, e, vt, u, c = NumRu::Lapack.cbdsqr( uplo, nru, d, e, vt, u, c, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 7)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 7)", argc);
  rblapack_uplo = argv[0];
  rblapack_nru = argv[1];
  rblapack_d = argv[2];
  rblapack_e = argv[3];
  rblapack_vt = argv[4];
  rblapack_u = argv[5];
  rblapack_c = argv[6];
  if (rb_options != Qnil) {
  }

  if (!NA_IsNArray(rblapack_c))
    rb_raise(rb_eArgError, "c (7th argument) must be NArray");
  if (NA_RANK(rblapack_c) != 2)
    rb_raise(rb_eArgError, "rank of c (7th argument) must be %d", 2);
  ncc = NA_SHAPE1(rblapack_c);
  ldc = NA_SHAPE0(rblapack_c);
  if (NA_TYPE(rblapack_c) != NA_SCOMPLEX)
    rblapack_c = na_change_type(rblapack_c, NA_SCOMPLEX);
  c = NA_PTR_TYPE(rblapack_c, complex*);
  if (!NA_IsNArray(rblapack_d))
    rb_raise(rb_eArgError, "d (3th argument) must be NArray");
  if (NA_RANK(rblapack_d) != 1)
    rb_raise(rb_eArgError, "rank of d (3th argument) must be %d", 1);
  n = NA_SHAPE0(rblapack_d);
  if (NA_TYPE(rblapack_d) != NA_SFLOAT)
    rblapack_d = na_change_type(rblapack_d, NA_SFLOAT);
  d = NA_PTR_TYPE(rblapack_d, real*);
  uplo = StringValueCStr(rblapack_uplo)[0];
  if (!NA_IsNArray(rblapack_u))
    rb_raise(rb_eArgError, "u (6th argument) must be NArray");
  if (NA_RANK(rblapack_u) != 2)
    rb_raise(rb_eArgError, "rank of u (6th argument) must be %d", 2);
  if (NA_SHAPE1(rblapack_u) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of u must be the same as shape 0 of d");
  ldu = NA_SHAPE0(rblapack_u);
  if (NA_TYPE(rblapack_u) != NA_SCOMPLEX)
    rblapack_u = na_change_type(rblapack_u, NA_SCOMPLEX);
  u = NA_PTR_TYPE(rblapack_u, complex*);
  if (!NA_IsNArray(rblapack_vt))
    rb_raise(rb_eArgError, "vt (5th argument) must be NArray");
  if (NA_RANK(rblapack_vt) != 2)
    rb_raise(rb_eArgError, "rank of vt (5th argument) must be %d", 2);
  ncvt = NA_SHAPE1(rblapack_vt);
  ldvt = NA_SHAPE0(rblapack_vt);
  if (NA_TYPE(rblapack_vt) != NA_SCOMPLEX)
    rblapack_vt = na_change_type(rblapack_vt, NA_SCOMPLEX);
  vt = NA_PTR_TYPE(rblapack_vt, complex*);
  nru = NUM2INT(rblapack_nru);
  if (!NA_IsNArray(rblapack_e))
    rb_raise(rb_eArgError, "e (4th argument) must be NArray");
  if (NA_RANK(rblapack_e) != 1)
    rb_raise(rb_eArgError, "rank of e (4th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_e) != (n-1))
    rb_raise(rb_eRuntimeError, "shape 0 of e must be %d", n-1);
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
    shape[0] = n-1;
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
    rblapack_vt_out__ = na_make_object(NA_SCOMPLEX, 2, shape, cNArray);
  }
  vt_out__ = NA_PTR_TYPE(rblapack_vt_out__, complex*);
  MEMCPY(vt_out__, vt, complex, NA_TOTAL(rblapack_vt));
  rblapack_vt = rblapack_vt_out__;
  vt = vt_out__;
  {
    int shape[2];
    shape[0] = ldu;
    shape[1] = n;
    rblapack_u_out__ = na_make_object(NA_SCOMPLEX, 2, shape, cNArray);
  }
  u_out__ = NA_PTR_TYPE(rblapack_u_out__, complex*);
  MEMCPY(u_out__, u, complex, NA_TOTAL(rblapack_u));
  rblapack_u = rblapack_u_out__;
  u = u_out__;
  {
    int shape[2];
    shape[0] = ldc;
    shape[1] = ncc;
    rblapack_c_out__ = na_make_object(NA_SCOMPLEX, 2, shape, cNArray);
  }
  c_out__ = NA_PTR_TYPE(rblapack_c_out__, complex*);
  MEMCPY(c_out__, c, complex, NA_TOTAL(rblapack_c));
  rblapack_c = rblapack_c_out__;
  c = c_out__;
  rwork = ALLOC_N(real, ((ncvt==nru)&&(nru==ncc)&&(ncc==0) ? 2*n : MAX(1, 4*n-4)));

  cbdsqr_(&uplo, &n, &ncvt, &nru, &ncc, d, e, vt, &ldvt, u, &ldu, c, &ldc, rwork, &info);

  free(rwork);
  rblapack_info = INT2NUM(info);
  return rb_ary_new3(6, rblapack_info, rblapack_d, rblapack_e, rblapack_vt, rblapack_u, rblapack_c);
}

void
init_lapack_cbdsqr(VALUE mLapack){
  rb_define_module_function(mLapack, "cbdsqr", rblapack_cbdsqr, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}

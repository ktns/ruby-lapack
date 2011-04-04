#include "rb_lapack.h"

extern VOID sorcsd_(char *jobu1, char *jobu2, char *jobv1t, char *jobv2t, char *trans, char *signs, integer *m, integer *p, integer *q, real *x11, integer *ldx11, real *x12, integer *ldx12, real *x21, integer *ldx21, real *x22, integer *ldx22, real *theta, real *u1, integer *ldu1, real *u2, integer *ldu2, real *v1t, integer *ldv1t, real *v2t, integer *ldv2t, real *work, integer *lwork, integer *iwork, integer *info);

static VALUE
rb_sorcsd(int argc, VALUE *argv, VALUE self){
  VALUE rb_jobu1;
  char jobu1; 
  VALUE rb_jobu2;
  char jobu2; 
  VALUE rb_jobv1t;
  char jobv1t; 
  VALUE rb_jobv2t;
  char jobv2t; 
  VALUE rb_trans;
  char trans; 
  VALUE rb_signs;
  char signs; 
  VALUE rb_m;
  integer m; 
  VALUE rb_x11;
  real *x11; 
  VALUE rb_x12;
  real *x12; 
  VALUE rb_x21;
  real *x21; 
  VALUE rb_x22;
  real *x22; 
  VALUE rb_ldu1;
  integer ldu1; 
  VALUE rb_ldu2;
  integer ldu2; 
  VALUE rb_ldv1t;
  integer ldv1t; 
  VALUE rb_ldv2t;
  integer ldv2t; 
  VALUE rb_lwork;
  integer lwork; 
  VALUE rb_theta;
  real *theta; 
  VALUE rb_u1;
  real *u1; 
  VALUE rb_u2;
  real *u2; 
  VALUE rb_v1t;
  real *v1t; 
  VALUE rb_v2t;
  real *v2t; 
  VALUE rb_info;
  integer info; 
  real *work;
  integer *iwork;

  integer ldx11;
  integer q;
  integer ldx12;
  integer ldx21;
  integer ldx22;
  integer p;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  theta, u1, u2, v1t, v2t, info = NumRu::Lapack.sorcsd( jobu1, jobu2, jobv1t, jobv2t, trans, signs, m, x11, x12, x21, x22, ldu1, ldu2, ldv1t, ldv2t, lwork)\n    or\n  NumRu::Lapack.sorcsd  # print help\n\n\nFORTRAN MANUAL\n      RECURSIVE SUBROUTINE SORCSD( JOBU1, JOBU2, JOBV1T, JOBV2T, TRANS, SIGNS, M, P, Q, X11, LDX11, X12, LDX12, X21, LDX21, X22, LDX22, THETA, U1, LDU1, U2, LDU2, V1T, LDV1T, V2T, LDV2T, WORK, LWORK, IWORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  SORCSD computes the CS decomposition of an M-by-M partitioned\n*  orthogonal matrix X:\n*\n*                                  [  I  0  0 |  0  0  0 ]\n*                                  [  0  C  0 |  0 -S  0 ]\n*      [ X11 | X12 ]   [ U1 |    ] [  0  0  0 |  0  0 -I ] [ V1 |    ]**T\n*  X = [-----------] = [---------] [---------------------] [---------]   .\n*      [ X21 | X22 ]   [    | U2 ] [  0  0  0 |  I  0  0 ] [    | V2 ]\n*                                  [  0  S  0 |  0  C  0 ]\n*                                  [  0  0  I |  0  0  0 ]\n*\n*  X11 is P-by-Q. The orthogonal matrices U1, U2, V1, and V2 are P-by-P,\n*  (M-P)-by-(M-P), Q-by-Q, and (M-Q)-by-(M-Q), respectively. C and S are\n*  R-by-R nonnegative diagonal matrices satisfying C^2 + S^2 = I, in\n*  which R = MIN(P,M-P,Q,M-Q).\n*\n\n*  Arguments\n*  =========\n*\n*  JOBU1   (input) CHARACTER\n*          = 'Y':      U1 is computed;\n*          otherwise:  U1 is not computed.\n*\n*  JOBU2   (input) CHARACTER\n*          = 'Y':      U2 is computed;\n*          otherwise:  U2 is not computed.\n*\n*  JOBV1T  (input) CHARACTER\n*          = 'Y':      V1T is computed;\n*          otherwise:  V1T is not computed.\n*\n*  JOBV2T  (input) CHARACTER\n*          = 'Y':      V2T is computed;\n*          otherwise:  V2T is not computed.\n*\n*  TRANS   (input) CHARACTER\n*          = 'T':      X, U1, U2, V1T, and V2T are stored in row-major\n*                      order;\n*          otherwise:  X, U1, U2, V1T, and V2T are stored in column-\n*                      major order.\n*\n*  SIGNS   (input) CHARACTER\n*          = 'O':      The lower-left block is made nonpositive (the\n*                      \"other\" convention);\n*          otherwise:  The upper-right block is made nonpositive (the\n*                      \"default\" convention).\n*\n*  M       (input) INTEGER\n*          The number of rows and columns in X.\n*\n*  P       (input) INTEGER\n*          The number of rows in X11 and X12. 0 <= P <= M.\n*\n*  Q       (input) INTEGER\n*          The number of columns in X11 and X21. 0 <= Q <= M.\n*\n*  X       (input/workspace) REAL array, dimension (LDX,M)\n*          On entry, the orthogonal matrix whose CSD is desired.\n*\n*  LDX     (input) INTEGER\n*          The leading dimension of X. LDX >= MAX(1,M).\n*\n*  THETA   (output) REAL array, dimension (R), in which R =\n*          MIN(P,M-P,Q,M-Q).\n*          C = DIAG( COS(THETA(1)), ... , COS(THETA(R)) ) and\n*          S = DIAG( SIN(THETA(1)), ... , SIN(THETA(R)) ).\n*\n*  U1      (output) REAL array, dimension (P)\n*          If JOBU1 = 'Y', U1 contains the P-by-P orthogonal matrix U1.\n*\n*  LDU1    (input) INTEGER\n*          The leading dimension of U1. If JOBU1 = 'Y', LDU1 >=\n*          MAX(1,P).\n*\n*  U2      (output) REAL array, dimension (M-P)\n*          If JOBU2 = 'Y', U2 contains the (M-P)-by-(M-P) orthogonal\n*          matrix U2.\n*\n*  LDU2    (input) INTEGER\n*          The leading dimension of U2. If JOBU2 = 'Y', LDU2 >=\n*          MAX(1,M-P).\n*\n*  V1T     (output) REAL array, dimension (Q)\n*          If JOBV1T = 'Y', V1T contains the Q-by-Q matrix orthogonal\n*          matrix V1**T.\n*\n*  LDV1T   (input) INTEGER\n*          The leading dimension of V1T. If JOBV1T = 'Y', LDV1T >=\n*          MAX(1,Q).\n*\n*  V2T     (output) REAL array, dimension (M-Q)\n*          If JOBV2T = 'Y', V2T contains the (M-Q)-by-(M-Q) orthogonal\n*          matrix V2**T.\n*\n*  LDV2T   (input) INTEGER\n*          The leading dimension of V2T. If JOBV2T = 'Y', LDV2T >=\n*          MAX(1,M-Q).\n*\n*  WORK    (workspace) REAL array, dimension (MAX(1,LWORK))\n*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.\n*          If INFO > 0 on exit, WORK(2:R) contains the values PHI(1),\n*          ..., PHI(R-1) that, together with THETA(1), ..., THETA(R),\n*          define the matrix in intermediate bidiagonal-block form\n*          remaining after nonconvergence. INFO specifies the number\n*          of nonzero PHI's.\n*\n*  LWORK   (input) INTEGER\n*          The dimension of the array WORK.\n*\n*          If LWORK = -1, then a workspace query is assumed; the routine\n*          only calculates the optimal size of the WORK array, returns\n*          this value as the first entry of the work array, and no error\n*          message related to LWORK is issued by XERBLA.\n*\n*  IWORK   (workspace) INTEGER array, dimension (M-Q)\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit.\n*          < 0:  if INFO = -i, the i-th argument had an illegal value.\n*          > 0:  SBBCSD did not converge. See the description of WORK\n*                above for details.\n*\n*  Reference\n*  =========\n*\n*  [1] Brian D. Sutton. Computing the complete CS decomposition. Numer.\n*      Algorithms, 50(1):33-65, 2009.\n*\n\n*  ===================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 16)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 16)", argc);
  rb_jobu1 = argv[0];
  rb_jobu2 = argv[1];
  rb_jobv1t = argv[2];
  rb_jobv2t = argv[3];
  rb_trans = argv[4];
  rb_signs = argv[5];
  rb_m = argv[6];
  rb_x11 = argv[7];
  rb_x12 = argv[8];
  rb_x21 = argv[9];
  rb_x22 = argv[10];
  rb_ldu1 = argv[11];
  rb_ldu2 = argv[12];
  rb_ldv1t = argv[13];
  rb_ldv2t = argv[14];
  rb_lwork = argv[15];

  trans = StringValueCStr(rb_trans)[0];
  jobv1t = StringValueCStr(rb_jobv1t)[0];
  jobv2t = StringValueCStr(rb_jobv2t)[0];
  lwork = NUM2INT(rb_lwork);
  signs = StringValueCStr(rb_signs)[0];
  if (!NA_IsNArray(rb_x21))
    rb_raise(rb_eArgError, "x21 (10th argument) must be NArray");
  if (NA_RANK(rb_x21) != 2)
    rb_raise(rb_eArgError, "rank of x21 (10th argument) must be %d", 2);
  q = NA_SHAPE1(rb_x21);
  ldx21 = NA_SHAPE0(rb_x21);
  if (ldx21 != (p))
    rb_raise(rb_eRuntimeError, "shape 0 of x21 must be %d", p);
  p = ldx21;
  if (NA_TYPE(rb_x21) != NA_SFLOAT)
    rb_x21 = na_change_type(rb_x21, NA_SFLOAT);
  x21 = NA_PTR_TYPE(rb_x21, real*);
  jobu1 = StringValueCStr(rb_jobu1)[0];
  jobu2 = StringValueCStr(rb_jobu2)[0];
  if (!NA_IsNArray(rb_x11))
    rb_raise(rb_eArgError, "x11 (8th argument) must be NArray");
  if (NA_RANK(rb_x11) != 2)
    rb_raise(rb_eArgError, "rank of x11 (8th argument) must be %d", 2);
  if (NA_SHAPE1(rb_x11) != q)
    rb_raise(rb_eRuntimeError, "shape 1 of x11 must be the same as shape 1 of x21");
  ldx11 = NA_SHAPE0(rb_x11);
  if (ldx11 != (p))
    rb_raise(rb_eRuntimeError, "shape 0 of x11 must be %d", p);
  p = ldx11;
  if (NA_TYPE(rb_x11) != NA_SFLOAT)
    rb_x11 = na_change_type(rb_x11, NA_SFLOAT);
  x11 = NA_PTR_TYPE(rb_x11, real*);
  m = NUM2INT(rb_m);
  if (!NA_IsNArray(rb_x22))
    rb_raise(rb_eArgError, "x22 (11th argument) must be NArray");
  if (NA_RANK(rb_x22) != 2)
    rb_raise(rb_eArgError, "rank of x22 (11th argument) must be %d", 2);
  if (NA_SHAPE1(rb_x22) != (m-q))
    rb_raise(rb_eRuntimeError, "shape 1 of x22 must be %d", m-q);
  ldx22 = NA_SHAPE0(rb_x22);
  if (ldx22 != (p))
    rb_raise(rb_eRuntimeError, "shape 0 of x22 must be %d", p);
  p = ldx22;
  if (NA_TYPE(rb_x22) != NA_SFLOAT)
    rb_x22 = na_change_type(rb_x22, NA_SFLOAT);
  x22 = NA_PTR_TYPE(rb_x22, real*);
  ldu1 = lsame_(&jobu1,"Y") ? MAX(1,p) : 0;
  ldu2 = lsame_(&jobu2,"Y") ? MAX(1,m-p) : 0;
  if (!NA_IsNArray(rb_x12))
    rb_raise(rb_eArgError, "x12 (9th argument) must be NArray");
  if (NA_RANK(rb_x12) != 2)
    rb_raise(rb_eArgError, "rank of x12 (9th argument) must be %d", 2);
  if (NA_SHAPE1(rb_x12) != (m-q))
    rb_raise(rb_eRuntimeError, "shape 1 of x12 must be %d", m-q);
  ldx12 = NA_SHAPE0(rb_x12);
  if (ldx12 != (p))
    rb_raise(rb_eRuntimeError, "shape 0 of x12 must be %d", p);
  p = ldx12;
  if (NA_TYPE(rb_x12) != NA_SFLOAT)
    rb_x12 = na_change_type(rb_x12, NA_SFLOAT);
  x12 = NA_PTR_TYPE(rb_x12, real*);
  ldv1t = lsame_(&jobv1t,"Y") ? MAX(1,q) : 0;
  ldx12 = p;
  ldv2t = lsame_(&jobv2t,"Y") ? MAX(1,m-q) : 0;
  ldx22 = p;
  ldx21 = p;
  ldx11 = p;
  {
    int shape[1];
    shape[0] = MIN(MIN(MIN(p,m-p),q),m-q);
    rb_theta = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  theta = NA_PTR_TYPE(rb_theta, real*);
  {
    int shape[1];
    shape[0] = p;
    rb_u1 = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  u1 = NA_PTR_TYPE(rb_u1, real*);
  {
    int shape[1];
    shape[0] = m-p;
    rb_u2 = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  u2 = NA_PTR_TYPE(rb_u2, real*);
  {
    int shape[1];
    shape[0] = q;
    rb_v1t = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  v1t = NA_PTR_TYPE(rb_v1t, real*);
  {
    int shape[1];
    shape[0] = m-q;
    rb_v2t = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  v2t = NA_PTR_TYPE(rb_v2t, real*);
  work = ALLOC_N(real, (MAX(1,lwork)));
  iwork = ALLOC_N(integer, (m-q));

  sorcsd_(&jobu1, &jobu2, &jobv1t, &jobv2t, &trans, &signs, &m, &p, &q, x11, &ldx11, x12, &ldx12, x21, &ldx21, x22, &ldx22, theta, u1, &ldu1, u2, &ldu2, v1t, &ldv1t, v2t, &ldv2t, work, &lwork, iwork, &info);

  free(work);
  free(iwork);
  rb_info = INT2NUM(info);
  return rb_ary_new3(6, rb_theta, rb_u1, rb_u2, rb_v1t, rb_v2t, rb_info);
}

void
init_lapack_sorcsd(VALUE mLapack){
  rb_define_module_function(mLapack, "sorcsd", rb_sorcsd, -1);
}

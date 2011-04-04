#include "rb_lapack.h"

extern VOID dorbdb_(char *trans, char *signs, integer *m, integer *p, integer *q, doublereal *x11, integer *ldx11, doublereal *x12, integer *ldx12, doublereal *x21, integer *ldx21, doublereal *x22, integer *ldx22, doublereal *theta, doublereal *phi, doublereal *taup1, doublereal *taup2, doublereal *tauq1, doublereal *tauq2, doublereal *work, integer *lwork, integer *info);

static VALUE
rb_dorbdb(int argc, VALUE *argv, VALUE self){
  VALUE rb_trans;
  char trans; 
  VALUE rb_signs;
  char signs; 
  VALUE rb_m;
  integer m; 
  VALUE rb_x11;
  doublereal *x11; 
  VALUE rb_x12;
  doublereal *x12; 
  VALUE rb_x21;
  doublereal *x21; 
  VALUE rb_x22;
  doublereal *x22; 
  VALUE rb_lwork;
  integer lwork; 
  VALUE rb_theta;
  doublereal *theta; 
  VALUE rb_phi;
  doublereal *phi; 
  VALUE rb_taup1;
  doublereal *taup1; 
  VALUE rb_taup2;
  doublereal *taup2; 
  VALUE rb_tauq1;
  doublereal *tauq1; 
  VALUE rb_tauq2;
  doublereal *tauq2; 
  VALUE rb_info;
  integer info; 
  VALUE rb_x11_out__;
  doublereal *x11_out__;
  VALUE rb_x12_out__;
  doublereal *x12_out__;
  VALUE rb_x21_out__;
  doublereal *x21_out__;
  VALUE rb_x22_out__;
  doublereal *x22_out__;
  doublereal *work;

  integer ldx11;
  integer q;
  integer ldx12;
  integer ldx21;
  integer ldx22;
  integer p;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  theta, phi, taup1, taup2, tauq1, tauq2, info, x11, x12, x21, x22 = NumRu::Lapack.dorbdb( trans, signs, m, x11, x12, x21, x22, lwork)\n    or\n  NumRu::Lapack.dorbdb  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE DORBDB( TRANS, SIGNS, M, P, Q, X11, LDX11, X12, LDX12, X21, LDX21, X22, LDX22, THETA, PHI, TAUP1, TAUP2, TAUQ1, TAUQ2, WORK, LWORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  DORBDB simultaneously bidiagonalizes the blocks of an M-by-M\n*  partitioned orthogonal matrix X:\n*\n*                                  [ B11 | B12 0  0 ]\n*      [ X11 | X12 ]   [ P1 |    ] [  0  |  0 -I  0 ] [ Q1 |    ]**T\n*  X = [-----------] = [---------] [----------------] [---------]   .\n*      [ X21 | X22 ]   [    | P2 ] [ B21 | B22 0  0 ] [    | Q2 ]\n*                                  [  0  |  0  0  I ]\n*\n*  X11 is P-by-Q. Q must be no larger than P, M-P, or M-Q. (If this is\n*  not the case, then X must be transposed and/or permuted. This can be\n*  done in constant time using the TRANS and SIGNS options. See DORCSD\n*  for details.)\n*\n*  The orthogonal matrices P1, P2, Q1, and Q2 are P-by-P, (M-P)-by-\n*  (M-P), Q-by-Q, and (M-Q)-by-(M-Q), respectively. They are\n*  represented implicitly by Householder vectors.\n*\n*  B11, B12, B21, and B22 are Q-by-Q bidiagonal matrices represented\n*  implicitly by angles THETA, PHI.\n*\n\n*  Arguments\n*  =========\n*\n*  TRANS   (input) CHARACTER\n*          = 'T':      X, U1, U2, V1T, and V2T are stored in row-major\n*                      order;\n*          otherwise:  X, U1, U2, V1T, and V2T are stored in column-\n*                      major order.\n*\n*  SIGNS   (input) CHARACTER\n*          = 'O':      The lower-left block is made nonpositive (the\n*                      \"other\" convention);\n*          otherwise:  The upper-right block is made nonpositive (the\n*                      \"default\" convention).\n*\n*  M       (input) INTEGER\n*          The number of rows and columns in X.\n*\n*  P       (input) INTEGER\n*          The number of rows in X11 and X12. 0 <= P <= M.\n*\n*  Q       (input) INTEGER\n*          The number of columns in X11 and X21. 0 <= Q <=\n*          MIN(P,M-P,M-Q).\n*\n*  X11     (input/output) DOUBLE PRECISION array, dimension (LDX11,Q)\n*          On entry, the top-left block of the orthogonal matrix to be\n*          reduced. On exit, the form depends on TRANS:\n*          If TRANS = 'N', then\n*             the columns of tril(X11) specify reflectors for P1,\n*             the rows of triu(X11,1) specify reflectors for Q1;\n*          else TRANS = 'T', and\n*             the rows of triu(X11) specify reflectors for P1,\n*             the columns of tril(X11,-1) specify reflectors for Q1.\n*\n*  LDX11   (input) INTEGER\n*          The leading dimension of X11. If TRANS = 'N', then LDX11 >=\n*          P; else LDX11 >= Q.\n*\n*  X12     (input/output) DOUBLE PRECISION array, dimension (LDX12,M-Q)\n*          On entry, the top-right block of the orthogonal matrix to\n*          be reduced. On exit, the form depends on TRANS:\n*          If TRANS = 'N', then\n*             the rows of triu(X12) specify the first P reflectors for\n*             Q2;\n*          else TRANS = 'T', and\n*             the columns of tril(X12) specify the first P reflectors\n*             for Q2.\n*\n*  LDX12   (input) INTEGER\n*          The leading dimension of X12. If TRANS = 'N', then LDX12 >=\n*          P; else LDX11 >= M-Q.\n*\n*  X21     (input/output) DOUBLE PRECISION array, dimension (LDX21,Q)\n*          On entry, the bottom-left block of the orthogonal matrix to\n*          be reduced. On exit, the form depends on TRANS:\n*          If TRANS = 'N', then\n*             the columns of tril(X21) specify reflectors for P2;\n*          else TRANS = 'T', and\n*             the rows of triu(X21) specify reflectors for P2.\n*\n*  LDX21   (input) INTEGER\n*          The leading dimension of X21. If TRANS = 'N', then LDX21 >=\n*          M-P; else LDX21 >= Q.\n*\n*  X22     (input/output) DOUBLE PRECISION array, dimension (LDX22,M-Q)\n*          On entry, the bottom-right block of the orthogonal matrix to\n*          be reduced. On exit, the form depends on TRANS:\n*          If TRANS = 'N', then\n*             the rows of triu(X22(Q+1:M-P,P+1:M-Q)) specify the last\n*             M-P-Q reflectors for Q2,\n*          else TRANS = 'T', and\n*             the columns of tril(X22(P+1:M-Q,Q+1:M-P)) specify the last\n*             M-P-Q reflectors for P2.\n*\n*  LDX22   (input) INTEGER\n*          The leading dimension of X22. If TRANS = 'N', then LDX22 >=\n*          M-P; else LDX22 >= M-Q.\n*\n*  THETA   (output) DOUBLE PRECISION array, dimension (Q)\n*          The entries of the bidiagonal blocks B11, B12, B21, B22 can\n*          be computed from the angles THETA and PHI. See Further\n*          Details.\n*\n*  PHI     (output) DOUBLE PRECISION array, dimension (Q-1)\n*          The entries of the bidiagonal blocks B11, B12, B21, B22 can\n*          be computed from the angles THETA and PHI. See Further\n*          Details.\n*\n*  TAUP1   (output) DOUBLE PRECISION array, dimension (P)\n*          The scalar factors of the elementary reflectors that define\n*          P1.\n*\n*  TAUP2   (output) DOUBLE PRECISION array, dimension (M-P)\n*          The scalar factors of the elementary reflectors that define\n*          P2.\n*\n*  TAUQ1   (output) DOUBLE PRECISION array, dimension (Q)\n*          The scalar factors of the elementary reflectors that define\n*          Q1.\n*\n*  TAUQ2   (output) DOUBLE PRECISION array, dimension (M-Q)\n*          The scalar factors of the elementary reflectors that define\n*          Q2.\n*\n*  WORK    (workspace) DOUBLE PRECISION array, dimension (LWORK)\n*\n*  LWORK   (input) INTEGER\n*          The dimension of the array WORK. LWORK >= M-Q.\n*\n*          If LWORK = -1, then a workspace query is assumed; the routine\n*          only calculates the optimal size of the WORK array, returns\n*          this value as the first entry of the WORK array, and no error\n*          message related to LWORK is issued by XERBLA.\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit.\n*          < 0:  if INFO = -i, the i-th argument had an illegal value.\n*\n\n*  Further Details\n*  ===============\n*\n*  The bidiagonal blocks B11, B12, B21, and B22 are represented\n*  implicitly by angles THETA(1), ..., THETA(Q) and PHI(1), ...,\n*  PHI(Q-1). B11 and B21 are upper bidiagonal, while B21 and B22 are\n*  lower bidiagonal. Every entry in each bidiagonal band is a product\n*  of a sine or cosine of a THETA with a sine or cosine of a PHI. See\n*  [1] or DORCSD for details.\n*\n*  P1, P2, Q1, and Q2 are represented as products of elementary\n*  reflectors. See DORCSD for details on generating P1, P2, Q1, and Q2\n*  using DORGQR and DORGLQ.\n*\n*  Reference\n*  =========\n*\n*  [1] Brian D. Sutton. Computing the complete CS decomposition. Numer.\n*      Algorithms, 50(1):33-65, 2009.\n*\n*  ====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 8)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 8)", argc);
  rb_trans = argv[0];
  rb_signs = argv[1];
  rb_m = argv[2];
  rb_x11 = argv[3];
  rb_x12 = argv[4];
  rb_x21 = argv[5];
  rb_x22 = argv[6];
  rb_lwork = argv[7];

  trans = StringValueCStr(rb_trans)[0];
  lwork = NUM2INT(rb_lwork);
  signs = StringValueCStr(rb_signs)[0];
  if (!NA_IsNArray(rb_x21))
    rb_raise(rb_eArgError, "x21 (6th argument) must be NArray");
  if (NA_RANK(rb_x21) != 2)
    rb_raise(rb_eArgError, "rank of x21 (6th argument) must be %d", 2);
  q = NA_SHAPE1(rb_x21);
  ldx21 = NA_SHAPE0(rb_x21);
  if (ldx21 != (p))
    rb_raise(rb_eRuntimeError, "shape 0 of x21 must be %d", p);
  p = ldx21;
  if (NA_TYPE(rb_x21) != NA_DFLOAT)
    rb_x21 = na_change_type(rb_x21, NA_DFLOAT);
  x21 = NA_PTR_TYPE(rb_x21, doublereal*);
  if (!NA_IsNArray(rb_x11))
    rb_raise(rb_eArgError, "x11 (4th argument) must be NArray");
  if (NA_RANK(rb_x11) != 2)
    rb_raise(rb_eArgError, "rank of x11 (4th argument) must be %d", 2);
  if (NA_SHAPE1(rb_x11) != q)
    rb_raise(rb_eRuntimeError, "shape 1 of x11 must be the same as shape 1 of x21");
  ldx11 = NA_SHAPE0(rb_x11);
  if (ldx11 != (p))
    rb_raise(rb_eRuntimeError, "shape 0 of x11 must be %d", p);
  p = ldx11;
  if (NA_TYPE(rb_x11) != NA_DFLOAT)
    rb_x11 = na_change_type(rb_x11, NA_DFLOAT);
  x11 = NA_PTR_TYPE(rb_x11, doublereal*);
  m = NUM2INT(rb_m);
  if (!NA_IsNArray(rb_x22))
    rb_raise(rb_eArgError, "x22 (7th argument) must be NArray");
  if (NA_RANK(rb_x22) != 2)
    rb_raise(rb_eArgError, "rank of x22 (7th argument) must be %d", 2);
  if (NA_SHAPE1(rb_x22) != (m-q))
    rb_raise(rb_eRuntimeError, "shape 1 of x22 must be %d", m-q);
  ldx22 = NA_SHAPE0(rb_x22);
  if (ldx22 != (p))
    rb_raise(rb_eRuntimeError, "shape 0 of x22 must be %d", p);
  p = ldx22;
  if (NA_TYPE(rb_x22) != NA_DFLOAT)
    rb_x22 = na_change_type(rb_x22, NA_DFLOAT);
  x22 = NA_PTR_TYPE(rb_x22, doublereal*);
  if (!NA_IsNArray(rb_x12))
    rb_raise(rb_eArgError, "x12 (5th argument) must be NArray");
  if (NA_RANK(rb_x12) != 2)
    rb_raise(rb_eArgError, "rank of x12 (5th argument) must be %d", 2);
  if (NA_SHAPE1(rb_x12) != (m-q))
    rb_raise(rb_eRuntimeError, "shape 1 of x12 must be %d", m-q);
  ldx12 = NA_SHAPE0(rb_x12);
  if (ldx12 != (p))
    rb_raise(rb_eRuntimeError, "shape 0 of x12 must be %d", p);
  p = ldx12;
  if (NA_TYPE(rb_x12) != NA_DFLOAT)
    rb_x12 = na_change_type(rb_x12, NA_DFLOAT);
  x12 = NA_PTR_TYPE(rb_x12, doublereal*);
  ldx12 = p;
  ldx22 = p;
  ldx21 = p;
  ldx11 = p;
  {
    int shape[1];
    shape[0] = q;
    rb_theta = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  theta = NA_PTR_TYPE(rb_theta, doublereal*);
  {
    int shape[1];
    shape[0] = q-1;
    rb_phi = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  phi = NA_PTR_TYPE(rb_phi, doublereal*);
  {
    int shape[1];
    shape[0] = p;
    rb_taup1 = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  taup1 = NA_PTR_TYPE(rb_taup1, doublereal*);
  {
    int shape[1];
    shape[0] = m-p;
    rb_taup2 = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  taup2 = NA_PTR_TYPE(rb_taup2, doublereal*);
  {
    int shape[1];
    shape[0] = q;
    rb_tauq1 = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  tauq1 = NA_PTR_TYPE(rb_tauq1, doublereal*);
  {
    int shape[1];
    shape[0] = m-q;
    rb_tauq2 = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  tauq2 = NA_PTR_TYPE(rb_tauq2, doublereal*);
  {
    int shape[2];
    shape[0] = ldx11;
    shape[1] = q;
    rb_x11_out__ = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  x11_out__ = NA_PTR_TYPE(rb_x11_out__, doublereal*);
  MEMCPY(x11_out__, x11, doublereal, NA_TOTAL(rb_x11));
  rb_x11 = rb_x11_out__;
  x11 = x11_out__;
  {
    int shape[2];
    shape[0] = ldx12;
    shape[1] = m-q;
    rb_x12_out__ = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  x12_out__ = NA_PTR_TYPE(rb_x12_out__, doublereal*);
  MEMCPY(x12_out__, x12, doublereal, NA_TOTAL(rb_x12));
  rb_x12 = rb_x12_out__;
  x12 = x12_out__;
  {
    int shape[2];
    shape[0] = ldx21;
    shape[1] = q;
    rb_x21_out__ = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  x21_out__ = NA_PTR_TYPE(rb_x21_out__, doublereal*);
  MEMCPY(x21_out__, x21, doublereal, NA_TOTAL(rb_x21));
  rb_x21 = rb_x21_out__;
  x21 = x21_out__;
  {
    int shape[2];
    shape[0] = ldx22;
    shape[1] = m-q;
    rb_x22_out__ = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  x22_out__ = NA_PTR_TYPE(rb_x22_out__, doublereal*);
  MEMCPY(x22_out__, x22, doublereal, NA_TOTAL(rb_x22));
  rb_x22 = rb_x22_out__;
  x22 = x22_out__;
  work = ALLOC_N(doublereal, (MAX(1,lwork)));

  dorbdb_(&trans, &signs, &m, &p, &q, x11, &ldx11, x12, &ldx12, x21, &ldx21, x22, &ldx22, theta, phi, taup1, taup2, tauq1, tauq2, work, &lwork, &info);

  free(work);
  rb_info = INT2NUM(info);
  return rb_ary_new3(11, rb_theta, rb_phi, rb_taup1, rb_taup2, rb_tauq1, rb_tauq2, rb_info, rb_x11, rb_x12, rb_x21, rb_x22);
}

void
init_lapack_dorbdb(VALUE mLapack){
  rb_define_module_function(mLapack, "dorbdb", rb_dorbdb, -1);
}

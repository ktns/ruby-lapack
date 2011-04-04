#include "rb_lapack.h"

extern VOID cbbcsd_(char *jobu1, char *jobu2, char *jobv1t, char *jobv2t, char *trans, integer *m, integer *p, integer *q, real *theta, real *phi, complex *u1, integer *ldu1, complex *u2, integer *ldu2, complex *v1t, integer *ldv1t, complex *v2t, integer *ldv2t, real *b11d, real *b11e, real *b12d, real *b12e, real *b21d, real *b21e, real *b22d, real *b22e, real *rwork, integer *lrwork, integer *info);

static VALUE
rb_cbbcsd(int argc, VALUE *argv, VALUE self){
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
  VALUE rb_m;
  integer m; 
  VALUE rb_theta;
  real *theta; 
  VALUE rb_phi;
  real *phi; 
  VALUE rb_u1;
  complex *u1; 
  VALUE rb_u2;
  complex *u2; 
  VALUE rb_v1t;
  complex *v1t; 
  VALUE rb_v2t;
  complex *v2t; 
  VALUE rb_lrwork;
  integer lrwork; 
  VALUE rb_b11d;
  real *b11d; 
  VALUE rb_b11e;
  real *b11e; 
  VALUE rb_b12d;
  real *b12d; 
  VALUE rb_b12e;
  real *b12e; 
  VALUE rb_b21d;
  real *b21d; 
  VALUE rb_b21e;
  real *b21e; 
  VALUE rb_b22d;
  real *b22d; 
  VALUE rb_b22e;
  real *b22e; 
  VALUE rb_info;
  integer info; 
  VALUE rb_theta_out__;
  real *theta_out__;
  VALUE rb_u1_out__;
  complex *u1_out__;
  VALUE rb_u2_out__;
  complex *u2_out__;
  VALUE rb_v1t_out__;
  complex *v1t_out__;
  VALUE rb_v2t_out__;
  complex *v2t_out__;
  real *rwork;

  integer q;
  integer ldu1;
  integer p;
  integer ldu2;
  integer ldv1t;
  integer ldv2t;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  b11d, b11e, b12d, b12e, b21d, b21e, b22d, b22e, info, theta, u1, u2, v1t, v2t = NumRu::Lapack.cbbcsd( jobu1, jobu2, jobv1t, jobv2t, trans, m, theta, phi, u1, u2, v1t, v2t, lrwork)\n    or\n  NumRu::Lapack.cbbcsd  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE CBBCSD( JOBU1, JOBU2, JOBV1T, JOBV2T, TRANS, M, P, Q, THETA, PHI, U1, LDU1, U2, LDU2, V1T, LDV1T, V2T, LDV2T, B11D, B11E, B12D, B12E, B21D, B21E, B22D, B22E, RWORK, LRWORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  CBBCSD computes the CS decomposition of a unitary matrix in\n*  bidiagonal-block form,\n*\n*\n*      [ B11 | B12 0  0 ]\n*      [  0  |  0 -I  0 ]\n*  X = [----------------]\n*      [ B21 | B22 0  0 ]\n*      [  0  |  0  0  I ]\n*\n*                                [  C | -S  0  0 ]\n*                    [ U1 |    ] [  0 |  0 -I  0 ] [ V1 |    ]**H\n*                  = [---------] [---------------] [---------]   .\n*                    [    | U2 ] [  S |  C  0  0 ] [    | V2 ]\n*                                [  0 |  0  0  I ]\n*\n*  X is M-by-M, its top-left block is P-by-Q, and Q must be no larger\n*  than P, M-P, or M-Q. (If Q is not the smallest index, then X must be\n*  transposed and/or permuted. This can be done in constant time using\n*  the TRANS and SIGNS options. See CUNCSD for details.)\n*\n*  The bidiagonal matrices B11, B12, B21, and B22 are represented\n*  implicitly by angles THETA(1:Q) and PHI(1:Q-1).\n*\n*  The unitary matrices U1, U2, V1T, and V2T are input/output.\n*  The input matrices are pre- or post-multiplied by the appropriate\n*  singular vector matrices.\n*\n\n*  Arguments\n*  =========\n*\n*  JOBU1   (input) CHARACTER\n*          = 'Y':      U1 is updated;\n*          otherwise:  U1 is not updated.\n*\n*  JOBU2   (input) CHARACTER\n*          = 'Y':      U2 is updated;\n*          otherwise:  U2 is not updated.\n*\n*  JOBV1T  (input) CHARACTER\n*          = 'Y':      V1T is updated;\n*          otherwise:  V1T is not updated.\n*\n*  JOBV2T  (input) CHARACTER\n*          = 'Y':      V2T is updated;\n*          otherwise:  V2T is not updated.\n*\n*  TRANS   (input) CHARACTER\n*          = 'T':      X, U1, U2, V1T, and V2T are stored in row-major\n*                      order;\n*          otherwise:  X, U1, U2, V1T, and V2T are stored in column-\n*                      major order.\n*\n*  M       (input) INTEGER\n*          The number of rows and columns in X, the unitary matrix in\n*          bidiagonal-block form.\n*\n*  P       (input) INTEGER\n*          The number of rows in the top-left block of X. 0 <= P <= M.\n*\n*  Q       (input) INTEGER\n*          The number of columns in the top-left block of X.\n*          0 <= Q <= MIN(P,M-P,M-Q).\n*\n*  THETA   (input/output) REAL array, dimension (Q)\n*          On entry, the angles THETA(1),...,THETA(Q) that, along with\n*          PHI(1), ...,PHI(Q-1), define the matrix in bidiagonal-block\n*          form. On exit, the angles whose cosines and sines define the\n*          diagonal blocks in the CS decomposition.\n*\n*  PHI     (input/workspace) REAL array, dimension (Q-1)\n*          The angles PHI(1),...,PHI(Q-1) that, along with THETA(1),...,\n*          THETA(Q), define the matrix in bidiagonal-block form.\n*\n*  U1      (input/output) COMPLEX array, dimension (LDU1,P)\n*          On entry, an LDU1-by-P matrix. On exit, U1 is postmultiplied\n*          by the left singular vector matrix common to [ B11 ; 0 ] and\n*          [ B12 0 0 ; 0 -I 0 0 ].\n*\n*  LDU1    (input) INTEGER\n*          The leading dimension of the array U1.\n*\n*  U2      (input/output) COMPLEX array, dimension (LDU2,M-P)\n*          On entry, an LDU2-by-(M-P) matrix. On exit, U2 is\n*          postmultiplied by the left singular vector matrix common to\n*          [ B21 ; 0 ] and [ B22 0 0 ; 0 0 I ].\n*\n*  LDU2    (input) INTEGER\n*          The leading dimension of the array U2.\n*\n*  V1T     (input/output) COMPLEX array, dimension (LDV1T,Q)\n*          On entry, a LDV1T-by-Q matrix. On exit, V1T is premultiplied\n*          by the conjugate transpose of the right singular vector\n*          matrix common to [ B11 ; 0 ] and [ B21 ; 0 ].\n*\n*  LDV1T   (input) INTEGER\n*          The leading dimension of the array V1T.\n*\n*  V2T     (input/output) COMPLEX array, dimenison (LDV2T,M-Q)\n*          On entry, a LDV2T-by-(M-Q) matrix. On exit, V2T is\n*          premultiplied by the conjugate transpose of the right\n*          singular vector matrix common to [ B12 0 0 ; 0 -I 0 ] and\n*          [ B22 0 0 ; 0 0 I ].\n*\n*  LDV2T   (input) INTEGER\n*          The leading dimension of the array V2T.\n*\n*  B11D    (output) REAL array, dimension (Q)\n*          When CBBCSD converges, B11D contains the cosines of THETA(1),\n*          ..., THETA(Q). If CBBCSD fails to converge, then B11D\n*          contains the diagonal of the partially reduced top-left\n*          block.\n*\n*  B11E    (output) REAL array, dimension (Q-1)\n*          When CBBCSD converges, B11E contains zeros. If CBBCSD fails\n*          to converge, then B11E contains the superdiagonal of the\n*          partially reduced top-left block.\n*\n*  B12D    (output) REAL array, dimension (Q)\n*          When CBBCSD converges, B12D contains the negative sines of\n*          THETA(1), ..., THETA(Q). If CBBCSD fails to converge, then\n*          B12D contains the diagonal of the partially reduced top-right\n*          block.\n*\n*  B12E    (output) REAL array, dimension (Q-1)\n*          When CBBCSD converges, B12E contains zeros. If CBBCSD fails\n*          to converge, then B12E contains the subdiagonal of the\n*          partially reduced top-right block.\n*\n*  RWORK   (workspace) REAL array, dimension (MAX(1,LWORK))\n*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.\n*\n*  LRWORK  (input) INTEGER\n*          The dimension of the array RWORK. LRWORK >= MAX(1,8*Q).\n*\n*          If LRWORK = -1, then a workspace query is assumed; the\n*          routine only calculates the optimal size of the RWORK array,\n*          returns this value as the first entry of the work array, and\n*          no error message related to LRWORK is issued by XERBLA.\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit.\n*          < 0:  if INFO = -i, the i-th argument had an illegal value.\n*          > 0:  if CBBCSD did not converge, INFO specifies the number\n*                of nonzero entries in PHI, and B11D, B11E, etc.,\n*                contain the partially reduced matrix.\n*\n*  Reference\n*  =========\n*\n*  [1] Brian D. Sutton. Computing the complete CS decomposition. Numer.\n*      Algorithms, 50(1):33-65, 2009.\n*\n*  Internal Parameters\n*  ===================\n*\n*  TOLMUL  REAL, default = MAX(10,MIN(100,EPS**(-1/8)))\n*          TOLMUL controls the convergence criterion of the QR loop.\n*          Angles THETA(i), PHI(i) are rounded to 0 or PI/2 when they\n*          are within TOLMUL*EPS of either bound.\n*\n\n*  ===================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 13)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 13)", argc);
  rb_jobu1 = argv[0];
  rb_jobu2 = argv[1];
  rb_jobv1t = argv[2];
  rb_jobv2t = argv[3];
  rb_trans = argv[4];
  rb_m = argv[5];
  rb_theta = argv[6];
  rb_phi = argv[7];
  rb_u1 = argv[8];
  rb_u2 = argv[9];
  rb_v1t = argv[10];
  rb_v2t = argv[11];
  rb_lrwork = argv[12];

  if (!NA_IsNArray(rb_theta))
    rb_raise(rb_eArgError, "theta (7th argument) must be NArray");
  if (NA_RANK(rb_theta) != 1)
    rb_raise(rb_eArgError, "rank of theta (7th argument) must be %d", 1);
  q = NA_SHAPE0(rb_theta);
  if (NA_TYPE(rb_theta) != NA_SFLOAT)
    rb_theta = na_change_type(rb_theta, NA_SFLOAT);
  theta = NA_PTR_TYPE(rb_theta, real*);
  jobu1 = StringValueCStr(rb_jobu1)[0];
  trans = StringValueCStr(rb_trans)[0];
  m = NUM2INT(rb_m);
  jobu2 = StringValueCStr(rb_jobu2)[0];
  jobv1t = StringValueCStr(rb_jobv1t)[0];
  jobv2t = StringValueCStr(rb_jobv2t)[0];
  if (!NA_IsNArray(rb_u1))
    rb_raise(rb_eArgError, "u1 (9th argument) must be NArray");
  if (NA_RANK(rb_u1) != 2)
    rb_raise(rb_eArgError, "rank of u1 (9th argument) must be %d", 2);
  p = NA_SHAPE1(rb_u1);
  ldu1 = NA_SHAPE0(rb_u1);
  if (NA_TYPE(rb_u1) != NA_SCOMPLEX)
    rb_u1 = na_change_type(rb_u1, NA_SCOMPLEX);
  u1 = NA_PTR_TYPE(rb_u1, complex*);
  if (!NA_IsNArray(rb_v1t))
    rb_raise(rb_eArgError, "v1t (11th argument) must be NArray");
  if (NA_RANK(rb_v1t) != 2)
    rb_raise(rb_eArgError, "rank of v1t (11th argument) must be %d", 2);
  if (NA_SHAPE1(rb_v1t) != q)
    rb_raise(rb_eRuntimeError, "shape 1 of v1t must be the same as shape 0 of theta");
  ldv1t = NA_SHAPE0(rb_v1t);
  if (NA_TYPE(rb_v1t) != NA_SCOMPLEX)
    rb_v1t = na_change_type(rb_v1t, NA_SCOMPLEX);
  v1t = NA_PTR_TYPE(rb_v1t, complex*);
  if (!NA_IsNArray(rb_u2))
    rb_raise(rb_eArgError, "u2 (10th argument) must be NArray");
  if (NA_RANK(rb_u2) != 2)
    rb_raise(rb_eArgError, "rank of u2 (10th argument) must be %d", 2);
  if (NA_SHAPE1(rb_u2) != (m-p))
    rb_raise(rb_eRuntimeError, "shape 1 of u2 must be %d", m-p);
  ldu2 = NA_SHAPE0(rb_u2);
  if (NA_TYPE(rb_u2) != NA_SCOMPLEX)
    rb_u2 = na_change_type(rb_u2, NA_SCOMPLEX);
  u2 = NA_PTR_TYPE(rb_u2, complex*);
  if (!NA_IsNArray(rb_phi))
    rb_raise(rb_eArgError, "phi (8th argument) must be NArray");
  if (NA_RANK(rb_phi) != 1)
    rb_raise(rb_eArgError, "rank of phi (8th argument) must be %d", 1);
  if (NA_SHAPE0(rb_phi) != (q-1))
    rb_raise(rb_eRuntimeError, "shape 0 of phi must be %d", q-1);
  if (NA_TYPE(rb_phi) != NA_SFLOAT)
    rb_phi = na_change_type(rb_phi, NA_SFLOAT);
  phi = NA_PTR_TYPE(rb_phi, real*);
  lrwork = MAX(1,8*q);
  if (!NA_IsNArray(rb_v2t))
    rb_raise(rb_eArgError, "v2t (12th argument) must be NArray");
  if (NA_RANK(rb_v2t) != 2)
    rb_raise(rb_eArgError, "rank of v2t (12th argument) must be %d", 2);
  if (NA_SHAPE1(rb_v2t) != (m-q))
    rb_raise(rb_eRuntimeError, "shape 1 of v2t must be %d", m-q);
  ldv2t = NA_SHAPE0(rb_v2t);
  if (NA_TYPE(rb_v2t) != NA_SCOMPLEX)
    rb_v2t = na_change_type(rb_v2t, NA_SCOMPLEX);
  v2t = NA_PTR_TYPE(rb_v2t, complex*);
  {
    int shape[1];
    shape[0] = q;
    rb_b11d = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  b11d = NA_PTR_TYPE(rb_b11d, real*);
  {
    int shape[1];
    shape[0] = q-1;
    rb_b11e = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  b11e = NA_PTR_TYPE(rb_b11e, real*);
  {
    int shape[1];
    shape[0] = q;
    rb_b12d = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  b12d = NA_PTR_TYPE(rb_b12d, real*);
  {
    int shape[1];
    shape[0] = q-1;
    rb_b12e = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  b12e = NA_PTR_TYPE(rb_b12e, real*);
  {
    int shape[1];
    shape[0] = q;
    rb_b21d = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  b21d = NA_PTR_TYPE(rb_b21d, real*);
  {
    int shape[1];
    shape[0] = q-1;
    rb_b21e = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  b21e = NA_PTR_TYPE(rb_b21e, real*);
  {
    int shape[1];
    shape[0] = q;
    rb_b22d = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  b22d = NA_PTR_TYPE(rb_b22d, real*);
  {
    int shape[1];
    shape[0] = q-1;
    rb_b22e = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  b22e = NA_PTR_TYPE(rb_b22e, real*);
  {
    int shape[1];
    shape[0] = q;
    rb_theta_out__ = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  theta_out__ = NA_PTR_TYPE(rb_theta_out__, real*);
  MEMCPY(theta_out__, theta, real, NA_TOTAL(rb_theta));
  rb_theta = rb_theta_out__;
  theta = theta_out__;
  {
    int shape[2];
    shape[0] = ldu1;
    shape[1] = p;
    rb_u1_out__ = na_make_object(NA_SCOMPLEX, 2, shape, cNArray);
  }
  u1_out__ = NA_PTR_TYPE(rb_u1_out__, complex*);
  MEMCPY(u1_out__, u1, complex, NA_TOTAL(rb_u1));
  rb_u1 = rb_u1_out__;
  u1 = u1_out__;
  {
    int shape[2];
    shape[0] = ldu2;
    shape[1] = m-p;
    rb_u2_out__ = na_make_object(NA_SCOMPLEX, 2, shape, cNArray);
  }
  u2_out__ = NA_PTR_TYPE(rb_u2_out__, complex*);
  MEMCPY(u2_out__, u2, complex, NA_TOTAL(rb_u2));
  rb_u2 = rb_u2_out__;
  u2 = u2_out__;
  {
    int shape[2];
    shape[0] = ldv1t;
    shape[1] = q;
    rb_v1t_out__ = na_make_object(NA_SCOMPLEX, 2, shape, cNArray);
  }
  v1t_out__ = NA_PTR_TYPE(rb_v1t_out__, complex*);
  MEMCPY(v1t_out__, v1t, complex, NA_TOTAL(rb_v1t));
  rb_v1t = rb_v1t_out__;
  v1t = v1t_out__;
  {
    int shape[2];
    shape[0] = ldv2t;
    shape[1] = m-q;
    rb_v2t_out__ = na_make_object(NA_SCOMPLEX, 2, shape, cNArray);
  }
  v2t_out__ = NA_PTR_TYPE(rb_v2t_out__, complex*);
  MEMCPY(v2t_out__, v2t, complex, NA_TOTAL(rb_v2t));
  rb_v2t = rb_v2t_out__;
  v2t = v2t_out__;
  rwork = ALLOC_N(real, (MAX(1,lrwork)));

  cbbcsd_(&jobu1, &jobu2, &jobv1t, &jobv2t, &trans, &m, &p, &q, theta, phi, u1, &ldu1, u2, &ldu2, v1t, &ldv1t, v2t, &ldv2t, b11d, b11e, b12d, b12e, b21d, b21e, b22d, b22e, rwork, &lrwork, &info);

  free(rwork);
  rb_info = INT2NUM(info);
  return rb_ary_new3(14, rb_b11d, rb_b11e, rb_b12d, rb_b12e, rb_b21d, rb_b21e, rb_b22d, rb_b22e, rb_info, rb_theta, rb_u1, rb_u2, rb_v1t, rb_v2t);
}

void
init_lapack_cbbcsd(VALUE mLapack){
  rb_define_module_function(mLapack, "cbbcsd", rb_cbbcsd, -1);
}

#include "rb_lapack.h"

extern VOID dbbcsd_(char *jobu1, char *jobu2, char *jobv1t, char *jobv2t, char *trans, integer *m, integer *p, integer *q, doublereal *theta, doublereal *phi, doublereal *u1, integer *ldu1, doublereal *u2, integer *ldu2, doublereal *v1t, integer *ldv1t, doublereal *v2t, integer *ldv2t, doublereal *b11d, doublereal *b11e, doublereal *b12d, doublereal *b12e, doublereal *b21d, doublereal *b21e, doublereal *b22d, doublereal *b22e, doublereal *work, integer *lwork, integer *info);

static VALUE sHelp, sUsage;

static VALUE
rblapack_dbbcsd(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_jobu1;
  char jobu1; 
  VALUE rblapack_jobu2;
  char jobu2; 
  VALUE rblapack_jobv1t;
  char jobv1t; 
  VALUE rblapack_jobv2t;
  char jobv2t; 
  VALUE rblapack_trans;
  char trans; 
  VALUE rblapack_m;
  integer m; 
  VALUE rblapack_theta;
  doublereal *theta; 
  VALUE rblapack_phi;
  doublereal *phi; 
  VALUE rblapack_u1;
  doublereal *u1; 
  VALUE rblapack_u2;
  doublereal *u2; 
  VALUE rblapack_v1t;
  doublereal *v1t; 
  VALUE rblapack_v2t;
  doublereal *v2t; 
  VALUE rblapack_lwork;
  integer lwork; 
  VALUE rblapack_b11d;
  doublereal *b11d; 
  VALUE rblapack_b11e;
  doublereal *b11e; 
  VALUE rblapack_b12d;
  doublereal *b12d; 
  VALUE rblapack_b12e;
  doublereal *b12e; 
  VALUE rblapack_b21d;
  doublereal *b21d; 
  VALUE rblapack_b21e;
  doublereal *b21e; 
  VALUE rblapack_b22d;
  doublereal *b22d; 
  VALUE rblapack_b22e;
  doublereal *b22e; 
  VALUE rblapack_info;
  integer info; 
  VALUE rblapack_theta_out__;
  doublereal *theta_out__;
  VALUE rblapack_u1_out__;
  doublereal *u1_out__;
  VALUE rblapack_u2_out__;
  doublereal *u2_out__;
  VALUE rblapack_v1t_out__;
  doublereal *v1t_out__;
  VALUE rblapack_v2t_out__;
  doublereal *v2t_out__;
  doublereal *work;

  integer q;
  integer ldu1;
  integer p;
  integer ldu2;
  integer ldv1t;
  integer ldv2t;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  b11d, b11e, b12d, b12e, b21d, b21e, b22d, b22e, info, theta, u1, u2, v1t, v2t = NumRu::Lapack.dbbcsd( jobu1, jobu2, jobv1t, jobv2t, trans, m, theta, phi, u1, u2, v1t, v2t, lwork, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE DBBCSD( JOBU1, JOBU2, JOBV1T, JOBV2T, TRANS, M, P, Q, THETA, PHI, U1, LDU1, U2, LDU2, V1T, LDV1T, V2T, LDV2T, B11D, B11E, B12D, B12E, B21D, B21E, B22D, B22E, WORK, LWORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  DBBCSD computes the CS decomposition of an orthogonal matrix in\n*  bidiagonal-block form,\n*\n*\n*      [ B11 | B12 0  0 ]\n*      [  0  |  0 -I  0 ]\n*  X = [----------------]\n*      [ B21 | B22 0  0 ]\n*      [  0  |  0  0  I ]\n*\n*                                [  C | -S  0  0 ]\n*                    [ U1 |    ] [  0 |  0 -I  0 ] [ V1 |    ]**T\n*                  = [---------] [---------------] [---------]   .\n*                    [    | U2 ] [  S |  C  0  0 ] [    | V2 ]\n*                                [  0 |  0  0  I ]\n*\n*  X is M-by-M, its top-left block is P-by-Q, and Q must be no larger\n*  than P, M-P, or M-Q. (If Q is not the smallest index, then X must be\n*  transposed and/or permuted. This can be done in constant time using\n*  the TRANS and SIGNS options. See DORCSD for details.)\n*\n*  The bidiagonal matrices B11, B12, B21, and B22 are represented\n*  implicitly by angles THETA(1:Q) and PHI(1:Q-1).\n*\n*  The orthogonal matrices U1, U2, V1T, and V2T are input/output.\n*  The input matrices are pre- or post-multiplied by the appropriate\n*  singular vector matrices.\n*\n\n*  Arguments\n*  =========\n*\n*  JOBU1   (input) CHARACTER\n*          = 'Y':      U1 is updated;\n*          otherwise:  U1 is not updated.\n*\n*  JOBU2   (input) CHARACTER\n*          = 'Y':      U2 is updated;\n*          otherwise:  U2 is not updated.\n*\n*  JOBV1T  (input) CHARACTER\n*          = 'Y':      V1T is updated;\n*          otherwise:  V1T is not updated.\n*\n*  JOBV2T  (input) CHARACTER\n*          = 'Y':      V2T is updated;\n*          otherwise:  V2T is not updated.\n*\n*  TRANS   (input) CHARACTER\n*          = 'T':      X, U1, U2, V1T, and V2T are stored in row-major\n*                      order;\n*          otherwise:  X, U1, U2, V1T, and V2T are stored in column-\n*                      major order.\n*\n*  M       (input) INTEGER\n*          The number of rows and columns in X, the orthogonal matrix in\n*          bidiagonal-block form.\n*\n*  P       (input) INTEGER\n*          The number of rows in the top-left block of X. 0 <= P <= M.\n*\n*  Q       (input) INTEGER\n*          The number of columns in the top-left block of X.\n*          0 <= Q <= MIN(P,M-P,M-Q).\n*\n*  THETA   (input/output) DOUBLE PRECISION array, dimension (Q)\n*          On entry, the angles THETA(1),...,THETA(Q) that, along with\n*          PHI(1), ...,PHI(Q-1), define the matrix in bidiagonal-block\n*          form. On exit, the angles whose cosines and sines define the\n*          diagonal blocks in the CS decomposition.\n*\n*  PHI     (input/workspace) DOUBLE PRECISION array, dimension (Q-1)\n*          The angles PHI(1),...,PHI(Q-1) that, along with THETA(1),...,\n*          THETA(Q), define the matrix in bidiagonal-block form.\n*\n*  U1      (input/output) DOUBLE PRECISION array, dimension (LDU1,P)\n*          On entry, an LDU1-by-P matrix. On exit, U1 is postmultiplied\n*          by the left singular vector matrix common to [ B11 ; 0 ] and\n*          [ B12 0 0 ; 0 -I 0 0 ].\n*\n*  LDU1    (input) INTEGER\n*          The leading dimension of the array U1.\n*\n*  U2      (input/output) DOUBLE PRECISION array, dimension (LDU2,M-P)\n*          On entry, an LDU2-by-(M-P) matrix. On exit, U2 is\n*          postmultiplied by the left singular vector matrix common to\n*          [ B21 ; 0 ] and [ B22 0 0 ; 0 0 I ].\n*\n*  LDU2    (input) INTEGER\n*          The leading dimension of the array U2.\n*\n*  V1T     (input/output) DOUBLE PRECISION array, dimension (LDV1T,Q)\n*          On entry, a LDV1T-by-Q matrix. On exit, V1T is premultiplied\n*          by the transpose of the right singular vector\n*          matrix common to [ B11 ; 0 ] and [ B21 ; 0 ].\n*\n*  LDV1T   (input) INTEGER\n*          The leading dimension of the array V1T.\n*\n*  V2T     (input/output) DOUBLE PRECISION array, dimenison (LDV2T,M-Q)\n*          On entry, a LDV2T-by-(M-Q) matrix. On exit, V2T is\n*          premultiplied by the transpose of the right\n*          singular vector matrix common to [ B12 0 0 ; 0 -I 0 ] and\n*          [ B22 0 0 ; 0 0 I ].\n*\n*  LDV2T   (input) INTEGER\n*          The leading dimension of the array V2T.\n*\n*  B11D    (output) DOUBLE PRECISION array, dimension (Q)\n*          When DBBCSD converges, B11D contains the cosines of THETA(1),\n*          ..., THETA(Q). If DBBCSD fails to converge, then B11D\n*          contains the diagonal of the partially reduced top-left\n*          block.\n*\n*  B11E    (output) DOUBLE PRECISION array, dimension (Q-1)\n*          When DBBCSD converges, B11E contains zeros. If DBBCSD fails\n*          to converge, then B11E contains the superdiagonal of the\n*          partially reduced top-left block.\n*\n*  B12D    (output) DOUBLE PRECISION array, dimension (Q)\n*          When DBBCSD converges, B12D contains the negative sines of\n*          THETA(1), ..., THETA(Q). If DBBCSD fails to converge, then\n*          B12D contains the diagonal of the partially reduced top-right\n*          block.\n*\n*  B12E    (output) DOUBLE PRECISION array, dimension (Q-1)\n*          When DBBCSD converges, B12E contains zeros. If DBBCSD fails\n*          to converge, then B12E contains the subdiagonal of the\n*          partially reduced top-right block.\n*\n*  WORK    (workspace) DOUBLE PRECISION array, dimension (MAX(1,LWORK))\n*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.\n*\n*  LWORK   (input) INTEGER\n*          The dimension of the array WORK. LWORK >= MAX(1,8*Q).\n*\n*          If LWORK = -1, then a workspace query is assumed; the\n*          routine only calculates the optimal size of the WORK array,\n*          returns this value as the first entry of the work array, and\n*          no error message related to LWORK is issued by XERBLA.\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit.\n*          < 0:  if INFO = -i, the i-th argument had an illegal value.\n*          > 0:  if DBBCSD did not converge, INFO specifies the number\n*                of nonzero entries in PHI, and B11D, B11E, etc.,\n*                contain the partially reduced matrix.\n*\n*  Reference\n*  =========\n*\n*  [1] Brian D. Sutton. Computing the complete CS decomposition. Numer.\n*      Algorithms, 50(1):33-65, 2009.\n*\n*  Internal Parameters\n*  ===================\n*\n*  TOLMUL  DOUBLE PRECISION, default = MAX(10,MIN(100,EPS**(-1/8)))\n*          TOLMUL controls the convergence criterion of the QR loop.\n*          Angles THETA(i), PHI(i) are rounded to 0 or PI/2 when they\n*          are within TOLMUL*EPS of either bound.\n*\n\n*  ===================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  b11d, b11e, b12d, b12e, b21d, b21e, b22d, b22e, info, theta, u1, u2, v1t, v2t = NumRu::Lapack.dbbcsd( jobu1, jobu2, jobv1t, jobv2t, trans, m, theta, phi, u1, u2, v1t, v2t, lwork, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 13)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 13)", argc);
  rblapack_jobu1 = argv[0];
  rblapack_jobu2 = argv[1];
  rblapack_jobv1t = argv[2];
  rblapack_jobv2t = argv[3];
  rblapack_trans = argv[4];
  rblapack_m = argv[5];
  rblapack_theta = argv[6];
  rblapack_phi = argv[7];
  rblapack_u1 = argv[8];
  rblapack_u2 = argv[9];
  rblapack_v1t = argv[10];
  rblapack_v2t = argv[11];
  rblapack_lwork = argv[12];
  if (rb_options != Qnil) {
  }

  if (!NA_IsNArray(rblapack_theta))
    rb_raise(rb_eArgError, "theta (7th argument) must be NArray");
  if (NA_RANK(rblapack_theta) != 1)
    rb_raise(rb_eArgError, "rank of theta (7th argument) must be %d", 1);
  q = NA_SHAPE0(rblapack_theta);
  if (NA_TYPE(rblapack_theta) != NA_DFLOAT)
    rblapack_theta = na_change_type(rblapack_theta, NA_DFLOAT);
  theta = NA_PTR_TYPE(rblapack_theta, doublereal*);
  jobu1 = StringValueCStr(rblapack_jobu1)[0];
  trans = StringValueCStr(rblapack_trans)[0];
  m = NUM2INT(rblapack_m);
  jobu2 = StringValueCStr(rblapack_jobu2)[0];
  if (!NA_IsNArray(rblapack_v1t))
    rb_raise(rb_eArgError, "v1t (11th argument) must be NArray");
  if (NA_RANK(rblapack_v1t) != 2)
    rb_raise(rb_eArgError, "rank of v1t (11th argument) must be %d", 2);
  if (NA_SHAPE1(rblapack_v1t) != q)
    rb_raise(rb_eRuntimeError, "shape 1 of v1t must be the same as shape 0 of theta");
  ldv1t = NA_SHAPE0(rblapack_v1t);
  if (NA_TYPE(rblapack_v1t) != NA_DFLOAT)
    rblapack_v1t = na_change_type(rblapack_v1t, NA_DFLOAT);
  v1t = NA_PTR_TYPE(rblapack_v1t, doublereal*);
  jobv2t = StringValueCStr(rblapack_jobv2t)[0];
  if (!NA_IsNArray(rblapack_u1))
    rb_raise(rb_eArgError, "u1 (9th argument) must be NArray");
  if (NA_RANK(rblapack_u1) != 2)
    rb_raise(rb_eArgError, "rank of u1 (9th argument) must be %d", 2);
  p = NA_SHAPE1(rblapack_u1);
  ldu1 = NA_SHAPE0(rblapack_u1);
  if (NA_TYPE(rblapack_u1) != NA_DFLOAT)
    rblapack_u1 = na_change_type(rblapack_u1, NA_DFLOAT);
  u1 = NA_PTR_TYPE(rblapack_u1, doublereal*);
  jobv1t = StringValueCStr(rblapack_jobv1t)[0];
  lwork = NUM2INT(rblapack_lwork);
  if (!NA_IsNArray(rblapack_u2))
    rb_raise(rb_eArgError, "u2 (10th argument) must be NArray");
  if (NA_RANK(rblapack_u2) != 2)
    rb_raise(rb_eArgError, "rank of u2 (10th argument) must be %d", 2);
  if (NA_SHAPE1(rblapack_u2) != (m-p))
    rb_raise(rb_eRuntimeError, "shape 1 of u2 must be %d", m-p);
  ldu2 = NA_SHAPE0(rblapack_u2);
  if (NA_TYPE(rblapack_u2) != NA_DFLOAT)
    rblapack_u2 = na_change_type(rblapack_u2, NA_DFLOAT);
  u2 = NA_PTR_TYPE(rblapack_u2, doublereal*);
  if (!NA_IsNArray(rblapack_phi))
    rb_raise(rb_eArgError, "phi (8th argument) must be NArray");
  if (NA_RANK(rblapack_phi) != 1)
    rb_raise(rb_eArgError, "rank of phi (8th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_phi) != (q-1))
    rb_raise(rb_eRuntimeError, "shape 0 of phi must be %d", q-1);
  if (NA_TYPE(rblapack_phi) != NA_DFLOAT)
    rblapack_phi = na_change_type(rblapack_phi, NA_DFLOAT);
  phi = NA_PTR_TYPE(rblapack_phi, doublereal*);
  if (!NA_IsNArray(rblapack_v2t))
    rb_raise(rb_eArgError, "v2t (12th argument) must be NArray");
  if (NA_RANK(rblapack_v2t) != 2)
    rb_raise(rb_eArgError, "rank of v2t (12th argument) must be %d", 2);
  if (NA_SHAPE1(rblapack_v2t) != (m-q))
    rb_raise(rb_eRuntimeError, "shape 1 of v2t must be %d", m-q);
  ldv2t = NA_SHAPE0(rblapack_v2t);
  if (NA_TYPE(rblapack_v2t) != NA_DFLOAT)
    rblapack_v2t = na_change_type(rblapack_v2t, NA_DFLOAT);
  v2t = NA_PTR_TYPE(rblapack_v2t, doublereal*);
  {
    int shape[1];
    shape[0] = q;
    rblapack_b11d = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  b11d = NA_PTR_TYPE(rblapack_b11d, doublereal*);
  {
    int shape[1];
    shape[0] = q-1;
    rblapack_b11e = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  b11e = NA_PTR_TYPE(rblapack_b11e, doublereal*);
  {
    int shape[1];
    shape[0] = q;
    rblapack_b12d = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  b12d = NA_PTR_TYPE(rblapack_b12d, doublereal*);
  {
    int shape[1];
    shape[0] = q-1;
    rblapack_b12e = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  b12e = NA_PTR_TYPE(rblapack_b12e, doublereal*);
  {
    int shape[1];
    shape[0] = q;
    rblapack_b21d = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  b21d = NA_PTR_TYPE(rblapack_b21d, doublereal*);
  {
    int shape[1];
    shape[0] = q-1;
    rblapack_b21e = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  b21e = NA_PTR_TYPE(rblapack_b21e, doublereal*);
  {
    int shape[1];
    shape[0] = q;
    rblapack_b22d = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  b22d = NA_PTR_TYPE(rblapack_b22d, doublereal*);
  {
    int shape[1];
    shape[0] = q-1;
    rblapack_b22e = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  b22e = NA_PTR_TYPE(rblapack_b22e, doublereal*);
  {
    int shape[1];
    shape[0] = q;
    rblapack_theta_out__ = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  theta_out__ = NA_PTR_TYPE(rblapack_theta_out__, doublereal*);
  MEMCPY(theta_out__, theta, doublereal, NA_TOTAL(rblapack_theta));
  rblapack_theta = rblapack_theta_out__;
  theta = theta_out__;
  {
    int shape[2];
    shape[0] = ldu1;
    shape[1] = p;
    rblapack_u1_out__ = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  u1_out__ = NA_PTR_TYPE(rblapack_u1_out__, doublereal*);
  MEMCPY(u1_out__, u1, doublereal, NA_TOTAL(rblapack_u1));
  rblapack_u1 = rblapack_u1_out__;
  u1 = u1_out__;
  {
    int shape[2];
    shape[0] = ldu2;
    shape[1] = m-p;
    rblapack_u2_out__ = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  u2_out__ = NA_PTR_TYPE(rblapack_u2_out__, doublereal*);
  MEMCPY(u2_out__, u2, doublereal, NA_TOTAL(rblapack_u2));
  rblapack_u2 = rblapack_u2_out__;
  u2 = u2_out__;
  {
    int shape[2];
    shape[0] = ldv1t;
    shape[1] = q;
    rblapack_v1t_out__ = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  v1t_out__ = NA_PTR_TYPE(rblapack_v1t_out__, doublereal*);
  MEMCPY(v1t_out__, v1t, doublereal, NA_TOTAL(rblapack_v1t));
  rblapack_v1t = rblapack_v1t_out__;
  v1t = v1t_out__;
  {
    int shape[2];
    shape[0] = ldv2t;
    shape[1] = m-q;
    rblapack_v2t_out__ = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  v2t_out__ = NA_PTR_TYPE(rblapack_v2t_out__, doublereal*);
  MEMCPY(v2t_out__, v2t, doublereal, NA_TOTAL(rblapack_v2t));
  rblapack_v2t = rblapack_v2t_out__;
  v2t = v2t_out__;
  work = ALLOC_N(doublereal, (MAX(1,lwork)));

  dbbcsd_(&jobu1, &jobu2, &jobv1t, &jobv2t, &trans, &m, &p, &q, theta, phi, u1, &ldu1, u2, &ldu2, v1t, &ldv1t, v2t, &ldv2t, b11d, b11e, b12d, b12e, b21d, b21e, b22d, b22e, work, &lwork, &info);

  free(work);
  rblapack_info = INT2NUM(info);
  return rb_ary_new3(14, rblapack_b11d, rblapack_b11e, rblapack_b12d, rblapack_b12e, rblapack_b21d, rblapack_b21e, rblapack_b22d, rblapack_b22e, rblapack_info, rblapack_theta, rblapack_u1, rblapack_u2, rblapack_v1t, rblapack_v2t);
}

void
init_lapack_dbbcsd(VALUE mLapack){
  rb_define_module_function(mLapack, "dbbcsd", rblapack_dbbcsd, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}

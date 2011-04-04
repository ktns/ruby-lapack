#include "rb_lapack.h"

extern VOID stgsy2_(char *trans, integer *ijob, integer *m, integer *n, real *a, integer *lda, real *b, integer *ldb, real *c, integer *ldc, real *d, integer *ldd, real *e, integer *lde, real *f, integer *ldf, real *scale, real *rdsum, real *rdscal, integer *iwork, integer *pq, integer *info);

static VALUE
rb_stgsy2(int argc, VALUE *argv, VALUE self){
  VALUE rb_trans;
  char trans; 
  VALUE rb_ijob;
  integer ijob; 
  VALUE rb_a;
  real *a; 
  VALUE rb_b;
  real *b; 
  VALUE rb_c;
  real *c; 
  VALUE rb_d;
  real *d; 
  VALUE rb_e;
  real *e; 
  VALUE rb_f;
  real *f; 
  VALUE rb_rdsum;
  real rdsum; 
  VALUE rb_rdscal;
  real rdscal; 
  VALUE rb_scale;
  real scale; 
  VALUE rb_pq;
  integer pq; 
  VALUE rb_info;
  integer info; 
  VALUE rb_c_out__;
  real *c_out__;
  VALUE rb_f_out__;
  real *f_out__;
  integer *iwork;

  integer lda;
  integer m;
  integer ldb;
  integer n;
  integer ldc;
  integer ldd;
  integer lde;
  integer ldf;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  scale, pq, info, c, f, rdsum, rdscal = NumRu::Lapack.stgsy2( trans, ijob, a, b, c, d, e, f, rdsum, rdscal)\n    or\n  NumRu::Lapack.stgsy2  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE STGSY2( TRANS, IJOB, M, N, A, LDA, B, LDB, C, LDC, D, LDD, E, LDE, F, LDF, SCALE, RDSUM, RDSCAL, IWORK, PQ, INFO )\n\n*  Purpose\n*  =======\n*\n*  STGSY2 solves the generalized Sylvester equation:\n*\n*              A * R - L * B = scale * C                (1)\n*              D * R - L * E = scale * F,\n*\n*  using Level 1 and 2 BLAS. where R and L are unknown M-by-N matrices,\n*  (A, D), (B, E) and (C, F) are given matrix pairs of size M-by-M,\n*  N-by-N and M-by-N, respectively, with real entries. (A, D) and (B, E)\n*  must be in generalized Schur canonical form, i.e. A, B are upper\n*  quasi triangular and D, E are upper triangular. The solution (R, L)\n*  overwrites (C, F). 0 <= SCALE <= 1 is an output scaling factor\n*  chosen to avoid overflow.\n*\n*  In matrix notation solving equation (1) corresponds to solve\n*  Z*x = scale*b, where Z is defined as\n*\n*         Z = [ kron(In, A)  -kron(B', Im) ]             (2)\n*             [ kron(In, D)  -kron(E', Im) ],\n*\n*  Ik is the identity matrix of size k and X' is the transpose of X.\n*  kron(X, Y) is the Kronecker product between the matrices X and Y.\n*  In the process of solving (1), we solve a number of such systems\n*  where Dim(In), Dim(In) = 1 or 2.\n*\n*  If TRANS = 'T', solve the transposed system Z'*y = scale*b for y,\n*  which is equivalent to solve for R and L in\n*\n*              A' * R  + D' * L   = scale *  C           (3)\n*              R  * B' + L  * E'  = scale * -F\n*\n*  This case is used to compute an estimate of Dif[(A, D), (B, E)] =\n*  sigma_min(Z) using reverse communicaton with SLACON.\n*\n*  STGSY2 also (IJOB >= 1) contributes to the computation in STGSYL\n*  of an upper bound on the separation between to matrix pairs. Then\n*  the input (A, D), (B, E) are sub-pencils of the matrix pair in\n*  STGSYL. See STGSYL for details.\n*\n\n*  Arguments\n*  =========\n*\n*  TRANS   (input) CHARACTER*1\n*          = 'N', solve the generalized Sylvester equation (1).\n*          = 'T': solve the 'transposed' system (3).\n*\n*  IJOB    (input) INTEGER\n*          Specifies what kind of functionality to be performed.\n*          = 0: solve (1) only.\n*          = 1: A contribution from this subsystem to a Frobenius\n*               norm-based estimate of the separation between two matrix\n*               pairs is computed. (look ahead strategy is used).\n*          = 2: A contribution from this subsystem to a Frobenius\n*               norm-based estimate of the separation between two matrix\n*               pairs is computed. (SGECON on sub-systems is used.)\n*          Not referenced if TRANS = 'T'.\n*\n*  M       (input) INTEGER\n*          On entry, M specifies the order of A and D, and the row\n*          dimension of C, F, R and L.\n*\n*  N       (input) INTEGER\n*          On entry, N specifies the order of B and E, and the column\n*          dimension of C, F, R and L.\n*\n*  A       (input) REAL array, dimension (LDA, M)\n*          On entry, A contains an upper quasi triangular matrix.\n*\n*  LDA     (input) INTEGER\n*          The leading dimension of the matrix A. LDA >= max(1, M).\n*\n*  B       (input) REAL array, dimension (LDB, N)\n*          On entry, B contains an upper quasi triangular matrix.\n*\n*  LDB     (input) INTEGER\n*          The leading dimension of the matrix B. LDB >= max(1, N).\n*\n*  C       (input/output) REAL array, dimension (LDC, N)\n*          On entry, C contains the right-hand-side of the first matrix\n*          equation in (1).\n*          On exit, if IJOB = 0, C has been overwritten by the\n*          solution R.\n*\n*  LDC     (input) INTEGER\n*          The leading dimension of the matrix C. LDC >= max(1, M).\n*\n*  D       (input) REAL array, dimension (LDD, M)\n*          On entry, D contains an upper triangular matrix.\n*\n*  LDD     (input) INTEGER\n*          The leading dimension of the matrix D. LDD >= max(1, M).\n*\n*  E       (input) REAL array, dimension (LDE, N)\n*          On entry, E contains an upper triangular matrix.\n*\n*  LDE     (input) INTEGER\n*          The leading dimension of the matrix E. LDE >= max(1, N).\n*\n*  F       (input/output) REAL array, dimension (LDF, N)\n*          On entry, F contains the right-hand-side of the second matrix\n*          equation in (1).\n*          On exit, if IJOB = 0, F has been overwritten by the\n*          solution L.\n*\n*  LDF     (input) INTEGER\n*          The leading dimension of the matrix F. LDF >= max(1, M).\n*\n*  SCALE   (output) REAL\n*          On exit, 0 <= SCALE <= 1. If 0 < SCALE < 1, the solutions\n*          R and L (C and F on entry) will hold the solutions to a\n*          slightly perturbed system but the input matrices A, B, D and\n*          E have not been changed. If SCALE = 0, R and L will hold the\n*          solutions to the homogeneous system with C = F = 0. Normally,\n*          SCALE = 1.\n*\n*  RDSUM   (input/output) REAL\n*          On entry, the sum of squares of computed contributions to\n*          the Dif-estimate under computation by STGSYL, where the\n*          scaling factor RDSCAL (see below) has been factored out.\n*          On exit, the corresponding sum of squares updated with the\n*          contributions from the current sub-system.\n*          If TRANS = 'T' RDSUM is not touched.\n*          NOTE: RDSUM only makes sense when STGSY2 is called by STGSYL.\n*\n*  RDSCAL  (input/output) REAL\n*          On entry, scaling factor used to prevent overflow in RDSUM.\n*          On exit, RDSCAL is updated w.r.t. the current contributions\n*          in RDSUM.\n*          If TRANS = 'T', RDSCAL is not touched.\n*          NOTE: RDSCAL only makes sense when STGSY2 is called by\n*                STGSYL.\n*\n*  IWORK   (workspace) INTEGER array, dimension (M+N+2)\n*\n*  PQ      (output) INTEGER\n*          On exit, the number of subsystems (of size 2-by-2, 4-by-4 and\n*          8-by-8) solved by this routine.\n*\n*  INFO    (output) INTEGER\n*          On exit, if INFO is set to\n*            =0: Successful exit\n*            <0: If INFO = -i, the i-th argument had an illegal value.\n*            >0: The matrix pairs (A, D) and (B, E) have common or very\n*                close eigenvalues.\n*\n\n*  Further Details\n*  ===============\n*\n*  Based on contributions by\n*     Bo Kagstrom and Peter Poromaa, Department of Computing Science,\n*     Umea University, S-901 87 Umea, Sweden.\n*\n*  =====================================================================\n*  Replaced various illegal calls to SCOPY by calls to SLASET.\n*  Sven Hammarling, 27/5/02.\n*\n\n");
    return Qnil;
  }
  if (argc != 10)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 10)", argc);
  rb_trans = argv[0];
  rb_ijob = argv[1];
  rb_a = argv[2];
  rb_b = argv[3];
  rb_c = argv[4];
  rb_d = argv[5];
  rb_e = argv[6];
  rb_f = argv[7];
  rb_rdsum = argv[8];
  rb_rdscal = argv[9];

  rdscal = (real)NUM2DBL(rb_rdscal);
  ijob = NUM2INT(rb_ijob);
  if (!NA_IsNArray(rb_a))
    rb_raise(rb_eArgError, "a (3th argument) must be NArray");
  if (NA_RANK(rb_a) != 2)
    rb_raise(rb_eArgError, "rank of a (3th argument) must be %d", 2);
  m = NA_SHAPE1(rb_a);
  lda = NA_SHAPE0(rb_a);
  if (NA_TYPE(rb_a) != NA_SFLOAT)
    rb_a = na_change_type(rb_a, NA_SFLOAT);
  a = NA_PTR_TYPE(rb_a, real*);
  if (!NA_IsNArray(rb_b))
    rb_raise(rb_eArgError, "b (4th argument) must be NArray");
  if (NA_RANK(rb_b) != 2)
    rb_raise(rb_eArgError, "rank of b (4th argument) must be %d", 2);
  n = NA_SHAPE1(rb_b);
  ldb = NA_SHAPE0(rb_b);
  if (NA_TYPE(rb_b) != NA_SFLOAT)
    rb_b = na_change_type(rb_b, NA_SFLOAT);
  b = NA_PTR_TYPE(rb_b, real*);
  if (!NA_IsNArray(rb_c))
    rb_raise(rb_eArgError, "c (5th argument) must be NArray");
  if (NA_RANK(rb_c) != 2)
    rb_raise(rb_eArgError, "rank of c (5th argument) must be %d", 2);
  if (NA_SHAPE1(rb_c) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of c must be the same as shape 1 of b");
  ldc = NA_SHAPE0(rb_c);
  if (NA_TYPE(rb_c) != NA_SFLOAT)
    rb_c = na_change_type(rb_c, NA_SFLOAT);
  c = NA_PTR_TYPE(rb_c, real*);
  if (!NA_IsNArray(rb_d))
    rb_raise(rb_eArgError, "d (6th argument) must be NArray");
  if (NA_RANK(rb_d) != 2)
    rb_raise(rb_eArgError, "rank of d (6th argument) must be %d", 2);
  if (NA_SHAPE1(rb_d) != m)
    rb_raise(rb_eRuntimeError, "shape 1 of d must be the same as shape 1 of a");
  ldd = NA_SHAPE0(rb_d);
  if (NA_TYPE(rb_d) != NA_SFLOAT)
    rb_d = na_change_type(rb_d, NA_SFLOAT);
  d = NA_PTR_TYPE(rb_d, real*);
  rdsum = (real)NUM2DBL(rb_rdsum);
  if (!NA_IsNArray(rb_e))
    rb_raise(rb_eArgError, "e (7th argument) must be NArray");
  if (NA_RANK(rb_e) != 2)
    rb_raise(rb_eArgError, "rank of e (7th argument) must be %d", 2);
  if (NA_SHAPE1(rb_e) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of e must be the same as shape 1 of b");
  lde = NA_SHAPE0(rb_e);
  if (NA_TYPE(rb_e) != NA_SFLOAT)
    rb_e = na_change_type(rb_e, NA_SFLOAT);
  e = NA_PTR_TYPE(rb_e, real*);
  if (!NA_IsNArray(rb_f))
    rb_raise(rb_eArgError, "f (8th argument) must be NArray");
  if (NA_RANK(rb_f) != 2)
    rb_raise(rb_eArgError, "rank of f (8th argument) must be %d", 2);
  if (NA_SHAPE1(rb_f) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of f must be the same as shape 1 of b");
  ldf = NA_SHAPE0(rb_f);
  if (NA_TYPE(rb_f) != NA_SFLOAT)
    rb_f = na_change_type(rb_f, NA_SFLOAT);
  f = NA_PTR_TYPE(rb_f, real*);
  trans = StringValueCStr(rb_trans)[0];
  {
    int shape[2];
    shape[0] = ldc;
    shape[1] = n;
    rb_c_out__ = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  c_out__ = NA_PTR_TYPE(rb_c_out__, real*);
  MEMCPY(c_out__, c, real, NA_TOTAL(rb_c));
  rb_c = rb_c_out__;
  c = c_out__;
  {
    int shape[2];
    shape[0] = ldf;
    shape[1] = n;
    rb_f_out__ = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  f_out__ = NA_PTR_TYPE(rb_f_out__, real*);
  MEMCPY(f_out__, f, real, NA_TOTAL(rb_f));
  rb_f = rb_f_out__;
  f = f_out__;
  iwork = ALLOC_N(integer, (m+n+2));

  stgsy2_(&trans, &ijob, &m, &n, a, &lda, b, &ldb, c, &ldc, d, &ldd, e, &lde, f, &ldf, &scale, &rdsum, &rdscal, iwork, &pq, &info);

  free(iwork);
  rb_scale = rb_float_new((double)scale);
  rb_pq = INT2NUM(pq);
  rb_info = INT2NUM(info);
  rb_rdsum = rb_float_new((double)rdsum);
  rb_rdscal = rb_float_new((double)rdscal);
  return rb_ary_new3(7, rb_scale, rb_pq, rb_info, rb_c, rb_f, rb_rdsum, rb_rdscal);
}

void
init_lapack_stgsy2(VALUE mLapack){
  rb_define_module_function(mLapack, "stgsy2", rb_stgsy2, -1);
}

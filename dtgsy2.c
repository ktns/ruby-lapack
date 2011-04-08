#include "rb_lapack.h"

extern VOID dtgsy2_(char *trans, integer *ijob, integer *m, integer *n, doublereal *a, integer *lda, doublereal *b, integer *ldb, doublereal *c, integer *ldc, doublereal *d, integer *ldd, doublereal *e, integer *lde, doublereal *f, integer *ldf, doublereal *scale, doublereal *rdsum, doublereal *rdscal, integer *iwork, integer *pq, integer *info);

static VALUE sHelp, sUsage;

static VALUE
rblapack_dtgsy2(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_trans;
  char trans; 
  VALUE rblapack_ijob;
  integer ijob; 
  VALUE rblapack_a;
  doublereal *a; 
  VALUE rblapack_b;
  doublereal *b; 
  VALUE rblapack_c;
  doublereal *c; 
  VALUE rblapack_d;
  doublereal *d; 
  VALUE rblapack_e;
  doublereal *e; 
  VALUE rblapack_f;
  doublereal *f; 
  VALUE rblapack_rdsum;
  doublereal rdsum; 
  VALUE rblapack_rdscal;
  doublereal rdscal; 
  VALUE rblapack_scale;
  doublereal scale; 
  VALUE rblapack_pq;
  integer pq; 
  VALUE rblapack_info;
  integer info; 
  VALUE rblapack_c_out__;
  doublereal *c_out__;
  VALUE rblapack_f_out__;
  doublereal *f_out__;
  integer *iwork;

  integer lda;
  integer m;
  integer ldb;
  integer n;
  integer ldc;
  integer ldd;
  integer lde;
  integer ldf;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  scale, pq, info, c, f, rdsum, rdscal = NumRu::Lapack.dtgsy2( trans, ijob, a, b, c, d, e, f, rdsum, rdscal, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE DTGSY2( TRANS, IJOB, M, N, A, LDA, B, LDB, C, LDC, D, LDD, E, LDE, F, LDF, SCALE, RDSUM, RDSCAL, IWORK, PQ, INFO )\n\n*  Purpose\n*  =======\n*\n*  DTGSY2 solves the generalized Sylvester equation:\n*\n*              A * R - L * B = scale * C                (1)\n*              D * R - L * E = scale * F,\n*\n*  using Level 1 and 2 BLAS. where R and L are unknown M-by-N matrices,\n*  (A, D), (B, E) and (C, F) are given matrix pairs of size M-by-M,\n*  N-by-N and M-by-N, respectively, with real entries. (A, D) and (B, E)\n*  must be in generalized Schur canonical form, i.e. A, B are upper\n*  quasi triangular and D, E are upper triangular. The solution (R, L)\n*  overwrites (C, F). 0 <= SCALE <= 1 is an output scaling factor\n*  chosen to avoid overflow.\n*\n*  In matrix notation solving equation (1) corresponds to solve\n*  Z*x = scale*b, where Z is defined as\n*\n*         Z = [ kron(In, A)  -kron(B', Im) ]             (2)\n*             [ kron(In, D)  -kron(E', Im) ],\n*\n*  Ik is the identity matrix of size k and X' is the transpose of X.\n*  kron(X, Y) is the Kronecker product between the matrices X and Y.\n*  In the process of solving (1), we solve a number of such systems\n*  where Dim(In), Dim(In) = 1 or 2.\n*\n*  If TRANS = 'T', solve the transposed system Z'*y = scale*b for y,\n*  which is equivalent to solve for R and L in\n*\n*              A' * R  + D' * L   = scale *  C           (3)\n*              R  * B' + L  * E'  = scale * -F\n*\n*  This case is used to compute an estimate of Dif[(A, D), (B, E)] =\n*  sigma_min(Z) using reverse communicaton with DLACON.\n*\n*  DTGSY2 also (IJOB >= 1) contributes to the computation in DTGSYL\n*  of an upper bound on the separation between to matrix pairs. Then\n*  the input (A, D), (B, E) are sub-pencils of the matrix pair in\n*  DTGSYL. See DTGSYL for details.\n*\n\n*  Arguments\n*  =========\n*\n*  TRANS   (input) CHARACTER*1\n*          = 'N', solve the generalized Sylvester equation (1).\n*          = 'T': solve the 'transposed' system (3).\n*\n*  IJOB    (input) INTEGER\n*          Specifies what kind of functionality to be performed.\n*          = 0: solve (1) only.\n*          = 1: A contribution from this subsystem to a Frobenius\n*               norm-based estimate of the separation between two matrix\n*               pairs is computed. (look ahead strategy is used).\n*          = 2: A contribution from this subsystem to a Frobenius\n*               norm-based estimate of the separation between two matrix\n*               pairs is computed. (DGECON on sub-systems is used.)\n*          Not referenced if TRANS = 'T'.\n*\n*  M       (input) INTEGER\n*          On entry, M specifies the order of A and D, and the row\n*          dimension of C, F, R and L.\n*\n*  N       (input) INTEGER\n*          On entry, N specifies the order of B and E, and the column\n*          dimension of C, F, R and L.\n*\n*  A       (input) DOUBLE PRECISION array, dimension (LDA, M)\n*          On entry, A contains an upper quasi triangular matrix.\n*\n*  LDA     (input) INTEGER\n*          The leading dimension of the matrix A. LDA >= max(1, M).\n*\n*  B       (input) DOUBLE PRECISION array, dimension (LDB, N)\n*          On entry, B contains an upper quasi triangular matrix.\n*\n*  LDB     (input) INTEGER\n*          The leading dimension of the matrix B. LDB >= max(1, N).\n*\n*  C       (input/output) DOUBLE PRECISION array, dimension (LDC, N)\n*          On entry, C contains the right-hand-side of the first matrix\n*          equation in (1).\n*          On exit, if IJOB = 0, C has been overwritten by the\n*          solution R.\n*\n*  LDC     (input) INTEGER\n*          The leading dimension of the matrix C. LDC >= max(1, M).\n*\n*  D       (input) DOUBLE PRECISION array, dimension (LDD, M)\n*          On entry, D contains an upper triangular matrix.\n*\n*  LDD     (input) INTEGER\n*          The leading dimension of the matrix D. LDD >= max(1, M).\n*\n*  E       (input) DOUBLE PRECISION array, dimension (LDE, N)\n*          On entry, E contains an upper triangular matrix.\n*\n*  LDE     (input) INTEGER\n*          The leading dimension of the matrix E. LDE >= max(1, N).\n*\n*  F       (input/output) DOUBLE PRECISION array, dimension (LDF, N)\n*          On entry, F contains the right-hand-side of the second matrix\n*          equation in (1).\n*          On exit, if IJOB = 0, F has been overwritten by the\n*          solution L.\n*\n*  LDF     (input) INTEGER\n*          The leading dimension of the matrix F. LDF >= max(1, M).\n*\n*  SCALE   (output) DOUBLE PRECISION\n*          On exit, 0 <= SCALE <= 1. If 0 < SCALE < 1, the solutions\n*          R and L (C and F on entry) will hold the solutions to a\n*          slightly perturbed system but the input matrices A, B, D and\n*          E have not been changed. If SCALE = 0, R and L will hold the\n*          solutions to the homogeneous system with C = F = 0. Normally,\n*          SCALE = 1.\n*\n*  RDSUM   (input/output) DOUBLE PRECISION\n*          On entry, the sum of squares of computed contributions to\n*          the Dif-estimate under computation by DTGSYL, where the\n*          scaling factor RDSCAL (see below) has been factored out.\n*          On exit, the corresponding sum of squares updated with the\n*          contributions from the current sub-system.\n*          If TRANS = 'T' RDSUM is not touched.\n*          NOTE: RDSUM only makes sense when DTGSY2 is called by DTGSYL.\n*\n*  RDSCAL  (input/output) DOUBLE PRECISION\n*          On entry, scaling factor used to prevent overflow in RDSUM.\n*          On exit, RDSCAL is updated w.r.t. the current contributions\n*          in RDSUM.\n*          If TRANS = 'T', RDSCAL is not touched.\n*          NOTE: RDSCAL only makes sense when DTGSY2 is called by\n*                DTGSYL.\n*\n*  IWORK   (workspace) INTEGER array, dimension (M+N+2)\n*\n*  PQ      (output) INTEGER\n*          On exit, the number of subsystems (of size 2-by-2, 4-by-4 and\n*          8-by-8) solved by this routine.\n*\n*  INFO    (output) INTEGER\n*          On exit, if INFO is set to\n*            =0: Successful exit\n*            <0: If INFO = -i, the i-th argument had an illegal value.\n*            >0: The matrix pairs (A, D) and (B, E) have common or very\n*                close eigenvalues.\n*\n\n*  Further Details\n*  ===============\n*\n*  Based on contributions by\n*     Bo Kagstrom and Peter Poromaa, Department of Computing Science,\n*     Umea University, S-901 87 Umea, Sweden.\n*\n*  =====================================================================\n*  Replaced various illegal calls to DCOPY by calls to DLASET.\n*  Sven Hammarling, 27/5/02.\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  scale, pq, info, c, f, rdsum, rdscal = NumRu::Lapack.dtgsy2( trans, ijob, a, b, c, d, e, f, rdsum, rdscal, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 10)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 10)", argc);
  rblapack_trans = argv[0];
  rblapack_ijob = argv[1];
  rblapack_a = argv[2];
  rblapack_b = argv[3];
  rblapack_c = argv[4];
  rblapack_d = argv[5];
  rblapack_e = argv[6];
  rblapack_f = argv[7];
  rblapack_rdsum = argv[8];
  rblapack_rdscal = argv[9];
  if (rb_options != Qnil) {
  }

  rdscal = NUM2DBL(rblapack_rdscal);
  ijob = NUM2INT(rblapack_ijob);
  if (!NA_IsNArray(rblapack_a))
    rb_raise(rb_eArgError, "a (3th argument) must be NArray");
  if (NA_RANK(rblapack_a) != 2)
    rb_raise(rb_eArgError, "rank of a (3th argument) must be %d", 2);
  m = NA_SHAPE1(rblapack_a);
  lda = NA_SHAPE0(rblapack_a);
  if (NA_TYPE(rblapack_a) != NA_DFLOAT)
    rblapack_a = na_change_type(rblapack_a, NA_DFLOAT);
  a = NA_PTR_TYPE(rblapack_a, doublereal*);
  if (!NA_IsNArray(rblapack_b))
    rb_raise(rb_eArgError, "b (4th argument) must be NArray");
  if (NA_RANK(rblapack_b) != 2)
    rb_raise(rb_eArgError, "rank of b (4th argument) must be %d", 2);
  n = NA_SHAPE1(rblapack_b);
  ldb = NA_SHAPE0(rblapack_b);
  if (NA_TYPE(rblapack_b) != NA_DFLOAT)
    rblapack_b = na_change_type(rblapack_b, NA_DFLOAT);
  b = NA_PTR_TYPE(rblapack_b, doublereal*);
  if (!NA_IsNArray(rblapack_c))
    rb_raise(rb_eArgError, "c (5th argument) must be NArray");
  if (NA_RANK(rblapack_c) != 2)
    rb_raise(rb_eArgError, "rank of c (5th argument) must be %d", 2);
  if (NA_SHAPE1(rblapack_c) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of c must be the same as shape 1 of b");
  ldc = NA_SHAPE0(rblapack_c);
  if (NA_TYPE(rblapack_c) != NA_DFLOAT)
    rblapack_c = na_change_type(rblapack_c, NA_DFLOAT);
  c = NA_PTR_TYPE(rblapack_c, doublereal*);
  if (!NA_IsNArray(rblapack_d))
    rb_raise(rb_eArgError, "d (6th argument) must be NArray");
  if (NA_RANK(rblapack_d) != 2)
    rb_raise(rb_eArgError, "rank of d (6th argument) must be %d", 2);
  if (NA_SHAPE1(rblapack_d) != m)
    rb_raise(rb_eRuntimeError, "shape 1 of d must be the same as shape 1 of a");
  ldd = NA_SHAPE0(rblapack_d);
  if (NA_TYPE(rblapack_d) != NA_DFLOAT)
    rblapack_d = na_change_type(rblapack_d, NA_DFLOAT);
  d = NA_PTR_TYPE(rblapack_d, doublereal*);
  rdsum = NUM2DBL(rblapack_rdsum);
  if (!NA_IsNArray(rblapack_e))
    rb_raise(rb_eArgError, "e (7th argument) must be NArray");
  if (NA_RANK(rblapack_e) != 2)
    rb_raise(rb_eArgError, "rank of e (7th argument) must be %d", 2);
  if (NA_SHAPE1(rblapack_e) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of e must be the same as shape 1 of b");
  lde = NA_SHAPE0(rblapack_e);
  if (NA_TYPE(rblapack_e) != NA_DFLOAT)
    rblapack_e = na_change_type(rblapack_e, NA_DFLOAT);
  e = NA_PTR_TYPE(rblapack_e, doublereal*);
  if (!NA_IsNArray(rblapack_f))
    rb_raise(rb_eArgError, "f (8th argument) must be NArray");
  if (NA_RANK(rblapack_f) != 2)
    rb_raise(rb_eArgError, "rank of f (8th argument) must be %d", 2);
  if (NA_SHAPE1(rblapack_f) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of f must be the same as shape 1 of b");
  ldf = NA_SHAPE0(rblapack_f);
  if (NA_TYPE(rblapack_f) != NA_DFLOAT)
    rblapack_f = na_change_type(rblapack_f, NA_DFLOAT);
  f = NA_PTR_TYPE(rblapack_f, doublereal*);
  trans = StringValueCStr(rblapack_trans)[0];
  {
    int shape[2];
    shape[0] = ldc;
    shape[1] = n;
    rblapack_c_out__ = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  c_out__ = NA_PTR_TYPE(rblapack_c_out__, doublereal*);
  MEMCPY(c_out__, c, doublereal, NA_TOTAL(rblapack_c));
  rblapack_c = rblapack_c_out__;
  c = c_out__;
  {
    int shape[2];
    shape[0] = ldf;
    shape[1] = n;
    rblapack_f_out__ = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  f_out__ = NA_PTR_TYPE(rblapack_f_out__, doublereal*);
  MEMCPY(f_out__, f, doublereal, NA_TOTAL(rblapack_f));
  rblapack_f = rblapack_f_out__;
  f = f_out__;
  iwork = ALLOC_N(integer, (m+n+2));

  dtgsy2_(&trans, &ijob, &m, &n, a, &lda, b, &ldb, c, &ldc, d, &ldd, e, &lde, f, &ldf, &scale, &rdsum, &rdscal, iwork, &pq, &info);

  free(iwork);
  rblapack_scale = rb_float_new((double)scale);
  rblapack_pq = INT2NUM(pq);
  rblapack_info = INT2NUM(info);
  rblapack_rdsum = rb_float_new((double)rdsum);
  rblapack_rdscal = rb_float_new((double)rdscal);
  return rb_ary_new3(7, rblapack_scale, rblapack_pq, rblapack_info, rblapack_c, rblapack_f, rblapack_rdsum, rblapack_rdscal);
}

void
init_lapack_dtgsy2(VALUE mLapack){
  rb_define_module_function(mLapack, "dtgsy2", rblapack_dtgsy2, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}

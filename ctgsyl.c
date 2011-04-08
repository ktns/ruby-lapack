#include "rb_lapack.h"

extern VOID ctgsyl_(char *trans, integer *ijob, integer *m, integer *n, complex *a, integer *lda, complex *b, integer *ldb, complex *c, integer *ldc, complex *d, integer *ldd, complex *e, integer *lde, complex *f, integer *ldf, real *scale, real *dif, complex *work, integer *lwork, integer *iwork, integer *info);

static VALUE sHelp, sUsage;

static VALUE
rblapack_ctgsyl(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_trans;
  char trans; 
  VALUE rblapack_ijob;
  integer ijob; 
  VALUE rblapack_a;
  complex *a; 
  VALUE rblapack_b;
  complex *b; 
  VALUE rblapack_c;
  complex *c; 
  VALUE rblapack_d;
  complex *d; 
  VALUE rblapack_e;
  complex *e; 
  VALUE rblapack_f;
  complex *f; 
  VALUE rblapack_lwork;
  integer lwork; 
  VALUE rblapack_scale;
  real scale; 
  VALUE rblapack_dif;
  real dif; 
  VALUE rblapack_work;
  complex *work; 
  VALUE rblapack_info;
  integer info; 
  VALUE rblapack_c_out__;
  complex *c_out__;
  VALUE rblapack_f_out__;
  complex *f_out__;
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
      printf("%s\n", "USAGE:\n  scale, dif, work, info, c, f = NumRu::Lapack.ctgsyl( trans, ijob, a, b, c, d, e, f, lwork, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE CTGSYL( TRANS, IJOB, M, N, A, LDA, B, LDB, C, LDC, D, LDD, E, LDE, F, LDF, SCALE, DIF, WORK, LWORK, IWORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  CTGSYL solves the generalized Sylvester equation:\n*\n*              A * R - L * B = scale * C            (1)\n*              D * R - L * E = scale * F\n*\n*  where R and L are unknown m-by-n matrices, (A, D), (B, E) and\n*  (C, F) are given matrix pairs of size m-by-m, n-by-n and m-by-n,\n*  respectively, with complex entries. A, B, D and E are upper\n*  triangular (i.e., (A,D) and (B,E) in generalized Schur form).\n*\n*  The solution (R, L) overwrites (C, F). 0 <= SCALE <= 1\n*  is an output scaling factor chosen to avoid overflow.\n*\n*  In matrix notation (1) is equivalent to solve Zx = scale*b, where Z\n*  is defined as\n*\n*         Z = [ kron(In, A)  -kron(B', Im) ]        (2)\n*             [ kron(In, D)  -kron(E', Im) ],\n*\n*  Here Ix is the identity matrix of size x and X' is the conjugate\n*  transpose of X. Kron(X, Y) is the Kronecker product between the\n*  matrices X and Y.\n*\n*  If TRANS = 'C', y in the conjugate transposed system Z'*y = scale*b\n*  is solved for, which is equivalent to solve for R and L in\n*\n*              A' * R + D' * L = scale * C           (3)\n*              R * B' + L * E' = scale * -F\n*\n*  This case (TRANS = 'C') is used to compute an one-norm-based estimate\n*  of Dif[(A,D), (B,E)], the separation between the matrix pairs (A,D)\n*  and (B,E), using CLACON.\n*\n*  If IJOB >= 1, CTGSYL computes a Frobenius norm-based estimate of\n*  Dif[(A,D),(B,E)]. That is, the reciprocal of a lower bound on the\n*  reciprocal of the smallest singular value of Z.\n*\n*  This is a level-3 BLAS algorithm.\n*\n\n*  Arguments\n*  =========\n*\n*  TRANS   (input) CHARACTER*1\n*          = 'N': solve the generalized sylvester equation (1).\n*          = 'C': solve the \"conjugate transposed\" system (3).\n*\n*  IJOB    (input) INTEGER\n*          Specifies what kind of functionality to be performed.\n*          =0: solve (1) only.\n*          =1: The functionality of 0 and 3.\n*          =2: The functionality of 0 and 4.\n*          =3: Only an estimate of Dif[(A,D), (B,E)] is computed.\n*              (look ahead strategy is used).\n*          =4: Only an estimate of Dif[(A,D), (B,E)] is computed.\n*              (CGECON on sub-systems is used).\n*          Not referenced if TRANS = 'C'.\n*\n*  M       (input) INTEGER\n*          The order of the matrices A and D, and the row dimension of\n*          the matrices C, F, R and L.\n*\n*  N       (input) INTEGER\n*          The order of the matrices B and E, and the column dimension\n*          of the matrices C, F, R and L.\n*\n*  A       (input) COMPLEX array, dimension (LDA, M)\n*          The upper triangular matrix A.\n*\n*  LDA     (input) INTEGER\n*          The leading dimension of the array A. LDA >= max(1, M).\n*\n*  B       (input) COMPLEX array, dimension (LDB, N)\n*          The upper triangular matrix B.\n*\n*  LDB     (input) INTEGER\n*          The leading dimension of the array B. LDB >= max(1, N).\n*\n*  C       (input/output) COMPLEX array, dimension (LDC, N)\n*          On entry, C contains the right-hand-side of the first matrix\n*          equation in (1) or (3).\n*          On exit, if IJOB = 0, 1 or 2, C has been overwritten by\n*          the solution R. If IJOB = 3 or 4 and TRANS = 'N', C holds R,\n*          the solution achieved during the computation of the\n*          Dif-estimate.\n*\n*  LDC     (input) INTEGER\n*          The leading dimension of the array C. LDC >= max(1, M).\n*\n*  D       (input) COMPLEX array, dimension (LDD, M)\n*          The upper triangular matrix D.\n*\n*  LDD     (input) INTEGER\n*          The leading dimension of the array D. LDD >= max(1, M).\n*\n*  E       (input) COMPLEX array, dimension (LDE, N)\n*          The upper triangular matrix E.\n*\n*  LDE     (input) INTEGER\n*          The leading dimension of the array E. LDE >= max(1, N).\n*\n*  F       (input/output) COMPLEX array, dimension (LDF, N)\n*          On entry, F contains the right-hand-side of the second matrix\n*          equation in (1) or (3).\n*          On exit, if IJOB = 0, 1 or 2, F has been overwritten by\n*          the solution L. If IJOB = 3 or 4 and TRANS = 'N', F holds L,\n*          the solution achieved during the computation of the\n*          Dif-estimate.\n*\n*  LDF     (input) INTEGER\n*          The leading dimension of the array F. LDF >= max(1, M).\n*\n*  DIF     (output) REAL\n*          On exit DIF is the reciprocal of a lower bound of the\n*          reciprocal of the Dif-function, i.e. DIF is an upper bound of\n*          Dif[(A,D), (B,E)] = sigma-min(Z), where Z as in (2).\n*          IF IJOB = 0 or TRANS = 'C', DIF is not referenced.\n*\n*  SCALE   (output) REAL\n*          On exit SCALE is the scaling factor in (1) or (3).\n*          If 0 < SCALE < 1, C and F hold the solutions R and L, resp.,\n*          to a slightly perturbed system but the input matrices A, B,\n*          D and E have not been changed. If SCALE = 0, R and L will\n*          hold the solutions to the homogenious system with C = F = 0.\n*\n*  WORK    (workspace/output) COMPLEX array, dimension (MAX(1,LWORK))\n*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.\n*\n*  LWORK   (input) INTEGER\n*          The dimension of the array WORK. LWORK > = 1.\n*          If IJOB = 1 or 2 and TRANS = 'N', LWORK >= max(1,2*M*N).\n*\n*          If LWORK = -1, then a workspace query is assumed; the routine\n*          only calculates the optimal size of the WORK array, returns\n*          this value as the first entry of the WORK array, and no error\n*          message related to LWORK is issued by XERBLA.\n*\n*  IWORK   (workspace) INTEGER array, dimension (M+N+2)\n*\n*  INFO    (output) INTEGER\n*            =0: successful exit\n*            <0: If INFO = -i, the i-th argument had an illegal value.\n*            >0: (A, D) and (B, E) have common or very close\n*                eigenvalues.\n*\n\n*  Further Details\n*  ===============\n*\n*  Based on contributions by\n*     Bo Kagstrom and Peter Poromaa, Department of Computing Science,\n*     Umea University, S-901 87 Umea, Sweden.\n*\n*  [1] B. Kagstrom and P. Poromaa, LAPACK-Style Algorithms and Software\n*      for Solving the Generalized Sylvester Equation and Estimating the\n*      Separation between Regular Matrix Pairs, Report UMINF - 93.23,\n*      Department of Computing Science, Umea University, S-901 87 Umea,\n*      Sweden, December 1993, Revised April 1994, Also as LAPACK Working\n*      Note 75.  To appear in ACM Trans. on Math. Software, Vol 22,\n*      No 1, 1996.\n*\n*  [2] B. Kagstrom, A Perturbation Analysis of the Generalized Sylvester\n*      Equation (AR - LB, DR - LE ) = (C, F), SIAM J. Matrix Anal.\n*      Appl., 15(4):1045-1060, 1994.\n*\n*  [3] B. Kagstrom and L. Westin, Generalized Schur Methods with\n*      Condition Estimators for Solving the Generalized Sylvester\n*      Equation, IEEE Transactions on Automatic Control, Vol. 34, No. 7,\n*      July 1989, pp 745-751.\n*\n*  =====================================================================\n*  Replaced various illegal calls to CCOPY by calls to CLASET.\n*  Sven Hammarling, 1/5/02.\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  scale, dif, work, info, c, f = NumRu::Lapack.ctgsyl( trans, ijob, a, b, c, d, e, f, lwork, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 9)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 9)", argc);
  rblapack_trans = argv[0];
  rblapack_ijob = argv[1];
  rblapack_a = argv[2];
  rblapack_b = argv[3];
  rblapack_c = argv[4];
  rblapack_d = argv[5];
  rblapack_e = argv[6];
  rblapack_f = argv[7];
  rblapack_lwork = argv[8];
  if (rb_options != Qnil) {
  }

  ijob = NUM2INT(rblapack_ijob);
  if (!NA_IsNArray(rblapack_a))
    rb_raise(rb_eArgError, "a (3th argument) must be NArray");
  if (NA_RANK(rblapack_a) != 2)
    rb_raise(rb_eArgError, "rank of a (3th argument) must be %d", 2);
  m = NA_SHAPE1(rblapack_a);
  lda = NA_SHAPE0(rblapack_a);
  if (NA_TYPE(rblapack_a) != NA_SCOMPLEX)
    rblapack_a = na_change_type(rblapack_a, NA_SCOMPLEX);
  a = NA_PTR_TYPE(rblapack_a, complex*);
  if (!NA_IsNArray(rblapack_b))
    rb_raise(rb_eArgError, "b (4th argument) must be NArray");
  if (NA_RANK(rblapack_b) != 2)
    rb_raise(rb_eArgError, "rank of b (4th argument) must be %d", 2);
  n = NA_SHAPE1(rblapack_b);
  ldb = NA_SHAPE0(rblapack_b);
  if (NA_TYPE(rblapack_b) != NA_SCOMPLEX)
    rblapack_b = na_change_type(rblapack_b, NA_SCOMPLEX);
  b = NA_PTR_TYPE(rblapack_b, complex*);
  if (!NA_IsNArray(rblapack_c))
    rb_raise(rb_eArgError, "c (5th argument) must be NArray");
  if (NA_RANK(rblapack_c) != 2)
    rb_raise(rb_eArgError, "rank of c (5th argument) must be %d", 2);
  if (NA_SHAPE1(rblapack_c) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of c must be the same as shape 1 of b");
  ldc = NA_SHAPE0(rblapack_c);
  if (NA_TYPE(rblapack_c) != NA_SCOMPLEX)
    rblapack_c = na_change_type(rblapack_c, NA_SCOMPLEX);
  c = NA_PTR_TYPE(rblapack_c, complex*);
  if (!NA_IsNArray(rblapack_d))
    rb_raise(rb_eArgError, "d (6th argument) must be NArray");
  if (NA_RANK(rblapack_d) != 2)
    rb_raise(rb_eArgError, "rank of d (6th argument) must be %d", 2);
  if (NA_SHAPE1(rblapack_d) != m)
    rb_raise(rb_eRuntimeError, "shape 1 of d must be the same as shape 1 of a");
  ldd = NA_SHAPE0(rblapack_d);
  if (NA_TYPE(rblapack_d) != NA_SCOMPLEX)
    rblapack_d = na_change_type(rblapack_d, NA_SCOMPLEX);
  d = NA_PTR_TYPE(rblapack_d, complex*);
  lwork = NUM2INT(rblapack_lwork);
  if (!NA_IsNArray(rblapack_e))
    rb_raise(rb_eArgError, "e (7th argument) must be NArray");
  if (NA_RANK(rblapack_e) != 2)
    rb_raise(rb_eArgError, "rank of e (7th argument) must be %d", 2);
  if (NA_SHAPE1(rblapack_e) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of e must be the same as shape 1 of b");
  lde = NA_SHAPE0(rblapack_e);
  if (NA_TYPE(rblapack_e) != NA_SCOMPLEX)
    rblapack_e = na_change_type(rblapack_e, NA_SCOMPLEX);
  e = NA_PTR_TYPE(rblapack_e, complex*);
  if (!NA_IsNArray(rblapack_f))
    rb_raise(rb_eArgError, "f (8th argument) must be NArray");
  if (NA_RANK(rblapack_f) != 2)
    rb_raise(rb_eArgError, "rank of f (8th argument) must be %d", 2);
  if (NA_SHAPE1(rblapack_f) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of f must be the same as shape 1 of b");
  ldf = NA_SHAPE0(rblapack_f);
  if (NA_TYPE(rblapack_f) != NA_SCOMPLEX)
    rblapack_f = na_change_type(rblapack_f, NA_SCOMPLEX);
  f = NA_PTR_TYPE(rblapack_f, complex*);
  trans = StringValueCStr(rblapack_trans)[0];
  {
    int shape[1];
    shape[0] = MAX(1,lwork);
    rblapack_work = na_make_object(NA_SCOMPLEX, 1, shape, cNArray);
  }
  work = NA_PTR_TYPE(rblapack_work, complex*);
  {
    int shape[2];
    shape[0] = ldc;
    shape[1] = n;
    rblapack_c_out__ = na_make_object(NA_SCOMPLEX, 2, shape, cNArray);
  }
  c_out__ = NA_PTR_TYPE(rblapack_c_out__, complex*);
  MEMCPY(c_out__, c, complex, NA_TOTAL(rblapack_c));
  rblapack_c = rblapack_c_out__;
  c = c_out__;
  {
    int shape[2];
    shape[0] = ldf;
    shape[1] = n;
    rblapack_f_out__ = na_make_object(NA_SCOMPLEX, 2, shape, cNArray);
  }
  f_out__ = NA_PTR_TYPE(rblapack_f_out__, complex*);
  MEMCPY(f_out__, f, complex, NA_TOTAL(rblapack_f));
  rblapack_f = rblapack_f_out__;
  f = f_out__;
  iwork = ALLOC_N(integer, (m+n+2));

  ctgsyl_(&trans, &ijob, &m, &n, a, &lda, b, &ldb, c, &ldc, d, &ldd, e, &lde, f, &ldf, &scale, &dif, work, &lwork, iwork, &info);

  free(iwork);
  rblapack_scale = rb_float_new((double)scale);
  rblapack_dif = rb_float_new((double)dif);
  rblapack_info = INT2NUM(info);
  return rb_ary_new3(6, rblapack_scale, rblapack_dif, rblapack_work, rblapack_info, rblapack_c, rblapack_f);
}

void
init_lapack_ctgsyl(VALUE mLapack){
  rb_define_module_function(mLapack, "ctgsyl", rblapack_ctgsyl, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}

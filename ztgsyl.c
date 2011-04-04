#include "rb_lapack.h"

extern VOID ztgsyl_(char *trans, integer *ijob, integer *m, integer *n, doublecomplex *a, integer *lda, doublecomplex *b, integer *ldb, doublecomplex *c, integer *ldc, doublecomplex *d, integer *ldd, doublecomplex *e, integer *lde, doublecomplex *f, integer *ldf, doublereal *scale, doublereal *dif, doublecomplex *work, integer *lwork, integer *iwork, integer *info);

static VALUE
rb_ztgsyl(int argc, VALUE *argv, VALUE self){
  VALUE rb_trans;
  char trans; 
  VALUE rb_ijob;
  integer ijob; 
  VALUE rb_a;
  doublecomplex *a; 
  VALUE rb_b;
  doublecomplex *b; 
  VALUE rb_c;
  doublecomplex *c; 
  VALUE rb_d;
  doublecomplex *d; 
  VALUE rb_e;
  doublecomplex *e; 
  VALUE rb_f;
  doublecomplex *f; 
  VALUE rb_lwork;
  integer lwork; 
  VALUE rb_scale;
  doublereal scale; 
  VALUE rb_dif;
  doublereal dif; 
  VALUE rb_work;
  doublecomplex *work; 
  VALUE rb_info;
  integer info; 
  VALUE rb_c_out__;
  doublecomplex *c_out__;
  VALUE rb_f_out__;
  doublecomplex *f_out__;
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
    printf("%s\n", "USAGE:\n  scale, dif, work, info, c, f = NumRu::Lapack.ztgsyl( trans, ijob, a, b, c, d, e, f, lwork)\n    or\n  NumRu::Lapack.ztgsyl  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE ZTGSYL( TRANS, IJOB, M, N, A, LDA, B, LDB, C, LDC, D, LDD, E, LDE, F, LDF, SCALE, DIF, WORK, LWORK, IWORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  ZTGSYL solves the generalized Sylvester equation:\n*\n*              A * R - L * B = scale * C            (1)\n*              D * R - L * E = scale * F\n*\n*  where R and L are unknown m-by-n matrices, (A, D), (B, E) and\n*  (C, F) are given matrix pairs of size m-by-m, n-by-n and m-by-n,\n*  respectively, with complex entries. A, B, D and E are upper\n*  triangular (i.e., (A,D) and (B,E) in generalized Schur form).\n*\n*  The solution (R, L) overwrites (C, F). 0 <= SCALE <= 1\n*  is an output scaling factor chosen to avoid overflow.\n*\n*  In matrix notation (1) is equivalent to solve Zx = scale*b, where Z\n*  is defined as\n*\n*         Z = [ kron(In, A)  -kron(B', Im) ]        (2)\n*             [ kron(In, D)  -kron(E', Im) ],\n*\n*  Here Ix is the identity matrix of size x and X' is the conjugate\n*  transpose of X. Kron(X, Y) is the Kronecker product between the\n*  matrices X and Y.\n*\n*  If TRANS = 'C', y in the conjugate transposed system Z'*y = scale*b\n*  is solved for, which is equivalent to solve for R and L in\n*\n*              A' * R + D' * L = scale * C           (3)\n*              R * B' + L * E' = scale * -F\n*\n*  This case (TRANS = 'C') is used to compute an one-norm-based estimate\n*  of Dif[(A,D), (B,E)], the separation between the matrix pairs (A,D)\n*  and (B,E), using ZLACON.\n*\n*  If IJOB >= 1, ZTGSYL computes a Frobenius norm-based estimate of\n*  Dif[(A,D),(B,E)]. That is, the reciprocal of a lower bound on the\n*  reciprocal of the smallest singular value of Z.\n*\n*  This is a level-3 BLAS algorithm.\n*\n\n*  Arguments\n*  =========\n*\n*  TRANS   (input) CHARACTER*1\n*          = 'N': solve the generalized sylvester equation (1).\n*          = 'C': solve the \"conjugate transposed\" system (3).\n*\n*  IJOB    (input) INTEGER\n*          Specifies what kind of functionality to be performed.\n*          =0: solve (1) only.\n*          =1: The functionality of 0 and 3.\n*          =2: The functionality of 0 and 4.\n*          =3: Only an estimate of Dif[(A,D), (B,E)] is computed.\n*              (look ahead strategy is used).\n*          =4: Only an estimate of Dif[(A,D), (B,E)] is computed.\n*              (ZGECON on sub-systems is used).\n*          Not referenced if TRANS = 'C'.\n*\n*  M       (input) INTEGER\n*          The order of the matrices A and D, and the row dimension of\n*          the matrices C, F, R and L.\n*\n*  N       (input) INTEGER\n*          The order of the matrices B and E, and the column dimension\n*          of the matrices C, F, R and L.\n*\n*  A       (input) COMPLEX*16 array, dimension (LDA, M)\n*          The upper triangular matrix A.\n*\n*  LDA     (input) INTEGER\n*          The leading dimension of the array A. LDA >= max(1, M).\n*\n*  B       (input) COMPLEX*16 array, dimension (LDB, N)\n*          The upper triangular matrix B.\n*\n*  LDB     (input) INTEGER\n*          The leading dimension of the array B. LDB >= max(1, N).\n*\n*  C       (input/output) COMPLEX*16 array, dimension (LDC, N)\n*          On entry, C contains the right-hand-side of the first matrix\n*          equation in (1) or (3).\n*          On exit, if IJOB = 0, 1 or 2, C has been overwritten by\n*          the solution R. If IJOB = 3 or 4 and TRANS = 'N', C holds R,\n*          the solution achieved during the computation of the\n*          Dif-estimate.\n*\n*  LDC     (input) INTEGER\n*          The leading dimension of the array C. LDC >= max(1, M).\n*\n*  D       (input) COMPLEX*16 array, dimension (LDD, M)\n*          The upper triangular matrix D.\n*\n*  LDD     (input) INTEGER\n*          The leading dimension of the array D. LDD >= max(1, M).\n*\n*  E       (input) COMPLEX*16 array, dimension (LDE, N)\n*          The upper triangular matrix E.\n*\n*  LDE     (input) INTEGER\n*          The leading dimension of the array E. LDE >= max(1, N).\n*\n*  F       (input/output) COMPLEX*16 array, dimension (LDF, N)\n*          On entry, F contains the right-hand-side of the second matrix\n*          equation in (1) or (3).\n*          On exit, if IJOB = 0, 1 or 2, F has been overwritten by\n*          the solution L. If IJOB = 3 or 4 and TRANS = 'N', F holds L,\n*          the solution achieved during the computation of the\n*          Dif-estimate.\n*\n*  LDF     (input) INTEGER\n*          The leading dimension of the array F. LDF >= max(1, M).\n*\n*  DIF     (output) DOUBLE PRECISION\n*          On exit DIF is the reciprocal of a lower bound of the\n*          reciprocal of the Dif-function, i.e. DIF is an upper bound of\n*          Dif[(A,D), (B,E)] = sigma-min(Z), where Z as in (2).\n*          IF IJOB = 0 or TRANS = 'C', DIF is not referenced.\n*\n*  SCALE   (output) DOUBLE PRECISION\n*          On exit SCALE is the scaling factor in (1) or (3).\n*          If 0 < SCALE < 1, C and F hold the solutions R and L, resp.,\n*          to a slightly perturbed system but the input matrices A, B,\n*          D and E have not been changed. If SCALE = 0, R and L will\n*          hold the solutions to the homogenious system with C = F = 0.\n*\n*  WORK    (workspace/output) COMPLEX*16 array, dimension (MAX(1,LWORK))\n*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.\n*\n*  LWORK   (input) INTEGER\n*          The dimension of the array WORK. LWORK > = 1.\n*          If IJOB = 1 or 2 and TRANS = 'N', LWORK >= max(1,2*M*N).\n*\n*          If LWORK = -1, then a workspace query is assumed; the routine\n*          only calculates the optimal size of the WORK array, returns\n*          this value as the first entry of the WORK array, and no error\n*          message related to LWORK is issued by XERBLA.\n*\n*  IWORK   (workspace) INTEGER array, dimension (M+N+2)\n*\n*  INFO    (output) INTEGER\n*            =0: successful exit\n*            <0: If INFO = -i, the i-th argument had an illegal value.\n*            >0: (A, D) and (B, E) have common or very close\n*                eigenvalues.\n*\n\n*  Further Details\n*  ===============\n*\n*  Based on contributions by\n*     Bo Kagstrom and Peter Poromaa, Department of Computing Science,\n*     Umea University, S-901 87 Umea, Sweden.\n*\n*  [1] B. Kagstrom and P. Poromaa, LAPACK-Style Algorithms and Software\n*      for Solving the Generalized Sylvester Equation and Estimating the\n*      Separation between Regular Matrix Pairs, Report UMINF - 93.23,\n*      Department of Computing Science, Umea University, S-901 87 Umea,\n*      Sweden, December 1993, Revised April 1994, Also as LAPACK Working\n*      Note 75.  To appear in ACM Trans. on Math. Software, Vol 22,\n*      No 1, 1996.\n*\n*  [2] B. Kagstrom, A Perturbation Analysis of the Generalized Sylvester\n*      Equation (AR - LB, DR - LE ) = (C, F), SIAM J. Matrix Anal.\n*      Appl., 15(4):1045-1060, 1994.\n*\n*  [3] B. Kagstrom and L. Westin, Generalized Schur Methods with\n*      Condition Estimators for Solving the Generalized Sylvester\n*      Equation, IEEE Transactions on Automatic Control, Vol. 34, No. 7,\n*      July 1989, pp 745-751.\n*\n*  =====================================================================\n*  Replaced various illegal calls to CCOPY by calls to CLASET.\n*  Sven Hammarling, 1/5/02.\n*\n\n");
    return Qnil;
  }
  if (argc != 9)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 9)", argc);
  rb_trans = argv[0];
  rb_ijob = argv[1];
  rb_a = argv[2];
  rb_b = argv[3];
  rb_c = argv[4];
  rb_d = argv[5];
  rb_e = argv[6];
  rb_f = argv[7];
  rb_lwork = argv[8];

  ijob = NUM2INT(rb_ijob);
  if (!NA_IsNArray(rb_a))
    rb_raise(rb_eArgError, "a (3th argument) must be NArray");
  if (NA_RANK(rb_a) != 2)
    rb_raise(rb_eArgError, "rank of a (3th argument) must be %d", 2);
  m = NA_SHAPE1(rb_a);
  lda = NA_SHAPE0(rb_a);
  if (NA_TYPE(rb_a) != NA_DCOMPLEX)
    rb_a = na_change_type(rb_a, NA_DCOMPLEX);
  a = NA_PTR_TYPE(rb_a, doublecomplex*);
  if (!NA_IsNArray(rb_b))
    rb_raise(rb_eArgError, "b (4th argument) must be NArray");
  if (NA_RANK(rb_b) != 2)
    rb_raise(rb_eArgError, "rank of b (4th argument) must be %d", 2);
  n = NA_SHAPE1(rb_b);
  ldb = NA_SHAPE0(rb_b);
  if (NA_TYPE(rb_b) != NA_DCOMPLEX)
    rb_b = na_change_type(rb_b, NA_DCOMPLEX);
  b = NA_PTR_TYPE(rb_b, doublecomplex*);
  if (!NA_IsNArray(rb_c))
    rb_raise(rb_eArgError, "c (5th argument) must be NArray");
  if (NA_RANK(rb_c) != 2)
    rb_raise(rb_eArgError, "rank of c (5th argument) must be %d", 2);
  if (NA_SHAPE1(rb_c) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of c must be the same as shape 1 of b");
  ldc = NA_SHAPE0(rb_c);
  if (NA_TYPE(rb_c) != NA_DCOMPLEX)
    rb_c = na_change_type(rb_c, NA_DCOMPLEX);
  c = NA_PTR_TYPE(rb_c, doublecomplex*);
  if (!NA_IsNArray(rb_d))
    rb_raise(rb_eArgError, "d (6th argument) must be NArray");
  if (NA_RANK(rb_d) != 2)
    rb_raise(rb_eArgError, "rank of d (6th argument) must be %d", 2);
  if (NA_SHAPE1(rb_d) != m)
    rb_raise(rb_eRuntimeError, "shape 1 of d must be the same as shape 1 of a");
  ldd = NA_SHAPE0(rb_d);
  if (NA_TYPE(rb_d) != NA_DCOMPLEX)
    rb_d = na_change_type(rb_d, NA_DCOMPLEX);
  d = NA_PTR_TYPE(rb_d, doublecomplex*);
  lwork = NUM2INT(rb_lwork);
  if (!NA_IsNArray(rb_e))
    rb_raise(rb_eArgError, "e (7th argument) must be NArray");
  if (NA_RANK(rb_e) != 2)
    rb_raise(rb_eArgError, "rank of e (7th argument) must be %d", 2);
  if (NA_SHAPE1(rb_e) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of e must be the same as shape 1 of b");
  lde = NA_SHAPE0(rb_e);
  if (NA_TYPE(rb_e) != NA_DCOMPLEX)
    rb_e = na_change_type(rb_e, NA_DCOMPLEX);
  e = NA_PTR_TYPE(rb_e, doublecomplex*);
  if (!NA_IsNArray(rb_f))
    rb_raise(rb_eArgError, "f (8th argument) must be NArray");
  if (NA_RANK(rb_f) != 2)
    rb_raise(rb_eArgError, "rank of f (8th argument) must be %d", 2);
  if (NA_SHAPE1(rb_f) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of f must be the same as shape 1 of b");
  ldf = NA_SHAPE0(rb_f);
  if (NA_TYPE(rb_f) != NA_DCOMPLEX)
    rb_f = na_change_type(rb_f, NA_DCOMPLEX);
  f = NA_PTR_TYPE(rb_f, doublecomplex*);
  trans = StringValueCStr(rb_trans)[0];
  {
    int shape[1];
    shape[0] = MAX(1,lwork);
    rb_work = na_make_object(NA_DCOMPLEX, 1, shape, cNArray);
  }
  work = NA_PTR_TYPE(rb_work, doublecomplex*);
  {
    int shape[2];
    shape[0] = ldc;
    shape[1] = n;
    rb_c_out__ = na_make_object(NA_DCOMPLEX, 2, shape, cNArray);
  }
  c_out__ = NA_PTR_TYPE(rb_c_out__, doublecomplex*);
  MEMCPY(c_out__, c, doublecomplex, NA_TOTAL(rb_c));
  rb_c = rb_c_out__;
  c = c_out__;
  {
    int shape[2];
    shape[0] = ldf;
    shape[1] = n;
    rb_f_out__ = na_make_object(NA_DCOMPLEX, 2, shape, cNArray);
  }
  f_out__ = NA_PTR_TYPE(rb_f_out__, doublecomplex*);
  MEMCPY(f_out__, f, doublecomplex, NA_TOTAL(rb_f));
  rb_f = rb_f_out__;
  f = f_out__;
  iwork = ALLOC_N(integer, (m+n+2));

  ztgsyl_(&trans, &ijob, &m, &n, a, &lda, b, &ldb, c, &ldc, d, &ldd, e, &lde, f, &ldf, &scale, &dif, work, &lwork, iwork, &info);

  free(iwork);
  rb_scale = rb_float_new((double)scale);
  rb_dif = rb_float_new((double)dif);
  rb_info = INT2NUM(info);
  return rb_ary_new3(6, rb_scale, rb_dif, rb_work, rb_info, rb_c, rb_f);
}

void
init_lapack_ztgsyl(VALUE mLapack){
  rb_define_module_function(mLapack, "ztgsyl", rb_ztgsyl, -1);
}

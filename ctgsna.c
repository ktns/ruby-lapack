#include "rb_lapack.h"

static VALUE
rb_ctgsna(int argc, VALUE *argv, VALUE self){
  VALUE rb_job;
  char job; 
  VALUE rb_howmny;
  char howmny; 
  VALUE rb_select;
  logical *select; 
  VALUE rb_a;
  complex *a; 
  VALUE rb_b;
  complex *b; 
  VALUE rb_vl;
  complex *vl; 
  VALUE rb_vr;
  complex *vr; 
  VALUE rb_lwork;
  integer lwork; 
  VALUE rb_s;
  real *s; 
  VALUE rb_dif;
  real *dif; 
  VALUE rb_m;
  integer m; 
  VALUE rb_work;
  complex *work; 
  VALUE rb_info;
  integer info; 
  integer *iwork;

  integer n;
  integer lda;
  integer ldb;
  integer ldvl;
  integer ldvr;
  integer mm;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  s, dif, m, work, info = NumRu::Lapack.ctgsna( job, howmny, select, a, b, vl, vr, lwork)\n    or\n  NumRu::Lapack.ctgsna  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE CTGSNA( JOB, HOWMNY, SELECT, N, A, LDA, B, LDB, VL, LDVL, VR, LDVR, S, DIF, MM, M, WORK, LWORK, IWORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  CTGSNA estimates reciprocal condition numbers for specified\n*  eigenvalues and/or eigenvectors of a matrix pair (A, B).\n*\n*  (A, B) must be in generalized Schur canonical form, that is, A and\n*  B are both upper triangular.\n*\n\n*  Arguments\n*  =========\n*\n*  JOB     (input) CHARACTER*1\n*          Specifies whether condition numbers are required for\n*          eigenvalues (S) or eigenvectors (DIF):\n*          = 'E': for eigenvalues only (S);\n*          = 'V': for eigenvectors only (DIF);\n*          = 'B': for both eigenvalues and eigenvectors (S and DIF).\n*\n*  HOWMNY  (input) CHARACTER*1\n*          = 'A': compute condition numbers for all eigenpairs;\n*          = 'S': compute condition numbers for selected eigenpairs\n*                 specified by the array SELECT.\n*\n*  SELECT  (input) LOGICAL array, dimension (N)\n*          If HOWMNY = 'S', SELECT specifies the eigenpairs for which\n*          condition numbers are required. To select condition numbers\n*          for the corresponding j-th eigenvalue and/or eigenvector,\n*          SELECT(j) must be set to .TRUE..\n*          If HOWMNY = 'A', SELECT is not referenced.\n*\n*  N       (input) INTEGER\n*          The order of the square matrix pair (A, B). N >= 0.\n*\n*  A       (input) COMPLEX array, dimension (LDA,N)\n*          The upper triangular matrix A in the pair (A,B).\n*\n*  LDA     (input) INTEGER\n*          The leading dimension of the array A. LDA >= max(1,N).\n*\n*  B       (input) COMPLEX array, dimension (LDB,N)\n*          The upper triangular matrix B in the pair (A, B).\n*\n*  LDB     (input) INTEGER\n*          The leading dimension of the array B. LDB >= max(1,N).\n*\n*  VL      (input) COMPLEX array, dimension (LDVL,M)\n*          IF JOB = 'E' or 'B', VL must contain left eigenvectors of\n*          (A, B), corresponding to the eigenpairs specified by HOWMNY\n*          and SELECT.  The eigenvectors must be stored in consecutive\n*          columns of VL, as returned by CTGEVC.\n*          If JOB = 'V', VL is not referenced.\n*\n*  LDVL    (input) INTEGER\n*          The leading dimension of the array VL. LDVL >= 1; and\n*          If JOB = 'E' or 'B', LDVL >= N.\n*\n*  VR      (input) COMPLEX array, dimension (LDVR,M)\n*          IF JOB = 'E' or 'B', VR must contain right eigenvectors of\n*          (A, B), corresponding to the eigenpairs specified by HOWMNY\n*          and SELECT.  The eigenvectors must be stored in consecutive\n*          columns of VR, as returned by CTGEVC.\n*          If JOB = 'V', VR is not referenced.\n*\n*  LDVR    (input) INTEGER\n*          The leading dimension of the array VR. LDVR >= 1;\n*          If JOB = 'E' or 'B', LDVR >= N.\n*\n*  S       (output) REAL array, dimension (MM)\n*          If JOB = 'E' or 'B', the reciprocal condition numbers of the\n*          selected eigenvalues, stored in consecutive elements of the\n*          array.\n*          If JOB = 'V', S is not referenced.\n*\n*  DIF     (output) REAL array, dimension (MM)\n*          If JOB = 'V' or 'B', the estimated reciprocal condition\n*          numbers of the selected eigenvectors, stored in consecutive\n*          elements of the array.\n*          If the eigenvalues cannot be reordered to compute DIF(j),\n*          DIF(j) is set to 0; this can only occur when the true value\n*          would be very small anyway.\n*          For each eigenvalue/vector specified by SELECT, DIF stores\n*          a Frobenius norm-based estimate of Difl.\n*          If JOB = 'E', DIF is not referenced.\n*\n*  MM      (input) INTEGER\n*          The number of elements in the arrays S and DIF. MM >= M.\n*\n*  M       (output) INTEGER\n*          The number of elements of the arrays S and DIF used to store\n*          the specified condition numbers; for each selected eigenvalue\n*          one element is used. If HOWMNY = 'A', M is set to N.\n*\n*  WORK    (workspace/output) COMPLEX array, dimension (MAX(1,LWORK))\n*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.\n*\n*  LWORK  (input) INTEGER\n*          The dimension of the array WORK. LWORK >= max(1,N).\n*          If JOB = 'V' or 'B', LWORK >= max(1,2*N*N).\n*\n*  IWORK   (workspace) INTEGER array, dimension (N+2)\n*          If JOB = 'E', IWORK is not referenced.\n*\n*  INFO    (output) INTEGER\n*          = 0: Successful exit\n*          < 0: If INFO = -i, the i-th argument had an illegal value\n*\n\n*  Further Details\n*  ===============\n*\n*  The reciprocal of the condition number of the i-th generalized\n*  eigenvalue w = (a, b) is defined as\n*\n*          S(I) = (|v'Au|**2 + |v'Bu|**2)**(1/2) / (norm(u)*norm(v))\n*\n*  where u and v are the right and left eigenvectors of (A, B)\n*  corresponding to w; |z| denotes the absolute value of the complex\n*  number, and norm(u) denotes the 2-norm of the vector u. The pair\n*  (a, b) corresponds to an eigenvalue w = a/b (= v'Au/v'Bu) of the\n*  matrix pair (A, B). If both a and b equal zero, then (A,B) is\n*  singular and S(I) = -1 is returned.\n*\n*  An approximate error bound on the chordal distance between the i-th\n*  computed generalized eigenvalue w and the corresponding exact\n*  eigenvalue lambda is\n*\n*          chord(w, lambda) <=   EPS * norm(A, B) / S(I),\n*\n*  where EPS is the machine precision.\n*\n*  The reciprocal of the condition number of the right eigenvector u\n*  and left eigenvector v corresponding to the generalized eigenvalue w\n*  is defined as follows. Suppose\n*\n*                   (A, B) = ( a   *  ) ( b  *  )  1\n*                            ( 0  A22 ),( 0 B22 )  n-1\n*                              1  n-1     1 n-1\n*\n*  Then the reciprocal condition number DIF(I) is\n*\n*          Difl[(a, b), (A22, B22)]  = sigma-min( Zl )\n*\n*  where sigma-min(Zl) denotes the smallest singular value of\n*\n*         Zl = [ kron(a, In-1) -kron(1, A22) ]\n*              [ kron(b, In-1) -kron(1, B22) ].\n*\n*  Here In-1 is the identity matrix of size n-1 and X' is the conjugate\n*  transpose of X. kron(X, Y) is the Kronecker product between the\n*  matrices X and Y.\n*\n*  We approximate the smallest singular value of Zl with an upper\n*  bound. This is done by CLATDF.\n*\n*  An approximate error bound for a computed eigenvector VL(i) or\n*  VR(i) is given by\n*\n*                      EPS * norm(A, B) / DIF(i).\n*\n*  See ref. [2-3] for more details and further references.\n*\n*  Based on contributions by\n*     Bo Kagstrom and Peter Poromaa, Department of Computing Science,\n*     Umea University, S-901 87 Umea, Sweden.\n*\n*  References\n*  ==========\n*\n*  [1] B. Kagstrom; A Direct Method for Reordering Eigenvalues in the\n*      Generalized Real Schur Form of a Regular Matrix Pair (A, B), in\n*      M.S. Moonen et al (eds), Linear Algebra for Large Scale and\n*      Real-Time Applications, Kluwer Academic Publ. 1993, pp 195-218.\n*\n*  [2] B. Kagstrom and P. Poromaa; Computing Eigenspaces with Specified\n*      Eigenvalues of a Regular Matrix Pair (A, B) and Condition\n*      Estimation: Theory, Algorithms and Software, Report\n*      UMINF - 94.04, Department of Computing Science, Umea University,\n*      S-901 87 Umea, Sweden, 1994. Also as LAPACK Working Note 87.\n*      To appear in Numerical Algorithms, 1996.\n*\n*  [3] B. Kagstrom and P. Poromaa, LAPACK-Style Algorithms and Software\n*      for Solving the Generalized Sylvester Equation and Estimating the\n*      Separation between Regular Matrix Pairs, Report UMINF - 93.23,\n*      Department of Computing Science, Umea University, S-901 87 Umea,\n*      Sweden, December 1993, Revised April 1994, Also as LAPACK Working\n*      Note 75.\n*      To appear in ACM Trans. on Math. Software, Vol 22, No 1, 1996.\n*\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 8)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 8)", argc);
  rb_job = argv[0];
  rb_howmny = argv[1];
  rb_select = argv[2];
  rb_a = argv[3];
  rb_b = argv[4];
  rb_vl = argv[5];
  rb_vr = argv[6];
  rb_lwork = argv[7];

  job = StringValueCStr(rb_job)[0];
  howmny = StringValueCStr(rb_howmny)[0];
  lwork = NUM2INT(rb_lwork);
  if (!NA_IsNArray(rb_select))
    rb_raise(rb_eArgError, "select (3th argument) must be NArray");
  if (NA_RANK(rb_select) != 1)
    rb_raise(rb_eArgError, "rank of select (3th argument) must be %d", 1);
  n = NA_SHAPE0(rb_select);
  if (NA_TYPE(rb_select) != NA_LINT)
    rb_select = na_change_type(rb_select, NA_LINT);
  select = NA_PTR_TYPE(rb_select, logical*);
  if (!NA_IsNArray(rb_a))
    rb_raise(rb_eArgError, "a (4th argument) must be NArray");
  if (NA_RANK(rb_a) != 2)
    rb_raise(rb_eArgError, "rank of a (4th argument) must be %d", 2);
  lda = NA_SHAPE0(rb_a);
  if (NA_SHAPE1(rb_a) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of a must be the same as shape 0 of select");
  if (NA_TYPE(rb_a) != NA_SCOMPLEX)
    rb_a = na_change_type(rb_a, NA_SCOMPLEX);
  a = NA_PTR_TYPE(rb_a, complex*);
  if (!NA_IsNArray(rb_b))
    rb_raise(rb_eArgError, "b (5th argument) must be NArray");
  if (NA_RANK(rb_b) != 2)
    rb_raise(rb_eArgError, "rank of b (5th argument) must be %d", 2);
  ldb = NA_SHAPE0(rb_b);
  if (NA_SHAPE1(rb_b) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of b must be the same as shape 0 of select");
  if (NA_TYPE(rb_b) != NA_SCOMPLEX)
    rb_b = na_change_type(rb_b, NA_SCOMPLEX);
  b = NA_PTR_TYPE(rb_b, complex*);
  if (!NA_IsNArray(rb_vl))
    rb_raise(rb_eArgError, "vl (6th argument) must be NArray");
  if (NA_RANK(rb_vl) != 2)
    rb_raise(rb_eArgError, "rank of vl (6th argument) must be %d", 2);
  ldvl = NA_SHAPE0(rb_vl);
  m = NA_SHAPE1(rb_vl);
  if (NA_TYPE(rb_vl) != NA_SCOMPLEX)
    rb_vl = na_change_type(rb_vl, NA_SCOMPLEX);
  vl = NA_PTR_TYPE(rb_vl, complex*);
  if (!NA_IsNArray(rb_vr))
    rb_raise(rb_eArgError, "vr (7th argument) must be NArray");
  if (NA_RANK(rb_vr) != 2)
    rb_raise(rb_eArgError, "rank of vr (7th argument) must be %d", 2);
  ldvr = NA_SHAPE0(rb_vr);
  if (NA_SHAPE1(rb_vr) != m)
    rb_raise(rb_eRuntimeError, "shape 1 of vr must be the same as shape 1 of vl");
  if (NA_TYPE(rb_vr) != NA_SCOMPLEX)
    rb_vr = na_change_type(rb_vr, NA_SCOMPLEX);
  vr = NA_PTR_TYPE(rb_vr, complex*);
  mm = m;
  {
    int shape[1];
    shape[0] = mm;
    rb_s = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  s = NA_PTR_TYPE(rb_s, real*);
  {
    int shape[1];
    shape[0] = mm;
    rb_dif = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  dif = NA_PTR_TYPE(rb_dif, real*);
  {
    int shape[1];
    shape[0] = MAX(1,lwork);
    rb_work = na_make_object(NA_SCOMPLEX, 1, shape, cNArray);
  }
  work = NA_PTR_TYPE(rb_work, complex*);
  iwork = ALLOC_N(integer, (lsame_(&job,"E") ? 0 : n+2));

  ctgsna_(&job, &howmny, select, &n, a, &lda, b, &ldb, vl, &ldvl, vr, &ldvr, s, dif, &mm, &m, work, &lwork, iwork, &info);

  free(iwork);
  rb_m = INT2NUM(m);
  rb_info = INT2NUM(info);
  return rb_ary_new3(5, rb_s, rb_dif, rb_m, rb_work, rb_info);
}

void
init_lapack_ctgsna(VALUE mLapack){
  rb_define_module_function(mLapack, "ctgsna", rb_ctgsna, -1);
}

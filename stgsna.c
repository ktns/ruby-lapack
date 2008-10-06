#include "rb_lapack.h"

static VALUE
rb_stgsna(int argc, VALUE *argv, VALUE self){
  VALUE rb_job;
  char job; 
  VALUE rb_howmny;
  char howmny; 
  VALUE rb_select;
  logical *select; 
  VALUE rb_a;
  real *a; 
  VALUE rb_b;
  real *b; 
  VALUE rb_vl;
  real *vl; 
  VALUE rb_vr;
  real *vr; 
  VALUE rb_lwork;
  integer lwork; 
  VALUE rb_s;
  real *s; 
  VALUE rb_dif;
  real *dif; 
  VALUE rb_m;
  integer m; 
  VALUE rb_work;
  real *work; 
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
    printf("%s\n", "USAGE:\n  s, dif, m, work, info = NumRu::Lapack.stgsna( job, howmny, select, a, b, vl, vr, lwork)\n    or\n  NumRu::Lapack.stgsna  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE STGSNA( JOB, HOWMNY, SELECT, N, A, LDA, B, LDB, VL, LDVL, VR, LDVR, S, DIF, MM, M, WORK, LWORK, IWORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  STGSNA estimates reciprocal condition numbers for specified\n*  eigenvalues and/or eigenvectors of a matrix pair (A, B) in\n*  generalized real Schur canonical form (or of any matrix pair\n*  (Q*A*Z', Q*B*Z') with orthogonal matrices Q and Z, where\n*  Z' denotes the transpose of Z.\n*\n*  (A, B) must be in generalized real Schur form (as returned by SGGES),\n*  i.e. A is block upper triangular with 1-by-1 and 2-by-2 diagonal\n*  blocks. B is upper triangular.\n*\n*\n\n*  Arguments\n*  =========\n*\n*  JOB     (input) CHARACTER*1\n*          Specifies whether condition numbers are required for\n*          eigenvalues (S) or eigenvectors (DIF):\n*          = 'E': for eigenvalues only (S);\n*          = 'V': for eigenvectors only (DIF);\n*          = 'B': for both eigenvalues and eigenvectors (S and DIF).\n*\n*  HOWMNY  (input) CHARACTER*1\n*          = 'A': compute condition numbers for all eigenpairs;\n*          = 'S': compute condition numbers for selected eigenpairs\n*                 specified by the array SELECT.\n*\n*  SELECT  (input) LOGICAL array, dimension (N)\n*          If HOWMNY = 'S', SELECT specifies the eigenpairs for which\n*          condition numbers are required. To select condition numbers\n*          for the eigenpair corresponding to a real eigenvalue w(j),\n*          SELECT(j) must be set to .TRUE.. To select condition numbers\n*          corresponding to a complex conjugate pair of eigenvalues w(j)\n*          and w(j+1), either SELECT(j) or SELECT(j+1) or both, must be\n*          set to .TRUE..\n*          If HOWMNY = 'A', SELECT is not referenced.\n*\n*  N       (input) INTEGER\n*          The order of the square matrix pair (A, B). N >= 0.\n*\n*  A       (input) REAL array, dimension (LDA,N)\n*          The upper quasi-triangular matrix A in the pair (A,B).\n*\n*  LDA     (input) INTEGER\n*          The leading dimension of the array A. LDA >= max(1,N).\n*\n*  B       (input) REAL array, dimension (LDB,N)\n*          The upper triangular matrix B in the pair (A,B).\n*\n*  LDB     (input) INTEGER\n*          The leading dimension of the array B. LDB >= max(1,N).\n*\n*  VL      (input) REAL array, dimension (LDVL,M)\n*          If JOB = 'E' or 'B', VL must contain left eigenvectors of\n*          (A, B), corresponding to the eigenpairs specified by HOWMNY\n*          and SELECT. The eigenvectors must be stored in consecutive\n*          columns of VL, as returned by STGEVC.\n*          If JOB = 'V', VL is not referenced.\n*\n*  LDVL    (input) INTEGER\n*          The leading dimension of the array VL. LDVL >= 1.\n*          If JOB = 'E' or 'B', LDVL >= N.\n*\n*  VR      (input) REAL array, dimension (LDVR,M)\n*          If JOB = 'E' or 'B', VR must contain right eigenvectors of\n*          (A, B), corresponding to the eigenpairs specified by HOWMNY\n*          and SELECT. The eigenvectors must be stored in consecutive\n*          columns ov VR, as returned by STGEVC.\n*          If JOB = 'V', VR is not referenced.\n*\n*  LDVR    (input) INTEGER\n*          The leading dimension of the array VR. LDVR >= 1.\n*          If JOB = 'E' or 'B', LDVR >= N.\n*\n*  S       (output) REAL array, dimension (MM)\n*          If JOB = 'E' or 'B', the reciprocal condition numbers of the\n*          selected eigenvalues, stored in consecutive elements of the\n*          array. For a complex conjugate pair of eigenvalues two\n*          consecutive elements of S are set to the same value. Thus\n*          S(j), DIF(j), and the j-th columns of VL and VR all\n*          correspond to the same eigenpair (but not in general the\n*          j-th eigenpair, unless all eigenpairs are selected).\n*          If JOB = 'V', S is not referenced.\n*\n*  DIF     (output) REAL array, dimension (MM)\n*          If JOB = 'V' or 'B', the estimated reciprocal condition\n*          numbers of the selected eigenvectors, stored in consecutive\n*          elements of the array. For a complex eigenvector two\n*          consecutive elements of DIF are set to the same value. If\n*          the eigenvalues cannot be reordered to compute DIF(j), DIF(j)\n*          is set to 0; this can only occur when the true value would be\n*          very small anyway.\n*          If JOB = 'E', DIF is not referenced.\n*\n*  MM      (input) INTEGER\n*          The number of elements in the arrays S and DIF. MM >= M.\n*\n*  M       (output) INTEGER\n*          The number of elements of the arrays S and DIF used to store\n*          the specified condition numbers; for each selected real\n*          eigenvalue one element is used, and for each selected complex\n*          conjugate pair of eigenvalues, two elements are used.\n*          If HOWMNY = 'A', M is set to N.\n*\n*  WORK    (workspace/output) REAL array, dimension (MAX(1,LWORK))\n*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.\n*\n*  LWORK   (input) INTEGER\n*          The dimension of the array WORK. LWORK >= max(1,N).\n*          If JOB = 'V' or 'B' LWORK >= 2*N*(N+2)+16.\n*\n*          If LWORK = -1, then a workspace query is assumed; the routine\n*          only calculates the optimal size of the WORK array, returns\n*          this value as the first entry of the WORK array, and no error\n*          message related to LWORK is issued by XERBLA.\n*\n*  IWORK   (workspace) INTEGER array, dimension (N + 6)\n*          If JOB = 'E', IWORK is not referenced.\n*\n*  INFO    (output) INTEGER\n*          =0: Successful exit\n*          <0: If INFO = -i, the i-th argument had an illegal value\n*\n*\n\n*  Further Details\n*  ===============\n*\n*  The reciprocal of the condition number of a generalized eigenvalue\n*  w = (a, b) is defined as\n*\n*       S(w) = (|u'Av|**2 + |u'Bv|**2)**(1/2) / (norm(u)*norm(v))\n*\n*  where u and v are the left and right eigenvectors of (A, B)\n*  corresponding to w; |z| denotes the absolute value of the complex\n*  number, and norm(u) denotes the 2-norm of the vector u.\n*  The pair (a, b) corresponds to an eigenvalue w = a/b (= u'Av/u'Bv)\n*  of the matrix pair (A, B). If both a and b equal zero, then (A B) is\n*  singular and S(I) = -1 is returned.\n*\n*  An approximate error bound on the chordal distance between the i-th\n*  computed generalized eigenvalue w and the corresponding exact\n*  eigenvalue lambda is\n*\n*       chord(w, lambda) <= EPS * norm(A, B) / S(I)\n*\n*  where EPS is the machine precision.\n*\n*  The reciprocal of the condition number DIF(i) of right eigenvector u\n*  and left eigenvector v corresponding to the generalized eigenvalue w\n*  is defined as follows:\n*\n*  a) If the i-th eigenvalue w = (a,b) is real\n*\n*     Suppose U and V are orthogonal transformations such that\n*\n*                U'*(A, B)*V  = (S, T) = ( a   *  ) ( b  *  )  1\n*                                        ( 0  S22 ),( 0 T22 )  n-1\n*                                          1  n-1     1 n-1\n*\n*     Then the reciprocal condition number DIF(i) is\n*\n*                Difl((a, b), (S22, T22)) = sigma-min( Zl ),\n*\n*     where sigma-min(Zl) denotes the smallest singular value of the\n*     2(n-1)-by-2(n-1) matrix\n*\n*         Zl = [ kron(a, In-1)  -kron(1, S22) ]\n*              [ kron(b, In-1)  -kron(1, T22) ] .\n*\n*     Here In-1 is the identity matrix of size n-1. kron(X, Y) is the\n*     Kronecker product between the matrices X and Y.\n*\n*     Note that if the default method for computing DIF(i) is wanted\n*     (see SLATDF), then the parameter DIFDRI (see below) should be\n*     changed from 3 to 4 (routine SLATDF(IJOB = 2 will be used)).\n*     See STGSYL for more details.\n*\n*  b) If the i-th and (i+1)-th eigenvalues are complex conjugate pair,\n*\n*     Suppose U and V are orthogonal transformations such that\n*\n*                U'*(A, B)*V = (S, T) = ( S11  *   ) ( T11  *  )  2\n*                                       ( 0    S22 ),( 0    T22) n-2\n*                                         2    n-2     2    n-2\n*\n*     and (S11, T11) corresponds to the complex conjugate eigenvalue\n*     pair (w, conjg(w)). There exist unitary matrices U1 and V1 such\n*     that\n*\n*         U1'*S11*V1 = ( s11 s12 )   and U1'*T11*V1 = ( t11 t12 )\n*                      (  0  s22 )                    (  0  t22 )\n*\n*     where the generalized eigenvalues w = s11/t11 and\n*     conjg(w) = s22/t22.\n*\n*     Then the reciprocal condition number DIF(i) is bounded by\n*\n*         min( d1, max( 1, |real(s11)/real(s22)| )*d2 )\n*\n*     where, d1 = Difl((s11, t11), (s22, t22)) = sigma-min(Z1), where\n*     Z1 is the complex 2-by-2 matrix\n*\n*              Z1 =  [ s11  -s22 ]\n*                    [ t11  -t22 ],\n*\n*     This is done by computing (using real arithmetic) the\n*     roots of the characteristical polynomial det(Z1' * Z1 - lambda I),\n*     where Z1' denotes the conjugate transpose of Z1 and det(X) denotes\n*     the determinant of X.\n*\n*     and d2 is an upper bound on Difl((S11, T11), (S22, T22)), i.e. an\n*     upper bound on sigma-min(Z2), where Z2 is (2n-2)-by-(2n-2)\n*\n*              Z2 = [ kron(S11', In-2)  -kron(I2, S22) ]\n*                   [ kron(T11', In-2)  -kron(I2, T22) ]\n*\n*     Note that if the default method for computing DIF is wanted (see\n*     SLATDF), then the parameter DIFDRI (see below) should be changed\n*     from 3 to 4 (routine SLATDF(IJOB = 2 will be used)). See STGSYL\n*     for more details.\n*\n*  For each eigenvalue/vector specified by SELECT, DIF stores a\n*  Frobenius norm-based estimate of Difl.\n*\n*  An approximate error bound for the i-th computed eigenvector VL(i) or\n*  VR(i) is given by\n*\n*             EPS * norm(A, B) / DIF(i).\n*\n*  See ref. [2-3] for more details and further references.\n*\n*  Based on contributions by\n*     Bo Kagstrom and Peter Poromaa, Department of Computing Science,\n*     Umea University, S-901 87 Umea, Sweden.\n*\n*  References\n*  ==========\n*\n*  [1] B. Kagstrom; A Direct Method for Reordering Eigenvalues in the\n*      Generalized Real Schur Form of a Regular Matrix Pair (A, B), in\n*      M.S. Moonen et al (eds), Linear Algebra for Large Scale and\n*      Real-Time Applications, Kluwer Academic Publ. 1993, pp 195-218.\n*\n*  [2] B. Kagstrom and P. Poromaa; Computing Eigenspaces with Specified\n*      Eigenvalues of a Regular Matrix Pair (A, B) and Condition\n*      Estimation: Theory, Algorithms and Software,\n*      Report UMINF - 94.04, Department of Computing Science, Umea\n*      University, S-901 87 Umea, Sweden, 1994. Also as LAPACK Working\n*      Note 87. To appear in Numerical Algorithms, 1996.\n*\n*  [3] B. Kagstrom and P. Poromaa, LAPACK-Style Algorithms and Software\n*      for Solving the Generalized Sylvester Equation and Estimating the\n*      Separation between Regular Matrix Pairs, Report UMINF - 93.23,\n*      Department of Computing Science, Umea University, S-901 87 Umea,\n*      Sweden, December 1993, Revised April 1994, Also as LAPACK Working\n*      Note 75.  To appear in ACM Trans. on Math. Software, Vol 22,\n*      No 1, 1996.\n*\n*  =====================================================================\n*\n\n");
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
  if (NA_TYPE(rb_a) != NA_SFLOAT)
    rb_a = na_change_type(rb_a, NA_SFLOAT);
  a = NA_PTR_TYPE(rb_a, real*);
  if (!NA_IsNArray(rb_b))
    rb_raise(rb_eArgError, "b (5th argument) must be NArray");
  if (NA_RANK(rb_b) != 2)
    rb_raise(rb_eArgError, "rank of b (5th argument) must be %d", 2);
  ldb = NA_SHAPE0(rb_b);
  if (NA_SHAPE1(rb_b) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of b must be the same as shape 0 of select");
  if (NA_TYPE(rb_b) != NA_SFLOAT)
    rb_b = na_change_type(rb_b, NA_SFLOAT);
  b = NA_PTR_TYPE(rb_b, real*);
  if (!NA_IsNArray(rb_vl))
    rb_raise(rb_eArgError, "vl (6th argument) must be NArray");
  if (NA_RANK(rb_vl) != 2)
    rb_raise(rb_eArgError, "rank of vl (6th argument) must be %d", 2);
  ldvl = NA_SHAPE0(rb_vl);
  m = NA_SHAPE1(rb_vl);
  if (NA_TYPE(rb_vl) != NA_SFLOAT)
    rb_vl = na_change_type(rb_vl, NA_SFLOAT);
  vl = NA_PTR_TYPE(rb_vl, real*);
  if (!NA_IsNArray(rb_vr))
    rb_raise(rb_eArgError, "vr (7th argument) must be NArray");
  if (NA_RANK(rb_vr) != 2)
    rb_raise(rb_eArgError, "rank of vr (7th argument) must be %d", 2);
  ldvr = NA_SHAPE0(rb_vr);
  if (NA_SHAPE1(rb_vr) != m)
    rb_raise(rb_eRuntimeError, "shape 1 of vr must be the same as shape 1 of vl");
  if (NA_TYPE(rb_vr) != NA_SFLOAT)
    rb_vr = na_change_type(rb_vr, NA_SFLOAT);
  vr = NA_PTR_TYPE(rb_vr, real*);
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
    rb_work = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  work = NA_PTR_TYPE(rb_work, real*);
  iwork = ALLOC_N(integer, (lsame_(&job,"E") ? 0 : n + 6));

  stgsna_(&job, &howmny, select, &n, a, &lda, b, &ldb, vl, &ldvl, vr, &ldvr, s, dif, &mm, &m, work, &lwork, iwork, &info);

  free(iwork);
  rb_m = INT2NUM(m);
  rb_info = INT2NUM(info);
  return rb_ary_new3(5, rb_s, rb_dif, rb_m, rb_work, rb_info);
}

void
init_lapack_stgsna(VALUE mLapack){
  rb_define_module_function(mLapack, "stgsna", rb_stgsna, -1);
}

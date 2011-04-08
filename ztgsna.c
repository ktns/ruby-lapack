#include "rb_lapack.h"

extern VOID ztgsna_(char *job, char *howmny, logical *select, integer *n, doublecomplex *a, integer *lda, doublecomplex *b, integer *ldb, doublecomplex *vl, integer *ldvl, doublecomplex *vr, integer *ldvr, doublereal *s, doublereal *dif, integer *mm, integer *m, doublecomplex *work, integer *lwork, integer *iwork, integer *info);

static VALUE sHelp, sUsage;

static VALUE
rblapack_ztgsna(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_job;
  char job; 
  VALUE rblapack_howmny;
  char howmny; 
  VALUE rblapack_select;
  logical *select; 
  VALUE rblapack_a;
  doublecomplex *a; 
  VALUE rblapack_b;
  doublecomplex *b; 
  VALUE rblapack_vl;
  doublecomplex *vl; 
  VALUE rblapack_vr;
  doublecomplex *vr; 
  VALUE rblapack_lwork;
  integer lwork; 
  VALUE rblapack_s;
  doublereal *s; 
  VALUE rblapack_dif;
  doublereal *dif; 
  VALUE rblapack_m;
  integer m; 
  VALUE rblapack_work;
  doublecomplex *work; 
  VALUE rblapack_info;
  integer info; 
  integer *iwork;

  integer n;
  integer lda;
  integer ldb;
  integer ldvl;
  integer ldvr;
  integer mm;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  s, dif, m, work, info = NumRu::Lapack.ztgsna( job, howmny, select, a, b, vl, vr, lwork, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE ZTGSNA( JOB, HOWMNY, SELECT, N, A, LDA, B, LDB, VL, LDVL, VR, LDVR, S, DIF, MM, M, WORK, LWORK, IWORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  ZTGSNA estimates reciprocal condition numbers for specified\n*  eigenvalues and/or eigenvectors of a matrix pair (A, B).\n*\n*  (A, B) must be in generalized Schur canonical form, that is, A and\n*  B are both upper triangular.\n*\n\n*  Arguments\n*  =========\n*\n*  JOB     (input) CHARACTER*1\n*          Specifies whether condition numbers are required for\n*          eigenvalues (S) or eigenvectors (DIF):\n*          = 'E': for eigenvalues only (S);\n*          = 'V': for eigenvectors only (DIF);\n*          = 'B': for both eigenvalues and eigenvectors (S and DIF).\n*\n*  HOWMNY  (input) CHARACTER*1\n*          = 'A': compute condition numbers for all eigenpairs;\n*          = 'S': compute condition numbers for selected eigenpairs\n*                 specified by the array SELECT.\n*\n*  SELECT  (input) LOGICAL array, dimension (N)\n*          If HOWMNY = 'S', SELECT specifies the eigenpairs for which\n*          condition numbers are required. To select condition numbers\n*          for the corresponding j-th eigenvalue and/or eigenvector,\n*          SELECT(j) must be set to .TRUE..\n*          If HOWMNY = 'A', SELECT is not referenced.\n*\n*  N       (input) INTEGER\n*          The order of the square matrix pair (A, B). N >= 0.\n*\n*  A       (input) COMPLEX*16 array, dimension (LDA,N)\n*          The upper triangular matrix A in the pair (A,B).\n*\n*  LDA     (input) INTEGER\n*          The leading dimension of the array A. LDA >= max(1,N).\n*\n*  B       (input) COMPLEX*16 array, dimension (LDB,N)\n*          The upper triangular matrix B in the pair (A, B).\n*\n*  LDB     (input) INTEGER\n*          The leading dimension of the array B. LDB >= max(1,N).\n*\n*  VL      (input) COMPLEX*16 array, dimension (LDVL,M)\n*          IF JOB = 'E' or 'B', VL must contain left eigenvectors of\n*          (A, B), corresponding to the eigenpairs specified by HOWMNY\n*          and SELECT.  The eigenvectors must be stored in consecutive\n*          columns of VL, as returned by ZTGEVC.\n*          If JOB = 'V', VL is not referenced.\n*\n*  LDVL    (input) INTEGER\n*          The leading dimension of the array VL. LDVL >= 1; and\n*          If JOB = 'E' or 'B', LDVL >= N.\n*\n*  VR      (input) COMPLEX*16 array, dimension (LDVR,M)\n*          IF JOB = 'E' or 'B', VR must contain right eigenvectors of\n*          (A, B), corresponding to the eigenpairs specified by HOWMNY\n*          and SELECT.  The eigenvectors must be stored in consecutive\n*          columns of VR, as returned by ZTGEVC.\n*          If JOB = 'V', VR is not referenced.\n*\n*  LDVR    (input) INTEGER\n*          The leading dimension of the array VR. LDVR >= 1;\n*          If JOB = 'E' or 'B', LDVR >= N.\n*\n*  S       (output) DOUBLE PRECISION array, dimension (MM)\n*          If JOB = 'E' or 'B', the reciprocal condition numbers of the\n*          selected eigenvalues, stored in consecutive elements of the\n*          array.\n*          If JOB = 'V', S is not referenced.\n*\n*  DIF     (output) DOUBLE PRECISION array, dimension (MM)\n*          If JOB = 'V' or 'B', the estimated reciprocal condition\n*          numbers of the selected eigenvectors, stored in consecutive\n*          elements of the array.\n*          If the eigenvalues cannot be reordered to compute DIF(j),\n*          DIF(j) is set to 0; this can only occur when the true value\n*          would be very small anyway.\n*          For each eigenvalue/vector specified by SELECT, DIF stores\n*          a Frobenius norm-based estimate of Difl.\n*          If JOB = 'E', DIF is not referenced.\n*\n*  MM      (input) INTEGER\n*          The number of elements in the arrays S and DIF. MM >= M.\n*\n*  M       (output) INTEGER\n*          The number of elements of the arrays S and DIF used to store\n*          the specified condition numbers; for each selected eigenvalue\n*          one element is used. If HOWMNY = 'A', M is set to N.\n*\n*  WORK    (workspace/output) COMPLEX*16 array, dimension (MAX(1,LWORK))\n*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.\n*\n*  LWORK  (input) INTEGER\n*          The dimension of the array WORK. LWORK >= max(1,N).\n*          If JOB = 'V' or 'B', LWORK >= max(1,2*N*N).\n*\n*  IWORK   (workspace) INTEGER array, dimension (N+2)\n*          If JOB = 'E', IWORK is not referenced.\n*\n*  INFO    (output) INTEGER\n*          = 0: Successful exit\n*          < 0: If INFO = -i, the i-th argument had an illegal value\n*\n\n*  Further Details\n*  ===============\n*\n*  The reciprocal of the condition number of the i-th generalized\n*  eigenvalue w = (a, b) is defined as\n*\n*          S(I) = (|v'Au|**2 + |v'Bu|**2)**(1/2) / (norm(u)*norm(v))\n*\n*  where u and v are the right and left eigenvectors of (A, B)\n*  corresponding to w; |z| denotes the absolute value of the complex\n*  number, and norm(u) denotes the 2-norm of the vector u. The pair\n*  (a, b) corresponds to an eigenvalue w = a/b (= v'Au/v'Bu) of the\n*  matrix pair (A, B). If both a and b equal zero, then (A,B) is\n*  singular and S(I) = -1 is returned.\n*\n*  An approximate error bound on the chordal distance between the i-th\n*  computed generalized eigenvalue w and the corresponding exact\n*  eigenvalue lambda is\n*\n*          chord(w, lambda) <=   EPS * norm(A, B) / S(I),\n*\n*  where EPS is the machine precision.\n*\n*  The reciprocal of the condition number of the right eigenvector u\n*  and left eigenvector v corresponding to the generalized eigenvalue w\n*  is defined as follows. Suppose\n*\n*                   (A, B) = ( a   *  ) ( b  *  )  1\n*                            ( 0  A22 ),( 0 B22 )  n-1\n*                              1  n-1     1 n-1\n*\n*  Then the reciprocal condition number DIF(I) is\n*\n*          Difl[(a, b), (A22, B22)]  = sigma-min( Zl )\n*\n*  where sigma-min(Zl) denotes the smallest singular value of\n*\n*         Zl = [ kron(a, In-1) -kron(1, A22) ]\n*              [ kron(b, In-1) -kron(1, B22) ].\n*\n*  Here In-1 is the identity matrix of size n-1 and X' is the conjugate\n*  transpose of X. kron(X, Y) is the Kronecker product between the\n*  matrices X and Y.\n*\n*  We approximate the smallest singular value of Zl with an upper\n*  bound. This is done by ZLATDF.\n*\n*  An approximate error bound for a computed eigenvector VL(i) or\n*  VR(i) is given by\n*\n*                      EPS * norm(A, B) / DIF(i).\n*\n*  See ref. [2-3] for more details and further references.\n*\n*  Based on contributions by\n*     Bo Kagstrom and Peter Poromaa, Department of Computing Science,\n*     Umea University, S-901 87 Umea, Sweden.\n*\n*  References\n*  ==========\n*\n*  [1] B. Kagstrom; A Direct Method for Reordering Eigenvalues in the\n*      Generalized Real Schur Form of a Regular Matrix Pair (A, B), in\n*      M.S. Moonen et al (eds), Linear Algebra for Large Scale and\n*      Real-Time Applications, Kluwer Academic Publ. 1993, pp 195-218.\n*\n*  [2] B. Kagstrom and P. Poromaa; Computing Eigenspaces with Specified\n*      Eigenvalues of a Regular Matrix Pair (A, B) and Condition\n*      Estimation: Theory, Algorithms and Software, Report\n*      UMINF - 94.04, Department of Computing Science, Umea University,\n*      S-901 87 Umea, Sweden, 1994. Also as LAPACK Working Note 87.\n*      To appear in Numerical Algorithms, 1996.\n*\n*  [3] B. Kagstrom and P. Poromaa, LAPACK-Style Algorithms and Software\n*      for Solving the Generalized Sylvester Equation and Estimating the\n*      Separation between Regular Matrix Pairs, Report UMINF - 93.23,\n*      Department of Computing Science, Umea University, S-901 87 Umea,\n*      Sweden, December 1993, Revised April 1994, Also as LAPACK Working\n*      Note 75.\n*      To appear in ACM Trans. on Math. Software, Vol 22, No 1, 1996.\n*\n*  =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  s, dif, m, work, info = NumRu::Lapack.ztgsna( job, howmny, select, a, b, vl, vr, lwork, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 8)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 8)", argc);
  rblapack_job = argv[0];
  rblapack_howmny = argv[1];
  rblapack_select = argv[2];
  rblapack_a = argv[3];
  rblapack_b = argv[4];
  rblapack_vl = argv[5];
  rblapack_vr = argv[6];
  rblapack_lwork = argv[7];
  if (rb_options != Qnil) {
  }

  howmny = StringValueCStr(rblapack_howmny)[0];
  if (!NA_IsNArray(rblapack_vl))
    rb_raise(rb_eArgError, "vl (6th argument) must be NArray");
  if (NA_RANK(rblapack_vl) != 2)
    rb_raise(rb_eArgError, "rank of vl (6th argument) must be %d", 2);
  m = NA_SHAPE1(rblapack_vl);
  ldvl = NA_SHAPE0(rblapack_vl);
  if (NA_TYPE(rblapack_vl) != NA_DCOMPLEX)
    rblapack_vl = na_change_type(rblapack_vl, NA_DCOMPLEX);
  vl = NA_PTR_TYPE(rblapack_vl, doublecomplex*);
  if (!NA_IsNArray(rblapack_a))
    rb_raise(rb_eArgError, "a (4th argument) must be NArray");
  if (NA_RANK(rblapack_a) != 2)
    rb_raise(rb_eArgError, "rank of a (4th argument) must be %d", 2);
  n = NA_SHAPE1(rblapack_a);
  lda = NA_SHAPE0(rblapack_a);
  if (NA_TYPE(rblapack_a) != NA_DCOMPLEX)
    rblapack_a = na_change_type(rblapack_a, NA_DCOMPLEX);
  a = NA_PTR_TYPE(rblapack_a, doublecomplex*);
  if (!NA_IsNArray(rblapack_b))
    rb_raise(rb_eArgError, "b (5th argument) must be NArray");
  if (NA_RANK(rblapack_b) != 2)
    rb_raise(rb_eArgError, "rank of b (5th argument) must be %d", 2);
  if (NA_SHAPE1(rblapack_b) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of b must be the same as shape 1 of a");
  ldb = NA_SHAPE0(rblapack_b);
  if (NA_TYPE(rblapack_b) != NA_DCOMPLEX)
    rblapack_b = na_change_type(rblapack_b, NA_DCOMPLEX);
  b = NA_PTR_TYPE(rblapack_b, doublecomplex*);
  if (!NA_IsNArray(rblapack_select))
    rb_raise(rb_eArgError, "select (3th argument) must be NArray");
  if (NA_RANK(rblapack_select) != 1)
    rb_raise(rb_eArgError, "rank of select (3th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_select) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of select must be the same as shape 1 of a");
  if (NA_TYPE(rblapack_select) != NA_LINT)
    rblapack_select = na_change_type(rblapack_select, NA_LINT);
  select = NA_PTR_TYPE(rblapack_select, logical*);
  if (!NA_IsNArray(rblapack_vr))
    rb_raise(rb_eArgError, "vr (7th argument) must be NArray");
  if (NA_RANK(rblapack_vr) != 2)
    rb_raise(rb_eArgError, "rank of vr (7th argument) must be %d", 2);
  if (NA_SHAPE1(rblapack_vr) != m)
    rb_raise(rb_eRuntimeError, "shape 1 of vr must be the same as shape 1 of vl");
  ldvr = NA_SHAPE0(rblapack_vr);
  if (NA_TYPE(rblapack_vr) != NA_DCOMPLEX)
    rblapack_vr = na_change_type(rblapack_vr, NA_DCOMPLEX);
  vr = NA_PTR_TYPE(rblapack_vr, doublecomplex*);
  job = StringValueCStr(rblapack_job)[0];
  lwork = NUM2INT(rblapack_lwork);
  mm = m;
  {
    int shape[1];
    shape[0] = mm;
    rblapack_s = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  s = NA_PTR_TYPE(rblapack_s, doublereal*);
  {
    int shape[1];
    shape[0] = mm;
    rblapack_dif = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  dif = NA_PTR_TYPE(rblapack_dif, doublereal*);
  {
    int shape[1];
    shape[0] = MAX(1,lwork);
    rblapack_work = na_make_object(NA_DCOMPLEX, 1, shape, cNArray);
  }
  work = NA_PTR_TYPE(rblapack_work, doublecomplex*);
  iwork = ALLOC_N(integer, (lsame_(&job,"E") ? 0 : n+2));

  ztgsna_(&job, &howmny, select, &n, a, &lda, b, &ldb, vl, &ldvl, vr, &ldvr, s, dif, &mm, &m, work, &lwork, iwork, &info);

  free(iwork);
  rblapack_m = INT2NUM(m);
  rblapack_info = INT2NUM(info);
  return rb_ary_new3(5, rblapack_s, rblapack_dif, rblapack_m, rblapack_work, rblapack_info);
}

void
init_lapack_ztgsna(VALUE mLapack){
  rb_define_module_function(mLapack, "ztgsna", rblapack_ztgsna, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}

#include "rb_lapack.h"

extern VOID stgex2_(logical *wantq, logical *wantz, integer *n, real *a, integer *lda, real *b, integer *ldb, real *q, integer *ldq, real *z, integer *ldz, integer *j1, integer *n1, integer *n2, real *work, integer *lwork, integer *info);

static VALUE
rb_stgex2(int argc, VALUE *argv, VALUE self){
  VALUE rb_wantq;
  logical wantq; 
  VALUE rb_wantz;
  logical wantz; 
  VALUE rb_a;
  real *a; 
  VALUE rb_b;
  real *b; 
  VALUE rb_q;
  real *q; 
  VALUE rb_ldq;
  integer ldq; 
  VALUE rb_z;
  real *z; 
  VALUE rb_j1;
  integer j1; 
  VALUE rb_n1;
  integer n1; 
  VALUE rb_n2;
  integer n2; 
  VALUE rb_info;
  integer info; 
  VALUE rb_a_out__;
  real *a_out__;
  VALUE rb_b_out__;
  real *b_out__;
  VALUE rb_q_out__;
  real *q_out__;
  VALUE rb_z_out__;
  real *z_out__;
  real *work;

  integer lda;
  integer n;
  integer ldb;
  integer ldz;
  integer lwork;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  info, a, b, q, z = NumRu::Lapack.stgex2( wantq, wantz, a, b, q, ldq, z, j1, n1, n2)\n    or\n  NumRu::Lapack.stgex2  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE STGEX2( WANTQ, WANTZ, N, A, LDA, B, LDB, Q, LDQ, Z, LDZ, J1, N1, N2, WORK, LWORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  STGEX2 swaps adjacent diagonal blocks (A11, B11) and (A22, B22)\n*  of size 1-by-1 or 2-by-2 in an upper (quasi) triangular matrix pair\n*  (A, B) by an orthogonal equivalence transformation.\n*\n*  (A, B) must be in generalized real Schur canonical form (as returned\n*  by SGGES), i.e. A is block upper triangular with 1-by-1 and 2-by-2\n*  diagonal blocks. B is upper triangular.\n*\n*  Optionally, the matrices Q and Z of generalized Schur vectors are\n*  updated.\n*\n*         Q(in) * A(in) * Z(in)' = Q(out) * A(out) * Z(out)'\n*         Q(in) * B(in) * Z(in)' = Q(out) * B(out) * Z(out)'\n*\n*\n\n*  Arguments\n*  =========\n*\n*  WANTQ   (input) LOGICAL\n*          .TRUE. : update the left transformation matrix Q;\n*          .FALSE.: do not update Q.\n*\n*  WANTZ   (input) LOGICAL\n*          .TRUE. : update the right transformation matrix Z;\n*          .FALSE.: do not update Z.\n*\n*  N       (input) INTEGER\n*          The order of the matrices A and B. N >= 0.\n*\n*  A      (input/output) REAL arrays, dimensions (LDA,N)\n*          On entry, the matrix A in the pair (A, B).\n*          On exit, the updated matrix A.\n*\n*  LDA     (input)  INTEGER\n*          The leading dimension of the array A. LDA >= max(1,N).\n*\n*  B      (input/output) REAL arrays, dimensions (LDB,N)\n*          On entry, the matrix B in the pair (A, B).\n*          On exit, the updated matrix B.\n*\n*  LDB     (input)  INTEGER\n*          The leading dimension of the array B. LDB >= max(1,N).\n*\n*  Q       (input/output) REAL array, dimension (LDZ,N)\n*          On entry, if WANTQ = .TRUE., the orthogonal matrix Q.\n*          On exit, the updated matrix Q.\n*          Not referenced if WANTQ = .FALSE..\n*\n*  LDQ     (input) INTEGER\n*          The leading dimension of the array Q. LDQ >= 1.\n*          If WANTQ = .TRUE., LDQ >= N.\n*\n*  Z       (input/output) REAL array, dimension (LDZ,N)\n*          On entry, if WANTZ =.TRUE., the orthogonal matrix Z.\n*          On exit, the updated matrix Z.\n*          Not referenced if WANTZ = .FALSE..\n*\n*  LDZ     (input) INTEGER\n*          The leading dimension of the array Z. LDZ >= 1.\n*          If WANTZ = .TRUE., LDZ >= N.\n*\n*  J1      (input) INTEGER\n*          The index to the first block (A11, B11). 1 <= J1 <= N.\n*\n*  N1      (input) INTEGER\n*          The order of the first block (A11, B11). N1 = 0, 1 or 2.\n*\n*  N2      (input) INTEGER\n*          The order of the second block (A22, B22). N2 = 0, 1 or 2.\n*\n*  WORK    (workspace) REAL array, dimension (MAX(1,LWORK)).\n*\n*  LWORK   (input) INTEGER\n*          The dimension of the array WORK.\n*          LWORK >=  MAX( N*(N2+N1), (N2+N1)*(N2+N1)*2 )\n*\n*  INFO    (output) INTEGER\n*            =0: Successful exit\n*            >0: If INFO = 1, the transformed matrix (A, B) would be\n*                too far from generalized Schur form; the blocks are\n*                not swapped and (A, B) and (Q, Z) are unchanged.\n*                The problem of swapping is too ill-conditioned.\n*            <0: If INFO = -16: LWORK is too small. Appropriate value\n*                for LWORK is returned in WORK(1).\n*\n\n*  Further Details\n*  ===============\n*\n*  Based on contributions by\n*     Bo Kagstrom and Peter Poromaa, Department of Computing Science,\n*     Umea University, S-901 87 Umea, Sweden.\n*\n*  In the current code both weak and strong stability tests are\n*  performed. The user can omit the strong stability test by changing\n*  the internal logical parameter WANDS to .FALSE.. See ref. [2] for\n*  details.\n*\n*  [1] B. Kagstrom; A Direct Method for Reordering Eigenvalues in the\n*      Generalized Real Schur Form of a Regular Matrix Pair (A, B), in\n*      M.S. Moonen et al (eds), Linear Algebra for Large Scale and\n*      Real-Time Applications, Kluwer Academic Publ. 1993, pp 195-218.\n*\n*  [2] B. Kagstrom and P. Poromaa; Computing Eigenspaces with Specified\n*      Eigenvalues of a Regular Matrix Pair (A, B) and Condition\n*      Estimation: Theory, Algorithms and Software,\n*      Report UMINF - 94.04, Department of Computing Science, Umea\n*      University, S-901 87 Umea, Sweden, 1994. Also as LAPACK Working\n*      Note 87. To appear in Numerical Algorithms, 1996.\n*\n*  =====================================================================\n*  Replaced various illegal calls to SCOPY by calls to SLASET, or by DO\n*  loops. Sven Hammarling, 1/5/02.\n*\n\n");
    return Qnil;
  }
  if (argc != 10)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 10)", argc);
  rb_wantq = argv[0];
  rb_wantz = argv[1];
  rb_a = argv[2];
  rb_b = argv[3];
  rb_q = argv[4];
  rb_ldq = argv[5];
  rb_z = argv[6];
  rb_j1 = argv[7];
  rb_n1 = argv[8];
  rb_n2 = argv[9];

  if (!NA_IsNArray(rb_a))
    rb_raise(rb_eArgError, "a (3th argument) must be NArray");
  if (NA_RANK(rb_a) != 2)
    rb_raise(rb_eArgError, "rank of a (3th argument) must be %d", 2);
  n = NA_SHAPE1(rb_a);
  lda = NA_SHAPE0(rb_a);
  if (NA_TYPE(rb_a) != NA_SFLOAT)
    rb_a = na_change_type(rb_a, NA_SFLOAT);
  a = NA_PTR_TYPE(rb_a, real*);
  wantz = (rb_wantz == Qtrue);
  if (!NA_IsNArray(rb_b))
    rb_raise(rb_eArgError, "b (4th argument) must be NArray");
  if (NA_RANK(rb_b) != 2)
    rb_raise(rb_eArgError, "rank of b (4th argument) must be %d", 2);
  if (NA_SHAPE1(rb_b) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of b must be the same as shape 1 of a");
  ldb = NA_SHAPE0(rb_b);
  if (NA_TYPE(rb_b) != NA_SFLOAT)
    rb_b = na_change_type(rb_b, NA_SFLOAT);
  b = NA_PTR_TYPE(rb_b, real*);
  n1 = NUM2INT(rb_n1);
  ldq = NUM2INT(rb_ldq);
  if (!NA_IsNArray(rb_q))
    rb_raise(rb_eArgError, "q (5th argument) must be NArray");
  if (NA_RANK(rb_q) != 2)
    rb_raise(rb_eArgError, "rank of q (5th argument) must be %d", 2);
  if (NA_SHAPE1(rb_q) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of q must be the same as shape 1 of a");
  ldz = NA_SHAPE0(rb_q);
  if (NA_TYPE(rb_q) != NA_SFLOAT)
    rb_q = na_change_type(rb_q, NA_SFLOAT);
  q = NA_PTR_TYPE(rb_q, real*);
  n2 = NUM2INT(rb_n2);
  if (!NA_IsNArray(rb_z))
    rb_raise(rb_eArgError, "z (7th argument) must be NArray");
  if (NA_RANK(rb_z) != 2)
    rb_raise(rb_eArgError, "rank of z (7th argument) must be %d", 2);
  if (NA_SHAPE1(rb_z) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of z must be the same as shape 1 of a");
  if (NA_SHAPE0(rb_z) != ldz)
    rb_raise(rb_eRuntimeError, "shape 0 of z must be the same as shape 0 of q");
  if (NA_TYPE(rb_z) != NA_SFLOAT)
    rb_z = na_change_type(rb_z, NA_SFLOAT);
  z = NA_PTR_TYPE(rb_z, real*);
  j1 = NUM2INT(rb_j1);
  wantq = (rb_wantq == Qtrue);
  lwork = MAX(1,(MAX(n*(n2+n1),(n2+n1)*(n2+n1)*2)));
  {
    int shape[2];
    shape[0] = lda;
    shape[1] = n;
    rb_a_out__ = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  a_out__ = NA_PTR_TYPE(rb_a_out__, real*);
  MEMCPY(a_out__, a, real, NA_TOTAL(rb_a));
  rb_a = rb_a_out__;
  a = a_out__;
  {
    int shape[2];
    shape[0] = ldb;
    shape[1] = n;
    rb_b_out__ = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  b_out__ = NA_PTR_TYPE(rb_b_out__, real*);
  MEMCPY(b_out__, b, real, NA_TOTAL(rb_b));
  rb_b = rb_b_out__;
  b = b_out__;
  {
    int shape[2];
    shape[0] = ldz;
    shape[1] = n;
    rb_q_out__ = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  q_out__ = NA_PTR_TYPE(rb_q_out__, real*);
  MEMCPY(q_out__, q, real, NA_TOTAL(rb_q));
  rb_q = rb_q_out__;
  q = q_out__;
  {
    int shape[2];
    shape[0] = ldz;
    shape[1] = n;
    rb_z_out__ = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  z_out__ = NA_PTR_TYPE(rb_z_out__, real*);
  MEMCPY(z_out__, z, real, NA_TOTAL(rb_z));
  rb_z = rb_z_out__;
  z = z_out__;
  work = ALLOC_N(real, (lwork));

  stgex2_(&wantq, &wantz, &n, a, &lda, b, &ldb, q, &ldq, z, &ldz, &j1, &n1, &n2, work, &lwork, &info);

  free(work);
  rb_info = INT2NUM(info);
  return rb_ary_new3(5, rb_info, rb_a, rb_b, rb_q, rb_z);
}

void
init_lapack_stgex2(VALUE mLapack){
  rb_define_module_function(mLapack, "stgex2", rb_stgex2, -1);
}

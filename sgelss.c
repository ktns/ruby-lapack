#include "rb_lapack.h"

extern VOID sgelss_(integer *m, integer *n, integer *nrhs, real *a, integer *lda, real *b, integer *ldb, real *s, real *rcond, integer *rank, real *work, integer *lwork, integer *info);

static VALUE
rb_sgelss(int argc, VALUE *argv, VALUE self){
  VALUE rb_m;
  integer m; 
  VALUE rb_a;
  real *a; 
  VALUE rb_b;
  real *b; 
  VALUE rb_rcond;
  real rcond; 
  VALUE rb_lwork;
  integer lwork; 
  VALUE rb_s;
  real *s; 
  VALUE rb_rank;
  integer rank; 
  VALUE rb_work;
  real *work; 
  VALUE rb_info;
  integer info; 
  VALUE rb_a_out__;
  real *a_out__;
  VALUE rb_b_out__;
  real *b_out__;

  integer lda;
  integer n;
  integer ldb;
  integer nrhs;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  s, rank, work, info, a, b = NumRu::Lapack.sgelss( m, a, b, rcond, lwork)\n    or\n  NumRu::Lapack.sgelss  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE SGELSS( M, N, NRHS, A, LDA, B, LDB, S, RCOND, RANK, WORK, LWORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  SGELSS computes the minimum norm solution to a real linear least\n*  squares problem:\n*\n*  Minimize 2-norm(| b - A*x |).\n*\n*  using the singular value decomposition (SVD) of A. A is an M-by-N\n*  matrix which may be rank-deficient.\n*\n*  Several right hand side vectors b and solution vectors x can be\n*  handled in a single call; they are stored as the columns of the\n*  M-by-NRHS right hand side matrix B and the N-by-NRHS solution matrix\n*  X.\n*\n*  The effective rank of A is determined by treating as zero those\n*  singular values which are less than RCOND times the largest singular\n*  value.\n*\n\n*  Arguments\n*  =========\n*\n*  M       (input) INTEGER\n*          The number of rows of the matrix A. M >= 0.\n*\n*  N       (input) INTEGER\n*          The number of columns of the matrix A. N >= 0.\n*\n*  NRHS    (input) INTEGER\n*          The number of right hand sides, i.e., the number of columns\n*          of the matrices B and X. NRHS >= 0.\n*\n*  A       (input/output) REAL array, dimension (LDA,N)\n*          On entry, the M-by-N matrix A.\n*          On exit, the first min(m,n) rows of A are overwritten with\n*          its right singular vectors, stored rowwise.\n*\n*  LDA     (input) INTEGER\n*          The leading dimension of the array A.  LDA >= max(1,M).\n*\n*  B       (input/output) REAL array, dimension (LDB,NRHS)\n*          On entry, the M-by-NRHS right hand side matrix B.\n*          On exit, B is overwritten by the N-by-NRHS solution\n*          matrix X.  If m >= n and RANK = n, the residual\n*          sum-of-squares for the solution in the i-th column is given\n*          by the sum of squares of elements n+1:m in that column.\n*\n*  LDB     (input) INTEGER\n*          The leading dimension of the array B. LDB >= max(1,max(M,N)).\n*\n*  S       (output) REAL array, dimension (min(M,N))\n*          The singular values of A in decreasing order.\n*          The condition number of A in the 2-norm = S(1)/S(min(m,n)).\n*\n*  RCOND   (input) REAL\n*          RCOND is used to determine the effective rank of A.\n*          Singular values S(i) <= RCOND*S(1) are treated as zero.\n*          If RCOND < 0, machine precision is used instead.\n*\n*  RANK    (output) INTEGER\n*          The effective rank of A, i.e., the number of singular values\n*          which are greater than RCOND*S(1).\n*\n*  WORK    (workspace/output) REAL array, dimension (MAX(1,LWORK))\n*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.\n*\n*  LWORK   (input) INTEGER\n*          The dimension of the array WORK. LWORK >= 1, and also:\n*          LWORK >= 3*min(M,N) + max( 2*min(M,N), max(M,N), NRHS )\n*          For good performance, LWORK should generally be larger.\n*\n*          If LWORK = -1, then a workspace query is assumed; the routine\n*          only calculates the optimal size of the WORK array, returns\n*          this value as the first entry of the WORK array, and no error\n*          message related to LWORK is issued by XERBLA.\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit\n*          < 0:  if INFO = -i, the i-th argument had an illegal value.\n*          > 0:  the algorithm for computing the SVD failed to converge;\n*                if INFO = i, i off-diagonal elements of an intermediate\n*                bidiagonal form did not converge to zero.\n*\n\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 5)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 5)", argc);
  rb_m = argv[0];
  rb_a = argv[1];
  rb_b = argv[2];
  rb_rcond = argv[3];
  rb_lwork = argv[4];

  rcond = (real)NUM2DBL(rb_rcond);
  if (!NA_IsNArray(rb_a))
    rb_raise(rb_eArgError, "a (2th argument) must be NArray");
  if (NA_RANK(rb_a) != 2)
    rb_raise(rb_eArgError, "rank of a (2th argument) must be %d", 2);
  n = NA_SHAPE1(rb_a);
  lda = NA_SHAPE0(rb_a);
  if (NA_TYPE(rb_a) != NA_SFLOAT)
    rb_a = na_change_type(rb_a, NA_SFLOAT);
  a = NA_PTR_TYPE(rb_a, real*);
  if (!NA_IsNArray(rb_b))
    rb_raise(rb_eArgError, "b (3th argument) must be NArray");
  if (NA_RANK(rb_b) != 2)
    rb_raise(rb_eArgError, "rank of b (3th argument) must be %d", 2);
  nrhs = NA_SHAPE1(rb_b);
  ldb = NA_SHAPE0(rb_b);
  if (NA_TYPE(rb_b) != NA_SFLOAT)
    rb_b = na_change_type(rb_b, NA_SFLOAT);
  b = NA_PTR_TYPE(rb_b, real*);
  m = NUM2INT(rb_m);
  lwork = NUM2INT(rb_lwork);
  {
    int shape[1];
    shape[0] = MIN(m,n);
    rb_s = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  s = NA_PTR_TYPE(rb_s, real*);
  {
    int shape[1];
    shape[0] = MAX(1,lwork);
    rb_work = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  work = NA_PTR_TYPE(rb_work, real*);
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
    shape[1] = nrhs;
    rb_b_out__ = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  b_out__ = NA_PTR_TYPE(rb_b_out__, real*);
  MEMCPY(b_out__, b, real, NA_TOTAL(rb_b));
  rb_b = rb_b_out__;
  b = b_out__;

  sgelss_(&m, &n, &nrhs, a, &lda, b, &ldb, s, &rcond, &rank, work, &lwork, &info);

  rb_rank = INT2NUM(rank);
  rb_info = INT2NUM(info);
  return rb_ary_new3(6, rb_s, rb_rank, rb_work, rb_info, rb_a, rb_b);
}

void
init_lapack_sgelss(VALUE mLapack){
  rb_define_module_function(mLapack, "sgelss", rb_sgelss, -1);
}

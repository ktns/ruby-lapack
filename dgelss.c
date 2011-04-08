#include "rb_lapack.h"

extern VOID dgelss_(integer *m, integer *n, integer *nrhs, doublereal *a, integer *lda, doublereal *b, integer *ldb, doublereal *s, doublereal *rcond, integer *rank, doublereal *work, integer *lwork, integer *info);

static VALUE sHelp, sUsage;

static VALUE
rblapack_dgelss(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_m;
  integer m; 
  VALUE rblapack_a;
  doublereal *a; 
  VALUE rblapack_b;
  doublereal *b; 
  VALUE rblapack_rcond;
  doublereal rcond; 
  VALUE rblapack_lwork;
  integer lwork; 
  VALUE rblapack_s;
  doublereal *s; 
  VALUE rblapack_rank;
  integer rank; 
  VALUE rblapack_work;
  doublereal *work; 
  VALUE rblapack_info;
  integer info; 
  VALUE rblapack_a_out__;
  doublereal *a_out__;
  VALUE rblapack_b_out__;
  doublereal *b_out__;

  integer lda;
  integer n;
  integer ldb;
  integer nrhs;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  s, rank, work, info, a, b = NumRu::Lapack.dgelss( m, a, b, rcond, lwork, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE DGELSS( M, N, NRHS, A, LDA, B, LDB, S, RCOND, RANK, WORK, LWORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  DGELSS computes the minimum norm solution to a real linear least\n*  squares problem:\n*\n*  Minimize 2-norm(| b - A*x |).\n*\n*  using the singular value decomposition (SVD) of A. A is an M-by-N\n*  matrix which may be rank-deficient.\n*\n*  Several right hand side vectors b and solution vectors x can be\n*  handled in a single call; they are stored as the columns of the\n*  M-by-NRHS right hand side matrix B and the N-by-NRHS solution matrix\n*  X.\n*\n*  The effective rank of A is determined by treating as zero those\n*  singular values which are less than RCOND times the largest singular\n*  value.\n*\n\n*  Arguments\n*  =========\n*\n*  M       (input) INTEGER\n*          The number of rows of the matrix A. M >= 0.\n*\n*  N       (input) INTEGER\n*          The number of columns of the matrix A. N >= 0.\n*\n*  NRHS    (input) INTEGER\n*          The number of right hand sides, i.e., the number of columns\n*          of the matrices B and X. NRHS >= 0.\n*\n*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)\n*          On entry, the M-by-N matrix A.\n*          On exit, the first min(m,n) rows of A are overwritten with\n*          its right singular vectors, stored rowwise.\n*\n*  LDA     (input) INTEGER\n*          The leading dimension of the array A.  LDA >= max(1,M).\n*\n*  B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)\n*          On entry, the M-by-NRHS right hand side matrix B.\n*          On exit, B is overwritten by the N-by-NRHS solution\n*          matrix X.  If m >= n and RANK = n, the residual\n*          sum-of-squares for the solution in the i-th column is given\n*          by the sum of squares of elements n+1:m in that column.\n*\n*  LDB     (input) INTEGER\n*          The leading dimension of the array B. LDB >= max(1,max(M,N)).\n*\n*  S       (output) DOUBLE PRECISION array, dimension (min(M,N))\n*          The singular values of A in decreasing order.\n*          The condition number of A in the 2-norm = S(1)/S(min(m,n)).\n*\n*  RCOND   (input) DOUBLE PRECISION\n*          RCOND is used to determine the effective rank of A.\n*          Singular values S(i) <= RCOND*S(1) are treated as zero.\n*          If RCOND < 0, machine precision is used instead.\n*\n*  RANK    (output) INTEGER\n*          The effective rank of A, i.e., the number of singular values\n*          which are greater than RCOND*S(1).\n*\n*  WORK    (workspace/output) DOUBLE PRECISION array, dimension (MAX(1,LWORK))\n*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.\n*\n*  LWORK   (input) INTEGER\n*          The dimension of the array WORK. LWORK >= 1, and also:\n*          LWORK >= 3*min(M,N) + max( 2*min(M,N), max(M,N), NRHS )\n*          For good performance, LWORK should generally be larger.\n*\n*          If LWORK = -1, then a workspace query is assumed; the routine\n*          only calculates the optimal size of the WORK array, returns\n*          this value as the first entry of the WORK array, and no error\n*          message related to LWORK is issued by XERBLA.\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit\n*          < 0:  if INFO = -i, the i-th argument had an illegal value.\n*          > 0:  the algorithm for computing the SVD failed to converge;\n*                if INFO = i, i off-diagonal elements of an intermediate\n*                bidiagonal form did not converge to zero.\n*\n\n*  =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  s, rank, work, info, a, b = NumRu::Lapack.dgelss( m, a, b, rcond, lwork, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 5)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 5)", argc);
  rblapack_m = argv[0];
  rblapack_a = argv[1];
  rblapack_b = argv[2];
  rblapack_rcond = argv[3];
  rblapack_lwork = argv[4];
  if (rb_options != Qnil) {
  }

  rcond = NUM2DBL(rblapack_rcond);
  if (!NA_IsNArray(rblapack_a))
    rb_raise(rb_eArgError, "a (2th argument) must be NArray");
  if (NA_RANK(rblapack_a) != 2)
    rb_raise(rb_eArgError, "rank of a (2th argument) must be %d", 2);
  n = NA_SHAPE1(rblapack_a);
  lda = NA_SHAPE0(rblapack_a);
  if (NA_TYPE(rblapack_a) != NA_DFLOAT)
    rblapack_a = na_change_type(rblapack_a, NA_DFLOAT);
  a = NA_PTR_TYPE(rblapack_a, doublereal*);
  if (!NA_IsNArray(rblapack_b))
    rb_raise(rb_eArgError, "b (3th argument) must be NArray");
  if (NA_RANK(rblapack_b) != 2)
    rb_raise(rb_eArgError, "rank of b (3th argument) must be %d", 2);
  nrhs = NA_SHAPE1(rblapack_b);
  ldb = NA_SHAPE0(rblapack_b);
  if (NA_TYPE(rblapack_b) != NA_DFLOAT)
    rblapack_b = na_change_type(rblapack_b, NA_DFLOAT);
  b = NA_PTR_TYPE(rblapack_b, doublereal*);
  m = NUM2INT(rblapack_m);
  lwork = NUM2INT(rblapack_lwork);
  {
    int shape[1];
    shape[0] = MIN(m,n);
    rblapack_s = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  s = NA_PTR_TYPE(rblapack_s, doublereal*);
  {
    int shape[1];
    shape[0] = MAX(1,lwork);
    rblapack_work = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  work = NA_PTR_TYPE(rblapack_work, doublereal*);
  {
    int shape[2];
    shape[0] = lda;
    shape[1] = n;
    rblapack_a_out__ = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  a_out__ = NA_PTR_TYPE(rblapack_a_out__, doublereal*);
  MEMCPY(a_out__, a, doublereal, NA_TOTAL(rblapack_a));
  rblapack_a = rblapack_a_out__;
  a = a_out__;
  {
    int shape[2];
    shape[0] = ldb;
    shape[1] = nrhs;
    rblapack_b_out__ = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  b_out__ = NA_PTR_TYPE(rblapack_b_out__, doublereal*);
  MEMCPY(b_out__, b, doublereal, NA_TOTAL(rblapack_b));
  rblapack_b = rblapack_b_out__;
  b = b_out__;

  dgelss_(&m, &n, &nrhs, a, &lda, b, &ldb, s, &rcond, &rank, work, &lwork, &info);

  rblapack_rank = INT2NUM(rank);
  rblapack_info = INT2NUM(info);
  return rb_ary_new3(6, rblapack_s, rblapack_rank, rblapack_work, rblapack_info, rblapack_a, rblapack_b);
}

void
init_lapack_dgelss(VALUE mLapack){
  rb_define_module_function(mLapack, "dgelss", rblapack_dgelss, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}

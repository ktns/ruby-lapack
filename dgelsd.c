#include "rb_lapack.h"

static VALUE
rb_dgelsd(int argc, VALUE *argv, VALUE self){
  VALUE rb_m;
  integer m; 
  VALUE rb_a;
  doublereal *a; 
  VALUE rb_b;
  doublereal *b; 
  VALUE rb_rcond;
  doublereal rcond; 
  VALUE rb_lwork;
  integer lwork; 
  VALUE rb_s;
  doublereal *s; 
  VALUE rb_rank;
  integer rank; 
  VALUE rb_work;
  doublereal *work; 
  VALUE rb_info;
  integer info; 
  VALUE rb_b_out__;
  doublereal *b_out__;
  integer *iwork;

  integer lda;
  integer n;
  integer ldb;
  integer nrhs;
  integer smlsiz;
  integer c__9;
  integer c__0;
  integer liwork;
  integer nlvl;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  s, rank, work, info, b = NumRu::Lapack.dgelsd( m, a, b, rcond, lwork)\n    or\n  NumRu::Lapack.dgelsd  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE DGELSD( M, N, NRHS, A, LDA, B, LDB, S, RCOND, RANK, WORK, LWORK, IWORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  DGELSD computes the minimum-norm solution to a real linear least\n*  squares problem:\n*      minimize 2-norm(| b - A*x |)\n*  using the singular value decomposition (SVD) of A. A is an M-by-N\n*  matrix which may be rank-deficient.\n*\n*  Several right hand side vectors b and solution vectors x can be\n*  handled in a single call; they are stored as the columns of the\n*  M-by-NRHS right hand side matrix B and the N-by-NRHS solution\n*  matrix X.\n*\n*  The problem is solved in three steps:\n*  (1) Reduce the coefficient matrix A to bidiagonal form with\n*      Householder transformations, reducing the original problem\n*      into a \"bidiagonal least squares problem\" (BLS)\n*  (2) Solve the BLS using a divide and conquer approach.\n*  (3) Apply back all the Householder tranformations to solve\n*      the original least squares problem.\n*\n*  The effective rank of A is determined by treating as zero those\n*  singular values which are less than RCOND times the largest singular\n*  value.\n*\n*  The divide and conquer algorithm makes very mild assumptions about\n*  floating point arithmetic. It will work on machines with a guard\n*  digit in add/subtract, or on those binary machines without guard\n*  digits which subtract like the Cray X-MP, Cray Y-MP, Cray C-90, or\n*  Cray-2. It could conceivably fail on hexadecimal or decimal machines\n*  without guard digits, but we know of none.\n*\n\n*  Arguments\n*  =========\n*\n*  M       (input) INTEGER\n*          The number of rows of A. M >= 0.\n*\n*  N       (input) INTEGER\n*          The number of columns of A. N >= 0.\n*\n*  NRHS    (input) INTEGER\n*          The number of right hand sides, i.e., the number of columns\n*          of the matrices B and X. NRHS >= 0.\n*\n*  A       (input) DOUBLE PRECISION array, dimension (LDA,N)\n*          On entry, the M-by-N matrix A.\n*          On exit, A has been destroyed.\n*\n*  LDA     (input) INTEGER\n*          The leading dimension of the array A.  LDA >= max(1,M).\n*\n*  B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)\n*          On entry, the M-by-NRHS right hand side matrix B.\n*          On exit, B is overwritten by the N-by-NRHS solution\n*          matrix X.  If m >= n and RANK = n, the residual\n*          sum-of-squares for the solution in the i-th column is given\n*          by the sum of squares of elements n+1:m in that column.\n*\n*  LDB     (input) INTEGER\n*          The leading dimension of the array B. LDB >= max(1,max(M,N)).\n*\n*  S       (output) DOUBLE PRECISION array, dimension (min(M,N))\n*          The singular values of A in decreasing order.\n*          The condition number of A in the 2-norm = S(1)/S(min(m,n)).\n*\n*  RCOND   (input) DOUBLE PRECISION\n*          RCOND is used to determine the effective rank of A.\n*          Singular values S(i) <= RCOND*S(1) are treated as zero.\n*          If RCOND < 0, machine precision is used instead.\n*\n*  RANK    (output) INTEGER\n*          The effective rank of A, i.e., the number of singular values\n*          which are greater than RCOND*S(1).\n*\n*  WORK    (workspace/output) DOUBLE PRECISION array, dimension (MAX(1,LWORK))\n*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.\n*\n*  LWORK   (input) INTEGER\n*          The dimension of the array WORK. LWORK must be at least 1.\n*          The exact minimum amount of workspace needed depends on M,\n*          N and NRHS. As long as LWORK is at least\n*              12*N + 2*N*SMLSIZ + 8*N*NLVL + N*NRHS + (SMLSIZ+1)**2,\n*          if M is greater than or equal to N or\n*              12*M + 2*M*SMLSIZ + 8*M*NLVL + M*NRHS + (SMLSIZ+1)**2,\n*          if M is less than N, the code will execute correctly.\n*          SMLSIZ is returned by ILAENV and is equal to the maximum\n*          size of the subproblems at the bottom of the computation\n*          tree (usually about 25), and\n*             NLVL = MAX( 0, INT( LOG_2( MIN( M,N )/(SMLSIZ+1) ) ) + 1 )\n*          For good performance, LWORK should generally be larger.\n*\n*          If LWORK = -1, then a workspace query is assumed; the routine\n*          only calculates the optimal size of the WORK array, returns\n*          this value as the first entry of the WORK array, and no error\n*          message related to LWORK is issued by XERBLA.\n*\n*  IWORK   (workspace) INTEGER array, dimension (MAX(1,LIWORK))\n*          LIWORK >= 3 * MINMN * NLVL + 11 * MINMN,\n*          where MINMN = MIN( M,N ).\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit\n*          < 0:  if INFO = -i, the i-th argument had an illegal value.\n*          > 0:  the algorithm for computing the SVD failed to converge;\n*                if INFO = i, i off-diagonal elements of an intermediate\n*                bidiagonal form did not converge to zero.\n*\n\n*  Further Details\n*  ===============\n*\n*  Based on contributions by\n*     Ming Gu and Ren-Cang Li, Computer Science Division, University of\n*       California at Berkeley, USA\n*     Osni Marques, LBNL/NERSC, USA\n*\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 5)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 5)", argc);
  rb_m = argv[0];
  rb_a = argv[1];
  rb_b = argv[2];
  rb_rcond = argv[3];
  rb_lwork = argv[4];

  m = NUM2INT(rb_m);
  rcond = NUM2DBL(rb_rcond);
  lwork = NUM2INT(rb_lwork);
  if (!NA_IsNArray(rb_a))
    rb_raise(rb_eArgError, "a (2th argument) must be NArray");
  if (NA_RANK(rb_a) != 2)
    rb_raise(rb_eArgError, "rank of a (2th argument) must be %d", 2);
  lda = NA_SHAPE0(rb_a);
  n = NA_SHAPE1(rb_a);
  if (NA_TYPE(rb_a) != NA_DFLOAT)
    rb_a = na_change_type(rb_a, NA_DFLOAT);
  a = NA_PTR_TYPE(rb_a, doublereal*);
  if (!NA_IsNArray(rb_b))
    rb_raise(rb_eArgError, "b (3th argument) must be NArray");
  if (NA_RANK(rb_b) != 2)
    rb_raise(rb_eArgError, "rank of b (3th argument) must be %d", 2);
  ldb = NA_SHAPE0(rb_b);
  nrhs = NA_SHAPE1(rb_b);
  if (NA_TYPE(rb_b) != NA_DFLOAT)
    rb_b = na_change_type(rb_b, NA_DFLOAT);
  b = NA_PTR_TYPE(rb_b, doublereal*);
  {
    int shape[1];
    shape[0] = MIN(m,n);
    rb_s = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  s = NA_PTR_TYPE(rb_s, doublereal*);
  {
    int shape[1];
    shape[0] = MAX(1,lwork);
    rb_work = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  work = NA_PTR_TYPE(rb_work, doublereal*);
  {
    int shape[2];
    shape[0] = ldb;
    shape[1] = nrhs;
    rb_b_out__ = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  b_out__ = NA_PTR_TYPE(rb_b_out__, doublereal*);
  MEMCPY(b_out__, b, doublereal, NA_TOTAL(rb_b));
  rb_b = rb_b_out__;
  b = b_out__;
  c__0 = 0;
  c__9 = 9;
  smlsiz = ilaenv_(&c__9,"DGELSD"," ",&c__0,&c__0,&c__0,&c__0);
  nlvl = MAX(0,((int)(log(((double)(MIN(m,n)))/(smlsiz+1))/log(2.0))+1));
  liwork = 3*(MIN(m,n))*nlvl+11*(MIN(m,n));
  iwork = ALLOC_N(integer, (MAX(1,liwork)));

  dgelsd_(&m, &n, &nrhs, a, &lda, b, &ldb, s, &rcond, &rank, work, &lwork, iwork, &info);

  free(iwork);
  rb_rank = INT2NUM(rank);
  rb_info = INT2NUM(info);
  return rb_ary_new3(5, rb_s, rb_rank, rb_work, rb_info, rb_b);
}

void
init_lapack_dgelsd(VALUE mLapack){
  rb_define_module_function(mLapack, "dgelsd", rb_dgelsd, -1);
}

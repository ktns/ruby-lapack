#include "rb_lapack.h"

extern VOID cgelsd_(integer *m, integer *n, integer *nrhs, complex *a, integer *lda, complex *b, integer *ldb, real *s, real *rcond, integer *rank, complex *work, integer *lwork, real *rwork, integer *iwork, integer *info);

static VALUE sHelp, sUsage;

static VALUE
rblapack_cgelsd(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_m;
  integer m; 
  VALUE rblapack_a;
  complex *a; 
  VALUE rblapack_b;
  complex *b; 
  VALUE rblapack_rcond;
  real rcond; 
  VALUE rblapack_lwork;
  integer lwork; 
  VALUE rblapack_s;
  real *s; 
  VALUE rblapack_rank;
  integer rank; 
  VALUE rblapack_work;
  complex *work; 
  VALUE rblapack_info;
  integer info; 
  VALUE rblapack_a_out__;
  complex *a_out__;
  VALUE rblapack_b_out__;
  complex *b_out__;
  real *rwork;
  integer *iwork;

  integer lda;
  integer n;
  integer ldb;
  integer nrhs;
  integer c__9;
  integer c__0;
  integer lrwork;
  integer liwork;
  integer smlsiz;
  integer nlvl;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  s, rank, work, info, a, b = NumRu::Lapack.cgelsd( m, a, b, rcond, lwork, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE CGELSD( M, N, NRHS, A, LDA, B, LDB, S, RCOND, RANK, WORK, LWORK, RWORK, IWORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  CGELSD computes the minimum-norm solution to a real linear least\n*  squares problem:\n*      minimize 2-norm(| b - A*x |)\n*  using the singular value decomposition (SVD) of A. A is an M-by-N\n*  matrix which may be rank-deficient.\n*\n*  Several right hand side vectors b and solution vectors x can be\n*  handled in a single call; they are stored as the columns of the\n*  M-by-NRHS right hand side matrix B and the N-by-NRHS solution\n*  matrix X.\n*\n*  The problem is solved in three steps:\n*  (1) Reduce the coefficient matrix A to bidiagonal form with\n*      Householder tranformations, reducing the original problem\n*      into a \"bidiagonal least squares problem\" (BLS)\n*  (2) Solve the BLS using a divide and conquer approach.\n*  (3) Apply back all the Householder tranformations to solve\n*      the original least squares problem.\n*\n*  The effective rank of A is determined by treating as zero those\n*  singular values which are less than RCOND times the largest singular\n*  value.\n*\n*  The divide and conquer algorithm makes very mild assumptions about\n*  floating point arithmetic. It will work on machines with a guard\n*  digit in add/subtract, or on those binary machines without guard\n*  digits which subtract like the Cray X-MP, Cray Y-MP, Cray C-90, or\n*  Cray-2. It could conceivably fail on hexadecimal or decimal machines\n*  without guard digits, but we know of none.\n*\n\n*  Arguments\n*  =========\n*\n*  M       (input) INTEGER\n*          The number of rows of the matrix A. M >= 0.\n*\n*  N       (input) INTEGER\n*          The number of columns of the matrix A. N >= 0.\n*\n*  NRHS    (input) INTEGER\n*          The number of right hand sides, i.e., the number of columns\n*          of the matrices B and X. NRHS >= 0.\n*\n*  A       (input/output) COMPLEX array, dimension (LDA,N)\n*          On entry, the M-by-N matrix A.\n*          On exit, A has been destroyed.\n*\n*  LDA     (input) INTEGER\n*          The leading dimension of the array A. LDA >= max(1,M).\n*\n*  B       (input/output) COMPLEX array, dimension (LDB,NRHS)\n*          On entry, the M-by-NRHS right hand side matrix B.\n*          On exit, B is overwritten by the N-by-NRHS solution matrix X.\n*          If m >= n and RANK = n, the residual sum-of-squares for\n*          the solution in the i-th column is given by the sum of\n*          squares of the modulus of elements n+1:m in that column.\n*\n*  LDB     (input) INTEGER\n*          The leading dimension of the array B.  LDB >= max(1,M,N).\n*\n*  S       (output) REAL array, dimension (min(M,N))\n*          The singular values of A in decreasing order.\n*          The condition number of A in the 2-norm = S(1)/S(min(m,n)).\n*\n*  RCOND   (input) REAL\n*          RCOND is used to determine the effective rank of A.\n*          Singular values S(i) <= RCOND*S(1) are treated as zero.\n*          If RCOND < 0, machine precision is used instead.\n*\n*  RANK    (output) INTEGER\n*          The effective rank of A, i.e., the number of singular values\n*          which are greater than RCOND*S(1).\n*\n*  WORK    (workspace/output) COMPLEX array, dimension (MAX(1,LWORK))\n*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.\n*\n*  LWORK   (input) INTEGER\n*          The dimension of the array WORK. LWORK must be at least 1.\n*          The exact minimum amount of workspace needed depends on M,\n*          N and NRHS. As long as LWORK is at least\n*              2 * N + N * NRHS\n*          if M is greater than or equal to N or\n*              2 * M + M * NRHS\n*          if M is less than N, the code will execute correctly.\n*          For good performance, LWORK should generally be larger.\n*\n*          If LWORK = -1, then a workspace query is assumed; the routine\n*          only calculates the optimal size of the array WORK and the\n*          minimum sizes of the arrays RWORK and IWORK, and returns\n*          these values as the first entries of the WORK, RWORK and\n*          IWORK arrays, and no error message related to LWORK is issued\n*          by XERBLA.\n*\n*  RWORK   (workspace) REAL array, dimension (MAX(1,LRWORK))\n*          LRWORK >=\n*             10*N + 2*N*SMLSIZ + 8*N*NLVL + 3*SMLSIZ*NRHS +\n*             MAX( (SMLSIZ+1)**2, N*(1+NRHS) + 2*NRHS )\n*          if M is greater than or equal to N or\n*             10*M + 2*M*SMLSIZ + 8*M*NLVL + 3*SMLSIZ*NRHS +\n*             MAX( (SMLSIZ+1)**2, N*(1+NRHS) + 2*NRHS )\n*          if M is less than N, the code will execute correctly.\n*          SMLSIZ is returned by ILAENV and is equal to the maximum\n*          size of the subproblems at the bottom of the computation\n*          tree (usually about 25), and\n*             NLVL = MAX( 0, INT( LOG_2( MIN( M,N )/(SMLSIZ+1) ) ) + 1 )\n*          On exit, if INFO = 0, RWORK(1) returns the minimum LRWORK.\n*\n*  IWORK   (workspace) INTEGER array, dimension (MAX(1,LIWORK))\n*          LIWORK >= max(1, 3*MINMN*NLVL + 11*MINMN),\n*          where MINMN = MIN( M,N ).\n*          On exit, if INFO = 0, IWORK(1) returns the minimum LIWORK.\n*\n*  INFO    (output) INTEGER\n*          = 0: successful exit\n*          < 0: if INFO = -i, the i-th argument had an illegal value.\n*          > 0:  the algorithm for computing the SVD failed to converge;\n*                if INFO = i, i off-diagonal elements of an intermediate\n*                bidiagonal form did not converge to zero.\n*\n\n*  Further Details\n*  ===============\n*\n*  Based on contributions by\n*     Ming Gu and Ren-Cang Li, Computer Science Division, University of\n*       California at Berkeley, USA\n*     Osni Marques, LBNL/NERSC, USA\n*\n*  =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  s, rank, work, info, a, b = NumRu::Lapack.cgelsd( m, a, b, rcond, lwork, [:usage => usage, :help => help])\n");
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

  c__9 = 9;
  rcond = (real)NUM2DBL(rblapack_rcond);
  if (!NA_IsNArray(rblapack_a))
    rb_raise(rb_eArgError, "a (2th argument) must be NArray");
  if (NA_RANK(rblapack_a) != 2)
    rb_raise(rb_eArgError, "rank of a (2th argument) must be %d", 2);
  n = NA_SHAPE1(rblapack_a);
  lda = NA_SHAPE0(rblapack_a);
  if (NA_TYPE(rblapack_a) != NA_SCOMPLEX)
    rblapack_a = na_change_type(rblapack_a, NA_SCOMPLEX);
  a = NA_PTR_TYPE(rblapack_a, complex*);
  c__0 = 0;
  if (!NA_IsNArray(rblapack_b))
    rb_raise(rb_eArgError, "b (3th argument) must be NArray");
  if (NA_RANK(rblapack_b) != 2)
    rb_raise(rb_eArgError, "rank of b (3th argument) must be %d", 2);
  nrhs = NA_SHAPE1(rblapack_b);
  ldb = NA_SHAPE0(rblapack_b);
  if (NA_TYPE(rblapack_b) != NA_SCOMPLEX)
    rblapack_b = na_change_type(rblapack_b, NA_SCOMPLEX);
  b = NA_PTR_TYPE(rblapack_b, complex*);
  lwork = NUM2INT(rblapack_lwork);
  m = NUM2INT(rblapack_m);
  smlsiz = ilaenv_(&c__9,"CGELSD"," ",&c__0,&c__0,&c__0,&c__0);
  nlvl = MAX(0,(int)(log(1.0*MIN(m,n)/(smlsiz+1))/log(2.0)));
  liwork = MAX(1,3*(MIN(m,n))*nlvl+11*(MIN(m,n)));
  lrwork = m>=n ? 10*n+2*n*smlsiz+8*n*nlvl+3*smlsiz*nrhs+(smlsiz+1)*(smlsiz+1) : 10*m+2*m*smlsiz+8*m*nlvl+2*smlsiz*nrhs+(smlsiz+1)*(smlsiz+1);
  {
    int shape[1];
    shape[0] = MIN(m,n);
    rblapack_s = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  s = NA_PTR_TYPE(rblapack_s, real*);
  {
    int shape[1];
    shape[0] = MAX(1,lwork);
    rblapack_work = na_make_object(NA_SCOMPLEX, 1, shape, cNArray);
  }
  work = NA_PTR_TYPE(rblapack_work, complex*);
  {
    int shape[2];
    shape[0] = lda;
    shape[1] = n;
    rblapack_a_out__ = na_make_object(NA_SCOMPLEX, 2, shape, cNArray);
  }
  a_out__ = NA_PTR_TYPE(rblapack_a_out__, complex*);
  MEMCPY(a_out__, a, complex, NA_TOTAL(rblapack_a));
  rblapack_a = rblapack_a_out__;
  a = a_out__;
  {
    int shape[2];
    shape[0] = ldb;
    shape[1] = nrhs;
    rblapack_b_out__ = na_make_object(NA_SCOMPLEX, 2, shape, cNArray);
  }
  b_out__ = NA_PTR_TYPE(rblapack_b_out__, complex*);
  MEMCPY(b_out__, b, complex, NA_TOTAL(rblapack_b));
  rblapack_b = rblapack_b_out__;
  b = b_out__;
  rwork = ALLOC_N(real, (MAX(1,lrwork)));
  iwork = ALLOC_N(integer, (MAX(1,liwork)));

  cgelsd_(&m, &n, &nrhs, a, &lda, b, &ldb, s, &rcond, &rank, work, &lwork, rwork, iwork, &info);

  free(rwork);
  free(iwork);
  rblapack_rank = INT2NUM(rank);
  rblapack_info = INT2NUM(info);
  return rb_ary_new3(6, rblapack_s, rblapack_rank, rblapack_work, rblapack_info, rblapack_a, rblapack_b);
}

void
init_lapack_cgelsd(VALUE mLapack){
  rb_define_module_function(mLapack, "cgelsd", rblapack_cgelsd, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}

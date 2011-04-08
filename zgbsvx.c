#include "rb_lapack.h"

extern VOID zgbsvx_(char *fact, char *trans, integer *n, integer *kl, integer *ku, integer *nrhs, doublecomplex *ab, integer *ldab, doublecomplex *afb, integer *ldafb, integer *ipiv, char *equed, doublereal *r, doublereal *c, doublecomplex *b, integer *ldb, doublecomplex *x, integer *ldx, doublereal *rcond, doublereal *ferr, doublereal *berr, doublecomplex *work, doublereal *rwork, integer *info);

static VALUE sHelp, sUsage;

static VALUE
rblapack_zgbsvx(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_fact;
  char fact; 
  VALUE rblapack_trans;
  char trans; 
  VALUE rblapack_kl;
  integer kl; 
  VALUE rblapack_ku;
  integer ku; 
  VALUE rblapack_ab;
  doublecomplex *ab; 
  VALUE rblapack_afb;
  doublecomplex *afb; 
  VALUE rblapack_ipiv;
  integer *ipiv; 
  VALUE rblapack_equed;
  char equed; 
  VALUE rblapack_r;
  doublereal *r; 
  VALUE rblapack_c;
  doublereal *c; 
  VALUE rblapack_b;
  doublecomplex *b; 
  VALUE rblapack_x;
  doublecomplex *x; 
  VALUE rblapack_rcond;
  doublereal rcond; 
  VALUE rblapack_ferr;
  doublereal *ferr; 
  VALUE rblapack_berr;
  doublereal *berr; 
  VALUE rblapack_rwork;
  doublereal *rwork; 
  VALUE rblapack_info;
  integer info; 
  VALUE rblapack_ab_out__;
  doublecomplex *ab_out__;
  VALUE rblapack_afb_out__;
  doublecomplex *afb_out__;
  VALUE rblapack_ipiv_out__;
  integer *ipiv_out__;
  VALUE rblapack_r_out__;
  doublereal *r_out__;
  VALUE rblapack_c_out__;
  doublereal *c_out__;
  VALUE rblapack_b_out__;
  doublecomplex *b_out__;
  doublecomplex *work;

  integer ldab;
  integer n;
  integer ldafb;
  integer ldb;
  integer nrhs;
  integer ldx;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  x, rcond, ferr, berr, rwork, info, ab, afb, ipiv, equed, r, c, b = NumRu::Lapack.zgbsvx( fact, trans, kl, ku, ab, afb, ipiv, equed, r, c, b, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE ZGBSVX( FACT, TRANS, N, KL, KU, NRHS, AB, LDAB, AFB, LDAFB, IPIV, EQUED, R, C, B, LDB, X, LDX, RCOND, FERR, BERR, WORK, RWORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  ZGBSVX uses the LU factorization to compute the solution to a complex\n*  system of linear equations A * X = B, A**T * X = B, or A**H * X = B,\n*  where A is a band matrix of order N with KL subdiagonals and KU\n*  superdiagonals, and X and B are N-by-NRHS matrices.\n*\n*  Error bounds on the solution and a condition estimate are also\n*  provided.\n*\n*  Description\n*  ===========\n*\n*  The following steps are performed by this subroutine:\n*\n*  1. If FACT = 'E', real scaling factors are computed to equilibrate\n*     the system:\n*        TRANS = 'N':  diag(R)*A*diag(C)     *inv(diag(C))*X = diag(R)*B\n*        TRANS = 'T': (diag(R)*A*diag(C))**T *inv(diag(R))*X = diag(C)*B\n*        TRANS = 'C': (diag(R)*A*diag(C))**H *inv(diag(R))*X = diag(C)*B\n*     Whether or not the system will be equilibrated depends on the\n*     scaling of the matrix A, but if equilibration is used, A is\n*     overwritten by diag(R)*A*diag(C) and B by diag(R)*B (if TRANS='N')\n*     or diag(C)*B (if TRANS = 'T' or 'C').\n*\n*  2. If FACT = 'N' or 'E', the LU decomposition is used to factor the\n*     matrix A (after equilibration if FACT = 'E') as\n*        A = L * U,\n*     where L is a product of permutation and unit lower triangular\n*     matrices with KL subdiagonals, and U is upper triangular with\n*     KL+KU superdiagonals.\n*\n*  3. If some U(i,i)=0, so that U is exactly singular, then the routine\n*     returns with INFO = i. Otherwise, the factored form of A is used\n*     to estimate the condition number of the matrix A.  If the\n*     reciprocal of the condition number is less than machine precision,\n*     INFO = N+1 is returned as a warning, but the routine still goes on\n*     to solve for X and compute error bounds as described below.\n*\n*  4. The system of equations is solved for X using the factored form\n*     of A.\n*\n*  5. Iterative refinement is applied to improve the computed solution\n*     matrix and calculate error bounds and backward error estimates\n*     for it.\n*\n*  6. If equilibration was used, the matrix X is premultiplied by\n*     diag(C) (if TRANS = 'N') or diag(R) (if TRANS = 'T' or 'C') so\n*     that it solves the original system before equilibration.\n*\n\n*  Arguments\n*  =========\n*\n*  FACT    (input) CHARACTER*1\n*          Specifies whether or not the factored form of the matrix A is\n*          supplied on entry, and if not, whether the matrix A should be\n*          equilibrated before it is factored.\n*          = 'F':  On entry, AFB and IPIV contain the factored form of\n*                  A.  If EQUED is not 'N', the matrix A has been\n*                  equilibrated with scaling factors given by R and C.\n*                  AB, AFB, and IPIV are not modified.\n*          = 'N':  The matrix A will be copied to AFB and factored.\n*          = 'E':  The matrix A will be equilibrated if necessary, then\n*                  copied to AFB and factored.\n*\n*  TRANS   (input) CHARACTER*1\n*          Specifies the form of the system of equations.\n*          = 'N':  A * X = B     (No transpose)\n*          = 'T':  A**T * X = B  (Transpose)\n*          = 'C':  A**H * X = B  (Conjugate transpose)\n*\n*  N       (input) INTEGER\n*          The number of linear equations, i.e., the order of the\n*          matrix A.  N >= 0.\n*\n*  KL      (input) INTEGER\n*          The number of subdiagonals within the band of A.  KL >= 0.\n*\n*  KU      (input) INTEGER\n*          The number of superdiagonals within the band of A.  KU >= 0.\n*\n*  NRHS    (input) INTEGER\n*          The number of right hand sides, i.e., the number of columns\n*          of the matrices B and X.  NRHS >= 0.\n*\n*  AB      (input/output) COMPLEX*16 array, dimension (LDAB,N)\n*          On entry, the matrix A in band storage, in rows 1 to KL+KU+1.\n*          The j-th column of A is stored in the j-th column of the\n*          array AB as follows:\n*          AB(KU+1+i-j,j) = A(i,j) for max(1,j-KU)<=i<=min(N,j+kl)\n*\n*          If FACT = 'F' and EQUED is not 'N', then A must have been\n*          equilibrated by the scaling factors in R and/or C.  AB is not\n*          modified if FACT = 'F' or 'N', or if FACT = 'E' and\n*          EQUED = 'N' on exit.\n*\n*          On exit, if EQUED .ne. 'N', A is scaled as follows:\n*          EQUED = 'R':  A := diag(R) * A\n*          EQUED = 'C':  A := A * diag(C)\n*          EQUED = 'B':  A := diag(R) * A * diag(C).\n*\n*  LDAB    (input) INTEGER\n*          The leading dimension of the array AB.  LDAB >= KL+KU+1.\n*\n*  AFB     (input or output) COMPLEX*16 array, dimension (LDAFB,N)\n*          If FACT = 'F', then AFB is an input argument and on entry\n*          contains details of the LU factorization of the band matrix\n*          A, as computed by ZGBTRF.  U is stored as an upper triangular\n*          band matrix with KL+KU superdiagonals in rows 1 to KL+KU+1,\n*          and the multipliers used during the factorization are stored\n*          in rows KL+KU+2 to 2*KL+KU+1.  If EQUED .ne. 'N', then AFB is\n*          the factored form of the equilibrated matrix A.\n*\n*          If FACT = 'N', then AFB is an output argument and on exit\n*          returns details of the LU factorization of A.\n*\n*          If FACT = 'E', then AFB is an output argument and on exit\n*          returns details of the LU factorization of the equilibrated\n*          matrix A (see the description of AB for the form of the\n*          equilibrated matrix).\n*\n*  LDAFB   (input) INTEGER\n*          The leading dimension of the array AFB.  LDAFB >= 2*KL+KU+1.\n*\n*  IPIV    (input or output) INTEGER array, dimension (N)\n*          If FACT = 'F', then IPIV is an input argument and on entry\n*          contains the pivot indices from the factorization A = L*U\n*          as computed by ZGBTRF; row i of the matrix was interchanged\n*          with row IPIV(i).\n*\n*          If FACT = 'N', then IPIV is an output argument and on exit\n*          contains the pivot indices from the factorization A = L*U\n*          of the original matrix A.\n*\n*          If FACT = 'E', then IPIV is an output argument and on exit\n*          contains the pivot indices from the factorization A = L*U\n*          of the equilibrated matrix A.\n*\n*  EQUED   (input or output) CHARACTER*1\n*          Specifies the form of equilibration that was done.\n*          = 'N':  No equilibration (always true if FACT = 'N').\n*          = 'R':  Row equilibration, i.e., A has been premultiplied by\n*                  diag(R).\n*          = 'C':  Column equilibration, i.e., A has been postmultiplied\n*                  by diag(C).\n*          = 'B':  Both row and column equilibration, i.e., A has been\n*                  replaced by diag(R) * A * diag(C).\n*          EQUED is an input argument if FACT = 'F'; otherwise, it is an\n*          output argument.\n*\n*  R       (input or output) DOUBLE PRECISION array, dimension (N)\n*          The row scale factors for A.  If EQUED = 'R' or 'B', A is\n*          multiplied on the left by diag(R); if EQUED = 'N' or 'C', R\n*          is not accessed.  R is an input argument if FACT = 'F';\n*          otherwise, R is an output argument.  If FACT = 'F' and\n*          EQUED = 'R' or 'B', each element of R must be positive.\n*\n*  C       (input or output) DOUBLE PRECISION array, dimension (N)\n*          The column scale factors for A.  If EQUED = 'C' or 'B', A is\n*          multiplied on the right by diag(C); if EQUED = 'N' or 'R', C\n*          is not accessed.  C is an input argument if FACT = 'F';\n*          otherwise, C is an output argument.  If FACT = 'F' and\n*          EQUED = 'C' or 'B', each element of C must be positive.\n*\n*  B       (input/output) COMPLEX*16 array, dimension (LDB,NRHS)\n*          On entry, the right hand side matrix B.\n*          On exit,\n*          if EQUED = 'N', B is not modified;\n*          if TRANS = 'N' and EQUED = 'R' or 'B', B is overwritten by\n*          diag(R)*B;\n*          if TRANS = 'T' or 'C' and EQUED = 'C' or 'B', B is\n*          overwritten by diag(C)*B.\n*\n*  LDB     (input) INTEGER\n*          The leading dimension of the array B.  LDB >= max(1,N).\n*\n*  X       (output) COMPLEX*16 array, dimension (LDX,NRHS)\n*          If INFO = 0 or INFO = N+1, the N-by-NRHS solution matrix X\n*          to the original system of equations.  Note that A and B are\n*          modified on exit if EQUED .ne. 'N', and the solution to the\n*          equilibrated system is inv(diag(C))*X if TRANS = 'N' and\n*          EQUED = 'C' or 'B', or inv(diag(R))*X if TRANS = 'T' or 'C'\n*          and EQUED = 'R' or 'B'.\n*\n*  LDX     (input) INTEGER\n*          The leading dimension of the array X.  LDX >= max(1,N).\n*\n*  RCOND   (output) DOUBLE PRECISION\n*          The estimate of the reciprocal condition number of the matrix\n*          A after equilibration (if done).  If RCOND is less than the\n*          machine precision (in particular, if RCOND = 0), the matrix\n*          is singular to working precision.  This condition is\n*          indicated by a return code of INFO > 0.\n*\n*  FERR    (output) DOUBLE PRECISION array, dimension (NRHS)\n*          The estimated forward error bound for each solution vector\n*          X(j) (the j-th column of the solution matrix X).\n*          If XTRUE is the true solution corresponding to X(j), FERR(j)\n*          is an estimated upper bound for the magnitude of the largest\n*          element in (X(j) - XTRUE) divided by the magnitude of the\n*          largest element in X(j).  The estimate is as reliable as\n*          the estimate for RCOND, and is almost always a slight\n*          overestimate of the true error.\n*\n*  BERR    (output) DOUBLE PRECISION array, dimension (NRHS)\n*          The componentwise relative backward error of each solution\n*          vector X(j) (i.e., the smallest relative change in\n*          any element of A or B that makes X(j) an exact solution).\n*\n*  WORK    (workspace) COMPLEX*16 array, dimension (2*N)\n*\n*  RWORK   (workspace/output) DOUBLE PRECISION array, dimension (N)\n*          On exit, RWORK(1) contains the reciprocal pivot growth\n*          factor norm(A)/norm(U). The \"max absolute element\" norm is\n*          used. If RWORK(1) is much less than 1, then the stability\n*          of the LU factorization of the (equilibrated) matrix A\n*          could be poor. This also means that the solution X, condition\n*          estimator RCOND, and forward error bound FERR could be\n*          unreliable. If factorization fails with 0<INFO<=N, then\n*          RWORK(1) contains the reciprocal pivot growth factor for the\n*          leading INFO columns of A.\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit\n*          < 0:  if INFO = -i, the i-th argument had an illegal value\n*          > 0:  if INFO = i, and i is\n*                <= N:  U(i,i) is exactly zero.  The factorization\n*                       has been completed, but the factor U is exactly\n*                       singular, so the solution and error bounds\n*                       could not be computed. RCOND = 0 is returned.\n*                = N+1: U is nonsingular, but RCOND is less than machine\n*                       precision, meaning that the matrix is singular\n*                       to working precision.  Nevertheless, the\n*                       solution and error bounds are computed because\n*                       there are a number of situations where the\n*                       computed solution can be more accurate than the\n*                       value of RCOND would suggest.\n*\n\n*  =====================================================================\n*  Moved setting of INFO = N+1 so INFO does not subsequently get\n*  overwritten.  Sven, 17 Mar 05. \n*  =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  x, rcond, ferr, berr, rwork, info, ab, afb, ipiv, equed, r, c, b = NumRu::Lapack.zgbsvx( fact, trans, kl, ku, ab, afb, ipiv, equed, r, c, b, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 11)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 11)", argc);
  rblapack_fact = argv[0];
  rblapack_trans = argv[1];
  rblapack_kl = argv[2];
  rblapack_ku = argv[3];
  rblapack_ab = argv[4];
  rblapack_afb = argv[5];
  rblapack_ipiv = argv[6];
  rblapack_equed = argv[7];
  rblapack_r = argv[8];
  rblapack_c = argv[9];
  rblapack_b = argv[10];
  if (rb_options != Qnil) {
  }

  trans = StringValueCStr(rblapack_trans)[0];
  if (!NA_IsNArray(rblapack_ipiv))
    rb_raise(rb_eArgError, "ipiv (7th argument) must be NArray");
  if (NA_RANK(rblapack_ipiv) != 1)
    rb_raise(rb_eArgError, "rank of ipiv (7th argument) must be %d", 1);
  n = NA_SHAPE0(rblapack_ipiv);
  if (NA_TYPE(rblapack_ipiv) != NA_LINT)
    rblapack_ipiv = na_change_type(rblapack_ipiv, NA_LINT);
  ipiv = NA_PTR_TYPE(rblapack_ipiv, integer*);
  if (!NA_IsNArray(rblapack_ab))
    rb_raise(rb_eArgError, "ab (5th argument) must be NArray");
  if (NA_RANK(rblapack_ab) != 2)
    rb_raise(rb_eArgError, "rank of ab (5th argument) must be %d", 2);
  if (NA_SHAPE1(rblapack_ab) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of ab must be the same as shape 0 of ipiv");
  ldab = NA_SHAPE0(rblapack_ab);
  if (NA_TYPE(rblapack_ab) != NA_DCOMPLEX)
    rblapack_ab = na_change_type(rblapack_ab, NA_DCOMPLEX);
  ab = NA_PTR_TYPE(rblapack_ab, doublecomplex*);
  kl = NUM2INT(rblapack_kl);
  if (!NA_IsNArray(rblapack_b))
    rb_raise(rb_eArgError, "b (11th argument) must be NArray");
  if (NA_RANK(rblapack_b) != 2)
    rb_raise(rb_eArgError, "rank of b (11th argument) must be %d", 2);
  nrhs = NA_SHAPE1(rblapack_b);
  ldb = NA_SHAPE0(rblapack_b);
  if (NA_TYPE(rblapack_b) != NA_DCOMPLEX)
    rblapack_b = na_change_type(rblapack_b, NA_DCOMPLEX);
  b = NA_PTR_TYPE(rblapack_b, doublecomplex*);
  if (!NA_IsNArray(rblapack_c))
    rb_raise(rb_eArgError, "c (10th argument) must be NArray");
  if (NA_RANK(rblapack_c) != 1)
    rb_raise(rb_eArgError, "rank of c (10th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_c) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of c must be the same as shape 0 of ipiv");
  if (NA_TYPE(rblapack_c) != NA_DFLOAT)
    rblapack_c = na_change_type(rblapack_c, NA_DFLOAT);
  c = NA_PTR_TYPE(rblapack_c, doublereal*);
  ku = NUM2INT(rblapack_ku);
  if (!NA_IsNArray(rblapack_r))
    rb_raise(rb_eArgError, "r (9th argument) must be NArray");
  if (NA_RANK(rblapack_r) != 1)
    rb_raise(rb_eArgError, "rank of r (9th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_r) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of r must be the same as shape 0 of ipiv");
  if (NA_TYPE(rblapack_r) != NA_DFLOAT)
    rblapack_r = na_change_type(rblapack_r, NA_DFLOAT);
  r = NA_PTR_TYPE(rblapack_r, doublereal*);
  if (!NA_IsNArray(rblapack_afb))
    rb_raise(rb_eArgError, "afb (6th argument) must be NArray");
  if (NA_RANK(rblapack_afb) != 2)
    rb_raise(rb_eArgError, "rank of afb (6th argument) must be %d", 2);
  if (NA_SHAPE1(rblapack_afb) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of afb must be the same as shape 0 of ipiv");
  ldafb = NA_SHAPE0(rblapack_afb);
  if (NA_TYPE(rblapack_afb) != NA_DCOMPLEX)
    rblapack_afb = na_change_type(rblapack_afb, NA_DCOMPLEX);
  afb = NA_PTR_TYPE(rblapack_afb, doublecomplex*);
  equed = StringValueCStr(rblapack_equed)[0];
  fact = StringValueCStr(rblapack_fact)[0];
  ldx = MAX(1,n);
  {
    int shape[2];
    shape[0] = ldx;
    shape[1] = nrhs;
    rblapack_x = na_make_object(NA_DCOMPLEX, 2, shape, cNArray);
  }
  x = NA_PTR_TYPE(rblapack_x, doublecomplex*);
  {
    int shape[1];
    shape[0] = nrhs;
    rblapack_ferr = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  ferr = NA_PTR_TYPE(rblapack_ferr, doublereal*);
  {
    int shape[1];
    shape[0] = nrhs;
    rblapack_berr = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  berr = NA_PTR_TYPE(rblapack_berr, doublereal*);
  {
    int shape[1];
    shape[0] = n;
    rblapack_rwork = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  rwork = NA_PTR_TYPE(rblapack_rwork, doublereal*);
  {
    int shape[2];
    shape[0] = ldab;
    shape[1] = n;
    rblapack_ab_out__ = na_make_object(NA_DCOMPLEX, 2, shape, cNArray);
  }
  ab_out__ = NA_PTR_TYPE(rblapack_ab_out__, doublecomplex*);
  MEMCPY(ab_out__, ab, doublecomplex, NA_TOTAL(rblapack_ab));
  rblapack_ab = rblapack_ab_out__;
  ab = ab_out__;
  {
    int shape[2];
    shape[0] = ldafb;
    shape[1] = n;
    rblapack_afb_out__ = na_make_object(NA_DCOMPLEX, 2, shape, cNArray);
  }
  afb_out__ = NA_PTR_TYPE(rblapack_afb_out__, doublecomplex*);
  MEMCPY(afb_out__, afb, doublecomplex, NA_TOTAL(rblapack_afb));
  rblapack_afb = rblapack_afb_out__;
  afb = afb_out__;
  {
    int shape[1];
    shape[0] = n;
    rblapack_ipiv_out__ = na_make_object(NA_LINT, 1, shape, cNArray);
  }
  ipiv_out__ = NA_PTR_TYPE(rblapack_ipiv_out__, integer*);
  MEMCPY(ipiv_out__, ipiv, integer, NA_TOTAL(rblapack_ipiv));
  rblapack_ipiv = rblapack_ipiv_out__;
  ipiv = ipiv_out__;
  {
    int shape[1];
    shape[0] = n;
    rblapack_r_out__ = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  r_out__ = NA_PTR_TYPE(rblapack_r_out__, doublereal*);
  MEMCPY(r_out__, r, doublereal, NA_TOTAL(rblapack_r));
  rblapack_r = rblapack_r_out__;
  r = r_out__;
  {
    int shape[1];
    shape[0] = n;
    rblapack_c_out__ = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  c_out__ = NA_PTR_TYPE(rblapack_c_out__, doublereal*);
  MEMCPY(c_out__, c, doublereal, NA_TOTAL(rblapack_c));
  rblapack_c = rblapack_c_out__;
  c = c_out__;
  {
    int shape[2];
    shape[0] = ldb;
    shape[1] = nrhs;
    rblapack_b_out__ = na_make_object(NA_DCOMPLEX, 2, shape, cNArray);
  }
  b_out__ = NA_PTR_TYPE(rblapack_b_out__, doublecomplex*);
  MEMCPY(b_out__, b, doublecomplex, NA_TOTAL(rblapack_b));
  rblapack_b = rblapack_b_out__;
  b = b_out__;
  work = ALLOC_N(doublecomplex, (2*n));

  zgbsvx_(&fact, &trans, &n, &kl, &ku, &nrhs, ab, &ldab, afb, &ldafb, ipiv, &equed, r, c, b, &ldb, x, &ldx, &rcond, ferr, berr, work, rwork, &info);

  free(work);
  rblapack_rcond = rb_float_new((double)rcond);
  rblapack_info = INT2NUM(info);
  rblapack_equed = rb_str_new(&equed,1);
  return rb_ary_new3(13, rblapack_x, rblapack_rcond, rblapack_ferr, rblapack_berr, rblapack_rwork, rblapack_info, rblapack_ab, rblapack_afb, rblapack_ipiv, rblapack_equed, rblapack_r, rblapack_c, rblapack_b);
}

void
init_lapack_zgbsvx(VALUE mLapack){
  rb_define_module_function(mLapack, "zgbsvx", rblapack_zgbsvx, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}

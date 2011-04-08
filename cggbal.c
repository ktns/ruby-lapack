#include "rb_lapack.h"

extern VOID cggbal_(char *job, integer *n, complex *a, integer *lda, complex *b, integer *ldb, integer *ilo, integer *ihi, real *lscale, real *rscale, real *work, integer *info);

static VALUE sHelp, sUsage;

static VALUE
rblapack_cggbal(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_job;
  char job; 
  VALUE rblapack_a;
  complex *a; 
  VALUE rblapack_b;
  complex *b; 
  VALUE rblapack_ilo;
  integer ilo; 
  VALUE rblapack_ihi;
  integer ihi; 
  VALUE rblapack_lscale;
  real *lscale; 
  VALUE rblapack_rscale;
  real *rscale; 
  VALUE rblapack_info;
  integer info; 
  VALUE rblapack_a_out__;
  complex *a_out__;
  VALUE rblapack_b_out__;
  complex *b_out__;
  real *work;

  integer lda;
  integer n;
  integer ldb;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  ilo, ihi, lscale, rscale, info, a, b = NumRu::Lapack.cggbal( job, a, b, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE CGGBAL( JOB, N, A, LDA, B, LDB, ILO, IHI, LSCALE, RSCALE, WORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  CGGBAL balances a pair of general complex matrices (A,B).  This\n*  involves, first, permuting A and B by similarity transformations to\n*  isolate eigenvalues in the first 1 to ILO$-$1 and last IHI+1 to N\n*  elements on the diagonal; and second, applying a diagonal similarity\n*  transformation to rows and columns ILO to IHI to make the rows\n*  and columns as close in norm as possible. Both steps are optional.\n*\n*  Balancing may reduce the 1-norm of the matrices, and improve the\n*  accuracy of the computed eigenvalues and/or eigenvectors in the\n*  generalized eigenvalue problem A*x = lambda*B*x.\n*\n\n*  Arguments\n*  =========\n*\n*  JOB     (input) CHARACTER*1\n*          Specifies the operations to be performed on A and B:\n*          = 'N':  none:  simply set ILO = 1, IHI = N, LSCALE(I) = 1.0\n*                  and RSCALE(I) = 1.0 for i=1,...,N;\n*          = 'P':  permute only;\n*          = 'S':  scale only;\n*          = 'B':  both permute and scale.\n*\n*  N       (input) INTEGER\n*          The order of the matrices A and B.  N >= 0.\n*\n*  A       (input/output) COMPLEX array, dimension (LDA,N)\n*          On entry, the input matrix A.\n*          On exit, A is overwritten by the balanced matrix.\n*          If JOB = 'N', A is not referenced.\n*\n*  LDA     (input) INTEGER\n*          The leading dimension of the array A. LDA >= max(1,N).\n*\n*  B       (input/output) COMPLEX array, dimension (LDB,N)\n*          On entry, the input matrix B.\n*          On exit, B is overwritten by the balanced matrix.\n*          If JOB = 'N', B is not referenced.\n*\n*  LDB     (input) INTEGER\n*          The leading dimension of the array B. LDB >= max(1,N).\n*\n*  ILO     (output) INTEGER\n*  IHI     (output) INTEGER\n*          ILO and IHI are set to integers such that on exit\n*          A(i,j) = 0 and B(i,j) = 0 if i > j and\n*          j = 1,...,ILO-1 or i = IHI+1,...,N.\n*          If JOB = 'N' or 'S', ILO = 1 and IHI = N.\n*\n*  LSCALE  (output) REAL array, dimension (N)\n*          Details of the permutations and scaling factors applied\n*          to the left side of A and B.  If P(j) is the index of the\n*          row interchanged with row j, and D(j) is the scaling factor\n*          applied to row j, then\n*            LSCALE(j) = P(j)    for J = 1,...,ILO-1\n*                      = D(j)    for J = ILO,...,IHI\n*                      = P(j)    for J = IHI+1,...,N.\n*          The order in which the interchanges are made is N to IHI+1,\n*          then 1 to ILO-1.\n*\n*  RSCALE  (output) REAL array, dimension (N)\n*          Details of the permutations and scaling factors applied\n*          to the right side of A and B.  If P(j) is the index of the\n*          column interchanged with column j, and D(j) is the scaling\n*          factor applied to column j, then\n*            RSCALE(j) = P(j)    for J = 1,...,ILO-1\n*                      = D(j)    for J = ILO,...,IHI\n*                      = P(j)    for J = IHI+1,...,N.\n*          The order in which the interchanges are made is N to IHI+1,\n*          then 1 to ILO-1.\n*\n*  WORK    (workspace) REAL array, dimension (lwork)\n*          lwork must be at least max(1,6*N) when JOB = 'S' or 'B', and\n*          at least 1 when JOB = 'N' or 'P'.\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit\n*          < 0:  if INFO = -i, the i-th argument had an illegal value.\n*\n\n*  Further Details\n*  ===============\n*\n*  See R.C. WARD, Balancing the generalized eigenvalue problem,\n*                 SIAM J. Sci. Stat. Comp. 2 (1981), 141-152.\n*\n*  =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  ilo, ihi, lscale, rscale, info, a, b = NumRu::Lapack.cggbal( job, a, b, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 3)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 3)", argc);
  rblapack_job = argv[0];
  rblapack_a = argv[1];
  rblapack_b = argv[2];
  if (rb_options != Qnil) {
  }

  if (!NA_IsNArray(rblapack_a))
    rb_raise(rb_eArgError, "a (2th argument) must be NArray");
  if (NA_RANK(rblapack_a) != 2)
    rb_raise(rb_eArgError, "rank of a (2th argument) must be %d", 2);
  n = NA_SHAPE1(rblapack_a);
  lda = NA_SHAPE0(rblapack_a);
  if (NA_TYPE(rblapack_a) != NA_SCOMPLEX)
    rblapack_a = na_change_type(rblapack_a, NA_SCOMPLEX);
  a = NA_PTR_TYPE(rblapack_a, complex*);
  if (!NA_IsNArray(rblapack_b))
    rb_raise(rb_eArgError, "b (3th argument) must be NArray");
  if (NA_RANK(rblapack_b) != 2)
    rb_raise(rb_eArgError, "rank of b (3th argument) must be %d", 2);
  if (NA_SHAPE1(rblapack_b) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of b must be the same as shape 1 of a");
  ldb = NA_SHAPE0(rblapack_b);
  if (NA_TYPE(rblapack_b) != NA_SCOMPLEX)
    rblapack_b = na_change_type(rblapack_b, NA_SCOMPLEX);
  b = NA_PTR_TYPE(rblapack_b, complex*);
  job = StringValueCStr(rblapack_job)[0];
  {
    int shape[1];
    shape[0] = n;
    rblapack_lscale = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  lscale = NA_PTR_TYPE(rblapack_lscale, real*);
  {
    int shape[1];
    shape[0] = n;
    rblapack_rscale = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  rscale = NA_PTR_TYPE(rblapack_rscale, real*);
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
    shape[1] = n;
    rblapack_b_out__ = na_make_object(NA_SCOMPLEX, 2, shape, cNArray);
  }
  b_out__ = NA_PTR_TYPE(rblapack_b_out__, complex*);
  MEMCPY(b_out__, b, complex, NA_TOTAL(rblapack_b));
  rblapack_b = rblapack_b_out__;
  b = b_out__;
  work = ALLOC_N(real, ((lsame_(&job,"S")||lsame_(&job,"B")) ? MAX(1,6*n) : (lsame_(&job,"N")||lsame_(&job,"P")) ? 1 : 0));

  cggbal_(&job, &n, a, &lda, b, &ldb, &ilo, &ihi, lscale, rscale, work, &info);

  free(work);
  rblapack_ilo = INT2NUM(ilo);
  rblapack_ihi = INT2NUM(ihi);
  rblapack_info = INT2NUM(info);
  return rb_ary_new3(7, rblapack_ilo, rblapack_ihi, rblapack_lscale, rblapack_rscale, rblapack_info, rblapack_a, rblapack_b);
}

void
init_lapack_cggbal(VALUE mLapack){
  rb_define_module_function(mLapack, "cggbal", rblapack_cggbal, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}

#include "rb_lapack.h"

extern VOID sggbal_(char *job, integer *n, real *a, integer *lda, real *b, integer *ldb, integer *ilo, integer *ihi, real *lscale, real *rscale, real *work, integer *info);

static VALUE
rb_sggbal(int argc, VALUE *argv, VALUE self){
  VALUE rb_job;
  char job; 
  VALUE rb_a;
  real *a; 
  VALUE rb_b;
  real *b; 
  VALUE rb_ilo;
  integer ilo; 
  VALUE rb_ihi;
  integer ihi; 
  VALUE rb_lscale;
  real *lscale; 
  VALUE rb_rscale;
  real *rscale; 
  VALUE rb_info;
  integer info; 
  VALUE rb_a_out__;
  real *a_out__;
  VALUE rb_b_out__;
  real *b_out__;
  real *work;

  integer lda;
  integer n;
  integer ldb;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  ilo, ihi, lscale, rscale, info, a, b = NumRu::Lapack.sggbal( job, a, b)\n    or\n  NumRu::Lapack.sggbal  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE SGGBAL( JOB, N, A, LDA, B, LDB, ILO, IHI, LSCALE, RSCALE, WORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  SGGBAL balances a pair of general real matrices (A,B).  This\n*  involves, first, permuting A and B by similarity transformations to\n*  isolate eigenvalues in the first 1 to ILO$-$1 and last IHI+1 to N\n*  elements on the diagonal; and second, applying a diagonal similarity\n*  transformation to rows and columns ILO to IHI to make the rows\n*  and columns as close in norm as possible. Both steps are optional.\n*\n*  Balancing may reduce the 1-norm of the matrices, and improve the\n*  accuracy of the computed eigenvalues and/or eigenvectors in the\n*  generalized eigenvalue problem A*x = lambda*B*x.\n*\n\n*  Arguments\n*  =========\n*\n*  JOB     (input) CHARACTER*1\n*          Specifies the operations to be performed on A and B:\n*          = 'N':  none:  simply set ILO = 1, IHI = N, LSCALE(I) = 1.0\n*                  and RSCALE(I) = 1.0 for i = 1,...,N.\n*          = 'P':  permute only;\n*          = 'S':  scale only;\n*          = 'B':  both permute and scale.\n*\n*  N       (input) INTEGER\n*          The order of the matrices A and B.  N >= 0.\n*\n*  A       (input/output) REAL array, dimension (LDA,N)\n*          On entry, the input matrix A.\n*          On exit,  A is overwritten by the balanced matrix.\n*          If JOB = 'N', A is not referenced.\n*\n*  LDA     (input) INTEGER\n*          The leading dimension of the array A. LDA >= max(1,N).\n*\n*  B       (input/output) REAL array, dimension (LDB,N)\n*          On entry, the input matrix B.\n*          On exit,  B is overwritten by the balanced matrix.\n*          If JOB = 'N', B is not referenced.\n*\n*  LDB     (input) INTEGER\n*          The leading dimension of the array B. LDB >= max(1,N).\n*\n*  ILO     (output) INTEGER\n*  IHI     (output) INTEGER\n*          ILO and IHI are set to integers such that on exit\n*          A(i,j) = 0 and B(i,j) = 0 if i > j and\n*          j = 1,...,ILO-1 or i = IHI+1,...,N.\n*          If JOB = 'N' or 'S', ILO = 1 and IHI = N.\n*\n*  LSCALE  (output) REAL array, dimension (N)\n*          Details of the permutations and scaling factors applied\n*          to the left side of A and B.  If P(j) is the index of the\n*          row interchanged with row j, and D(j)\n*          is the scaling factor applied to row j, then\n*            LSCALE(j) = P(j)    for J = 1,...,ILO-1\n*                      = D(j)    for J = ILO,...,IHI\n*                      = P(j)    for J = IHI+1,...,N.\n*          The order in which the interchanges are made is N to IHI+1,\n*          then 1 to ILO-1.\n*\n*  RSCALE  (output) REAL array, dimension (N)\n*          Details of the permutations and scaling factors applied\n*          to the right side of A and B.  If P(j) is the index of the\n*          column interchanged with column j, and D(j)\n*          is the scaling factor applied to column j, then\n*            LSCALE(j) = P(j)    for J = 1,...,ILO-1\n*                      = D(j)    for J = ILO,...,IHI\n*                      = P(j)    for J = IHI+1,...,N.\n*          The order in which the interchanges are made is N to IHI+1,\n*          then 1 to ILO-1.\n*\n*  WORK    (workspace) REAL array, dimension (lwork)\n*          lwork must be at least max(1,6*N) when JOB = 'S' or 'B', and\n*          at least 1 when JOB = 'N' or 'P'.\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit\n*          < 0:  if INFO = -i, the i-th argument had an illegal value.\n*\n\n*  Further Details\n*  ===============\n*\n*  See R.C. WARD, Balancing the generalized eigenvalue problem,\n*                 SIAM J. Sci. Stat. Comp. 2 (1981), 141-152.\n*\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 3)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 3)", argc);
  rb_job = argv[0];
  rb_a = argv[1];
  rb_b = argv[2];

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
  if (NA_SHAPE1(rb_b) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of b must be the same as shape 1 of a");
  ldb = NA_SHAPE0(rb_b);
  if (NA_TYPE(rb_b) != NA_SFLOAT)
    rb_b = na_change_type(rb_b, NA_SFLOAT);
  b = NA_PTR_TYPE(rb_b, real*);
  job = StringValueCStr(rb_job)[0];
  {
    int shape[1];
    shape[0] = n;
    rb_lscale = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  lscale = NA_PTR_TYPE(rb_lscale, real*);
  {
    int shape[1];
    shape[0] = n;
    rb_rscale = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  rscale = NA_PTR_TYPE(rb_rscale, real*);
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
  work = ALLOC_N(real, ((lsame_(&job,"S")||lsame_(&job,"B")) ? MAX(1,6*n) : (lsame_(&job,"N")||lsame_(&job,"P")) ? 1 : 0));

  sggbal_(&job, &n, a, &lda, b, &ldb, &ilo, &ihi, lscale, rscale, work, &info);

  free(work);
  rb_ilo = INT2NUM(ilo);
  rb_ihi = INT2NUM(ihi);
  rb_info = INT2NUM(info);
  return rb_ary_new3(7, rb_ilo, rb_ihi, rb_lscale, rb_rscale, rb_info, rb_a, rb_b);
}

void
init_lapack_sggbal(VALUE mLapack){
  rb_define_module_function(mLapack, "sggbal", rb_sggbal, -1);
}

#include "rb_lapack.h"

extern VOID sgebal_(char *job, integer *n, real *a, integer *lda, integer *ilo, integer *ihi, real *scale, integer *info);

static VALUE
rb_sgebal(int argc, VALUE *argv, VALUE self){
  VALUE rb_job;
  char job; 
  VALUE rb_a;
  real *a; 
  VALUE rb_ilo;
  integer ilo; 
  VALUE rb_ihi;
  integer ihi; 
  VALUE rb_scale;
  real *scale; 
  VALUE rb_info;
  integer info; 
  VALUE rb_a_out__;
  real *a_out__;

  integer lda;
  integer n;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  ilo, ihi, scale, info, a = NumRu::Lapack.sgebal( job, a)\n    or\n  NumRu::Lapack.sgebal  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE SGEBAL( JOB, N, A, LDA, ILO, IHI, SCALE, INFO )\n\n*  Purpose\n*  =======\n*\n*  SGEBAL balances a general real matrix A.  This involves, first,\n*  permuting A by a similarity transformation to isolate eigenvalues\n*  in the first 1 to ILO-1 and last IHI+1 to N elements on the\n*  diagonal; and second, applying a diagonal similarity transformation\n*  to rows and columns ILO to IHI to make the rows and columns as\n*  close in norm as possible.  Both steps are optional.\n*\n*  Balancing may reduce the 1-norm of the matrix, and improve the\n*  accuracy of the computed eigenvalues and/or eigenvectors.\n*\n\n*  Arguments\n*  =========\n*\n*  JOB     (input) CHARACTER*1\n*          Specifies the operations to be performed on A:\n*          = 'N':  none:  simply set ILO = 1, IHI = N, SCALE(I) = 1.0\n*                  for i = 1,...,N;\n*          = 'P':  permute only;\n*          = 'S':  scale only;\n*          = 'B':  both permute and scale.\n*\n*  N       (input) INTEGER\n*          The order of the matrix A.  N >= 0.\n*\n*  A       (input/output) REAL array, dimension (LDA,N)\n*          On entry, the input matrix A.\n*          On exit,  A is overwritten by the balanced matrix.\n*          If JOB = 'N', A is not referenced.\n*          See Further Details.\n*\n*  LDA     (input) INTEGER\n*          The leading dimension of the array A.  LDA >= max(1,N).\n*\n*  ILO     (output) INTEGER\n*  IHI     (output) INTEGER\n*          ILO and IHI are set to integers such that on exit\n*          A(i,j) = 0 if i > j and j = 1,...,ILO-1 or I = IHI+1,...,N.\n*          If JOB = 'N' or 'S', ILO = 1 and IHI = N.\n*\n*  SCALE   (output) REAL array, dimension (N)\n*          Details of the permutations and scaling factors applied to\n*          A.  If P(j) is the index of the row and column interchanged\n*          with row and column j and D(j) is the scaling factor\n*          applied to row and column j, then\n*          SCALE(j) = P(j)    for j = 1,...,ILO-1\n*                   = D(j)    for j = ILO,...,IHI\n*                   = P(j)    for j = IHI+1,...,N.\n*          The order in which the interchanges are made is N to IHI+1,\n*          then 1 to ILO-1.\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit.\n*          < 0:  if INFO = -i, the i-th argument had an illegal value.\n*\n\n*  Further Details\n*  ===============\n*\n*  The permutations consist of row and column interchanges which put\n*  the matrix in the form\n*\n*             ( T1   X   Y  )\n*     P A P = (  0   B   Z  )\n*             (  0   0   T2 )\n*\n*  where T1 and T2 are upper triangular matrices whose eigenvalues lie\n*  along the diagonal.  The column indices ILO and IHI mark the starting\n*  and ending columns of the submatrix B. Balancing consists of applying\n*  a diagonal similarity transformation inv(D) * B * D to make the\n*  1-norms of each row of B and its corresponding column nearly equal.\n*  The output matrix is\n*\n*     ( T1     X*D          Y    )\n*     (  0  inv(D)*B*D  inv(D)*Z ).\n*     (  0      0           T2   )\n*\n*  Information about the permutations P and the diagonal matrix D is\n*  returned in the vector SCALE.\n*\n*  This subroutine is based on the EISPACK routine BALANC.\n*\n*  Modified by Tzu-Yi Chen, Computer Science Division, University of\n*    California at Berkeley, USA\n*\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 2)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 2)", argc);
  rb_job = argv[0];
  rb_a = argv[1];

  if (!NA_IsNArray(rb_a))
    rb_raise(rb_eArgError, "a (2th argument) must be NArray");
  if (NA_RANK(rb_a) != 2)
    rb_raise(rb_eArgError, "rank of a (2th argument) must be %d", 2);
  n = NA_SHAPE1(rb_a);
  lda = NA_SHAPE0(rb_a);
  if (NA_TYPE(rb_a) != NA_SFLOAT)
    rb_a = na_change_type(rb_a, NA_SFLOAT);
  a = NA_PTR_TYPE(rb_a, real*);
  job = StringValueCStr(rb_job)[0];
  {
    int shape[1];
    shape[0] = n;
    rb_scale = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  scale = NA_PTR_TYPE(rb_scale, real*);
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

  sgebal_(&job, &n, a, &lda, &ilo, &ihi, scale, &info);

  rb_ilo = INT2NUM(ilo);
  rb_ihi = INT2NUM(ihi);
  rb_info = INT2NUM(info);
  return rb_ary_new3(5, rb_ilo, rb_ihi, rb_scale, rb_info, rb_a);
}

void
init_lapack_sgebal(VALUE mLapack){
  rb_define_module_function(mLapack, "sgebal", rb_sgebal, -1);
}

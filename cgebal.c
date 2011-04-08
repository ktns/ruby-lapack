#include "rb_lapack.h"

extern VOID cgebal_(char *job, integer *n, complex *a, integer *lda, integer *ilo, integer *ihi, real *scale, integer *info);

static VALUE sHelp, sUsage;

static VALUE
rblapack_cgebal(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_job;
  char job; 
  VALUE rblapack_a;
  complex *a; 
  VALUE rblapack_ilo;
  integer ilo; 
  VALUE rblapack_ihi;
  integer ihi; 
  VALUE rblapack_scale;
  real *scale; 
  VALUE rblapack_info;
  integer info; 
  VALUE rblapack_a_out__;
  complex *a_out__;

  integer lda;
  integer n;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  ilo, ihi, scale, info, a = NumRu::Lapack.cgebal( job, a, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE CGEBAL( JOB, N, A, LDA, ILO, IHI, SCALE, INFO )\n\n*  Purpose\n*  =======\n*\n*  CGEBAL balances a general complex matrix A.  This involves, first,\n*  permuting A by a similarity transformation to isolate eigenvalues\n*  in the first 1 to ILO-1 and last IHI+1 to N elements on the\n*  diagonal; and second, applying a diagonal similarity transformation\n*  to rows and columns ILO to IHI to make the rows and columns as\n*  close in norm as possible.  Both steps are optional.\n*\n*  Balancing may reduce the 1-norm of the matrix, and improve the\n*  accuracy of the computed eigenvalues and/or eigenvectors.\n*\n\n*  Arguments\n*  =========\n*\n*  JOB     (input) CHARACTER*1\n*          Specifies the operations to be performed on A:\n*          = 'N':  none:  simply set ILO = 1, IHI = N, SCALE(I) = 1.0\n*                  for i = 1,...,N;\n*          = 'P':  permute only;\n*          = 'S':  scale only;\n*          = 'B':  both permute and scale.\n*\n*  N       (input) INTEGER\n*          The order of the matrix A.  N >= 0.\n*\n*  A       (input/output) COMPLEX array, dimension (LDA,N)\n*          On entry, the input matrix A.\n*          On exit,  A is overwritten by the balanced matrix.\n*          If JOB = 'N', A is not referenced.\n*          See Further Details.\n*\n*  LDA     (input) INTEGER\n*          The leading dimension of the array A.  LDA >= max(1,N).\n*\n*  ILO     (output) INTEGER\n*  IHI     (output) INTEGER\n*          ILO and IHI are set to integers such that on exit\n*          A(i,j) = 0 if i > j and j = 1,...,ILO-1 or I = IHI+1,...,N.\n*          If JOB = 'N' or 'S', ILO = 1 and IHI = N.\n*\n*  SCALE   (output) REAL array, dimension (N)\n*          Details of the permutations and scaling factors applied to\n*          A.  If P(j) is the index of the row and column interchanged\n*          with row and column j and D(j) is the scaling factor\n*          applied to row and column j, then\n*          SCALE(j) = P(j)    for j = 1,...,ILO-1\n*                   = D(j)    for j = ILO,...,IHI\n*                   = P(j)    for j = IHI+1,...,N.\n*          The order in which the interchanges are made is N to IHI+1,\n*          then 1 to ILO-1.\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit.\n*          < 0:  if INFO = -i, the i-th argument had an illegal value.\n*\n\n*  Further Details\n*  ===============\n*\n*  The permutations consist of row and column interchanges which put\n*  the matrix in the form\n*\n*             ( T1   X   Y  )\n*     P A P = (  0   B   Z  )\n*             (  0   0   T2 )\n*\n*  where T1 and T2 are upper triangular matrices whose eigenvalues lie\n*  along the diagonal.  The column indices ILO and IHI mark the starting\n*  and ending columns of the submatrix B. Balancing consists of applying\n*  a diagonal similarity transformation inv(D) * B * D to make the\n*  1-norms of each row of B and its corresponding column nearly equal.\n*  The output matrix is\n*\n*     ( T1     X*D          Y    )\n*     (  0  inv(D)*B*D  inv(D)*Z ).\n*     (  0      0           T2   )\n*\n*  Information about the permutations P and the diagonal matrix D is\n*  returned in the vector SCALE.\n*\n*  This subroutine is based on the EISPACK routine CBAL.\n*\n*  Modified by Tzu-Yi Chen, Computer Science Division, University of\n*    California at Berkeley, USA\n*\n*  =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  ilo, ihi, scale, info, a = NumRu::Lapack.cgebal( job, a, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 2)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 2)", argc);
  rblapack_job = argv[0];
  rblapack_a = argv[1];
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
  job = StringValueCStr(rblapack_job)[0];
  {
    int shape[1];
    shape[0] = n;
    rblapack_scale = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  scale = NA_PTR_TYPE(rblapack_scale, real*);
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

  cgebal_(&job, &n, a, &lda, &ilo, &ihi, scale, &info);

  rblapack_ilo = INT2NUM(ilo);
  rblapack_ihi = INT2NUM(ihi);
  rblapack_info = INT2NUM(info);
  return rb_ary_new3(5, rblapack_ilo, rblapack_ihi, rblapack_scale, rblapack_info, rblapack_a);
}

void
init_lapack_cgebal(VALUE mLapack){
  rb_define_module_function(mLapack, "cgebal", rblapack_cgebal, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}

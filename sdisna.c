#include "rb_lapack.h"

static VALUE
rb_sdisna(int argc, VALUE *argv, VALUE self){
  VALUE rb_job;
  char job; 
  VALUE rb_n;
  integer n; 
  VALUE rb_d;
  real *d; 
  VALUE rb_sep;
  real *sep; 
  VALUE rb_info;
  integer info; 

  integer m;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  sep, info = NumRu::Lapack.sdisna( job, n, d)\n    or\n  NumRu::Lapack.sdisna  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE SDISNA( JOB, M, N, D, SEP, INFO )\n\n*  Purpose\n*  =======\n*\n*  SDISNA computes the reciprocal condition numbers for the eigenvectors\n*  of a real symmetric or complex Hermitian matrix or for the left or\n*  right singular vectors of a general m-by-n matrix. The reciprocal\n*  condition number is the 'gap' between the corresponding eigenvalue or\n*  singular value and the nearest other one.\n*\n*  The bound on the error, measured by angle in radians, in the I-th\n*  computed vector is given by\n*\n*         SLAMCH( 'E' ) * ( ANORM / SEP( I ) )\n*\n*  where ANORM = 2-norm(A) = max( abs( D(j) ) ).  SEP(I) is not allowed\n*  to be smaller than SLAMCH( 'E' )*ANORM in order to limit the size of\n*  the error bound.\n*\n*  SDISNA may also be used to compute error bounds for eigenvectors of\n*  the generalized symmetric definite eigenproblem.\n*\n\n*  Arguments\n*  =========\n*\n*  JOB     (input) CHARACTER*1\n*          Specifies for which problem the reciprocal condition numbers\n*          should be computed:\n*          = 'E':  the eigenvectors of a symmetric/Hermitian matrix;\n*          = 'L':  the left singular vectors of a general matrix;\n*          = 'R':  the right singular vectors of a general matrix.\n*\n*  M       (input) INTEGER\n*          The number of rows of the matrix. M >= 0.\n*\n*  N       (input) INTEGER\n*          If JOB = 'L' or 'R', the number of columns of the matrix,\n*          in which case N >= 0. Ignored if JOB = 'E'.\n*\n*  D       (input) REAL array, dimension (M) if JOB = 'E'\n*                              dimension (min(M,N)) if JOB = 'L' or 'R'\n*          The eigenvalues (if JOB = 'E') or singular values (if JOB =\n*          'L' or 'R') of the matrix, in either increasing or decreasing\n*          order. If singular values, they must be non-negative.\n*\n*  SEP     (output) REAL array, dimension (M) if JOB = 'E'\n*                               dimension (min(M,N)) if JOB = 'L' or 'R'\n*          The reciprocal condition numbers of the vectors.\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit.\n*          < 0:  if INFO = -i, the i-th argument had an illegal value.\n*\n\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 3)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 3)", argc);
  rb_job = argv[0];
  rb_n = argv[1];
  rb_d = argv[2];

  job = StringValueCStr(rb_job)[0];
  n = NUM2INT(rb_n);
  if (!NA_IsNArray(rb_d))
    rb_raise(rb_eArgError, "d (3th argument) must be NArray");
  if (NA_RANK(rb_d) != 1)
    rb_raise(rb_eArgError, "rank of d (3th argument) must be %d", 1);
  m = NA_SHAPE0(rb_d);
  if (NA_TYPE(rb_d) != NA_SFLOAT)
    rb_d = na_change_type(rb_d, NA_SFLOAT);
  d = NA_PTR_TYPE(rb_d, real*);
  {
    int shape[1];
    shape[0] = lsame_(&job,"E") ? m : ((lsame_(&job,"L")) || (lsame_(&job,"R"))) ? MIN(m,n) : 0;
    rb_sep = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  sep = NA_PTR_TYPE(rb_sep, real*);

  sdisna_(&job, &m, &n, d, sep, &info);

  rb_info = INT2NUM(info);
  return rb_ary_new3(2, rb_sep, rb_info);
}

void
init_lapack_sdisna(VALUE mLapack){
  rb_define_module_function(mLapack, "sdisna", rb_sdisna, -1);
}

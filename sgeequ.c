#include "rb_lapack.h"

extern VOID sgeequ_(integer *m, integer *n, real *a, integer *lda, real *r, real *c, real *rowcnd, real *colcnd, real *amax, integer *info);

static VALUE
rb_sgeequ(int argc, VALUE *argv, VALUE self){
  VALUE rb_a;
  real *a; 
  VALUE rb_r;
  real *r; 
  VALUE rb_c;
  real *c; 
  VALUE rb_rowcnd;
  real rowcnd; 
  VALUE rb_colcnd;
  real colcnd; 
  VALUE rb_amax;
  real amax; 
  VALUE rb_info;
  integer info; 

  integer lda;
  integer n;
  integer m;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  r, c, rowcnd, colcnd, amax, info = NumRu::Lapack.sgeequ( a)\n    or\n  NumRu::Lapack.sgeequ  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE SGEEQU( M, N, A, LDA, R, C, ROWCND, COLCND, AMAX, INFO )\n\n*  Purpose\n*  =======\n*\n*  SGEEQU computes row and column scalings intended to equilibrate an\n*  M-by-N matrix A and reduce its condition number.  R returns the row\n*  scale factors and C the column scale factors, chosen to try to make\n*  the largest element in each row and column of the matrix B with\n*  elements B(i,j)=R(i)*A(i,j)*C(j) have absolute value 1.\n*\n*  R(i) and C(j) are restricted to be between SMLNUM = smallest safe\n*  number and BIGNUM = largest safe number.  Use of these scaling\n*  factors is not guaranteed to reduce the condition number of A but\n*  works well in practice.\n*\n\n*  Arguments\n*  =========\n*\n*  M       (input) INTEGER\n*          The number of rows of the matrix A.  M >= 0.\n*\n*  N       (input) INTEGER\n*          The number of columns of the matrix A.  N >= 0.\n*\n*  A       (input) REAL array, dimension (LDA,N)\n*          The M-by-N matrix whose equilibration factors are\n*          to be computed.\n*\n*  LDA     (input) INTEGER\n*          The leading dimension of the array A.  LDA >= max(1,M).\n*\n*  R       (output) REAL array, dimension (M)\n*          If INFO = 0 or INFO > M, R contains the row scale factors\n*          for A.\n*\n*  C       (output) REAL array, dimension (N)\n*          If INFO = 0,  C contains the column scale factors for A.\n*\n*  ROWCND  (output) REAL\n*          If INFO = 0 or INFO > M, ROWCND contains the ratio of the\n*          smallest R(i) to the largest R(i).  If ROWCND >= 0.1 and\n*          AMAX is neither too large nor too small, it is not worth\n*          scaling by R.\n*\n*  COLCND  (output) REAL\n*          If INFO = 0, COLCND contains the ratio of the smallest\n*          C(i) to the largest C(i).  If COLCND >= 0.1, it is not\n*          worth scaling by C.\n*\n*  AMAX    (output) REAL\n*          Absolute value of largest matrix element.  If AMAX is very\n*          close to overflow or very close to underflow, the matrix\n*          should be scaled.\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit\n*          < 0:  if INFO = -i, the i-th argument had an illegal value\n*          > 0:  if INFO = i,  and i is\n*                <= M:  the i-th row of A is exactly zero\n*                >  M:  the (i-M)-th column of A is exactly zero\n*\n\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 1)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 1)", argc);
  rb_a = argv[0];

  if (!NA_IsNArray(rb_a))
    rb_raise(rb_eArgError, "a (1th argument) must be NArray");
  if (NA_RANK(rb_a) != 2)
    rb_raise(rb_eArgError, "rank of a (1th argument) must be %d", 2);
  n = NA_SHAPE1(rb_a);
  lda = NA_SHAPE0(rb_a);
  if (NA_TYPE(rb_a) != NA_SFLOAT)
    rb_a = na_change_type(rb_a, NA_SFLOAT);
  a = NA_PTR_TYPE(rb_a, real*);
  m = lda;
  {
    int shape[1];
    shape[0] = m;
    rb_r = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  r = NA_PTR_TYPE(rb_r, real*);
  {
    int shape[1];
    shape[0] = n;
    rb_c = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  c = NA_PTR_TYPE(rb_c, real*);

  sgeequ_(&m, &n, a, &lda, r, c, &rowcnd, &colcnd, &amax, &info);

  rb_rowcnd = rb_float_new((double)rowcnd);
  rb_colcnd = rb_float_new((double)colcnd);
  rb_amax = rb_float_new((double)amax);
  rb_info = INT2NUM(info);
  return rb_ary_new3(6, rb_r, rb_c, rb_rowcnd, rb_colcnd, rb_amax, rb_info);
}

void
init_lapack_sgeequ(VALUE mLapack){
  rb_define_module_function(mLapack, "sgeequ", rb_sgeequ, -1);
}

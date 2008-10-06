#include "rb_lapack.h"

static VALUE
rb_clarcm(int argc, VALUE *argv, VALUE self){
  VALUE rb_a;
  real *a; 
  VALUE rb_b;
  complex *b; 
  VALUE rb_c;
  complex *c; 
  real *rwork;

  integer lda;
  integer m;
  integer ldb;
  integer n;
  integer ldc;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  c = NumRu::Lapack.clarcm( a, b)\n    or\n  NumRu::Lapack.clarcm  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE CLARCM( M, N, A, LDA, B, LDB, C, LDC, RWORK )\n\n*  Purpose\n*  =======\n*\n*  CLARCM performs a very simple matrix-matrix multiplication:\n*           C := A * B,\n*  where A is M by M and real; B is M by N and complex;\n*  C is M by N and complex.\n*\n\n*  Arguments\n*  =========\n*\n*  M       (input) INTEGER\n*          The number of rows of the matrix A and of the matrix C.\n*          M >= 0.\n*\n*  N       (input) INTEGER\n*          The number of columns and rows of the matrix B and\n*          the number of columns of the matrix C.\n*          N >= 0.\n*\n*  A       (input) REAL array, dimension (LDA, M)\n*          A contains the M by M matrix A.\n*\n*  LDA     (input) INTEGER\n*          The leading dimension of the array A. LDA >=max(1,M).\n*\n*  B       (input) REAL array, dimension (LDB, N)\n*          B contains the M by N matrix B.\n*\n*  LDB     (input) INTEGER\n*          The leading dimension of the array B. LDB >=max(1,M).\n*\n*  C       (input) COMPLEX array, dimension (LDC, N)\n*          C contains the M by N matrix C.\n*\n*  LDC     (input) INTEGER\n*          The leading dimension of the array C. LDC >=max(1,M).\n*\n*  RWORK   (workspace) REAL array, dimension (2*M*N)\n*\n\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 2)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 2)", argc);
  rb_a = argv[0];
  rb_b = argv[1];

  if (!NA_IsNArray(rb_a))
    rb_raise(rb_eArgError, "a (1th argument) must be NArray");
  if (NA_RANK(rb_a) != 2)
    rb_raise(rb_eArgError, "rank of a (1th argument) must be %d", 2);
  lda = NA_SHAPE0(rb_a);
  m = NA_SHAPE1(rb_a);
  if (NA_TYPE(rb_a) != NA_SFLOAT)
    rb_a = na_change_type(rb_a, NA_SFLOAT);
  a = NA_PTR_TYPE(rb_a, real*);
  if (!NA_IsNArray(rb_b))
    rb_raise(rb_eArgError, "b (2th argument) must be NArray");
  if (NA_RANK(rb_b) != 2)
    rb_raise(rb_eArgError, "rank of b (2th argument) must be %d", 2);
  ldb = NA_SHAPE0(rb_b);
  n = NA_SHAPE1(rb_b);
  if (NA_TYPE(rb_b) != NA_SCOMPLEX)
    rb_b = na_change_type(rb_b, NA_SCOMPLEX);
  b = NA_PTR_TYPE(rb_b, complex*);
  ldc = MAX(1,m);
  {
    int shape[2];
    shape[0] = ldc;
    shape[1] = n;
    rb_c = na_make_object(NA_SCOMPLEX, 2, shape, cNArray);
  }
  c = NA_PTR_TYPE(rb_c, complex*);
  rwork = ALLOC_N(real, (2*m*n));

  clarcm_(&m, &n, a, &lda, b, &ldb, c, &ldc, rwork);

  free(rwork);
  return rb_c;
}

void
init_lapack_clarcm(VALUE mLapack){
  rb_define_module_function(mLapack, "clarcm", rb_clarcm, -1);
}

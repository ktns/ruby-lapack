#include "rb_lapack.h"

extern VOID zlarcm_(integer *m, integer *n, doublereal *a, integer *lda, doublecomplex *b, integer *ldb, doublecomplex *c, integer *ldc, doublereal *rwork);

static VALUE
rb_zlarcm(int argc, VALUE *argv, VALUE self){
  VALUE rb_a;
  doublereal *a; 
  VALUE rb_b;
  doublecomplex *b; 
  VALUE rb_c;
  doublecomplex *c; 
  doublereal *rwork;

  integer lda;
  integer m;
  integer ldb;
  integer n;
  integer ldc;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  c = NumRu::Lapack.zlarcm( a, b)\n    or\n  NumRu::Lapack.zlarcm  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE ZLARCM( M, N, A, LDA, B, LDB, C, LDC, RWORK )\n\n*  Purpose\n*  =======\n*\n*  ZLARCM performs a very simple matrix-matrix multiplication:\n*           C := A * B,\n*  where A is M by M and real; B is M by N and complex;\n*  C is M by N and complex.\n*\n\n*  Arguments\n*  =========\n*\n*  M       (input) INTEGER\n*          The number of rows of the matrix A and of the matrix C.\n*          M >= 0.\n*\n*  N       (input) INTEGER\n*          The number of columns and rows of the matrix B and\n*          the number of columns of the matrix C.\n*          N >= 0.\n*\n*  A       (input) DOUBLE PRECISION array, dimension (LDA, M)\n*          A contains the M by M matrix A.\n*\n*  LDA     (input) INTEGER\n*          The leading dimension of the array A. LDA >=max(1,M).\n*\n*  B       (input) DOUBLE PRECISION array, dimension (LDB, N)\n*          B contains the M by N matrix B.\n*\n*  LDB     (input) INTEGER\n*          The leading dimension of the array B. LDB >=max(1,M).\n*\n*  C       (input) COMPLEX*16 array, dimension (LDC, N)\n*          C contains the M by N matrix C.\n*\n*  LDC     (input) INTEGER\n*          The leading dimension of the array C. LDC >=max(1,M).\n*\n*  RWORK   (workspace) DOUBLE PRECISION array, dimension (2*M*N)\n*\n\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 2)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 2)", argc);
  rb_a = argv[0];
  rb_b = argv[1];

  if (!NA_IsNArray(rb_b))
    rb_raise(rb_eArgError, "b (2th argument) must be NArray");
  if (NA_RANK(rb_b) != 2)
    rb_raise(rb_eArgError, "rank of b (2th argument) must be %d", 2);
  n = NA_SHAPE1(rb_b);
  ldb = NA_SHAPE0(rb_b);
  if (NA_TYPE(rb_b) != NA_DCOMPLEX)
    rb_b = na_change_type(rb_b, NA_DCOMPLEX);
  b = NA_PTR_TYPE(rb_b, doublecomplex*);
  if (!NA_IsNArray(rb_a))
    rb_raise(rb_eArgError, "a (1th argument) must be NArray");
  if (NA_RANK(rb_a) != 2)
    rb_raise(rb_eArgError, "rank of a (1th argument) must be %d", 2);
  m = NA_SHAPE1(rb_a);
  lda = NA_SHAPE0(rb_a);
  if (NA_TYPE(rb_a) != NA_DFLOAT)
    rb_a = na_change_type(rb_a, NA_DFLOAT);
  a = NA_PTR_TYPE(rb_a, doublereal*);
  ldc = MAX(1,m);
  {
    int shape[2];
    shape[0] = ldc;
    shape[1] = n;
    rb_c = na_make_object(NA_DCOMPLEX, 2, shape, cNArray);
  }
  c = NA_PTR_TYPE(rb_c, doublecomplex*);
  rwork = ALLOC_N(doublereal, (2*m*n));

  zlarcm_(&m, &n, a, &lda, b, &ldb, c, &ldc, rwork);

  free(rwork);
  return rb_c;
}

void
init_lapack_zlarcm(VALUE mLapack){
  rb_define_module_function(mLapack, "zlarcm", rb_zlarcm, -1);
}

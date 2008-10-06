#include "rb_lapack.h"

static VALUE
rb_zlacrm(int argc, VALUE *argv, VALUE self){
  VALUE rb_m;
  integer m; 
  VALUE rb_a;
  doublecomplex *a; 
  VALUE rb_b;
  doublereal *b; 
  VALUE rb_c;
  doublecomplex *c; 
  doublereal *rwork;

  integer lda;
  integer n;
  integer ldb;
  integer ldc;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  c = NumRu::Lapack.zlacrm( m, a, b)\n    or\n  NumRu::Lapack.zlacrm  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE ZLACRM( M, N, A, LDA, B, LDB, C, LDC, RWORK )\n\n*  Purpose\n*  =======\n*\n*  ZLACRM performs a very simple matrix-matrix multiplication:\n*           C := A * B,\n*  where A is M by N and complex; B is N by N and real;\n*  C is M by N and complex.\n*\n\n*  Arguments\n*  =========\n*\n*  M       (input) INTEGER\n*          The number of rows of the matrix A and of the matrix C.\n*          M >= 0.\n*\n*  N       (input) INTEGER\n*          The number of columns and rows of the matrix B and\n*          the number of columns of the matrix C.\n*          N >= 0.\n*\n*  A       (input) COMPLEX*16 array, dimension (LDA, N)\n*          A contains the M by N matrix A.\n*\n*  LDA     (input) INTEGER\n*          The leading dimension of the array A. LDA >=max(1,M).\n*\n*  B       (input) DOUBLE PRECISION array, dimension (LDB, N)\n*          B contains the N by N matrix B.\n*\n*  LDB     (input) INTEGER\n*          The leading dimension of the array B. LDB >=max(1,N).\n*\n*  C       (input) COMPLEX*16 array, dimension (LDC, N)\n*          C contains the M by N matrix C.\n*\n*  LDC     (input) INTEGER\n*          The leading dimension of the array C. LDC >=max(1,N).\n*\n*  RWORK   (workspace) DOUBLE PRECISION array, dimension (2*M*N)\n*\n\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 3)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 3)", argc);
  rb_m = argv[0];
  rb_a = argv[1];
  rb_b = argv[2];

  m = NUM2INT(rb_m);
  if (!NA_IsNArray(rb_a))
    rb_raise(rb_eArgError, "a (2th argument) must be NArray");
  if (NA_RANK(rb_a) != 2)
    rb_raise(rb_eArgError, "rank of a (2th argument) must be %d", 2);
  lda = NA_SHAPE0(rb_a);
  n = NA_SHAPE1(rb_a);
  if (NA_TYPE(rb_a) != NA_DCOMPLEX)
    rb_a = na_change_type(rb_a, NA_DCOMPLEX);
  a = NA_PTR_TYPE(rb_a, doublecomplex*);
  if (!NA_IsNArray(rb_b))
    rb_raise(rb_eArgError, "b (3th argument) must be NArray");
  if (NA_RANK(rb_b) != 2)
    rb_raise(rb_eArgError, "rank of b (3th argument) must be %d", 2);
  ldb = NA_SHAPE0(rb_b);
  if (NA_SHAPE1(rb_b) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of b must be the same as shape 1 of a");
  if (NA_TYPE(rb_b) != NA_DFLOAT)
    rb_b = na_change_type(rb_b, NA_DFLOAT);
  b = NA_PTR_TYPE(rb_b, doublereal*);
  ldc = MAX(1,n);
  {
    int shape[2];
    shape[0] = ldc;
    shape[1] = n;
    rb_c = na_make_object(NA_DCOMPLEX, 2, shape, cNArray);
  }
  c = NA_PTR_TYPE(rb_c, doublecomplex*);
  rwork = ALLOC_N(doublereal, (2*m*n));

  zlacrm_(&m, &n, a, &lda, b, &ldb, c, &ldc, rwork);

  free(rwork);
  return rb_c;
}

void
init_lapack_zlacrm(VALUE mLapack){
  rb_define_module_function(mLapack, "zlacrm", rb_zlacrm, -1);
}

#include "rb_lapack.h"

extern VOID zlarcm_(integer *m, integer *n, doublereal *a, integer *lda, doublecomplex *b, integer *ldb, doublecomplex *c, integer *ldc, doublereal *rwork);

static VALUE sHelp, sUsage;

static VALUE
rblapack_zlarcm(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_a;
  doublereal *a; 
  VALUE rblapack_b;
  doublecomplex *b; 
  VALUE rblapack_c;
  doublecomplex *c; 
  doublereal *rwork;

  integer lda;
  integer m;
  integer ldb;
  integer n;
  integer ldc;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  c = NumRu::Lapack.zlarcm( a, b, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE ZLARCM( M, N, A, LDA, B, LDB, C, LDC, RWORK )\n\n*  Purpose\n*  =======\n*\n*  ZLARCM performs a very simple matrix-matrix multiplication:\n*           C := A * B,\n*  where A is M by M and real; B is M by N and complex;\n*  C is M by N and complex.\n*\n\n*  Arguments\n*  =========\n*\n*  M       (input) INTEGER\n*          The number of rows of the matrix A and of the matrix C.\n*          M >= 0.\n*\n*  N       (input) INTEGER\n*          The number of columns and rows of the matrix B and\n*          the number of columns of the matrix C.\n*          N >= 0.\n*\n*  A       (input) DOUBLE PRECISION array, dimension (LDA, M)\n*          A contains the M by M matrix A.\n*\n*  LDA     (input) INTEGER\n*          The leading dimension of the array A. LDA >=max(1,M).\n*\n*  B       (input) DOUBLE PRECISION array, dimension (LDB, N)\n*          B contains the M by N matrix B.\n*\n*  LDB     (input) INTEGER\n*          The leading dimension of the array B. LDB >=max(1,M).\n*\n*  C       (input) COMPLEX*16 array, dimension (LDC, N)\n*          C contains the M by N matrix C.\n*\n*  LDC     (input) INTEGER\n*          The leading dimension of the array C. LDC >=max(1,M).\n*\n*  RWORK   (workspace) DOUBLE PRECISION array, dimension (2*M*N)\n*\n\n*  =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  c = NumRu::Lapack.zlarcm( a, b, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 2)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 2)", argc);
  rblapack_a = argv[0];
  rblapack_b = argv[1];
  if (rb_options != Qnil) {
  }

  if (!NA_IsNArray(rblapack_b))
    rb_raise(rb_eArgError, "b (2th argument) must be NArray");
  if (NA_RANK(rblapack_b) != 2)
    rb_raise(rb_eArgError, "rank of b (2th argument) must be %d", 2);
  n = NA_SHAPE1(rblapack_b);
  ldb = NA_SHAPE0(rblapack_b);
  if (NA_TYPE(rblapack_b) != NA_DCOMPLEX)
    rblapack_b = na_change_type(rblapack_b, NA_DCOMPLEX);
  b = NA_PTR_TYPE(rblapack_b, doublecomplex*);
  if (!NA_IsNArray(rblapack_a))
    rb_raise(rb_eArgError, "a (1th argument) must be NArray");
  if (NA_RANK(rblapack_a) != 2)
    rb_raise(rb_eArgError, "rank of a (1th argument) must be %d", 2);
  m = NA_SHAPE1(rblapack_a);
  lda = NA_SHAPE0(rblapack_a);
  if (NA_TYPE(rblapack_a) != NA_DFLOAT)
    rblapack_a = na_change_type(rblapack_a, NA_DFLOAT);
  a = NA_PTR_TYPE(rblapack_a, doublereal*);
  ldc = MAX(1,m);
  {
    int shape[2];
    shape[0] = ldc;
    shape[1] = n;
    rblapack_c = na_make_object(NA_DCOMPLEX, 2, shape, cNArray);
  }
  c = NA_PTR_TYPE(rblapack_c, doublecomplex*);
  rwork = ALLOC_N(doublereal, (2*m*n));

  zlarcm_(&m, &n, a, &lda, b, &ldb, c, &ldc, rwork);

  free(rwork);
  return rblapack_c;
}

void
init_lapack_zlarcm(VALUE mLapack){
  rb_define_module_function(mLapack, "zlarcm", rblapack_zlarcm, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}

#include "rb_lapack.h"

static VALUE
rb_clacpy(int argc, VALUE *argv, VALUE self){
  VALUE rb_uplo;
  char uplo; 
  VALUE rb_m;
  integer m; 
  VALUE rb_a;
  complex *a; 
  VALUE rb_b;
  complex *b; 

  integer lda;
  integer n;
  integer ldb;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  b = NumRu::Lapack.clacpy( uplo, m, a)\n    or\n  NumRu::Lapack.clacpy  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE CLACPY( UPLO, M, N, A, LDA, B, LDB )\n\n*  Purpose\n*  =======\n*\n*  CLACPY copies all or part of a two-dimensional matrix A to another\n*  matrix B.\n*\n\n*  Arguments\n*  =========\n*\n*  UPLO    (input) CHARACTER*1\n*          Specifies the part of the matrix A to be copied to B.\n*          = 'U':      Upper triangular part\n*          = 'L':      Lower triangular part\n*          Otherwise:  All of the matrix A\n*\n*  M       (input) INTEGER\n*          The number of rows of the matrix A.  M >= 0.\n*\n*  N       (input) INTEGER\n*          The number of columns of the matrix A.  N >= 0.\n*\n*  A       (input) COMPLEX array, dimension (LDA,N)\n*          The m by n matrix A.  If UPLO = 'U', only the upper trapezium\n*          is accessed; if UPLO = 'L', only the lower trapezium is\n*          accessed.\n*\n*  LDA     (input) INTEGER\n*          The leading dimension of the array A.  LDA >= max(1,M).\n*\n*  B       (output) COMPLEX array, dimension (LDB,N)\n*          On exit, B = A in the locations specified by UPLO.\n*\n*  LDB     (input) INTEGER\n*          The leading dimension of the array B.  LDB >= max(1,M).\n*\n\n*  =====================================================================\n*\n*     .. Local Scalars ..\n      INTEGER            I, J\n*     ..\n*     .. External Functions ..\n      LOGICAL            LSAME\n      EXTERNAL           LSAME\n*     ..\n*     .. Intrinsic Functions ..\n      INTRINSIC          MIN\n*     ..\n\n");
    return Qnil;
  }
  if (argc != 3)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 3)", argc);
  rb_uplo = argv[0];
  rb_m = argv[1];
  rb_a = argv[2];

  uplo = StringValueCStr(rb_uplo)[0];
  m = NUM2INT(rb_m);
  if (!NA_IsNArray(rb_a))
    rb_raise(rb_eArgError, "a (3th argument) must be NArray");
  if (NA_RANK(rb_a) != 2)
    rb_raise(rb_eArgError, "rank of a (3th argument) must be %d", 2);
  lda = NA_SHAPE0(rb_a);
  n = NA_SHAPE1(rb_a);
  if (NA_TYPE(rb_a) != NA_SCOMPLEX)
    rb_a = na_change_type(rb_a, NA_SCOMPLEX);
  a = NA_PTR_TYPE(rb_a, complex*);
  ldb = MAX(1,m);
  {
    int shape[2];
    shape[0] = ldb;
    shape[1] = n;
    rb_b = na_make_object(NA_SCOMPLEX, 2, shape, cNArray);
  }
  b = NA_PTR_TYPE(rb_b, complex*);

  clacpy_(&uplo, &m, &n, a, &lda, b, &ldb);

  return rb_b;
}

void
init_lapack_clacpy(VALUE mLapack){
  rb_define_module_function(mLapack, "clacpy", rb_clacpy, -1);
}

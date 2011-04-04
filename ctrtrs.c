#include "rb_lapack.h"

extern VOID ctrtrs_(char *uplo, char *trans, char *diag, integer *n, integer *nrhs, complex *a, integer *lda, complex *b, integer *ldb, integer *info);

static VALUE
rb_ctrtrs(int argc, VALUE *argv, VALUE self){
  VALUE rb_uplo;
  char uplo; 
  VALUE rb_trans;
  char trans; 
  VALUE rb_diag;
  char diag; 
  VALUE rb_a;
  complex *a; 
  VALUE rb_b;
  complex *b; 
  VALUE rb_info;
  integer info; 
  VALUE rb_b_out__;
  complex *b_out__;

  integer lda;
  integer n;
  integer ldb;
  integer nrhs;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  info, b = NumRu::Lapack.ctrtrs( uplo, trans, diag, a, b)\n    or\n  NumRu::Lapack.ctrtrs  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE CTRTRS( UPLO, TRANS, DIAG, N, NRHS, A, LDA, B, LDB, INFO )\n\n*  Purpose\n*  =======\n*\n*  CTRTRS solves a triangular system of the form\n*\n*     A * X = B,  A**T * X = B,  or  A**H * X = B,\n*\n*  where A is a triangular matrix of order N, and B is an N-by-NRHS\n*  matrix.  A check is made to verify that A is nonsingular.\n*\n\n*  Arguments\n*  =========\n*\n*  UPLO    (input) CHARACTER*1\n*          = 'U':  A is upper triangular;\n*          = 'L':  A is lower triangular.\n*\n*  TRANS   (input) CHARACTER*1\n*          Specifies the form of the system of equations:\n*          = 'N':  A * X = B     (No transpose)\n*          = 'T':  A**T * X = B  (Transpose)\n*          = 'C':  A**H * X = B  (Conjugate transpose)\n*\n*  DIAG    (input) CHARACTER*1\n*          = 'N':  A is non-unit triangular;\n*          = 'U':  A is unit triangular.\n*\n*  N       (input) INTEGER\n*          The order of the matrix A.  N >= 0.\n*\n*  NRHS    (input) INTEGER\n*          The number of right hand sides, i.e., the number of columns\n*          of the matrix B.  NRHS >= 0.\n*\n*  A       (input) COMPLEX array, dimension (LDA,N)\n*          The triangular matrix A.  If UPLO = 'U', the leading N-by-N\n*          upper triangular part of the array A contains the upper\n*          triangular matrix, and the strictly lower triangular part of\n*          A is not referenced.  If UPLO = 'L', the leading N-by-N lower\n*          triangular part of the array A contains the lower triangular\n*          matrix, and the strictly upper triangular part of A is not\n*          referenced.  If DIAG = 'U', the diagonal elements of A are\n*          also not referenced and are assumed to be 1.\n*\n*  LDA     (input) INTEGER\n*          The leading dimension of the array A.  LDA >= max(1,N).\n*\n*  B       (input/output) COMPLEX array, dimension (LDB,NRHS)\n*          On entry, the right hand side matrix B.\n*          On exit, if INFO = 0, the solution matrix X.\n*\n*  LDB     (input) INTEGER\n*          The leading dimension of the array B.  LDB >= max(1,N).\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit\n*          < 0: if INFO = -i, the i-th argument had an illegal value\n*          > 0: if INFO = i, the i-th diagonal element of A is zero,\n*               indicating that the matrix is singular and the solutions\n*               X have not been computed.\n*\n\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 5)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 5)", argc);
  rb_uplo = argv[0];
  rb_trans = argv[1];
  rb_diag = argv[2];
  rb_a = argv[3];
  rb_b = argv[4];

  if (!NA_IsNArray(rb_a))
    rb_raise(rb_eArgError, "a (4th argument) must be NArray");
  if (NA_RANK(rb_a) != 2)
    rb_raise(rb_eArgError, "rank of a (4th argument) must be %d", 2);
  n = NA_SHAPE1(rb_a);
  lda = NA_SHAPE0(rb_a);
  if (NA_TYPE(rb_a) != NA_SCOMPLEX)
    rb_a = na_change_type(rb_a, NA_SCOMPLEX);
  a = NA_PTR_TYPE(rb_a, complex*);
  if (!NA_IsNArray(rb_b))
    rb_raise(rb_eArgError, "b (5th argument) must be NArray");
  if (NA_RANK(rb_b) != 2)
    rb_raise(rb_eArgError, "rank of b (5th argument) must be %d", 2);
  nrhs = NA_SHAPE1(rb_b);
  ldb = NA_SHAPE0(rb_b);
  if (NA_TYPE(rb_b) != NA_SCOMPLEX)
    rb_b = na_change_type(rb_b, NA_SCOMPLEX);
  b = NA_PTR_TYPE(rb_b, complex*);
  diag = StringValueCStr(rb_diag)[0];
  trans = StringValueCStr(rb_trans)[0];
  uplo = StringValueCStr(rb_uplo)[0];
  {
    int shape[2];
    shape[0] = ldb;
    shape[1] = nrhs;
    rb_b_out__ = na_make_object(NA_SCOMPLEX, 2, shape, cNArray);
  }
  b_out__ = NA_PTR_TYPE(rb_b_out__, complex*);
  MEMCPY(b_out__, b, complex, NA_TOTAL(rb_b));
  rb_b = rb_b_out__;
  b = b_out__;

  ctrtrs_(&uplo, &trans, &diag, &n, &nrhs, a, &lda, b, &ldb, &info);

  rb_info = INT2NUM(info);
  return rb_ary_new3(2, rb_info, rb_b);
}

void
init_lapack_ctrtrs(VALUE mLapack){
  rb_define_module_function(mLapack, "ctrtrs", rb_ctrtrs, -1);
}

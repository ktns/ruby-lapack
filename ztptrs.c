#include "rb_lapack.h"

extern VOID ztptrs_(char *uplo, char *trans, char *diag, integer *n, integer *nrhs, doublecomplex *ap, doublecomplex *b, integer *ldb, integer *info);

static VALUE
rb_ztptrs(int argc, VALUE *argv, VALUE self){
  VALUE rb_uplo;
  char uplo; 
  VALUE rb_trans;
  char trans; 
  VALUE rb_diag;
  char diag; 
  VALUE rb_n;
  integer n; 
  VALUE rb_ap;
  doublecomplex *ap; 
  VALUE rb_b;
  doublecomplex *b; 
  VALUE rb_info;
  integer info; 
  VALUE rb_b_out__;
  doublecomplex *b_out__;

  integer ldb;
  integer nrhs;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  info, b = NumRu::Lapack.ztptrs( uplo, trans, diag, n, ap, b)\n    or\n  NumRu::Lapack.ztptrs  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE ZTPTRS( UPLO, TRANS, DIAG, N, NRHS, AP, B, LDB, INFO )\n\n*  Purpose\n*  =======\n*\n*  ZTPTRS solves a triangular system of the form\n*\n*     A * X = B,  A**T * X = B,  or  A**H * X = B,\n*\n*  where A is a triangular matrix of order N stored in packed format,\n*  and B is an N-by-NRHS matrix.  A check is made to verify that A is\n*  nonsingular.\n*\n\n*  Arguments\n*  =========\n*\n*  UPLO    (input) CHARACTER*1\n*          = 'U':  A is upper triangular;\n*          = 'L':  A is lower triangular.\n*\n*  TRANS   (input) CHARACTER*1\n*          Specifies the form of the system of equations:\n*          = 'N':  A * X = B     (No transpose)\n*          = 'T':  A**T * X = B  (Transpose)\n*          = 'C':  A**H * X = B  (Conjugate transpose)\n*\n*  DIAG    (input) CHARACTER*1\n*          = 'N':  A is non-unit triangular;\n*          = 'U':  A is unit triangular.\n*\n*  N       (input) INTEGER\n*          The order of the matrix A.  N >= 0.\n*\n*  NRHS    (input) INTEGER\n*          The number of right hand sides, i.e., the number of columns\n*          of the matrix B.  NRHS >= 0.\n*\n*  AP      (input) COMPLEX*16 array, dimension (N*(N+1)/2)\n*          The upper or lower triangular matrix A, packed columnwise in\n*          a linear array.  The j-th column of A is stored in the array\n*          AP as follows:\n*          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;\n*          if UPLO = 'L', AP(i + (j-1)*(2*n-j)/2) = A(i,j) for j<=i<=n.\n*\n*  B       (input/output) COMPLEX*16 array, dimension (LDB,NRHS)\n*          On entry, the right hand side matrix B.\n*          On exit, if INFO = 0, the solution matrix X.\n*\n*  LDB     (input) INTEGER\n*          The leading dimension of the array B.  LDB >= max(1,N).\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit\n*          < 0:  if INFO = -i, the i-th argument had an illegal value\n*          > 0:  if INFO = i, the i-th diagonal element of A is zero,\n*                indicating that the matrix is singular and the\n*                solutions X have not been computed.\n*\n\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 6)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 6)", argc);
  rb_uplo = argv[0];
  rb_trans = argv[1];
  rb_diag = argv[2];
  rb_n = argv[3];
  rb_ap = argv[4];
  rb_b = argv[5];

  if (!NA_IsNArray(rb_b))
    rb_raise(rb_eArgError, "b (6th argument) must be NArray");
  if (NA_RANK(rb_b) != 2)
    rb_raise(rb_eArgError, "rank of b (6th argument) must be %d", 2);
  nrhs = NA_SHAPE1(rb_b);
  ldb = NA_SHAPE0(rb_b);
  if (NA_TYPE(rb_b) != NA_DCOMPLEX)
    rb_b = na_change_type(rb_b, NA_DCOMPLEX);
  b = NA_PTR_TYPE(rb_b, doublecomplex*);
  diag = StringValueCStr(rb_diag)[0];
  n = NUM2INT(rb_n);
  trans = StringValueCStr(rb_trans)[0];
  uplo = StringValueCStr(rb_uplo)[0];
  if (!NA_IsNArray(rb_ap))
    rb_raise(rb_eArgError, "ap (5th argument) must be NArray");
  if (NA_RANK(rb_ap) != 1)
    rb_raise(rb_eArgError, "rank of ap (5th argument) must be %d", 1);
  if (NA_SHAPE0(rb_ap) != (n*(n+1)/2))
    rb_raise(rb_eRuntimeError, "shape 0 of ap must be %d", n*(n+1)/2);
  if (NA_TYPE(rb_ap) != NA_DCOMPLEX)
    rb_ap = na_change_type(rb_ap, NA_DCOMPLEX);
  ap = NA_PTR_TYPE(rb_ap, doublecomplex*);
  {
    int shape[2];
    shape[0] = ldb;
    shape[1] = nrhs;
    rb_b_out__ = na_make_object(NA_DCOMPLEX, 2, shape, cNArray);
  }
  b_out__ = NA_PTR_TYPE(rb_b_out__, doublecomplex*);
  MEMCPY(b_out__, b, doublecomplex, NA_TOTAL(rb_b));
  rb_b = rb_b_out__;
  b = b_out__;

  ztptrs_(&uplo, &trans, &diag, &n, &nrhs, ap, b, &ldb, &info);

  rb_info = INT2NUM(info);
  return rb_ary_new3(2, rb_info, rb_b);
}

void
init_lapack_ztptrs(VALUE mLapack){
  rb_define_module_function(mLapack, "ztptrs", rb_ztptrs, -1);
}

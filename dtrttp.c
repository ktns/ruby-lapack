#include "rb_lapack.h"

static VALUE
rb_dtrttp(int argc, VALUE *argv, VALUE self){
  VALUE rb_uplo;
  char uplo; 
  VALUE rb_a;
  doublereal *a; 
  VALUE rb_ap;
  doublereal *ap; 
  VALUE rb_info;
  integer info; 

  integer lda;
  integer n;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  ap, info = NumRu::Lapack.dtrttp( uplo, a)\n    or\n  NumRu::Lapack.dtrttp  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE DTRTTP( UPLO, N, A, LDA, AP, INFO )\n\n*  Purpose\n*  =======\n*\n*  DTRTTP copies a triangular matrix A from full format (TR) to standard\n*  packed format (TP).\n*\n\n*  Arguments\n*  =========\n*\n*  UPLO    (input) CHARACTER\n*          = 'U':  A is upper triangular.\n*          = 'L':  A is lower triangular.\n*\n*  N       (input) INTEGER\n*          The order of the matrices AP and A.  N >= 0.\n*\n*  A       (input) DOUBLE PRECISION array, dimension (LDA,N)\n*          On exit, the triangular matrix A.  If UPLO = 'U', the leading\n*          N-by-N upper triangular part of A contains the upper\n*          triangular part of the matrix A, and the strictly lower\n*          triangular part of A is not referenced.  If UPLO = 'L', the\n*          leading N-by-N lower triangular part of A contains the lower\n*          triangular part of the matrix A, and the strictly upper\n*          triangular part of A is not referenced.\n*\n*  LDA     (input) INTEGER\n*          The leading dimension of the array A.  LDA >= max(1,N).\n*\n*  AP      (output) DOUBLE PRECISION array, dimension (N*(N+1)/2\n*          On exit, the upper or lower triangular matrix A, packed\n*          columnwise in a linear array. The j-th column of A is stored\n*          in the array AP as follows:\n*          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;\n*          if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = A(i,j) for j<=i<=n.\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit\n*          < 0:  if INFO = -i, the i-th argument had an illegal value\n*\n\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 2)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 2)", argc);
  rb_uplo = argv[0];
  rb_a = argv[1];

  uplo = StringValueCStr(rb_uplo)[0];
  if (!NA_IsNArray(rb_a))
    rb_raise(rb_eArgError, "a (2th argument) must be NArray");
  if (NA_RANK(rb_a) != 2)
    rb_raise(rb_eArgError, "rank of a (2th argument) must be %d", 2);
  lda = NA_SHAPE0(rb_a);
  n = NA_SHAPE1(rb_a);
  if (NA_TYPE(rb_a) != NA_DFLOAT)
    rb_a = na_change_type(rb_a, NA_DFLOAT);
  a = NA_PTR_TYPE(rb_a, doublereal*);
  {
    int shape[1];
    shape[0] = n*(n+1)/2;
    rb_ap = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  ap = NA_PTR_TYPE(rb_ap, doublereal*);

  dtrttp_(&uplo, &n, a, &lda, ap, &info);

  rb_info = INT2NUM(info);
  return rb_ary_new3(2, rb_ap, rb_info);
}

void
init_lapack_dtrttp(VALUE mLapack){
  rb_define_module_function(mLapack, "dtrttp", rb_dtrttp, -1);
}

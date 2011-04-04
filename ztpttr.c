#include "rb_lapack.h"

extern VOID ztpttr_(char *uplo, integer *n, doublecomplex *ap, doublecomplex *a, integer *lda, integer *info);

static VALUE
rb_ztpttr(int argc, VALUE *argv, VALUE self){
  VALUE rb_uplo;
  char uplo; 
  VALUE rb_ap;
  doublecomplex *ap; 
  VALUE rb_a;
  doublecomplex *a; 
  VALUE rb_info;
  integer info; 

  integer ldap;
  integer lda;
  integer n;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  a, info = NumRu::Lapack.ztpttr( uplo, ap)\n    or\n  NumRu::Lapack.ztpttr  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE ZTPTTR( UPLO, N, AP, A, LDA, INFO )\n\n*  Purpose\n*  =======\n*\n*  ZTPTTR copies a triangular matrix A from standard packed format (TP)\n*  to standard full format (TR).\n*\n\n*  Arguments\n*  =========\n*\n*  UPLO    (input) CHARACTER*1\n*          = 'U':  A is upper triangular.\n*          = 'L':  A is lower triangular.\n*\n*  N       (input) INTEGER\n*          The order of the matrix A. N >= 0.\n*\n*  AP      (input) COMPLEX*16 array, dimension ( N*(N+1)/2 ),\n*          On entry, the upper or lower triangular matrix A, packed\n*          columnwise in a linear array. The j-th column of A is stored\n*          in the array AP as follows:\n*          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;\n*          if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = A(i,j) for j<=i<=n.\n*\n*  A       (output) COMPLEX*16 array, dimension ( LDA, N )\n*          On exit, the triangular matrix A.  If UPLO = 'U', the leading\n*          N-by-N upper triangular part of A contains the upper\n*          triangular part of the matrix A, and the strictly lower\n*          triangular part of A is not referenced.  If UPLO = 'L', the\n*          leading N-by-N lower triangular part of A contains the lower\n*          triangular part of the matrix A, and the strictly upper\n*          triangular part of A is not referenced.\n*\n*  LDA     (input) INTEGER\n*          The leading dimension of the array A.  LDA >= max(1,N).\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit\n*          < 0:  if INFO = -i, the i-th argument had an illegal value\n*\n\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 2)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 2)", argc);
  rb_uplo = argv[0];
  rb_ap = argv[1];

  uplo = StringValueCStr(rb_uplo)[0];
  if (!NA_IsNArray(rb_ap))
    rb_raise(rb_eArgError, "ap (2th argument) must be NArray");
  if (NA_RANK(rb_ap) != 1)
    rb_raise(rb_eArgError, "rank of ap (2th argument) must be %d", 1);
  ldap = NA_SHAPE0(rb_ap);
  if (NA_TYPE(rb_ap) != NA_DCOMPLEX)
    rb_ap = na_change_type(rb_ap, NA_DCOMPLEX);
  ap = NA_PTR_TYPE(rb_ap, doublecomplex*);
  n = ((int)sqrtf(ldap*8-1.0f)-1)/2;
  lda = MAX(1,n);
  {
    int shape[2];
    shape[0] = lda;
    shape[1] = n;
    rb_a = na_make_object(NA_DCOMPLEX, 2, shape, cNArray);
  }
  a = NA_PTR_TYPE(rb_a, doublecomplex*);

  ztpttr_(&uplo, &n, ap, a, &lda, &info);

  rb_info = INT2NUM(info);
  return rb_ary_new3(2, rb_a, rb_info);
}

void
init_lapack_ztpttr(VALUE mLapack){
  rb_define_module_function(mLapack, "ztpttr", rb_ztpttr, -1);
}

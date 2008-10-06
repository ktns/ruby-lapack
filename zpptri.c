#include "rb_lapack.h"

static VALUE
rb_zpptri(int argc, VALUE *argv, VALUE self){
  VALUE rb_uplo;
  char uplo; 
  VALUE rb_n;
  integer n; 
  VALUE rb_ap;
  doublecomplex *ap; 
  VALUE rb_info;
  integer info; 
  VALUE rb_ap_out__;
  doublecomplex *ap_out__;


  if (argc == 0) {
    printf("%s\n", "USAGE:\n  info, ap = NumRu::Lapack.zpptri( uplo, n, ap)\n    or\n  NumRu::Lapack.zpptri  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE ZPPTRI( UPLO, N, AP, INFO )\n\n*  Purpose\n*  =======\n*\n*  ZPPTRI computes the inverse of a complex Hermitian positive definite\n*  matrix A using the Cholesky factorization A = U**H*U or A = L*L**H\n*  computed by ZPPTRF.\n*\n\n*  Arguments\n*  =========\n*\n*  UPLO    (input) CHARACTER*1\n*          = 'U':  Upper triangular factor is stored in AP;\n*          = 'L':  Lower triangular factor is stored in AP.\n*\n*  N       (input) INTEGER\n*          The order of the matrix A.  N >= 0.\n*\n*  AP      (input/output) COMPLEX*16 array, dimension (N*(N+1)/2)\n*          On entry, the triangular factor U or L from the Cholesky\n*          factorization A = U**H*U or A = L*L**H, packed columnwise as\n*          a linear array.  The j-th column of U or L is stored in the\n*          array AP as follows:\n*          if UPLO = 'U', AP(i + (j-1)*j/2) = U(i,j) for 1<=i<=j;\n*          if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = L(i,j) for j<=i<=n.\n*\n*          On exit, the upper or lower triangle of the (Hermitian)\n*          inverse of A, overwriting the input factor U or L.\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit\n*          < 0:  if INFO = -i, the i-th argument had an illegal value\n*          > 0:  if INFO = i, the (i,i) element of the factor U or L is\n*                zero, and the inverse could not be computed.\n*\n\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 3)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 3)", argc);
  rb_uplo = argv[0];
  rb_n = argv[1];
  rb_ap = argv[2];

  uplo = StringValueCStr(rb_uplo)[0];
  n = NUM2INT(rb_n);
  if (!NA_IsNArray(rb_ap))
    rb_raise(rb_eArgError, "ap (3th argument) must be NArray");
  if (NA_RANK(rb_ap) != 1)
    rb_raise(rb_eArgError, "rank of ap (3th argument) must be %d", 1);
  if (NA_SHAPE0(rb_ap) != (n*(n+1)/2))
    rb_raise(rb_eRuntimeError, "shape 0 of ap must be %d", n*(n+1)/2);
  if (NA_TYPE(rb_ap) != NA_DCOMPLEX)
    rb_ap = na_change_type(rb_ap, NA_DCOMPLEX);
  ap = NA_PTR_TYPE(rb_ap, doublecomplex*);
  {
    int shape[1];
    shape[0] = n*(n+1)/2;
    rb_ap_out__ = na_make_object(NA_DCOMPLEX, 1, shape, cNArray);
  }
  ap_out__ = NA_PTR_TYPE(rb_ap_out__, doublecomplex*);
  MEMCPY(ap_out__, ap, doublecomplex, NA_TOTAL(rb_ap));
  rb_ap = rb_ap_out__;
  ap = ap_out__;

  zpptri_(&uplo, &n, ap, &info);

  rb_info = INT2NUM(info);
  return rb_ary_new3(2, rb_info, rb_ap);
}

void
init_lapack_zpptri(VALUE mLapack){
  rb_define_module_function(mLapack, "zpptri", rb_zpptri, -1);
}

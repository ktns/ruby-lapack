#include "rb_lapack.h"

extern VOID zpptrf_(char *uplo, integer *n, doublecomplex *ap, integer *info);

static VALUE
rb_zpptrf(int argc, VALUE *argv, VALUE self){
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
    printf("%s\n", "USAGE:\n  info, ap = NumRu::Lapack.zpptrf( uplo, n, ap)\n    or\n  NumRu::Lapack.zpptrf  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE ZPPTRF( UPLO, N, AP, INFO )\n\n*  Purpose\n*  =======\n*\n*  ZPPTRF computes the Cholesky factorization of a complex Hermitian\n*  positive definite matrix A stored in packed format.\n*\n*  The factorization has the form\n*     A = U**H * U,  if UPLO = 'U', or\n*     A = L  * L**H,  if UPLO = 'L',\n*  where U is an upper triangular matrix and L is lower triangular.\n*\n\n*  Arguments\n*  =========\n*\n*  UPLO    (input) CHARACTER*1\n*          = 'U':  Upper triangle of A is stored;\n*          = 'L':  Lower triangle of A is stored.\n*\n*  N       (input) INTEGER\n*          The order of the matrix A.  N >= 0.\n*\n*  AP      (input/output) COMPLEX*16 array, dimension (N*(N+1)/2)\n*          On entry, the upper or lower triangle of the Hermitian matrix\n*          A, packed columnwise in a linear array.  The j-th column of A\n*          is stored in the array AP as follows:\n*          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;\n*          if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = A(i,j) for j<=i<=n.\n*          See below for further details.\n*\n*          On exit, if INFO = 0, the triangular factor U or L from the\n*          Cholesky factorization A = U**H*U or A = L*L**H, in the same\n*          storage format as A.\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit\n*          < 0:  if INFO = -i, the i-th argument had an illegal value\n*          > 0:  if INFO = i, the leading minor of order i is not\n*                positive definite, and the factorization could not be\n*                completed.\n*\n\n*  Further Details\n*  ===============\n*\n*  The packed storage scheme is illustrated by the following example\n*  when N = 4, UPLO = 'U':\n*\n*  Two-dimensional storage of the Hermitian matrix A:\n*\n*     a11 a12 a13 a14\n*         a22 a23 a24\n*             a33 a34     (aij = conjg(aji))\n*                 a44\n*\n*  Packed storage of the upper triangle of A:\n*\n*  AP = [ a11, a12, a22, a13, a23, a33, a14, a24, a34, a44 ]\n*\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 3)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 3)", argc);
  rb_uplo = argv[0];
  rb_n = argv[1];
  rb_ap = argv[2];

  n = NUM2INT(rb_n);
  uplo = StringValueCStr(rb_uplo)[0];
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

  zpptrf_(&uplo, &n, ap, &info);

  rb_info = INT2NUM(info);
  return rb_ary_new3(2, rb_info, rb_ap);
}

void
init_lapack_zpptrf(VALUE mLapack){
  rb_define_module_function(mLapack, "zpptrf", rb_zpptrf, -1);
}

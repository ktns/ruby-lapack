#include "rb_lapack.h"

static VALUE
rb_stptri(int argc, VALUE *argv, VALUE self){
  VALUE rb_uplo;
  char uplo; 
  VALUE rb_diag;
  char diag; 
  VALUE rb_n;
  integer n; 
  VALUE rb_ap;
  real *ap; 
  VALUE rb_info;
  integer info; 
  VALUE rb_ap_out__;
  real *ap_out__;


  if (argc == 0) {
    printf("%s\n", "USAGE:\n  info, ap = NumRu::Lapack.stptri( uplo, diag, n, ap)\n    or\n  NumRu::Lapack.stptri  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE STPTRI( UPLO, DIAG, N, AP, INFO )\n\n*  Purpose\n*  =======\n*\n*  STPTRI computes the inverse of a real upper or lower triangular\n*  matrix A stored in packed format.\n*\n\n*  Arguments\n*  =========\n*\n*  UPLO    (input) CHARACTER*1\n*          = 'U':  A is upper triangular;\n*          = 'L':  A is lower triangular.\n*\n*  DIAG    (input) CHARACTER*1\n*          = 'N':  A is non-unit triangular;\n*          = 'U':  A is unit triangular.\n*\n*  N       (input) INTEGER\n*          The order of the matrix A.  N >= 0.\n*\n*  AP      (input/output) REAL array, dimension (N*(N+1)/2)\n*          On entry, the upper or lower triangular matrix A, stored\n*          columnwise in a linear array.  The j-th column of A is stored\n*          in the array AP as follows:\n*          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;\n*          if UPLO = 'L', AP(i + (j-1)*((2*n-j)/2) = A(i,j) for j<=i<=n.\n*          See below for further details.\n*          On exit, the (triangular) inverse of the original matrix, in\n*          the same packed storage format.\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit\n*          < 0:  if INFO = -i, the i-th argument had an illegal value\n*          > 0:  if INFO = i, A(i,i) is exactly zero.  The triangular\n*                matrix is singular and its inverse can not be computed.\n*\n\n*  Further Details\n*  ===============\n*\n*  A triangular matrix A can be transferred to packed storage using one\n*  of the following program segments:\n*\n*  UPLO = 'U':                      UPLO = 'L':\n*\n*        JC = 1                           JC = 1\n*        DO 2 J = 1, N                    DO 2 J = 1, N\n*           DO 1 I = 1, J                    DO 1 I = J, N\n*              AP(JC+I-1) = A(I,J)              AP(JC+I-J) = A(I,J)\n*      1    CONTINUE                    1    CONTINUE\n*           JC = JC + J                      JC = JC + N - J + 1\n*      2 CONTINUE                       2 CONTINUE\n*\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 4)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 4)", argc);
  rb_uplo = argv[0];
  rb_diag = argv[1];
  rb_n = argv[2];
  rb_ap = argv[3];

  uplo = StringValueCStr(rb_uplo)[0];
  diag = StringValueCStr(rb_diag)[0];
  n = NUM2INT(rb_n);
  if (!NA_IsNArray(rb_ap))
    rb_raise(rb_eArgError, "ap (4th argument) must be NArray");
  if (NA_RANK(rb_ap) != 1)
    rb_raise(rb_eArgError, "rank of ap (4th argument) must be %d", 1);
  if (NA_SHAPE0(rb_ap) != (n*(n+1)/2))
    rb_raise(rb_eRuntimeError, "shape 0 of ap must be %d", n*(n+1)/2);
  if (NA_TYPE(rb_ap) != NA_SFLOAT)
    rb_ap = na_change_type(rb_ap, NA_SFLOAT);
  ap = NA_PTR_TYPE(rb_ap, real*);
  {
    int shape[1];
    shape[0] = n*(n+1)/2;
    rb_ap_out__ = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  ap_out__ = NA_PTR_TYPE(rb_ap_out__, real*);
  MEMCPY(ap_out__, ap, real, NA_TOTAL(rb_ap));
  rb_ap = rb_ap_out__;
  ap = ap_out__;

  stptri_(&uplo, &diag, &n, ap, &info);

  rb_info = INT2NUM(info);
  return rb_ary_new3(2, rb_info, rb_ap);
}

void
init_lapack_stptri(VALUE mLapack){
  rb_define_module_function(mLapack, "stptri", rb_stptri, -1);
}

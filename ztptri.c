#include "rb_lapack.h"

extern VOID ztptri_(char *uplo, char *diag, integer *n, doublecomplex *ap, integer *info);

static VALUE sHelp, sUsage;

static VALUE
rblapack_ztptri(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_uplo;
  char uplo; 
  VALUE rblapack_diag;
  char diag; 
  VALUE rblapack_n;
  integer n; 
  VALUE rblapack_ap;
  doublecomplex *ap; 
  VALUE rblapack_info;
  integer info; 
  VALUE rblapack_ap_out__;
  doublecomplex *ap_out__;


  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  info, ap = NumRu::Lapack.ztptri( uplo, diag, n, ap, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE ZTPTRI( UPLO, DIAG, N, AP, INFO )\n\n*  Purpose\n*  =======\n*\n*  ZTPTRI computes the inverse of a complex upper or lower triangular\n*  matrix A stored in packed format.\n*\n\n*  Arguments\n*  =========\n*\n*  UPLO    (input) CHARACTER*1\n*          = 'U':  A is upper triangular;\n*          = 'L':  A is lower triangular.\n*\n*  DIAG    (input) CHARACTER*1\n*          = 'N':  A is non-unit triangular;\n*          = 'U':  A is unit triangular.\n*\n*  N       (input) INTEGER\n*          The order of the matrix A.  N >= 0.\n*\n*  AP      (input/output) COMPLEX*16 array, dimension (N*(N+1)/2)\n*          On entry, the upper or lower triangular matrix A, stored\n*          columnwise in a linear array.  The j-th column of A is stored\n*          in the array AP as follows:\n*          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;\n*          if UPLO = 'L', AP(i + (j-1)*((2*n-j)/2) = A(i,j) for j<=i<=n.\n*          See below for further details.\n*          On exit, the (triangular) inverse of the original matrix, in\n*          the same packed storage format.\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit\n*          < 0:  if INFO = -i, the i-th argument had an illegal value\n*          > 0:  if INFO = i, A(i,i) is exactly zero.  The triangular\n*                matrix is singular and its inverse can not be computed.\n*\n\n*  Further Details\n*  ===============\n*\n*  A triangular matrix A can be transferred to packed storage using one\n*  of the following program segments:\n*\n*  UPLO = 'U':                      UPLO = 'L':\n*\n*        JC = 1                           JC = 1\n*        DO 2 J = 1, N                    DO 2 J = 1, N\n*           DO 1 I = 1, J                    DO 1 I = J, N\n*              AP(JC+I-1) = A(I,J)              AP(JC+I-J) = A(I,J)\n*      1    CONTINUE                    1    CONTINUE\n*           JC = JC + J                      JC = JC + N - J + 1\n*      2 CONTINUE                       2 CONTINUE\n*\n*  =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  info, ap = NumRu::Lapack.ztptri( uplo, diag, n, ap, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 4)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 4)", argc);
  rblapack_uplo = argv[0];
  rblapack_diag = argv[1];
  rblapack_n = argv[2];
  rblapack_ap = argv[3];
  if (rb_options != Qnil) {
  }

  diag = StringValueCStr(rblapack_diag)[0];
  n = NUM2INT(rblapack_n);
  uplo = StringValueCStr(rblapack_uplo)[0];
  if (!NA_IsNArray(rblapack_ap))
    rb_raise(rb_eArgError, "ap (4th argument) must be NArray");
  if (NA_RANK(rblapack_ap) != 1)
    rb_raise(rb_eArgError, "rank of ap (4th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_ap) != (n*(n+1)/2))
    rb_raise(rb_eRuntimeError, "shape 0 of ap must be %d", n*(n+1)/2);
  if (NA_TYPE(rblapack_ap) != NA_DCOMPLEX)
    rblapack_ap = na_change_type(rblapack_ap, NA_DCOMPLEX);
  ap = NA_PTR_TYPE(rblapack_ap, doublecomplex*);
  {
    int shape[1];
    shape[0] = n*(n+1)/2;
    rblapack_ap_out__ = na_make_object(NA_DCOMPLEX, 1, shape, cNArray);
  }
  ap_out__ = NA_PTR_TYPE(rblapack_ap_out__, doublecomplex*);
  MEMCPY(ap_out__, ap, doublecomplex, NA_TOTAL(rblapack_ap));
  rblapack_ap = rblapack_ap_out__;
  ap = ap_out__;

  ztptri_(&uplo, &diag, &n, ap, &info);

  rblapack_info = INT2NUM(info);
  return rb_ary_new3(2, rblapack_info, rblapack_ap);
}

void
init_lapack_ztptri(VALUE mLapack){
  rb_define_module_function(mLapack, "ztptri", rblapack_ztptri, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}

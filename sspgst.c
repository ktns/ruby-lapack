#include "rb_lapack.h"

static VALUE
rb_sspgst(int argc, VALUE *argv, VALUE self){
  VALUE rb_itype;
  integer itype; 
  VALUE rb_uplo;
  char uplo; 
  VALUE rb_n;
  integer n; 
  VALUE rb_ap;
  real *ap; 
  VALUE rb_bp;
  real *bp; 
  VALUE rb_info;
  integer info; 
  VALUE rb_ap_out__;
  real *ap_out__;


  if (argc == 0) {
    printf("%s\n", "USAGE:\n  info, ap = NumRu::Lapack.sspgst( itype, uplo, n, ap, bp)\n    or\n  NumRu::Lapack.sspgst  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE SSPGST( ITYPE, UPLO, N, AP, BP, INFO )\n\n*  Purpose\n*  =======\n*\n*  SSPGST reduces a real symmetric-definite generalized eigenproblem\n*  to standard form, using packed storage.\n*\n*  If ITYPE = 1, the problem is A*x = lambda*B*x,\n*  and A is overwritten by inv(U**T)*A*inv(U) or inv(L)*A*inv(L**T)\n*\n*  If ITYPE = 2 or 3, the problem is A*B*x = lambda*x or\n*  B*A*x = lambda*x, and A is overwritten by U*A*U**T or L**T*A*L.\n*\n*  B must have been previously factorized as U**T*U or L*L**T by SPPTRF.\n*\n\n*  Arguments\n*  =========\n*\n*  ITYPE   (input) INTEGER\n*          = 1: compute inv(U**T)*A*inv(U) or inv(L)*A*inv(L**T);\n*          = 2 or 3: compute U*A*U**T or L**T*A*L.\n*\n*  UPLO    (input) CHARACTER*1\n*          = 'U':  Upper triangle of A is stored and B is factored as\n*                  U**T*U;\n*          = 'L':  Lower triangle of A is stored and B is factored as\n*                  L*L**T.\n*\n*  N       (input) INTEGER\n*          The order of the matrices A and B.  N >= 0.\n*\n*  AP      (input/output) REAL array, dimension (N*(N+1)/2)\n*          On entry, the upper or lower triangle of the symmetric matrix\n*          A, packed columnwise in a linear array.  The j-th column of A\n*          is stored in the array AP as follows:\n*          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;\n*          if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = A(i,j) for j<=i<=n.\n*\n*          On exit, if INFO = 0, the transformed matrix, stored in the\n*          same format as A.\n*\n*  BP      (input) REAL array, dimension (N*(N+1)/2)\n*          The triangular factor from the Cholesky factorization of B,\n*          stored in the same format as A, as returned by SPPTRF.\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit\n*          < 0:  if INFO = -i, the i-th argument had an illegal value\n*\n\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 5)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 5)", argc);
  rb_itype = argv[0];
  rb_uplo = argv[1];
  rb_n = argv[2];
  rb_ap = argv[3];
  rb_bp = argv[4];

  itype = NUM2INT(rb_itype);
  uplo = StringValueCStr(rb_uplo)[0];
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
  if (!NA_IsNArray(rb_bp))
    rb_raise(rb_eArgError, "bp (5th argument) must be NArray");
  if (NA_RANK(rb_bp) != 1)
    rb_raise(rb_eArgError, "rank of bp (5th argument) must be %d", 1);
  if (NA_SHAPE0(rb_bp) != (n*(n+1)/2))
    rb_raise(rb_eRuntimeError, "shape 0 of bp must be %d", n*(n+1)/2);
  if (NA_TYPE(rb_bp) != NA_SFLOAT)
    rb_bp = na_change_type(rb_bp, NA_SFLOAT);
  bp = NA_PTR_TYPE(rb_bp, real*);
  {
    int shape[1];
    shape[0] = n*(n+1)/2;
    rb_ap_out__ = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  ap_out__ = NA_PTR_TYPE(rb_ap_out__, real*);
  MEMCPY(ap_out__, ap, real, NA_TOTAL(rb_ap));
  rb_ap = rb_ap_out__;
  ap = ap_out__;

  sspgst_(&itype, &uplo, &n, ap, bp, &info);

  rb_info = INT2NUM(info);
  return rb_ary_new3(2, rb_info, rb_ap);
}

void
init_lapack_sspgst(VALUE mLapack){
  rb_define_module_function(mLapack, "sspgst", rb_sspgst, -1);
}

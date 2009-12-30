#include "rb_lapack.h"

static VALUE
rb_zspsv(int argc, VALUE *argv, VALUE self){
  VALUE rb_uplo;
  char uplo; 
  VALUE rb_ap;
  doublecomplex *ap; 
  VALUE rb_b;
  doublecomplex *b; 
  VALUE rb_ipiv;
  integer *ipiv; 
  VALUE rb_info;
  integer info; 
  VALUE rb_ap_out__;
  doublecomplex *ap_out__;
  VALUE rb_b_out__;
  doublecomplex *b_out__;

  integer n;
  integer ldb;
  integer nrhs;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  ipiv, info, ap, b = NumRu::Lapack.zspsv( uplo, ap, b)\n    or\n  NumRu::Lapack.zspsv  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE ZSPSV( UPLO, N, NRHS, AP, IPIV, B, LDB, INFO )\n\n*  Purpose\n*  =======\n*\n*  ZSPSV computes the solution to a complex system of linear equations\n*     A * X = B,\n*  where A is an N-by-N symmetric matrix stored in packed format and X\n*  and B are N-by-NRHS matrices.\n*\n*  The diagonal pivoting method is used to factor A as\n*     A = U * D * U**T,  if UPLO = 'U', or\n*     A = L * D * L**T,  if UPLO = 'L',\n*  where U (or L) is a product of permutation and unit upper (lower)\n*  triangular matrices, D is symmetric and block diagonal with 1-by-1\n*  and 2-by-2 diagonal blocks.  The factored form of A is then used to\n*  solve the system of equations A * X = B.\n*\n\n*  Arguments\n*  =========\n*\n*  UPLO    (input) CHARACTER*1\n*          = 'U':  Upper triangle of A is stored;\n*          = 'L':  Lower triangle of A is stored.\n*\n*  N       (input) INTEGER\n*          The number of linear equations, i.e., the order of the\n*          matrix A.  N >= 0.\n*\n*  NRHS    (input) INTEGER\n*          The number of right hand sides, i.e., the number of columns\n*          of the matrix B.  NRHS >= 0.\n*\n*  AP      (input/output) COMPLEX*16 array, dimension (N*(N+1)/2)\n*          On entry, the upper or lower triangle of the symmetric matrix\n*          A, packed columnwise in a linear array.  The j-th column of A\n*          is stored in the array AP as follows:\n*          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;\n*          if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = A(i,j) for j<=i<=n.\n*          See below for further details.\n*\n*          On exit, the block diagonal matrix D and the multipliers used\n*          to obtain the factor U or L from the factorization\n*          A = U*D*U**T or A = L*D*L**T as computed by ZSPTRF, stored as\n*          a packed triangular matrix in the same storage format as A.\n*\n*  IPIV    (output) INTEGER array, dimension (N)\n*          Details of the interchanges and the block structure of D, as\n*          determined by ZSPTRF.  If IPIV(k) > 0, then rows and columns\n*          k and IPIV(k) were interchanged, and D(k,k) is a 1-by-1\n*          diagonal block.  If UPLO = 'U' and IPIV(k) = IPIV(k-1) < 0,\n*          then rows and columns k-1 and -IPIV(k) were interchanged and\n*          D(k-1:k,k-1:k) is a 2-by-2 diagonal block.  If UPLO = 'L' and\n*          IPIV(k) = IPIV(k+1) < 0, then rows and columns k+1 and\n*          -IPIV(k) were interchanged and D(k:k+1,k:k+1) is a 2-by-2\n*          diagonal block.\n*\n*  B       (input/output) COMPLEX*16 array, dimension (LDB,NRHS)\n*          On entry, the N-by-NRHS right hand side matrix B.\n*          On exit, if INFO = 0, the N-by-NRHS solution matrix X.\n*\n*  LDB     (input) INTEGER\n*          The leading dimension of the array B.  LDB >= max(1,N).\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit\n*          < 0:  if INFO = -i, the i-th argument had an illegal value\n*          > 0:  if INFO = i, D(i,i) is exactly zero.  The factorization\n*                has been completed, but the block diagonal matrix D is\n*                exactly singular, so the solution could not be\n*                computed.\n*\n\n*  Further Details\n*  ===============\n*\n*  The packed storage scheme is illustrated by the following example\n*  when N = 4, UPLO = 'U':\n*\n*  Two-dimensional storage of the symmetric matrix A:\n*\n*     a11 a12 a13 a14\n*         a22 a23 a24\n*             a33 a34     (aij = aji)\n*                 a44\n*\n*  Packed storage of the upper triangle of A:\n*\n*  AP = [ a11, a12, a22, a13, a23, a33, a14, a24, a34, a44 ]\n*\n*  =====================================================================\n*\n*     .. External Functions ..\n      LOGICAL            LSAME\n      EXTERNAL           LSAME\n*     ..\n*     .. External Subroutines ..\n      EXTERNAL           XERBLA, ZSPTRF, ZSPTRS\n*     ..\n*     .. Intrinsic Functions ..\n      INTRINSIC          MAX\n*     ..\n\n");
    return Qnil;
  }
  if (argc != 3)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 3)", argc);
  rb_uplo = argv[0];
  rb_ap = argv[1];
  rb_b = argv[2];

  uplo = StringValueCStr(rb_uplo)[0];
  if (!NA_IsNArray(rb_b))
    rb_raise(rb_eArgError, "b (2th argument) must be NArray");
  if (NA_RANK(rb_b) != 2)
    rb_raise(rb_eArgError, "rank of b (2th argument) must be %d", 2);
  ldb = NA_SHAPE0(rb_b);
  nrhs = NA_SHAPE1(rb_b);
  if (NA_TYPE(rb_b) != NA_DCOMPLEX)
    rb_b = na_change_type(rb_b, NA_DCOMPLEX);
  b = NA_PTR_TYPE(rb_b, doublecomplex*);
  if (!NA_IsNArray(rb_ap))
    rb_raise(rb_eArgError, "ap (3th argument) must be NArray");
  if (NA_RANK(rb_ap) != 1)
    rb_raise(rb_eArgError, "rank of ap (3th argument) must be %d", 1);
  n = ldb;
  if (NA_SHAPE0(rb_ap) != (n*(n+1)/2))
    rb_raise(rb_eRuntimeError, "shape 0 of ap must be %d", n*(n+1)/2);
  if (NA_TYPE(rb_ap) != NA_DCOMPLEX)
    rb_ap = na_change_type(rb_ap, NA_DCOMPLEX);
  ap = NA_PTR_TYPE(rb_ap, doublecomplex*);
  {
    int shape[1];
    shape[0] = DIM_LEN(n);
    rb_ipiv = na_make_object(NA_LINT, 1, shape, cNArray);
  }
  ipiv = NA_PTR_TYPE(rb_ipiv, integer*);
  {
    int shape[1];
    shape[0] = n*(n+1)/2;
    rb_ap_out__ = na_make_object(NA_DCOMPLEX, 1, shape, cNArray);
  }
  ap_out__ = NA_PTR_TYPE(rb_ap_out__, doublecomplex*);
  MEMCPY(ap_out__, ap, doublecomplex, NA_TOTAL(rb_ap));
  rb_ap = rb_ap_out__;
  ap = ap_out__;
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

  zspsv_(&uplo, &n, &nrhs, ap, ipiv, b, &ldb, &info);

  rb_info = INT2NUM(info);
  return rb_ary_new3(4, rb_ipiv, rb_info, rb_ap, rb_b);
}

void
init_lapack_zspsv(VALUE mLapack){
  rb_define_module_function(mLapack, "zspsv", rb_zspsv, -1);
}

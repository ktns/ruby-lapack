#include "rb_lapack.h"

static VALUE
rb_zpbtf2(int argc, VALUE *argv, VALUE self){
  VALUE rb_uplo;
  char uplo; 
  VALUE rb_kd;
  integer kd; 
  VALUE rb_ab;
  doublecomplex *ab; 
  VALUE rb_info;
  integer info; 
  VALUE rb_ab_out__;
  doublecomplex *ab_out__;

  integer ldab;
  integer n;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  info, ab = NumRu::Lapack.zpbtf2( uplo, kd, ab)\n    or\n  NumRu::Lapack.zpbtf2  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE ZPBTF2( UPLO, N, KD, AB, LDAB, INFO )\n\n*  Purpose\n*  =======\n*\n*  ZPBTF2 computes the Cholesky factorization of a complex Hermitian\n*  positive definite band matrix A.\n*\n*  The factorization has the form\n*     A = U' * U ,  if UPLO = 'U', or\n*     A = L  * L',  if UPLO = 'L',\n*  where U is an upper triangular matrix, U' is the conjugate transpose\n*  of U, and L is lower triangular.\n*\n*  This is the unblocked version of the algorithm, calling Level 2 BLAS.\n*\n\n*  Arguments\n*  =========\n*\n*  UPLO    (input) CHARACTER*1\n*          Specifies whether the upper or lower triangular part of the\n*          Hermitian matrix A is stored:\n*          = 'U':  Upper triangular\n*          = 'L':  Lower triangular\n*\n*  N       (input) INTEGER\n*          The order of the matrix A.  N >= 0.\n*\n*  KD      (input) INTEGER\n*          The number of super-diagonals of the matrix A if UPLO = 'U',\n*          or the number of sub-diagonals if UPLO = 'L'.  KD >= 0.\n*\n*  AB      (input/output) COMPLEX*16 array, dimension (LDAB,N)\n*          On entry, the upper or lower triangle of the Hermitian band\n*          matrix A, stored in the first KD+1 rows of the array.  The\n*          j-th column of A is stored in the j-th column of the array AB\n*          as follows:\n*          if UPLO = 'U', AB(kd+1+i-j,j) = A(i,j) for max(1,j-kd)<=i<=j;\n*          if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=min(n,j+kd).\n*\n*          On exit, if INFO = 0, the triangular factor U or L from the\n*          Cholesky factorization A = U'*U or A = L*L' of the band\n*          matrix A, in the same storage format as A.\n*\n*  LDAB    (input) INTEGER\n*          The leading dimension of the array AB.  LDAB >= KD+1.\n*\n*  INFO    (output) INTEGER\n*          = 0: successful exit\n*          < 0: if INFO = -k, the k-th argument had an illegal value\n*          > 0: if INFO = k, the leading minor of order k is not\n*               positive definite, and the factorization could not be\n*               completed.\n*\n\n*  Further Details\n*  ===============\n*\n*  The band storage scheme is illustrated by the following example, when\n*  N = 6, KD = 2, and UPLO = 'U':\n*\n*  On entry:                       On exit:\n*\n*      *    *   a13  a24  a35  a46      *    *   u13  u24  u35  u46\n*      *   a12  a23  a34  a45  a56      *   u12  u23  u34  u45  u56\n*     a11  a22  a33  a44  a55  a66     u11  u22  u33  u44  u55  u66\n*\n*  Similarly, if UPLO = 'L' the format of A is as follows:\n*\n*  On entry:                       On exit:\n*\n*     a11  a22  a33  a44  a55  a66     l11  l22  l33  l44  l55  l66\n*     a21  a32  a43  a54  a65   *      l21  l32  l43  l54  l65   *\n*     a31  a42  a53  a64   *    *      l31  l42  l53  l64   *    *\n*\n*  Array elements marked * are not used by the routine.\n*\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 3)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 3)", argc);
  rb_uplo = argv[0];
  rb_kd = argv[1];
  rb_ab = argv[2];

  uplo = StringValueCStr(rb_uplo)[0];
  kd = NUM2INT(rb_kd);
  if (!NA_IsNArray(rb_ab))
    rb_raise(rb_eArgError, "ab (3th argument) must be NArray");
  if (NA_RANK(rb_ab) != 2)
    rb_raise(rb_eArgError, "rank of ab (3th argument) must be %d", 2);
  ldab = NA_SHAPE0(rb_ab);
  n = NA_SHAPE1(rb_ab);
  if (NA_TYPE(rb_ab) != NA_DCOMPLEX)
    rb_ab = na_change_type(rb_ab, NA_DCOMPLEX);
  ab = NA_PTR_TYPE(rb_ab, doublecomplex*);
  {
    int shape[2];
    shape[0] = ldab;
    shape[1] = n;
    rb_ab_out__ = na_make_object(NA_DCOMPLEX, 2, shape, cNArray);
  }
  ab_out__ = NA_PTR_TYPE(rb_ab_out__, doublecomplex*);
  MEMCPY(ab_out__, ab, doublecomplex, NA_TOTAL(rb_ab));
  rb_ab = rb_ab_out__;
  ab = ab_out__;

  zpbtf2_(&uplo, &n, &kd, ab, &ldab, &info);

  rb_info = INT2NUM(info);
  return rb_ary_new3(2, rb_info, rb_ab);
}

void
init_lapack_zpbtf2(VALUE mLapack){
  rb_define_module_function(mLapack, "zpbtf2", rb_zpbtf2, -1);
}

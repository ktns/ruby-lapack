#include "rb_lapack.h"

extern doublereal zlanhb_(char *norm, char *uplo, integer *n, integer *k, doublecomplex *ab, integer *ldab, doublereal *work);

static VALUE
rb_zlanhb(int argc, VALUE *argv, VALUE self){
  VALUE rb_norm;
  char norm; 
  VALUE rb_uplo;
  char uplo; 
  VALUE rb_k;
  integer k; 
  VALUE rb_ab;
  doublecomplex *ab; 
  VALUE rb___out__;
  doublereal __out__; 
  doublereal *work;

  integer ldab;
  integer n;
  integer lwork;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  __out__ = NumRu::Lapack.zlanhb( norm, uplo, k, ab)\n    or\n  NumRu::Lapack.zlanhb  # print help\n\n\nFORTRAN MANUAL\n      DOUBLE PRECISION FUNCTION ZLANHB( NORM, UPLO, N, K, AB, LDAB, WORK )\n\n*  Purpose\n*  =======\n*\n*  ZLANHB  returns the value of the one norm,  or the Frobenius norm, or\n*  the  infinity norm,  or the element of  largest absolute value  of an\n*  n by n hermitian band matrix A,  with k super-diagonals.\n*\n*  Description\n*  ===========\n*\n*  ZLANHB returns the value\n*\n*     ZLANHB = ( max(abs(A(i,j))), NORM = 'M' or 'm'\n*              (\n*              ( norm1(A),         NORM = '1', 'O' or 'o'\n*              (\n*              ( normI(A),         NORM = 'I' or 'i'\n*              (\n*              ( normF(A),         NORM = 'F', 'f', 'E' or 'e'\n*\n*  where  norm1  denotes the  one norm of a matrix (maximum column sum),\n*  normI  denotes the  infinity norm  of a matrix  (maximum row sum) and\n*  normF  denotes the  Frobenius norm of a matrix (square root of sum of\n*  squares).  Note that  max(abs(A(i,j)))  is not a consistent matrix norm.\n*\n\n*  Arguments\n*  =========\n*\n*  NORM    (input) CHARACTER*1\n*          Specifies the value to be returned in ZLANHB as described\n*          above.\n*\n*  UPLO    (input) CHARACTER*1\n*          Specifies whether the upper or lower triangular part of the\n*          band matrix A is supplied.\n*          = 'U':  Upper triangular\n*          = 'L':  Lower triangular\n*\n*  N       (input) INTEGER\n*          The order of the matrix A.  N >= 0.  When N = 0, ZLANHB is\n*          set to zero.\n*\n*  K       (input) INTEGER\n*          The number of super-diagonals or sub-diagonals of the\n*          band matrix A.  K >= 0.\n*\n*  AB      (input) COMPLEX*16 array, dimension (LDAB,N)\n*          The upper or lower triangle of the hermitian band matrix A,\n*          stored in the first K+1 rows of AB.  The j-th column of A is\n*          stored in the j-th column of the array AB as follows:\n*          if UPLO = 'U', AB(k+1+i-j,j) = A(i,j) for max(1,j-k)<=i<=j;\n*          if UPLO = 'L', AB(1+i-j,j)   = A(i,j) for j<=i<=min(n,j+k).\n*          Note that the imaginary parts of the diagonal elements need\n*          not be set and are assumed to be zero.\n*\n*  LDAB    (input) INTEGER\n*          The leading dimension of the array AB.  LDAB >= K+1.\n*\n*  WORK    (workspace) DOUBLE PRECISION array, dimension (MAX(1,LWORK)),\n*          where LWORK >= N when NORM = 'I' or '1' or 'O'; otherwise,\n*          WORK is not referenced.\n*\n\n* =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 4)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 4)", argc);
  rb_norm = argv[0];
  rb_uplo = argv[1];
  rb_k = argv[2];
  rb_ab = argv[3];

  k = NUM2INT(rb_k);
  if (!NA_IsNArray(rb_ab))
    rb_raise(rb_eArgError, "ab (4th argument) must be NArray");
  if (NA_RANK(rb_ab) != 2)
    rb_raise(rb_eArgError, "rank of ab (4th argument) must be %d", 2);
  n = NA_SHAPE1(rb_ab);
  ldab = NA_SHAPE0(rb_ab);
  if (NA_TYPE(rb_ab) != NA_DCOMPLEX)
    rb_ab = na_change_type(rb_ab, NA_DCOMPLEX);
  ab = NA_PTR_TYPE(rb_ab, doublecomplex*);
  norm = StringValueCStr(rb_norm)[0];
  uplo = StringValueCStr(rb_uplo)[0];
  lwork = ((lsame_(&norm,"I")) || ((('1') || ('o')))) ? n : 0;
  work = ALLOC_N(doublereal, (MAX(1,lwork)));

  __out__ = zlanhb_(&norm, &uplo, &n, &k, ab, &ldab, work);

  free(work);
  rb___out__ = rb_float_new((double)__out__);
  return rb___out__;
}

void
init_lapack_zlanhb(VALUE mLapack){
  rb_define_module_function(mLapack, "zlanhb", rb_zlanhb, -1);
}

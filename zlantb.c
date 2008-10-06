#include "rb_lapack.h"

extern VOID zlantb_(doublereal *__out__, char *norm, char *uplo, char *diag, integer *n, integer *k, doublecomplex *ab, integer *ldab, doublereal *work);
static VALUE
rb_zlantb(int argc, VALUE *argv, VALUE self){
  VALUE rb_norm;
  char norm; 
  VALUE rb_uplo;
  char uplo; 
  VALUE rb_diag;
  char diag; 
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
    printf("%s\n", "USAGE:\n  __out__ = NumRu::Lapack.zlantb( norm, uplo, diag, k, ab)\n    or\n  NumRu::Lapack.zlantb  # print help\n\n\nFORTRAN MANUAL\n      DOUBLE PRECISION FUNCTION ZLANTB( NORM, UPLO, DIAG, N, K, AB, LDAB, WORK )\n\n*  Purpose\n*  =======\n*\n*  ZLANTB  returns the value of the one norm,  or the Frobenius norm, or\n*  the  infinity norm,  or the element of  largest absolute value  of an\n*  n by n triangular band matrix A,  with ( k + 1 ) diagonals.\n*\n*  Description\n*  ===========\n*\n*  ZLANTB returns the value\n*\n*     ZLANTB = ( max(abs(A(i,j))), NORM = 'M' or 'm'\n*              (\n*              ( norm1(A),         NORM = '1', 'O' or 'o'\n*              (\n*              ( normI(A),         NORM = 'I' or 'i'\n*              (\n*              ( normF(A),         NORM = 'F', 'f', 'E' or 'e'\n*\n*  where  norm1  denotes the  one norm of a matrix (maximum column sum),\n*  normI  denotes the  infinity norm  of a matrix  (maximum row sum) and\n*  normF  denotes the  Frobenius norm of a matrix (square root of sum of\n*  squares).  Note that  max(abs(A(i,j)))  is not a consistent matrix norm.\n*\n\n*  Arguments\n*  =========\n*\n*  NORM    (input) CHARACTER*1\n*          Specifies the value to be returned in ZLANTB as described\n*          above.\n*\n*  UPLO    (input) CHARACTER*1\n*          Specifies whether the matrix A is upper or lower triangular.\n*          = 'U':  Upper triangular\n*          = 'L':  Lower triangular\n*\n*  DIAG    (input) CHARACTER*1\n*          Specifies whether or not the matrix A is unit triangular.\n*          = 'N':  Non-unit triangular\n*          = 'U':  Unit triangular\n*\n*  N       (input) INTEGER\n*          The order of the matrix A.  N >= 0.  When N = 0, ZLANTB is\n*          set to zero.\n*\n*  K       (input) INTEGER\n*          The number of super-diagonals of the matrix A if UPLO = 'U',\n*          or the number of sub-diagonals of the matrix A if UPLO = 'L'.\n*          K >= 0.\n*\n*  AB      (input) COMPLEX*16 array, dimension (LDAB,N)\n*          The upper or lower triangular band matrix A, stored in the\n*          first k+1 rows of AB.  The j-th column of A is stored\n*          in the j-th column of the array AB as follows:\n*          if UPLO = 'U', AB(k+1+i-j,j) = A(i,j) for max(1,j-k)<=i<=j;\n*          if UPLO = 'L', AB(1+i-j,j)   = A(i,j) for j<=i<=min(n,j+k).\n*          Note that when DIAG = 'U', the elements of the array AB\n*          corresponding to the diagonal elements of the matrix A are\n*          not referenced, but are assumed to be one.\n*\n*  LDAB    (input) INTEGER\n*          The leading dimension of the array AB.  LDAB >= K+1.\n*\n*  WORK    (workspace) DOUBLE PRECISION array, dimension (MAX(1,LWORK)),\n*          where LWORK >= N when NORM = 'I'; otherwise, WORK is not\n*          referenced.\n*\n\n* =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 5)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 5)", argc);
  rb_norm = argv[0];
  rb_uplo = argv[1];
  rb_diag = argv[2];
  rb_k = argv[3];
  rb_ab = argv[4];

  norm = StringValueCStr(rb_norm)[0];
  uplo = StringValueCStr(rb_uplo)[0];
  diag = StringValueCStr(rb_diag)[0];
  k = NUM2INT(rb_k);
  if (!NA_IsNArray(rb_ab))
    rb_raise(rb_eArgError, "ab (5th argument) must be NArray");
  if (NA_RANK(rb_ab) != 2)
    rb_raise(rb_eArgError, "rank of ab (5th argument) must be %d", 2);
  ldab = NA_SHAPE0(rb_ab);
  n = NA_SHAPE1(rb_ab);
  if (NA_TYPE(rb_ab) != NA_DCOMPLEX)
    rb_ab = na_change_type(rb_ab, NA_DCOMPLEX);
  ab = NA_PTR_TYPE(rb_ab, doublecomplex*);
  lwork = lsame_(&norm,"I") ? n : 0;
  work = ALLOC_N(doublereal, (MAX(1,lwork)));

  zlantb_(&__out__, &norm, &uplo, &diag, &n, &k, ab, &ldab, work);

  free(work);
  rb___out__ = rb_float_new((double)__out__);
  return rb___out__;
}

void
init_lapack_zlantb(VALUE mLapack){
  rb_define_module_function(mLapack, "zlantb", rb_zlantb, -1);
}

#include "rb_lapack.h"

static VALUE
rb_zlantr(int argc, VALUE *argv, VALUE self){
  VALUE rb_norm;
  char norm; 
  VALUE rb_uplo;
  char uplo; 
  VALUE rb_diag;
  char diag; 
  VALUE rb_m;
  integer m; 
  VALUE rb_a;
  doublecomplex *a; 
  VALUE rb___out__;
  doublereal __out__; 
  doublereal *work;

  integer lda;
  integer n;
  integer lwork;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  __out__ = NumRu::Lapack.zlantr( norm, uplo, diag, m, a)\n    or\n  NumRu::Lapack.zlantr  # print help\n\n\nFORTRAN MANUAL\n      DOUBLE PRECISION FUNCTION ZLANTR( NORM, UPLO, DIAG, M, N, A, LDA, WORK )\n\n*  Purpose\n*  =======\n*\n*  ZLANTR  returns the value of the one norm,  or the Frobenius norm, or\n*  the  infinity norm,  or the  element of  largest absolute value  of a\n*  trapezoidal or triangular matrix A.\n*\n*  Description\n*  ===========\n*\n*  ZLANTR returns the value\n*\n*     ZLANTR = ( max(abs(A(i,j))), NORM = 'M' or 'm'\n*              (\n*              ( norm1(A),         NORM = '1', 'O' or 'o'\n*              (\n*              ( normI(A),         NORM = 'I' or 'i'\n*              (\n*              ( normF(A),         NORM = 'F', 'f', 'E' or 'e'\n*\n*  where  norm1  denotes the  one norm of a matrix (maximum column sum),\n*  normI  denotes the  infinity norm  of a matrix  (maximum row sum) and\n*  normF  denotes the  Frobenius norm of a matrix (square root of sum of\n*  squares).  Note that  max(abs(A(i,j)))  is not a consistent matrix norm.\n*\n\n*  Arguments\n*  =========\n*\n*  NORM    (input) CHARACTER*1\n*          Specifies the value to be returned in ZLANTR as described\n*          above.\n*\n*  UPLO    (input) CHARACTER*1\n*          Specifies whether the matrix A is upper or lower trapezoidal.\n*          = 'U':  Upper trapezoidal\n*          = 'L':  Lower trapezoidal\n*          Note that A is triangular instead of trapezoidal if M = N.\n*\n*  DIAG    (input) CHARACTER*1\n*          Specifies whether or not the matrix A has unit diagonal.\n*          = 'N':  Non-unit diagonal\n*          = 'U':  Unit diagonal\n*\n*  M       (input) INTEGER\n*          The number of rows of the matrix A.  M >= 0, and if\n*          UPLO = 'U', M <= N.  When M = 0, ZLANTR is set to zero.\n*\n*  N       (input) INTEGER\n*          The number of columns of the matrix A.  N >= 0, and if\n*          UPLO = 'L', N <= M.  When N = 0, ZLANTR is set to zero.\n*\n*  A       (input) COMPLEX*16 array, dimension (LDA,N)\n*          The trapezoidal matrix A (A is triangular if M = N).\n*          If UPLO = 'U', the leading m by n upper trapezoidal part of\n*          the array A contains the upper trapezoidal matrix, and the\n*          strictly lower triangular part of A is not referenced.\n*          If UPLO = 'L', the leading m by n lower trapezoidal part of\n*          the array A contains the lower trapezoidal matrix, and the\n*          strictly upper triangular part of A is not referenced.  Note\n*          that when DIAG = 'U', the diagonal elements of A are not\n*          referenced and are assumed to be one.\n*\n*  LDA     (input) INTEGER\n*          The leading dimension of the array A.  LDA >= max(M,1).\n*\n*  WORK    (workspace) DOUBLE PRECISION array, dimension (MAX(1,LWORK)),\n*          where LWORK >= M when NORM = 'I'; otherwise, WORK is not\n*          referenced.\n*\n\n* =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 5)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 5)", argc);
  rb_norm = argv[0];
  rb_uplo = argv[1];
  rb_diag = argv[2];
  rb_m = argv[3];
  rb_a = argv[4];

  norm = StringValueCStr(rb_norm)[0];
  uplo = StringValueCStr(rb_uplo)[0];
  diag = StringValueCStr(rb_diag)[0];
  m = NUM2INT(rb_m);
  if (!NA_IsNArray(rb_a))
    rb_raise(rb_eArgError, "a (5th argument) must be NArray");
  if (NA_RANK(rb_a) != 2)
    rb_raise(rb_eArgError, "rank of a (5th argument) must be %d", 2);
  lda = NA_SHAPE0(rb_a);
  n = NA_SHAPE1(rb_a);
  if (NA_TYPE(rb_a) != NA_DCOMPLEX)
    rb_a = na_change_type(rb_a, NA_DCOMPLEX);
  a = NA_PTR_TYPE(rb_a, doublecomplex*);
  lwork = lsame_(&norm,"I") ? m : 0;
  work = ALLOC_N(doublereal, (MAX(1,lwork)));

  __out__ = zlantr_(&norm, &uplo, &diag, &m, &n, a, &lda, work);

  free(work);
  rb___out__ = rb_float_new((double)__out__);
  return rb___out__;
}

void
init_lapack_zlantr(VALUE mLapack){
  rb_define_module_function(mLapack, "zlantr", rb_zlantr, -1);
}

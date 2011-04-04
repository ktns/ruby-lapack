#include "rb_lapack.h"

extern doublereal dlantp_(char *norm, char *uplo, char *diag, integer *n, doublereal *ap, doublereal *work);

static VALUE
rb_dlantp(int argc, VALUE *argv, VALUE self){
  VALUE rb_norm;
  char norm; 
  VALUE rb_uplo;
  char uplo; 
  VALUE rb_diag;
  char diag; 
  VALUE rb_n;
  integer n; 
  VALUE rb_ap;
  doublereal *ap; 
  VALUE rb___out__;
  doublereal __out__; 
  doublereal *work;

  integer lwork;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  __out__ = NumRu::Lapack.dlantp( norm, uplo, diag, n, ap)\n    or\n  NumRu::Lapack.dlantp  # print help\n\n\nFORTRAN MANUAL\n      DOUBLE PRECISION FUNCTION DLANTP( NORM, UPLO, DIAG, N, AP, WORK )\n\n*  Purpose\n*  =======\n*\n*  DLANTP  returns the value of the one norm,  or the Frobenius norm, or\n*  the  infinity norm,  or the  element of  largest absolute value  of a\n*  triangular matrix A, supplied in packed form.\n*\n*  Description\n*  ===========\n*\n*  DLANTP returns the value\n*\n*     DLANTP = ( max(abs(A(i,j))), NORM = 'M' or 'm'\n*              (\n*              ( norm1(A),         NORM = '1', 'O' or 'o'\n*              (\n*              ( normI(A),         NORM = 'I' or 'i'\n*              (\n*              ( normF(A),         NORM = 'F', 'f', 'E' or 'e'\n*\n*  where  norm1  denotes the  one norm of a matrix (maximum column sum),\n*  normI  denotes the  infinity norm  of a matrix  (maximum row sum) and\n*  normF  denotes the  Frobenius norm of a matrix (square root of sum of\n*  squares).  Note that  max(abs(A(i,j)))  is not a consistent matrix norm.\n*\n\n*  Arguments\n*  =========\n*\n*  NORM    (input) CHARACTER*1\n*          Specifies the value to be returned in DLANTP as described\n*          above.\n*\n*  UPLO    (input) CHARACTER*1\n*          Specifies whether the matrix A is upper or lower triangular.\n*          = 'U':  Upper triangular\n*          = 'L':  Lower triangular\n*\n*  DIAG    (input) CHARACTER*1\n*          Specifies whether or not the matrix A is unit triangular.\n*          = 'N':  Non-unit triangular\n*          = 'U':  Unit triangular\n*\n*  N       (input) INTEGER\n*          The order of the matrix A.  N >= 0.  When N = 0, DLANTP is\n*          set to zero.\n*\n*  AP      (input) DOUBLE PRECISION array, dimension (N*(N+1)/2)\n*          The upper or lower triangular matrix A, packed columnwise in\n*          a linear array.  The j-th column of A is stored in the array\n*          AP as follows:\n*          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;\n*          if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = A(i,j) for j<=i<=n.\n*          Note that when DIAG = 'U', the elements of the array AP\n*          corresponding to the diagonal elements of the matrix A are\n*          not referenced, but are assumed to be one.\n*\n*  WORK    (workspace) DOUBLE PRECISION array, dimension (MAX(1,LWORK)),\n*          where LWORK >= N when NORM = 'I'; otherwise, WORK is not\n*          referenced.\n*\n\n* =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 5)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 5)", argc);
  rb_norm = argv[0];
  rb_uplo = argv[1];
  rb_diag = argv[2];
  rb_n = argv[3];
  rb_ap = argv[4];

  diag = StringValueCStr(rb_diag)[0];
  n = NUM2INT(rb_n);
  norm = StringValueCStr(rb_norm)[0];
  uplo = StringValueCStr(rb_uplo)[0];
  if (!NA_IsNArray(rb_ap))
    rb_raise(rb_eArgError, "ap (5th argument) must be NArray");
  if (NA_RANK(rb_ap) != 1)
    rb_raise(rb_eArgError, "rank of ap (5th argument) must be %d", 1);
  if (NA_SHAPE0(rb_ap) != (n*(n+1)/2))
    rb_raise(rb_eRuntimeError, "shape 0 of ap must be %d", n*(n+1)/2);
  if (NA_TYPE(rb_ap) != NA_DFLOAT)
    rb_ap = na_change_type(rb_ap, NA_DFLOAT);
  ap = NA_PTR_TYPE(rb_ap, doublereal*);
  lwork = lsame_(&norm,"I") ? n : 0;
  work = ALLOC_N(doublereal, (MAX(1,lwork)));

  __out__ = dlantp_(&norm, &uplo, &diag, &n, ap, work);

  free(work);
  rb___out__ = rb_float_new((double)__out__);
  return rb___out__;
}

void
init_lapack_dlantp(VALUE mLapack){
  rb_define_module_function(mLapack, "dlantp", rb_dlantp, -1);
}

#include "rb_lapack.h"

extern real clanhp_(char *norm, char *uplo, integer *n, complex *ap, real *work);

static VALUE
rb_clanhp(int argc, VALUE *argv, VALUE self){
  VALUE rb_norm;
  char norm; 
  VALUE rb_uplo;
  char uplo; 
  VALUE rb_n;
  integer n; 
  VALUE rb_ap;
  complex *ap; 
  VALUE rb___out__;
  real __out__; 
  real *work;


  if (argc == 0) {
    printf("%s\n", "USAGE:\n  __out__ = NumRu::Lapack.clanhp( norm, uplo, n, ap)\n    or\n  NumRu::Lapack.clanhp  # print help\n\n\nFORTRAN MANUAL\n      REAL             FUNCTION CLANHP( NORM, UPLO, N, AP, WORK )\n\n*  Purpose\n*  =======\n*\n*  CLANHP  returns the value of the one norm,  or the Frobenius norm, or\n*  the  infinity norm,  or the  element of  largest absolute value  of a\n*  complex hermitian matrix A,  supplied in packed form.\n*\n*  Description\n*  ===========\n*\n*  CLANHP returns the value\n*\n*     CLANHP = ( max(abs(A(i,j))), NORM = 'M' or 'm'\n*              (\n*              ( norm1(A),         NORM = '1', 'O' or 'o'\n*              (\n*              ( normI(A),         NORM = 'I' or 'i'\n*              (\n*              ( normF(A),         NORM = 'F', 'f', 'E' or 'e'\n*\n*  where  norm1  denotes the  one norm of a matrix (maximum column sum),\n*  normI  denotes the  infinity norm  of a matrix  (maximum row sum) and\n*  normF  denotes the  Frobenius norm of a matrix (square root of sum of\n*  squares).  Note that  max(abs(A(i,j)))  is not a consistent matrix norm.\n*\n\n*  Arguments\n*  =========\n*\n*  NORM    (input) CHARACTER*1\n*          Specifies the value to be returned in CLANHP as described\n*          above.\n*\n*  UPLO    (input) CHARACTER*1\n*          Specifies whether the upper or lower triangular part of the\n*          hermitian matrix A is supplied.\n*          = 'U':  Upper triangular part of A is supplied\n*          = 'L':  Lower triangular part of A is supplied\n*\n*  N       (input) INTEGER\n*          The order of the matrix A.  N >= 0.  When N = 0, CLANHP is\n*          set to zero.\n*\n*  AP      (input) COMPLEX array, dimension (N*(N+1)/2)\n*          The upper or lower triangle of the hermitian matrix A, packed\n*          columnwise in a linear array.  The j-th column of A is stored\n*          in the array AP as follows:\n*          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;\n*          if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = A(i,j) for j<=i<=n.\n*          Note that the  imaginary parts of the diagonal elements need\n*          not be set and are assumed to be zero.\n*\n*  WORK    (workspace) REAL array, dimension (MAX(1,LWORK)),\n*          where LWORK >= N when NORM = 'I' or '1' or 'O'; otherwise,\n*          WORK is not referenced.\n*\n\n* =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 4)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 4)", argc);
  rb_norm = argv[0];
  rb_uplo = argv[1];
  rb_n = argv[2];
  rb_ap = argv[3];

  n = NUM2INT(rb_n);
  uplo = StringValueCStr(rb_uplo)[0];
  norm = StringValueCStr(rb_norm)[0];
  if (!NA_IsNArray(rb_ap))
    rb_raise(rb_eArgError, "ap (4th argument) must be NArray");
  if (NA_RANK(rb_ap) != 1)
    rb_raise(rb_eArgError, "rank of ap (4th argument) must be %d", 1);
  if (NA_SHAPE0(rb_ap) != (n*(n+1)/2))
    rb_raise(rb_eRuntimeError, "shape 0 of ap must be %d", n*(n+1)/2);
  if (NA_TYPE(rb_ap) != NA_SCOMPLEX)
    rb_ap = na_change_type(rb_ap, NA_SCOMPLEX);
  ap = NA_PTR_TYPE(rb_ap, complex*);
  work = ALLOC_N(real, (MAX(1,(lsame_(&norm,"I")||lsame_(&norm,"1")||lsame_(&norm,"O")) ? n : 0)));

  __out__ = clanhp_(&norm, &uplo, &n, ap, work);

  free(work);
  rb___out__ = rb_float_new((double)__out__);
  return rb___out__;
}

void
init_lapack_clanhp(VALUE mLapack){
  rb_define_module_function(mLapack, "clanhp", rb_clanhp, -1);
}

#include "rb_lapack.h"

extern doublereal dlansp_(char *norm, char *uplo, integer *n, doublereal *ap, doublereal *work);

static VALUE
rb_dlansp(int argc, VALUE *argv, VALUE self){
  VALUE rb_norm;
  char norm; 
  VALUE rb_uplo;
  char uplo; 
  VALUE rb_n;
  integer n; 
  VALUE rb_ap;
  doublereal *ap; 
  VALUE rb___out__;
  doublereal __out__; 
  doublereal *work;

  integer lwork;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  __out__ = NumRu::Lapack.dlansp( norm, uplo, n, ap)\n    or\n  NumRu::Lapack.dlansp  # print help\n\n\nFORTRAN MANUAL\n      DOUBLE PRECISION FUNCTION DLANSP( NORM, UPLO, N, AP, WORK )\n\n*  Purpose\n*  =======\n*\n*  DLANSP  returns the value of the one norm,  or the Frobenius norm, or\n*  the  infinity norm,  or the  element of  largest absolute value  of a\n*  real symmetric matrix A,  supplied in packed form.\n*\n*  Description\n*  ===========\n*\n*  DLANSP returns the value\n*\n*     DLANSP = ( max(abs(A(i,j))), NORM = 'M' or 'm'\n*              (\n*              ( norm1(A),         NORM = '1', 'O' or 'o'\n*              (\n*              ( normI(A),         NORM = 'I' or 'i'\n*              (\n*              ( normF(A),         NORM = 'F', 'f', 'E' or 'e'\n*\n*  where  norm1  denotes the  one norm of a matrix (maximum column sum),\n*  normI  denotes the  infinity norm  of a matrix  (maximum row sum) and\n*  normF  denotes the  Frobenius norm of a matrix (square root of sum of\n*  squares).  Note that  max(abs(A(i,j)))  is not a consistent matrix norm.\n*\n\n*  Arguments\n*  =========\n*\n*  NORM    (input) CHARACTER*1\n*          Specifies the value to be returned in DLANSP as described\n*          above.\n*\n*  UPLO    (input) CHARACTER*1\n*          Specifies whether the upper or lower triangular part of the\n*          symmetric matrix A is supplied.\n*          = 'U':  Upper triangular part of A is supplied\n*          = 'L':  Lower triangular part of A is supplied\n*\n*  N       (input) INTEGER\n*          The order of the matrix A.  N >= 0.  When N = 0, DLANSP is\n*          set to zero.\n*\n*  AP      (input) DOUBLE PRECISION array, dimension (N*(N+1)/2)\n*          The upper or lower triangle of the symmetric matrix A, packed\n*          columnwise in a linear array.  The j-th column of A is stored\n*          in the array AP as follows:\n*          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;\n*          if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = A(i,j) for j<=i<=n.\n*\n*  WORK    (workspace) DOUBLE PRECISION array, dimension (MAX(1,LWORK)),\n*          where LWORK >= N when NORM = 'I' or '1' or 'O'; otherwise,\n*          WORK is not referenced.\n*\n\n* =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 4)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 4)", argc);
  rb_norm = argv[0];
  rb_uplo = argv[1];
  rb_n = argv[2];
  rb_ap = argv[3];

  n = NUM2INT(rb_n);
  norm = StringValueCStr(rb_norm)[0];
  uplo = StringValueCStr(rb_uplo)[0];
  if (!NA_IsNArray(rb_ap))
    rb_raise(rb_eArgError, "ap (4th argument) must be NArray");
  if (NA_RANK(rb_ap) != 1)
    rb_raise(rb_eArgError, "rank of ap (4th argument) must be %d", 1);
  if (NA_SHAPE0(rb_ap) != (n*(n+1)/2))
    rb_raise(rb_eRuntimeError, "shape 0 of ap must be %d", n*(n+1)/2);
  if (NA_TYPE(rb_ap) != NA_DFLOAT)
    rb_ap = na_change_type(rb_ap, NA_DFLOAT);
  ap = NA_PTR_TYPE(rb_ap, doublereal*);
  lwork = ((lsame_(&norm,"I")) || ((('1') || ('o')))) ? n : 0;
  work = ALLOC_N(doublereal, (MAX(1,lwork)));

  __out__ = dlansp_(&norm, &uplo, &n, ap, work);

  free(work);
  rb___out__ = rb_float_new((double)__out__);
  return rb___out__;
}

void
init_lapack_dlansp(VALUE mLapack){
  rb_define_module_function(mLapack, "dlansp", rb_dlansp, -1);
}

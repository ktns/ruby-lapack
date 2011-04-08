#include "rb_lapack.h"

extern doublereal zlantp_(char *norm, char *uplo, char *diag, integer *n, doublecomplex *ap, doublereal *work);

static VALUE sHelp, sUsage;

static VALUE
rblapack_zlantp(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_norm;
  char norm; 
  VALUE rblapack_uplo;
  char uplo; 
  VALUE rblapack_diag;
  char diag; 
  VALUE rblapack_n;
  integer n; 
  VALUE rblapack_ap;
  doublecomplex *ap; 
  VALUE rblapack___out__;
  doublereal __out__; 
  doublereal *work;

  integer lwork;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  __out__ = NumRu::Lapack.zlantp( norm, uplo, diag, n, ap, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      DOUBLE PRECISION FUNCTION ZLANTP( NORM, UPLO, DIAG, N, AP, WORK )\n\n*  Purpose\n*  =======\n*\n*  ZLANTP  returns the value of the one norm,  or the Frobenius norm, or\n*  the  infinity norm,  or the  element of  largest absolute value  of a\n*  triangular matrix A, supplied in packed form.\n*\n*  Description\n*  ===========\n*\n*  ZLANTP returns the value\n*\n*     ZLANTP = ( max(abs(A(i,j))), NORM = 'M' or 'm'\n*              (\n*              ( norm1(A),         NORM = '1', 'O' or 'o'\n*              (\n*              ( normI(A),         NORM = 'I' or 'i'\n*              (\n*              ( normF(A),         NORM = 'F', 'f', 'E' or 'e'\n*\n*  where  norm1  denotes the  one norm of a matrix (maximum column sum),\n*  normI  denotes the  infinity norm  of a matrix  (maximum row sum) and\n*  normF  denotes the  Frobenius norm of a matrix (square root of sum of\n*  squares).  Note that  max(abs(A(i,j)))  is not a consistent matrix norm.\n*\n\n*  Arguments\n*  =========\n*\n*  NORM    (input) CHARACTER*1\n*          Specifies the value to be returned in ZLANTP as described\n*          above.\n*\n*  UPLO    (input) CHARACTER*1\n*          Specifies whether the matrix A is upper or lower triangular.\n*          = 'U':  Upper triangular\n*          = 'L':  Lower triangular\n*\n*  DIAG    (input) CHARACTER*1\n*          Specifies whether or not the matrix A is unit triangular.\n*          = 'N':  Non-unit triangular\n*          = 'U':  Unit triangular\n*\n*  N       (input) INTEGER\n*          The order of the matrix A.  N >= 0.  When N = 0, ZLANTP is\n*          set to zero.\n*\n*  AP      (input) COMPLEX*16 array, dimension (N*(N+1)/2)\n*          The upper or lower triangular matrix A, packed columnwise in\n*          a linear array.  The j-th column of A is stored in the array\n*          AP as follows:\n*          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;\n*          if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = A(i,j) for j<=i<=n.\n*          Note that when DIAG = 'U', the elements of the array AP\n*          corresponding to the diagonal elements of the matrix A are\n*          not referenced, but are assumed to be one.\n*\n*  WORK    (workspace) DOUBLE PRECISION array, dimension (MAX(1,LWORK)),\n*          where LWORK >= N when NORM = 'I'; otherwise, WORK is not\n*          referenced.\n*\n\n* =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  __out__ = NumRu::Lapack.zlantp( norm, uplo, diag, n, ap, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 5)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 5)", argc);
  rblapack_norm = argv[0];
  rblapack_uplo = argv[1];
  rblapack_diag = argv[2];
  rblapack_n = argv[3];
  rblapack_ap = argv[4];
  if (rb_options != Qnil) {
  }

  diag = StringValueCStr(rblapack_diag)[0];
  n = NUM2INT(rblapack_n);
  norm = StringValueCStr(rblapack_norm)[0];
  uplo = StringValueCStr(rblapack_uplo)[0];
  if (!NA_IsNArray(rblapack_ap))
    rb_raise(rb_eArgError, "ap (5th argument) must be NArray");
  if (NA_RANK(rblapack_ap) != 1)
    rb_raise(rb_eArgError, "rank of ap (5th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_ap) != (n*(n+1)/2))
    rb_raise(rb_eRuntimeError, "shape 0 of ap must be %d", n*(n+1)/2);
  if (NA_TYPE(rblapack_ap) != NA_DCOMPLEX)
    rblapack_ap = na_change_type(rblapack_ap, NA_DCOMPLEX);
  ap = NA_PTR_TYPE(rblapack_ap, doublecomplex*);
  lwork = lsame_(&norm,"I") ? n : 0;
  work = ALLOC_N(doublereal, (MAX(1,lwork)));

  __out__ = zlantp_(&norm, &uplo, &diag, &n, ap, work);

  free(work);
  rblapack___out__ = rb_float_new((double)__out__);
  return rblapack___out__;
}

void
init_lapack_zlantp(VALUE mLapack){
  rb_define_module_function(mLapack, "zlantp", rblapack_zlantp, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}

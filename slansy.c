#include "rb_lapack.h"

extern real slansy_(char *norm, char *uplo, integer *n, real *a, integer *lda, real *work);

static VALUE sHelp, sUsage;

static VALUE
rblapack_slansy(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_norm;
  char norm; 
  VALUE rblapack_uplo;
  char uplo; 
  VALUE rblapack_a;
  real *a; 
  VALUE rblapack___out__;
  real __out__; 
  real *work;

  integer lda;
  integer n;
  integer lwork;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  __out__ = NumRu::Lapack.slansy( norm, uplo, a, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      REAL             FUNCTION SLANSY( NORM, UPLO, N, A, LDA, WORK )\n\n*  Purpose\n*  =======\n*\n*  SLANSY  returns the value of the one norm,  or the Frobenius norm, or\n*  the  infinity norm,  or the  element of  largest absolute value  of a\n*  real symmetric matrix A.\n*\n*  Description\n*  ===========\n*\n*  SLANSY returns the value\n*\n*     SLANSY = ( max(abs(A(i,j))), NORM = 'M' or 'm'\n*              (\n*              ( norm1(A),         NORM = '1', 'O' or 'o'\n*              (\n*              ( normI(A),         NORM = 'I' or 'i'\n*              (\n*              ( normF(A),         NORM = 'F', 'f', 'E' or 'e'\n*\n*  where  norm1  denotes the  one norm of a matrix (maximum column sum),\n*  normI  denotes the  infinity norm  of a matrix  (maximum row sum) and\n*  normF  denotes the  Frobenius norm of a matrix (square root of sum of\n*  squares).  Note that  max(abs(A(i,j)))  is not a consistent matrix norm.\n*\n\n*  Arguments\n*  =========\n*\n*  NORM    (input) CHARACTER*1\n*          Specifies the value to be returned in SLANSY as described\n*          above.\n*\n*  UPLO    (input) CHARACTER*1\n*          Specifies whether the upper or lower triangular part of the\n*          symmetric matrix A is to be referenced.\n*          = 'U':  Upper triangular part of A is referenced\n*          = 'L':  Lower triangular part of A is referenced\n*\n*  N       (input) INTEGER\n*          The order of the matrix A.  N >= 0.  When N = 0, SLANSY is\n*          set to zero.\n*\n*  A       (input) REAL array, dimension (LDA,N)\n*          The symmetric matrix A.  If UPLO = 'U', the leading n by n\n*          upper triangular part of A contains the upper triangular part\n*          of the matrix A, and the strictly lower triangular part of A\n*          is not referenced.  If UPLO = 'L', the leading n by n lower\n*          triangular part of A contains the lower triangular part of\n*          the matrix A, and the strictly upper triangular part of A is\n*          not referenced.\n*\n*  LDA     (input) INTEGER\n*          The leading dimension of the array A.  LDA >= max(N,1).\n*\n*  WORK    (workspace) REAL array, dimension (MAX(1,LWORK)),\n*          where LWORK >= N when NORM = 'I' or '1' or 'O'; otherwise,\n*          WORK is not referenced.\n*\n\n* =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  __out__ = NumRu::Lapack.slansy( norm, uplo, a, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 3)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 3)", argc);
  rblapack_norm = argv[0];
  rblapack_uplo = argv[1];
  rblapack_a = argv[2];
  if (rb_options != Qnil) {
  }

  if (!NA_IsNArray(rblapack_a))
    rb_raise(rb_eArgError, "a (3th argument) must be NArray");
  if (NA_RANK(rblapack_a) != 2)
    rb_raise(rb_eArgError, "rank of a (3th argument) must be %d", 2);
  n = NA_SHAPE1(rblapack_a);
  lda = NA_SHAPE0(rblapack_a);
  if (NA_TYPE(rblapack_a) != NA_SFLOAT)
    rblapack_a = na_change_type(rblapack_a, NA_SFLOAT);
  a = NA_PTR_TYPE(rblapack_a, real*);
  uplo = StringValueCStr(rblapack_uplo)[0];
  norm = StringValueCStr(rblapack_norm)[0];
  lwork = ((lsame_(&norm,"I")) || ((('1') || ('o')))) ? n : 0;
  work = ALLOC_N(real, (MAX(1,lwork)));

  __out__ = slansy_(&norm, &uplo, &n, a, &lda, work);

  free(work);
  rblapack___out__ = rb_float_new((double)__out__);
  return rblapack___out__;
}

void
init_lapack_slansy(VALUE mLapack){
  rb_define_module_function(mLapack, "slansy", rblapack_slansy, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}

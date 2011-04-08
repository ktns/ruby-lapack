#include "rb_lapack.h"

extern doublereal dlantr_(char *norm, char *uplo, char *diag, integer *m, integer *n, doublereal *a, integer *lda, doublereal *work);

static VALUE sHelp, sUsage;

static VALUE
rblapack_dlantr(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_norm;
  char norm; 
  VALUE rblapack_uplo;
  char uplo; 
  VALUE rblapack_diag;
  char diag; 
  VALUE rblapack_m;
  integer m; 
  VALUE rblapack_a;
  doublereal *a; 
  VALUE rblapack___out__;
  doublereal __out__; 
  doublereal *work;

  integer lda;
  integer n;
  integer lwork;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  __out__ = NumRu::Lapack.dlantr( norm, uplo, diag, m, a, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      DOUBLE PRECISION FUNCTION DLANTR( NORM, UPLO, DIAG, M, N, A, LDA, WORK )\n\n*  Purpose\n*  =======\n*\n*  DLANTR  returns the value of the one norm,  or the Frobenius norm, or\n*  the  infinity norm,  or the  element of  largest absolute value  of a\n*  trapezoidal or triangular matrix A.\n*\n*  Description\n*  ===========\n*\n*  DLANTR returns the value\n*\n*     DLANTR = ( max(abs(A(i,j))), NORM = 'M' or 'm'\n*              (\n*              ( norm1(A),         NORM = '1', 'O' or 'o'\n*              (\n*              ( normI(A),         NORM = 'I' or 'i'\n*              (\n*              ( normF(A),         NORM = 'F', 'f', 'E' or 'e'\n*\n*  where  norm1  denotes the  one norm of a matrix (maximum column sum),\n*  normI  denotes the  infinity norm  of a matrix  (maximum row sum) and\n*  normF  denotes the  Frobenius norm of a matrix (square root of sum of\n*  squares).  Note that  max(abs(A(i,j)))  is not a consistent matrix norm.\n*\n\n*  Arguments\n*  =========\n*\n*  NORM    (input) CHARACTER*1\n*          Specifies the value to be returned in DLANTR as described\n*          above.\n*\n*  UPLO    (input) CHARACTER*1\n*          Specifies whether the matrix A is upper or lower trapezoidal.\n*          = 'U':  Upper trapezoidal\n*          = 'L':  Lower trapezoidal\n*          Note that A is triangular instead of trapezoidal if M = N.\n*\n*  DIAG    (input) CHARACTER*1\n*          Specifies whether or not the matrix A has unit diagonal.\n*          = 'N':  Non-unit diagonal\n*          = 'U':  Unit diagonal\n*\n*  M       (input) INTEGER\n*          The number of rows of the matrix A.  M >= 0, and if\n*          UPLO = 'U', M <= N.  When M = 0, DLANTR is set to zero.\n*\n*  N       (input) INTEGER\n*          The number of columns of the matrix A.  N >= 0, and if\n*          UPLO = 'L', N <= M.  When N = 0, DLANTR is set to zero.\n*\n*  A       (input) DOUBLE PRECISION array, dimension (LDA,N)\n*          The trapezoidal matrix A (A is triangular if M = N).\n*          If UPLO = 'U', the leading m by n upper trapezoidal part of\n*          the array A contains the upper trapezoidal matrix, and the\n*          strictly lower triangular part of A is not referenced.\n*          If UPLO = 'L', the leading m by n lower trapezoidal part of\n*          the array A contains the lower trapezoidal matrix, and the\n*          strictly upper triangular part of A is not referenced.  Note\n*          that when DIAG = 'U', the diagonal elements of A are not\n*          referenced and are assumed to be one.\n*\n*  LDA     (input) INTEGER\n*          The leading dimension of the array A.  LDA >= max(M,1).\n*\n*  WORK    (workspace) DOUBLE PRECISION array, dimension (MAX(1,LWORK)),\n*          where LWORK >= M when NORM = 'I'; otherwise, WORK is not\n*          referenced.\n*\n\n* =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  __out__ = NumRu::Lapack.dlantr( norm, uplo, diag, m, a, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 5)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 5)", argc);
  rblapack_norm = argv[0];
  rblapack_uplo = argv[1];
  rblapack_diag = argv[2];
  rblapack_m = argv[3];
  rblapack_a = argv[4];
  if (rb_options != Qnil) {
  }

  if (!NA_IsNArray(rblapack_a))
    rb_raise(rb_eArgError, "a (5th argument) must be NArray");
  if (NA_RANK(rblapack_a) != 2)
    rb_raise(rb_eArgError, "rank of a (5th argument) must be %d", 2);
  n = NA_SHAPE1(rblapack_a);
  lda = NA_SHAPE0(rblapack_a);
  if (NA_TYPE(rblapack_a) != NA_DFLOAT)
    rblapack_a = na_change_type(rblapack_a, NA_DFLOAT);
  a = NA_PTR_TYPE(rblapack_a, doublereal*);
  m = NUM2INT(rblapack_m);
  diag = StringValueCStr(rblapack_diag)[0];
  norm = StringValueCStr(rblapack_norm)[0];
  uplo = StringValueCStr(rblapack_uplo)[0];
  lwork = lsame_(&norm,"I") ? m : 0;
  work = ALLOC_N(doublereal, (MAX(1,lwork)));

  __out__ = dlantr_(&norm, &uplo, &diag, &m, &n, a, &lda, work);

  free(work);
  rblapack___out__ = rb_float_new((double)__out__);
  return rblapack___out__;
}

void
init_lapack_dlantr(VALUE mLapack){
  rb_define_module_function(mLapack, "dlantr", rblapack_dlantr, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}

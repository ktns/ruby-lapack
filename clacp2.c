#include "rb_lapack.h"

extern VOID clacp2_(char *uplo, integer *m, integer *n, real *a, integer *lda, complex *b, integer *ldb);

static VALUE sHelp, sUsage;

static VALUE
rblapack_clacp2(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_uplo;
  char uplo; 
  VALUE rblapack_m;
  integer m; 
  VALUE rblapack_a;
  real *a; 
  VALUE rblapack_b;
  complex *b; 

  integer lda;
  integer n;
  integer ldb;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  b = NumRu::Lapack.clacp2( uplo, m, a, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE CLACP2( UPLO, M, N, A, LDA, B, LDB )\n\n*  Purpose\n*  =======\n*\n*  CLACP2 copies all or part of a real two-dimensional matrix A to a\n*  complex matrix B.\n*\n\n*  Arguments\n*  =========\n*\n*  UPLO    (input) CHARACTER*1\n*          Specifies the part of the matrix A to be copied to B.\n*          = 'U':      Upper triangular part\n*          = 'L':      Lower triangular part\n*          Otherwise:  All of the matrix A\n*\n*  M       (input) INTEGER\n*          The number of rows of the matrix A.  M >= 0.\n*\n*  N       (input) INTEGER\n*          The number of columns of the matrix A.  N >= 0.\n*\n*  A       (input) REAL array, dimension (LDA,N)\n*          The m by n matrix A.  If UPLO = 'U', only the upper trapezium\n*          is accessed; if UPLO = 'L', only the lower trapezium is\n*          accessed.\n*\n*  LDA     (input) INTEGER\n*          The leading dimension of the array A.  LDA >= max(1,M).\n*\n*  B       (output) COMPLEX array, dimension (LDB,N)\n*          On exit, B = A in the locations specified by UPLO.\n*\n*  LDB     (input) INTEGER\n*          The leading dimension of the array B.  LDB >= max(1,M).\n*\n\n*  =====================================================================\n*\n*     .. Local Scalars ..\n      INTEGER            I, J\n*     ..\n*     .. External Functions ..\n      LOGICAL            LSAME\n      EXTERNAL           LSAME\n*     ..\n*     .. Intrinsic Functions ..\n      INTRINSIC          MIN\n*     ..\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  b = NumRu::Lapack.clacp2( uplo, m, a, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 3)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 3)", argc);
  rblapack_uplo = argv[0];
  rblapack_m = argv[1];
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
  m = NUM2INT(rblapack_m);
  uplo = StringValueCStr(rblapack_uplo)[0];
  ldb = MAX(1,m);
  {
    int shape[2];
    shape[0] = ldb;
    shape[1] = n;
    rblapack_b = na_make_object(NA_SCOMPLEX, 2, shape, cNArray);
  }
  b = NA_PTR_TYPE(rblapack_b, complex*);

  clacp2_(&uplo, &m, &n, a, &lda, b, &ldb);

  return rblapack_b;
}

void
init_lapack_clacp2(VALUE mLapack){
  rb_define_module_function(mLapack, "clacp2", rblapack_clacp2, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}

#include "rb_lapack.h"

extern VOID zlat2c_(char *uplo, integer *n, doublecomplex *a, integer *lda, complex *sa, integer *ldsa, integer *info);

static VALUE sHelp, sUsage;

static VALUE
rblapack_zlat2c(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_uplo;
  char uplo; 
  VALUE rblapack_a;
  doublecomplex *a; 
  VALUE rblapack_sa;
  complex *sa; 
  VALUE rblapack_info;
  integer info; 

  integer lda;
  integer n;
  integer ldsa;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  sa, info = NumRu::Lapack.zlat2c( uplo, a, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE ZLAT2C( UPLO, N, A, LDA, SA, LDSA, INFO )\n\n*  Purpose\n*  =======\n*\n*  ZLAT2C converts a COMPLEX*16 triangular matrix, SA, to a COMPLEX\n*  triangular matrix, A.\n*\n*  RMAX is the overflow for the SINGLE PRECISION arithmetic\n*  ZLAT2C checks that all the entries of A are between -RMAX and\n*  RMAX. If not the convertion is aborted and a flag is raised.\n*\n*  This is an auxiliary routine so there is no argument checking.\n*\n\n*  Arguments\n*  =========\n*\n*  UPLO    (input) CHARACTER*1\n*          = 'U':  A is upper triangular;\n*          = 'L':  A is lower triangular.\n*\n*  N       (input) INTEGER\n*          The number of rows and columns of the matrix A.  N >= 0.\n*\n*  A       (input) COMPLEX*16 array, dimension (LDA,N)\n*          On entry, the N-by-N triangular coefficient matrix A.\n*\n*  LDA     (input) INTEGER\n*          The leading dimension of the array A.  LDA >= max(1,N).\n*\n*  SA      (output) COMPLEX array, dimension (LDSA,N)\n*          Only the UPLO part of SA is referenced.  On exit, if INFO=0,\n*          the N-by-N coefficient matrix SA; if INFO>0, the content of\n*          the UPLO part of SA is unspecified.\n*\n*  LDSA    (input) INTEGER\n*          The leading dimension of the array SA.  LDSA >= max(1,M).\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit.\n*          = 1:  an entry of the matrix A is greater than the SINGLE\n*                PRECISION overflow threshold, in this case, the content\n*                of the UPLO part of SA in exit is unspecified.\n*\n*  =========\n*\n*     .. Local Scalars ..\n      INTEGER            I, J\n      DOUBLE PRECISION   RMAX\n      LOGICAL            UPPER\n*     ..\n*     .. Intrinsic Functions ..\n      INTRINSIC          DBLE, DIMAG\n*     ..\n*     .. External Functions ..\n      REAL               SLAMCH\n      LOGICAL            LSAME\n      EXTERNAL           SLAMCH, LSAME\n*     ..\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  sa, info = NumRu::Lapack.zlat2c( uplo, a, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 2)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 2)", argc);
  rblapack_uplo = argv[0];
  rblapack_a = argv[1];
  if (rb_options != Qnil) {
  }

  if (!NA_IsNArray(rblapack_a))
    rb_raise(rb_eArgError, "a (2th argument) must be NArray");
  if (NA_RANK(rblapack_a) != 2)
    rb_raise(rb_eArgError, "rank of a (2th argument) must be %d", 2);
  n = NA_SHAPE1(rblapack_a);
  lda = NA_SHAPE0(rblapack_a);
  if (NA_TYPE(rblapack_a) != NA_DCOMPLEX)
    rblapack_a = na_change_type(rblapack_a, NA_DCOMPLEX);
  a = NA_PTR_TYPE(rblapack_a, doublecomplex*);
  uplo = StringValueCStr(rblapack_uplo)[0];
  ldsa = MAX(1,n);
  {
    int shape[2];
    shape[0] = ldsa;
    shape[1] = n;
    rblapack_sa = na_make_object(NA_SCOMPLEX, 2, shape, cNArray);
  }
  sa = NA_PTR_TYPE(rblapack_sa, complex*);

  zlat2c_(&uplo, &n, a, &lda, sa, &ldsa, &info);

  rblapack_info = INT2NUM(info);
  return rb_ary_new3(2, rblapack_sa, rblapack_info);
}

void
init_lapack_zlat2c(VALUE mLapack){
  rb_define_module_function(mLapack, "zlat2c", rblapack_zlat2c, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}

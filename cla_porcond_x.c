#include "rb_lapack.h"

extern real cla_porcond_x_(char *uplo, integer *n, complex *a, integer *lda, complex *af, integer *ldaf, complex *x, integer *info, complex *work, real *rwork);

static VALUE sHelp, sUsage;

static VALUE
rblapack_cla_porcond_x(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_uplo;
  char uplo; 
  VALUE rblapack_a;
  complex *a; 
  VALUE rblapack_af;
  complex *af; 
  VALUE rblapack_x;
  complex *x; 
  VALUE rblapack_work;
  complex *work; 
  VALUE rblapack_rwork;
  real *rwork; 
  VALUE rblapack_info;
  integer info; 
  VALUE rblapack___out__;
  real __out__; 

  integer lda;
  integer n;
  integer ldaf;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  info, __out__ = NumRu::Lapack.cla_porcond_x( uplo, a, af, x, work, rwork, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      REAL FUNCTION CLA_PORCOND_X( UPLO, N, A, LDA, AF, LDAF, X, INFO, WORK, RWORK )\n\n*  Purpose\n*  =======\n*\n*     CLA_PORCOND_X Computes the infinity norm condition number of\n*     op(A) * diag(X) where X is a COMPLEX vector.\n*\n\n*  Arguments\n*  =========\n*\n*     UPLO    (input) CHARACTER*1\n*       = 'U':  Upper triangle of A is stored;\n*       = 'L':  Lower triangle of A is stored.\n*\n*     N       (input) INTEGER\n*     The number of linear equations, i.e., the order of the\n*     matrix A.  N >= 0.\n*\n*     A       (input) COMPLEX array, dimension (LDA,N)\n*     On entry, the N-by-N matrix A.\n*\n*     LDA     (input) INTEGER\n*     The leading dimension of the array A.  LDA >= max(1,N).\n*\n*     AF      (input) COMPLEX array, dimension (LDAF,N)\n*     The triangular factor U or L from the Cholesky factorization\n*     A = U**T*U or A = L*L**T, as computed by CPOTRF.\n*\n*     LDAF    (input) INTEGER\n*     The leading dimension of the array AF.  LDAF >= max(1,N).\n*\n*     X       (input) COMPLEX array, dimension (N)\n*     The vector X in the formula op(A) * diag(X).\n*\n*     INFO    (output) INTEGER\n*       = 0:  Successful exit.\n*     i > 0:  The ith argument is invalid.\n*\n*     WORK    (input) COMPLEX array, dimension (2*N).\n*     Workspace.\n*\n*     RWORK   (input) REAL array, dimension (N).\n*     Workspace.\n*\n\n*  =====================================================================\n*\n*     .. Local Scalars ..\n      INTEGER            KASE, I, J\n      REAL               AINVNM, ANORM, TMP\n      LOGICAL            UP\n      COMPLEX            ZDUM\n*     ..\n*     .. Local Arrays ..\n      INTEGER            ISAVE( 3 )\n*     ..\n*     .. External Functions ..\n      LOGICAL            LSAME\n      EXTERNAL           LSAME\n*     ..\n*     .. External Subroutines ..\n      EXTERNAL           CLACN2, CPOTRS, XERBLA\n*     ..\n*     .. Intrinsic Functions ..\n      INTRINSIC          ABS, MAX, REAL, AIMAG\n*     ..\n*     .. Statement Functions ..\n      REAL CABS1\n*     ..\n*     .. Statement Function Definitions ..\n      CABS1( ZDUM ) = ABS( REAL( ZDUM ) ) + ABS( AIMAG( ZDUM ) )\n*     ..\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  info, __out__ = NumRu::Lapack.cla_porcond_x( uplo, a, af, x, work, rwork, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 6)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 6)", argc);
  rblapack_uplo = argv[0];
  rblapack_a = argv[1];
  rblapack_af = argv[2];
  rblapack_x = argv[3];
  rblapack_work = argv[4];
  rblapack_rwork = argv[5];
  if (rb_options != Qnil) {
  }

  if (!NA_IsNArray(rblapack_a))
    rb_raise(rb_eArgError, "a (2th argument) must be NArray");
  if (NA_RANK(rblapack_a) != 2)
    rb_raise(rb_eArgError, "rank of a (2th argument) must be %d", 2);
  n = NA_SHAPE1(rblapack_a);
  lda = NA_SHAPE0(rblapack_a);
  if (NA_TYPE(rblapack_a) != NA_SCOMPLEX)
    rblapack_a = na_change_type(rblapack_a, NA_SCOMPLEX);
  a = NA_PTR_TYPE(rblapack_a, complex*);
  if (!NA_IsNArray(rblapack_x))
    rb_raise(rb_eArgError, "x (4th argument) must be NArray");
  if (NA_RANK(rblapack_x) != 1)
    rb_raise(rb_eArgError, "rank of x (4th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_x) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of x must be the same as shape 1 of a");
  if (NA_TYPE(rblapack_x) != NA_SCOMPLEX)
    rblapack_x = na_change_type(rblapack_x, NA_SCOMPLEX);
  x = NA_PTR_TYPE(rblapack_x, complex*);
  if (!NA_IsNArray(rblapack_rwork))
    rb_raise(rb_eArgError, "rwork (6th argument) must be NArray");
  if (NA_RANK(rblapack_rwork) != 1)
    rb_raise(rb_eArgError, "rank of rwork (6th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_rwork) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of rwork must be the same as shape 1 of a");
  if (NA_TYPE(rblapack_rwork) != NA_SFLOAT)
    rblapack_rwork = na_change_type(rblapack_rwork, NA_SFLOAT);
  rwork = NA_PTR_TYPE(rblapack_rwork, real*);
  uplo = StringValueCStr(rblapack_uplo)[0];
  if (!NA_IsNArray(rblapack_af))
    rb_raise(rb_eArgError, "af (3th argument) must be NArray");
  if (NA_RANK(rblapack_af) != 2)
    rb_raise(rb_eArgError, "rank of af (3th argument) must be %d", 2);
  if (NA_SHAPE1(rblapack_af) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of af must be the same as shape 1 of a");
  ldaf = NA_SHAPE0(rblapack_af);
  if (NA_TYPE(rblapack_af) != NA_SCOMPLEX)
    rblapack_af = na_change_type(rblapack_af, NA_SCOMPLEX);
  af = NA_PTR_TYPE(rblapack_af, complex*);
  if (!NA_IsNArray(rblapack_work))
    rb_raise(rb_eArgError, "work (5th argument) must be NArray");
  if (NA_RANK(rblapack_work) != 1)
    rb_raise(rb_eArgError, "rank of work (5th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_work) != (2*n))
    rb_raise(rb_eRuntimeError, "shape 0 of work must be %d", 2*n);
  if (NA_TYPE(rblapack_work) != NA_SCOMPLEX)
    rblapack_work = na_change_type(rblapack_work, NA_SCOMPLEX);
  work = NA_PTR_TYPE(rblapack_work, complex*);

  __out__ = cla_porcond_x_(&uplo, &n, a, &lda, af, &ldaf, x, &info, work, rwork);

  rblapack_info = INT2NUM(info);
  rblapack___out__ = rb_float_new((double)__out__);
  return rb_ary_new3(2, rblapack_info, rblapack___out__);
}

void
init_lapack_cla_porcond_x(VALUE mLapack){
  rb_define_module_function(mLapack, "cla_porcond_x", rblapack_cla_porcond_x, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}

#include "rb_lapack.h"

extern doublereal zla_syrpvgrw_(char *uplo, integer *n, integer *info, doublecomplex *a, integer *lda, doublecomplex *af, integer *ldaf, integer *ipiv, doublecomplex *work);

static VALUE sHelp, sUsage;

static VALUE
rblapack_zla_syrpvgrw(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_uplo;
  char uplo; 
  VALUE rblapack_info;
  integer info; 
  VALUE rblapack_a;
  doublecomplex *a; 
  VALUE rblapack_af;
  doublecomplex *af; 
  VALUE rblapack_ipiv;
  integer *ipiv; 
  VALUE rblapack_work;
  doublecomplex *work; 
  VALUE rblapack___out__;
  doublereal __out__; 

  integer lda;
  integer n;
  integer ldaf;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  __out__ = NumRu::Lapack.zla_syrpvgrw( uplo, info, a, af, ipiv, work, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      DOUBLE PRECISION FUNCTION ZLA_SYRPVGRW( UPLO, N, INFO, A, LDA, AF, LDAF, IPIV, WORK )\n\n*  Purpose\n*  =======\n* \n*  ZLA_SYRPVGRW computes the reciprocal pivot growth factor\n*  norm(A)/norm(U). The \"max absolute element\" norm is used. If this is\n*  much less than 1, the stability of the LU factorization of the\n*  (equilibrated) matrix A could be poor. This also means that the\n*  solution X, estimated condition numbers, and error bounds could be\n*  unreliable.\n*\n\n*  Arguments\n*  =========\n*\n*     UPLO    (input) CHARACTER*1\n*       = 'U':  Upper triangle of A is stored;\n*       = 'L':  Lower triangle of A is stored.\n*\n*     N       (input) INTEGER\n*     The number of linear equations, i.e., the order of the\n*     matrix A.  N >= 0.\n*\n*     INFO    (input) INTEGER\n*     The value of INFO returned from ZSYTRF, .i.e., the pivot in\n*     column INFO is exactly 0.\n*\n*     NCOLS   (input) INTEGER\n*     The number of columns of the matrix A. NCOLS >= 0.\n*\n*     A       (input) COMPLEX*16 array, dimension (LDA,N)\n*     On entry, the N-by-N matrix A.\n*\n*     LDA     (input) INTEGER\n*     The leading dimension of the array A.  LDA >= max(1,N).\n*\n*     AF      (input) COMPLEX*16 array, dimension (LDAF,N)\n*     The block diagonal matrix D and the multipliers used to\n*     obtain the factor U or L as computed by ZSYTRF.\n*\n*     LDAF    (input) INTEGER\n*     The leading dimension of the array AF.  LDAF >= max(1,N).\n*\n*     IPIV    (input) INTEGER array, dimension (N)\n*     Details of the interchanges and the block structure of D\n*     as determined by ZSYTRF.\n*\n*     WORK    (input) COMPLEX*16 array, dimension (2*N)\n*\n\n*  =====================================================================\n*\n*     .. Local Scalars ..\n      INTEGER            NCOLS, I, J, K, KP\n      DOUBLE PRECISION   AMAX, UMAX, RPVGRW, TMP\n      LOGICAL            UPPER\n      COMPLEX*16         ZDUM\n*     ..\n*     .. Intrinsic Functions ..\n      INTRINSIC          ABS, REAL, DIMAG, MAX, MIN\n*     ..\n*     .. External Subroutines ..\n      EXTERNAL           LSAME, ZLASET\n      LOGICAL            LSAME\n*     ..\n*     .. Statement Functions ..\n      DOUBLE PRECISION   CABS1\n*     ..\n*     .. Statement Function Definitions ..\n      CABS1( ZDUM ) = ABS( DBLE ( ZDUM ) ) + ABS( DIMAG ( ZDUM ) )\n*     ..\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  __out__ = NumRu::Lapack.zla_syrpvgrw( uplo, info, a, af, ipiv, work, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 6)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 6)", argc);
  rblapack_uplo = argv[0];
  rblapack_info = argv[1];
  rblapack_a = argv[2];
  rblapack_af = argv[3];
  rblapack_ipiv = argv[4];
  rblapack_work = argv[5];
  if (rb_options != Qnil) {
  }

  if (!NA_IsNArray(rblapack_ipiv))
    rb_raise(rb_eArgError, "ipiv (5th argument) must be NArray");
  if (NA_RANK(rblapack_ipiv) != 1)
    rb_raise(rb_eArgError, "rank of ipiv (5th argument) must be %d", 1);
  n = NA_SHAPE0(rblapack_ipiv);
  if (NA_TYPE(rblapack_ipiv) != NA_LINT)
    rblapack_ipiv = na_change_type(rblapack_ipiv, NA_LINT);
  ipiv = NA_PTR_TYPE(rblapack_ipiv, integer*);
  if (!NA_IsNArray(rblapack_a))
    rb_raise(rb_eArgError, "a (3th argument) must be NArray");
  if (NA_RANK(rblapack_a) != 2)
    rb_raise(rb_eArgError, "rank of a (3th argument) must be %d", 2);
  if (NA_SHAPE1(rblapack_a) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of a must be the same as shape 0 of ipiv");
  lda = NA_SHAPE0(rblapack_a);
  if (NA_TYPE(rblapack_a) != NA_DCOMPLEX)
    rblapack_a = na_change_type(rblapack_a, NA_DCOMPLEX);
  a = NA_PTR_TYPE(rblapack_a, doublecomplex*);
  if (!NA_IsNArray(rblapack_af))
    rb_raise(rb_eArgError, "af (4th argument) must be NArray");
  if (NA_RANK(rblapack_af) != 2)
    rb_raise(rb_eArgError, "rank of af (4th argument) must be %d", 2);
  if (NA_SHAPE1(rblapack_af) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of af must be the same as shape 0 of ipiv");
  ldaf = NA_SHAPE0(rblapack_af);
  if (NA_TYPE(rblapack_af) != NA_DCOMPLEX)
    rblapack_af = na_change_type(rblapack_af, NA_DCOMPLEX);
  af = NA_PTR_TYPE(rblapack_af, doublecomplex*);
  uplo = StringValueCStr(rblapack_uplo)[0];
  info = NUM2INT(rblapack_info);
  if (!NA_IsNArray(rblapack_work))
    rb_raise(rb_eArgError, "work (6th argument) must be NArray");
  if (NA_RANK(rblapack_work) != 1)
    rb_raise(rb_eArgError, "rank of work (6th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_work) != (2*n))
    rb_raise(rb_eRuntimeError, "shape 0 of work must be %d", 2*n);
  if (NA_TYPE(rblapack_work) != NA_DCOMPLEX)
    rblapack_work = na_change_type(rblapack_work, NA_DCOMPLEX);
  work = NA_PTR_TYPE(rblapack_work, doublecomplex*);

  __out__ = zla_syrpvgrw_(&uplo, &n, &info, a, &lda, af, &ldaf, ipiv, work);

  rblapack___out__ = rb_float_new((double)__out__);
  return rblapack___out__;
}

void
init_lapack_zla_syrpvgrw(VALUE mLapack){
  rb_define_module_function(mLapack, "zla_syrpvgrw", rblapack_zla_syrpvgrw, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}

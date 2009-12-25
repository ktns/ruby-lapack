#include "rb_lapack.h"

extern VOID zla_herpvgrw_(doublereal *__out__, char *uplo, integer *n, integer *info, doublecomplex *a, integer *lda, doublecomplex *af, integer *ldaf, integer *ipiv, doublecomplex *work);
static VALUE
rb_zla_herpvgrw(int argc, VALUE *argv, VALUE self){
  VALUE rb_uplo;
  char uplo; 
  VALUE rb_info;
  integer info; 
  VALUE rb_a;
  doublecomplex *a; 
  VALUE rb_af;
  doublecomplex *af; 
  VALUE rb_ipiv;
  integer *ipiv; 
  VALUE rb_work;
  doublecomplex *work; 
  VALUE rb___out__;
  doublereal __out__; 

  integer lda;
  integer n;
  integer ldaf;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  __out__ = NumRu::Lapack.zla_herpvgrw( uplo, info, a, af, ipiv, work)\n    or\n  NumRu::Lapack.zla_herpvgrw  # print help\n\n\nFORTRAN MANUAL\n      DOUBLE PRECISION FUNCTION ZLA_HERPVGRW( UPLO, N, INFO, A, LDA, AF, LDAF, IPIV, WORK )\n\n*  Purpose\n*  =======\n* \n*  ZLA_HERPVGRW computes the reciprocal pivot growth factor\n*  norm(A)/norm(U). The \"max absolute element\" norm is used. If this is\n*  much less than 1, the stability of the LU factorization of the\n*  (equilibrated) matrix A could be poor. This also means that the\n*  solution X, estimated condition numbers, and error bounds could be\n*  unreliable.\n*\n\n*  Arguments\n*  =========\n*\n*     UPLO    (input) CHARACTER*1\n*       = 'U':  Upper triangle of A is stored;\n*       = 'L':  Lower triangle of A is stored.\n*\n*     N       (input) INTEGER\n*     The number of linear equations, i.e., the order of the\n*     matrix A.  N >= 0.\n*\n*     INFO    (input) INTEGER\n*     The value of INFO returned from ZHETRF, .i.e., the pivot in\n*     column INFO is exactly 0.\n*\n*     NCOLS   (input) INTEGER\n*     The number of columns of the matrix A. NCOLS >= 0.\n*\n*     A       (input) COMPLEX*16 array, dimension (LDA,N)\n*     On entry, the N-by-N matrix A.\n*\n*     LDA     (input) INTEGER\n*     The leading dimension of the array A.  LDA >= max(1,N).\n*\n*     AF      (input) COMPLEX*16 array, dimension (LDAF,N)\n*     The block diagonal matrix D and the multipliers used to\n*     obtain the factor U or L as computed by ZHETRF.\n*\n*     LDAF    (input) INTEGER\n*     The leading dimension of the array AF.  LDAF >= max(1,N).\n*\n*     IPIV    (input) INTEGER array, dimension (N)\n*     Details of the interchanges and the block structure of D\n*     as determined by ZHETRF.\n*\n*     WORK    (input) COMPLEX*16 array, dimension (2*N)\n*\n\n*  =====================================================================\n*\n*     .. Local Scalars ..\n      INTEGER            NCOLS, I, J, K, KP\n      DOUBLE PRECISION   AMAX, UMAX, RPVGRW, TMP\n      LOGICAL            UPPER, LSAME\n      COMPLEX*16         ZDUM\n*     ..\n*     .. External Functions ..\n      EXTERNAL           LSAME, ZLASET\n*     ..\n*     .. Intrinsic Functions ..\n      INTRINSIC          ABS, REAL, DIMAG, MAX, MIN\n*     ..\n*     .. Statement Functions ..\n      DOUBLE PRECISION   CABS1\n*     ..\n*     .. Statement Function Definitions ..\n      CABS1( ZDUM ) = ABS( DBLE ( ZDUM ) ) + ABS( DIMAG ( ZDUM ) )\n*     ..\n\n");
    return Qnil;
  }
  if (argc != 6)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 6)", argc);
  rb_uplo = argv[0];
  rb_info = argv[1];
  rb_a = argv[2];
  rb_af = argv[3];
  rb_ipiv = argv[4];
  rb_work = argv[5];

  uplo = StringValueCStr(rb_uplo)[0];
  info = NUM2INT(rb_info);
  if (!NA_IsNArray(rb_a))
    rb_raise(rb_eArgError, "a (3th argument) must be NArray");
  if (NA_RANK(rb_a) != 2)
    rb_raise(rb_eArgError, "rank of a (3th argument) must be %d", 2);
  lda = NA_SHAPE0(rb_a);
  n = NA_SHAPE1(rb_a);
  if (NA_TYPE(rb_a) != NA_DCOMPLEX)
    rb_a = na_change_type(rb_a, NA_DCOMPLEX);
  a = NA_PTR_TYPE(rb_a, doublecomplex*);
  if (!NA_IsNArray(rb_af))
    rb_raise(rb_eArgError, "af (4th argument) must be NArray");
  if (NA_RANK(rb_af) != 2)
    rb_raise(rb_eArgError, "rank of af (4th argument) must be %d", 2);
  ldaf = NA_SHAPE0(rb_af);
  if (NA_SHAPE1(rb_af) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of af must be the same as shape 1 of a");
  if (NA_TYPE(rb_af) != NA_DCOMPLEX)
    rb_af = na_change_type(rb_af, NA_DCOMPLEX);
  af = NA_PTR_TYPE(rb_af, doublecomplex*);
  if (!NA_IsNArray(rb_ipiv))
    rb_raise(rb_eArgError, "ipiv (5th argument) must be NArray");
  if (NA_RANK(rb_ipiv) != 1)
    rb_raise(rb_eArgError, "rank of ipiv (5th argument) must be %d", 1);
  if (NA_SHAPE0(rb_ipiv) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of ipiv must be the same as shape 1 of a");
  if (NA_TYPE(rb_ipiv) != NA_LINT)
    rb_ipiv = na_change_type(rb_ipiv, NA_LINT);
  ipiv = NA_PTR_TYPE(rb_ipiv, integer*);
  if (!NA_IsNArray(rb_work))
    rb_raise(rb_eArgError, "work (6th argument) must be NArray");
  if (NA_RANK(rb_work) != 1)
    rb_raise(rb_eArgError, "rank of work (6th argument) must be %d", 1);
  if (NA_SHAPE0(rb_work) != (2*n))
    rb_raise(rb_eRuntimeError, "shape 0 of work must be %d", 2*n);
  if (NA_TYPE(rb_work) != NA_DCOMPLEX)
    rb_work = na_change_type(rb_work, NA_DCOMPLEX);
  work = NA_PTR_TYPE(rb_work, doublecomplex*);

  zla_herpvgrw_(&__out__, &uplo, &n, &info, a, &lda, af, &ldaf, ipiv, work);

  rb___out__ = rb_float_new((double)__out__);
  return rb___out__;
}

void
init_lapack_zla_herpvgrw(VALUE mLapack){
  rb_define_module_function(mLapack, "zla_herpvgrw", rb_zla_herpvgrw, -1);
}

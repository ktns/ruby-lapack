#include "rb_lapack.h"

extern VOID zla_gercond_x_(doublereal *__out__, char *trans, integer *n, doublecomplex *a, integer *lda, doublecomplex *af, integer *ldaf, integer *ipiv, doublecomplex *x, integer *info, doublecomplex *work, doublereal *rwork);
static VALUE
rb_zla_gercond_x(int argc, VALUE *argv, VALUE self){
  VALUE rb_trans;
  char trans; 
  VALUE rb_a;
  doublecomplex *a; 
  VALUE rb_af;
  doublecomplex *af; 
  VALUE rb_ipiv;
  integer *ipiv; 
  VALUE rb_x;
  doublecomplex *x; 
  VALUE rb_work;
  doublecomplex *work; 
  VALUE rb_rwork;
  doublereal *rwork; 
  VALUE rb_info;
  integer info; 
  VALUE rb___out__;
  doublereal __out__; 

  integer lda;
  integer n;
  integer ldaf;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  info, __out__ = NumRu::Lapack.zla_gercond_x( trans, a, af, ipiv, x, work, rwork)\n    or\n  NumRu::Lapack.zla_gercond_x  # print help\n\n\nFORTRAN MANUAL\n      DOUBLE PRECISION FUNCTION ZLA_GERCOND_X( TRANS, N, A, LDA, AF, LDAF, IPIV, X, INFO, WORK, RWORK )\n\n*  Purpose\n*  =======\n*\n*     ZLA_GERCOND_X computes the infinity norm condition number of\n*     op(A) * diag(X) where X is a COMPLEX*16 vector.\n*\n\n*  Arguments\n*  =========\n*\n*     TRANS   (input) CHARACTER*1\n*     Specifies the form of the system of equations:\n*       = 'N':  A * X = B     (No transpose)\n*       = 'T':  A**T * X = B  (Transpose)\n*       = 'C':  A**H * X = B  (Conjugate Transpose = Transpose)\n*\n*     N       (input) INTEGER\n*     The number of linear equations, i.e., the order of the\n*     matrix A.  N >= 0.\n*\n*     A       (input) COMPLEX*16 array, dimension (LDA,N)\n*     On entry, the N-by-N matrix A.\n*\n*     LDA     (input) INTEGER\n*     The leading dimension of the array A.  LDA >= max(1,N).\n*\n*     AF      (input) COMPLEX*16 array, dimension (LDAF,N)\n*     The factors L and U from the factorization\n*     A = P*L*U as computed by ZGETRF.\n*\n*     LDAF    (input) INTEGER\n*     The leading dimension of the array AF.  LDAF >= max(1,N).\n*\n*     IPIV    (input) INTEGER array, dimension (N)\n*     The pivot indices from the factorization A = P*L*U\n*     as computed by ZGETRF; row i of the matrix was interchanged\n*     with row IPIV(i).\n*\n*     X       (input) COMPLEX*16 array, dimension (N)\n*     The vector X in the formula op(A) * diag(X).\n*\n*     INFO    (output) INTEGER\n*       = 0:  Successful exit.\n*     i > 0:  The ith argument is invalid.\n*\n*     WORK    (input) COMPLEX*16 array, dimension (2*N).\n*     Workspace.\n*\n*     RWORK   (input) DOUBLE PRECISION array, dimension (N).\n*     Workspace.\n*\n\n*  =====================================================================\n*\n*     .. Local Scalars ..\n      LOGICAL            NOTRANS\n      INTEGER            KASE\n      DOUBLE PRECISION   AINVNM, ANORM, TMP\n      INTEGER            I, J\n      COMPLEX*16         ZDUM\n*     ..\n*     .. Local Arrays ..\n      INTEGER            ISAVE( 3 )\n*     ..\n*     .. External Functions ..\n      LOGICAL            LSAME\n      EXTERNAL           LSAME\n*     ..\n*     .. External Subroutines ..\n      EXTERNAL           ZLACN2, ZGETRS, XERBLA\n*     ..\n*     .. Intrinsic Functions ..\n      INTRINSIC          ABS, MAX, REAL, DIMAG\n*     ..\n*     .. Statement Functions ..\n      DOUBLE PRECISION   CABS1\n*     ..\n*     .. Statement Function Definitions ..\n      CABS1( ZDUM ) = ABS( DBLE( ZDUM ) ) + ABS( DIMAG( ZDUM ) )\n*     ..\n\n");
    return Qnil;
  }
  if (argc != 7)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 7)", argc);
  rb_trans = argv[0];
  rb_a = argv[1];
  rb_af = argv[2];
  rb_ipiv = argv[3];
  rb_x = argv[4];
  rb_work = argv[5];
  rb_rwork = argv[6];

  trans = StringValueCStr(rb_trans)[0];
  if (!NA_IsNArray(rb_a))
    rb_raise(rb_eArgError, "a (2th argument) must be NArray");
  if (NA_RANK(rb_a) != 2)
    rb_raise(rb_eArgError, "rank of a (2th argument) must be %d", 2);
  lda = NA_SHAPE0(rb_a);
  n = NA_SHAPE1(rb_a);
  if (NA_TYPE(rb_a) != NA_DCOMPLEX)
    rb_a = na_change_type(rb_a, NA_DCOMPLEX);
  a = NA_PTR_TYPE(rb_a, doublecomplex*);
  if (!NA_IsNArray(rb_af))
    rb_raise(rb_eArgError, "af (3th argument) must be NArray");
  if (NA_RANK(rb_af) != 2)
    rb_raise(rb_eArgError, "rank of af (3th argument) must be %d", 2);
  ldaf = NA_SHAPE0(rb_af);
  if (NA_SHAPE1(rb_af) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of af must be the same as shape 1 of a");
  if (NA_TYPE(rb_af) != NA_DCOMPLEX)
    rb_af = na_change_type(rb_af, NA_DCOMPLEX);
  af = NA_PTR_TYPE(rb_af, doublecomplex*);
  if (!NA_IsNArray(rb_ipiv))
    rb_raise(rb_eArgError, "ipiv (4th argument) must be NArray");
  if (NA_RANK(rb_ipiv) != 1)
    rb_raise(rb_eArgError, "rank of ipiv (4th argument) must be %d", 1);
  if (NA_SHAPE0(rb_ipiv) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of ipiv must be the same as shape 1 of a");
  if (NA_TYPE(rb_ipiv) != NA_LINT)
    rb_ipiv = na_change_type(rb_ipiv, NA_LINT);
  ipiv = NA_PTR_TYPE(rb_ipiv, integer*);
  if (!NA_IsNArray(rb_x))
    rb_raise(rb_eArgError, "x (5th argument) must be NArray");
  if (NA_RANK(rb_x) != 1)
    rb_raise(rb_eArgError, "rank of x (5th argument) must be %d", 1);
  if (NA_SHAPE0(rb_x) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of x must be the same as shape 1 of a");
  if (NA_TYPE(rb_x) != NA_DCOMPLEX)
    rb_x = na_change_type(rb_x, NA_DCOMPLEX);
  x = NA_PTR_TYPE(rb_x, doublecomplex*);
  if (!NA_IsNArray(rb_work))
    rb_raise(rb_eArgError, "work (6th argument) must be NArray");
  if (NA_RANK(rb_work) != 1)
    rb_raise(rb_eArgError, "rank of work (6th argument) must be %d", 1);
  if (NA_SHAPE0(rb_work) != (2*n))
    rb_raise(rb_eRuntimeError, "shape 0 of work must be %d", 2*n);
  if (NA_TYPE(rb_work) != NA_DCOMPLEX)
    rb_work = na_change_type(rb_work, NA_DCOMPLEX);
  work = NA_PTR_TYPE(rb_work, doublecomplex*);
  if (!NA_IsNArray(rb_rwork))
    rb_raise(rb_eArgError, "rwork (7th argument) must be NArray");
  if (NA_RANK(rb_rwork) != 1)
    rb_raise(rb_eArgError, "rank of rwork (7th argument) must be %d", 1);
  if (NA_SHAPE0(rb_rwork) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of rwork must be the same as shape 1 of a");
  if (NA_TYPE(rb_rwork) != NA_DFLOAT)
    rb_rwork = na_change_type(rb_rwork, NA_DFLOAT);
  rwork = NA_PTR_TYPE(rb_rwork, doublereal*);

  zla_gercond_x_(&__out__, &trans, &n, a, &lda, af, &ldaf, ipiv, x, &info, work, rwork);

  rb_info = INT2NUM(info);
  rb___out__ = rb_float_new((double)__out__);
  return rb_ary_new3(2, rb_info, rb___out__);
}

void
init_lapack_zla_gercond_x(VALUE mLapack){
  rb_define_module_function(mLapack, "zla_gercond_x", rb_zla_gercond_x, -1);
}

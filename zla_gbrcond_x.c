#include "rb_lapack.h"

extern doublereal zla_gbrcond_x_(char *trans, integer *n, integer *kl, integer *ku, doublecomplex *ab, integer *ldab, doublecomplex *afb, integer *ldafb, integer *ipiv, doublecomplex *x, integer *info, doublecomplex *work, doublereal *rwork);

static VALUE
rb_zla_gbrcond_x(int argc, VALUE *argv, VALUE self){
  VALUE rb_trans;
  char trans; 
  VALUE rb_kl;
  integer kl; 
  VALUE rb_ku;
  integer ku; 
  VALUE rb_ab;
  doublecomplex *ab; 
  VALUE rb_afb;
  doublecomplex *afb; 
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

  integer ldab;
  integer n;
  integer ldafb;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  info, __out__ = NumRu::Lapack.zla_gbrcond_x( trans, kl, ku, ab, afb, ipiv, x, work, rwork)\n    or\n  NumRu::Lapack.zla_gbrcond_x  # print help\n\n\nFORTRAN MANUAL\n      DOUBLE PRECISION FUNCTION ZLA_GBRCOND_X( TRANS, N, KL, KU, AB, LDAB, AFB, LDAFB, IPIV, X, INFO, WORK, RWORK )\n\n*  Purpose\n*  =======\n*\n*     ZLA_GBRCOND_X Computes the infinity norm condition number of\n*     op(A) * diag(X) where X is a COMPLEX*16 vector.\n*\n\n*  Arguments\n*  =========\n*\n*     TRANS   (input) CHARACTER*1\n*     Specifies the form of the system of equations:\n*       = 'N':  A * X = B     (No transpose)\n*       = 'T':  A**T * X = B  (Transpose)\n*       = 'C':  A**H * X = B  (Conjugate Transpose = Transpose)\n*\n*     N       (input) INTEGER\n*     The number of linear equations, i.e., the order of the\n*     matrix A.  N >= 0.\n*\n*     KL      (input) INTEGER\n*     The number of subdiagonals within the band of A.  KL >= 0.\n*\n*     KU      (input) INTEGER\n*     The number of superdiagonals within the band of A.  KU >= 0.\n*\n*     AB      (input) COMPLEX*16 array, dimension (LDAB,N)\n*     On entry, the matrix A in band storage, in rows 1 to KL+KU+1.\n*     The j-th column of A is stored in the j-th column of the\n*     array AB as follows:\n*     AB(KU+1+i-j,j) = A(i,j) for max(1,j-KU)<=i<=min(N,j+kl)\n*\n*     LDAB    (input) INTEGER\n*     The leading dimension of the array AB.  LDAB >= KL+KU+1.\n*\n*     AFB     (input) COMPLEX*16 array, dimension (LDAFB,N)\n*     Details of the LU factorization of the band matrix A, as\n*     computed by ZGBTRF.  U is stored as an upper triangular\n*     band matrix with KL+KU superdiagonals in rows 1 to KL+KU+1,\n*     and the multipliers used during the factorization are stored\n*     in rows KL+KU+2 to 2*KL+KU+1.\n*\n*     LDAFB   (input) INTEGER\n*     The leading dimension of the array AFB.  LDAFB >= 2*KL+KU+1.\n*\n*     IPIV    (input) INTEGER array, dimension (N)\n*     The pivot indices from the factorization A = P*L*U\n*     as computed by ZGBTRF; row i of the matrix was interchanged\n*     with row IPIV(i).\n*\n*     X       (input) COMPLEX*16 array, dimension (N)\n*     The vector X in the formula op(A) * diag(X).\n*\n*     INFO    (output) INTEGER\n*       = 0:  Successful exit.\n*     i > 0:  The ith argument is invalid.\n*\n*     WORK    (input) COMPLEX*16 array, dimension (2*N).\n*     Workspace.\n*\n*     RWORK   (input) DOUBLE PRECISION array, dimension (N).\n*     Workspace.\n*\n\n*  =====================================================================\n*\n*     .. Local Scalars ..\n      LOGICAL            NOTRANS\n      INTEGER            KASE, I, J\n      DOUBLE PRECISION   AINVNM, ANORM, TMP\n      COMPLEX*16         ZDUM\n*     ..\n*     .. Local Arrays ..\n      INTEGER            ISAVE( 3 )\n*     ..\n*     .. External Functions ..\n      LOGICAL            LSAME\n      EXTERNAL           LSAME\n*     ..\n*     .. External Subroutines ..\n      EXTERNAL           ZLACN2, ZGBTRS, XERBLA\n*     ..\n*     .. Intrinsic Functions ..\n      INTRINSIC          ABS, MAX\n*     ..\n*     .. Statement Functions ..\n      DOUBLE PRECISION   CABS1\n*     ..\n*     .. Statement Function Definitions ..\n      CABS1( ZDUM ) = ABS( DBLE( ZDUM ) ) + ABS( DIMAG( ZDUM ) )\n*     ..\n\n");
    return Qnil;
  }
  if (argc != 9)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 9)", argc);
  rb_trans = argv[0];
  rb_kl = argv[1];
  rb_ku = argv[2];
  rb_ab = argv[3];
  rb_afb = argv[4];
  rb_ipiv = argv[5];
  rb_x = argv[6];
  rb_work = argv[7];
  rb_rwork = argv[8];

  if (!NA_IsNArray(rb_ipiv))
    rb_raise(rb_eArgError, "ipiv (6th argument) must be NArray");
  if (NA_RANK(rb_ipiv) != 1)
    rb_raise(rb_eArgError, "rank of ipiv (6th argument) must be %d", 1);
  n = NA_SHAPE0(rb_ipiv);
  if (NA_TYPE(rb_ipiv) != NA_LINT)
    rb_ipiv = na_change_type(rb_ipiv, NA_LINT);
  ipiv = NA_PTR_TYPE(rb_ipiv, integer*);
  if (!NA_IsNArray(rb_ab))
    rb_raise(rb_eArgError, "ab (4th argument) must be NArray");
  if (NA_RANK(rb_ab) != 2)
    rb_raise(rb_eArgError, "rank of ab (4th argument) must be %d", 2);
  if (NA_SHAPE1(rb_ab) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of ab must be the same as shape 0 of ipiv");
  ldab = NA_SHAPE0(rb_ab);
  if (NA_TYPE(rb_ab) != NA_DCOMPLEX)
    rb_ab = na_change_type(rb_ab, NA_DCOMPLEX);
  ab = NA_PTR_TYPE(rb_ab, doublecomplex*);
  kl = NUM2INT(rb_kl);
  if (!NA_IsNArray(rb_x))
    rb_raise(rb_eArgError, "x (7th argument) must be NArray");
  if (NA_RANK(rb_x) != 1)
    rb_raise(rb_eArgError, "rank of x (7th argument) must be %d", 1);
  if (NA_SHAPE0(rb_x) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of x must be the same as shape 0 of ipiv");
  if (NA_TYPE(rb_x) != NA_DCOMPLEX)
    rb_x = na_change_type(rb_x, NA_DCOMPLEX);
  x = NA_PTR_TYPE(rb_x, doublecomplex*);
  ku = NUM2INT(rb_ku);
  if (!NA_IsNArray(rb_afb))
    rb_raise(rb_eArgError, "afb (5th argument) must be NArray");
  if (NA_RANK(rb_afb) != 2)
    rb_raise(rb_eArgError, "rank of afb (5th argument) must be %d", 2);
  if (NA_SHAPE1(rb_afb) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of afb must be the same as shape 0 of ipiv");
  ldafb = NA_SHAPE0(rb_afb);
  if (NA_TYPE(rb_afb) != NA_DCOMPLEX)
    rb_afb = na_change_type(rb_afb, NA_DCOMPLEX);
  afb = NA_PTR_TYPE(rb_afb, doublecomplex*);
  trans = StringValueCStr(rb_trans)[0];
  if (!NA_IsNArray(rb_rwork))
    rb_raise(rb_eArgError, "rwork (9th argument) must be NArray");
  if (NA_RANK(rb_rwork) != 1)
    rb_raise(rb_eArgError, "rank of rwork (9th argument) must be %d", 1);
  if (NA_SHAPE0(rb_rwork) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of rwork must be the same as shape 0 of ipiv");
  if (NA_TYPE(rb_rwork) != NA_DFLOAT)
    rb_rwork = na_change_type(rb_rwork, NA_DFLOAT);
  rwork = NA_PTR_TYPE(rb_rwork, doublereal*);
  if (!NA_IsNArray(rb_work))
    rb_raise(rb_eArgError, "work (8th argument) must be NArray");
  if (NA_RANK(rb_work) != 1)
    rb_raise(rb_eArgError, "rank of work (8th argument) must be %d", 1);
  if (NA_SHAPE0(rb_work) != (2*n))
    rb_raise(rb_eRuntimeError, "shape 0 of work must be %d", 2*n);
  if (NA_TYPE(rb_work) != NA_DCOMPLEX)
    rb_work = na_change_type(rb_work, NA_DCOMPLEX);
  work = NA_PTR_TYPE(rb_work, doublecomplex*);

  __out__ = zla_gbrcond_x_(&trans, &n, &kl, &ku, ab, &ldab, afb, &ldafb, ipiv, x, &info, work, rwork);

  rb_info = INT2NUM(info);
  rb___out__ = rb_float_new((double)__out__);
  return rb_ary_new3(2, rb_info, rb___out__);
}

void
init_lapack_zla_gbrcond_x(VALUE mLapack){
  rb_define_module_function(mLapack, "zla_gbrcond_x", rb_zla_gbrcond_x, -1);
}

#include "rb_lapack.h"

extern real cla_gbrcond_c_(char *trans, integer *n, integer *kl, integer *ku, complex *ab, integer *ldab, complex *afb, integer *ldafb, integer *ipiv, real *c, logical *capply, integer *info, complex *work, real *rwork);

static VALUE sHelp, sUsage;

static VALUE
rblapack_cla_gbrcond_c(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_trans;
  char trans; 
  VALUE rblapack_kl;
  integer kl; 
  VALUE rblapack_ku;
  integer ku; 
  VALUE rblapack_ab;
  complex *ab; 
  VALUE rblapack_afb;
  complex *afb; 
  VALUE rblapack_ipiv;
  integer *ipiv; 
  VALUE rblapack_c;
  real *c; 
  VALUE rblapack_capply;
  logical capply; 
  VALUE rblapack_work;
  complex *work; 
  VALUE rblapack_rwork;
  real *rwork; 
  VALUE rblapack_info;
  integer info; 
  VALUE rblapack___out__;
  real __out__; 

  integer ldab;
  integer n;
  integer ldafb;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  info, __out__ = NumRu::Lapack.cla_gbrcond_c( trans, kl, ku, ab, afb, ipiv, c, capply, work, rwork, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      REAL FUNCTION CLA_GBRCOND_C( TRANS, N, KL, KU, AB, LDAB, AFB, LDAFB, IPIV, C, CAPPLY, INFO, WORK, RWORK )\n\n*  Purpose\n*  =======\n*\n*     CLA_GBRCOND_C Computes the infinity norm condition number of\n*     op(A) * inv(diag(C)) where C is a REAL vector.\n*\n\n*  Arguments\n*  =========\n*\n*     TRANS   (input) CHARACTER*1\n*     Specifies the form of the system of equations:\n*       = 'N':  A * X = B     (No transpose)\n*       = 'T':  A**T * X = B  (Transpose)\n*       = 'C':  A**H * X = B  (Conjugate Transpose = Transpose)\n*\n*     N       (input) INTEGER\n*     The number of linear equations, i.e., the order of the\n*     matrix A.  N >= 0.\n*\n*     KL      (input) INTEGER\n*     The number of subdiagonals within the band of A.  KL >= 0.\n*\n*     KU      (input) INTEGER\n*     The number of superdiagonals within the band of A.  KU >= 0.\n*\n*     AB      (input) COMPLEX array, dimension (LDAB,N)\n*     On entry, the matrix A in band storage, in rows 1 to KL+KU+1.\n*     The j-th column of A is stored in the j-th column of the\n*     array AB as follows:\n*     AB(KU+1+i-j,j) = A(i,j) for max(1,j-KU)<=i<=min(N,j+kl)\n*\n*     LDAB    (input) INTEGER\n*     The leading dimension of the array AB.  LDAB >= KL+KU+1.\n*\n*     AFB     (input) COMPLEX array, dimension (LDAFB,N)\n*     Details of the LU factorization of the band matrix A, as\n*     computed by CGBTRF.  U is stored as an upper triangular\n*     band matrix with KL+KU superdiagonals in rows 1 to KL+KU+1,\n*     and the multipliers used during the factorization are stored\n*     in rows KL+KU+2 to 2*KL+KU+1.\n*\n*     LDAFB   (input) INTEGER\n*     The leading dimension of the array AFB.  LDAFB >= 2*KL+KU+1.\n*\n*     IPIV    (input) INTEGER array, dimension (N)\n*     The pivot indices from the factorization A = P*L*U\n*     as computed by CGBTRF; row i of the matrix was interchanged\n*     with row IPIV(i).\n*\n*     C       (input) REAL array, dimension (N)\n*     The vector C in the formula op(A) * inv(diag(C)).\n*\n*     CAPPLY  (input) LOGICAL\n*     If .TRUE. then access the vector C in the formula above.\n*\n*     INFO    (output) INTEGER\n*       = 0:  Successful exit.\n*     i > 0:  The ith argument is invalid.\n*\n*     WORK    (input) COMPLEX array, dimension (2*N).\n*     Workspace.\n*\n*     RWORK   (input) REAL array, dimension (N).\n*     Workspace.\n*\n\n*  =====================================================================\n*\n*     .. Local Scalars ..\n      LOGICAL            NOTRANS\n      INTEGER            KASE, I, J\n      REAL               AINVNM, ANORM, TMP\n      COMPLEX            ZDUM\n*     ..\n*     .. Local Arrays ..\n      INTEGER            ISAVE( 3 )\n*     ..\n*     .. External Functions ..\n      LOGICAL            LSAME\n      EXTERNAL           LSAME\n*     ..\n*     .. External Subroutines ..\n      EXTERNAL           CLACN2, CGBTRS, XERBLA\n*     ..\n*     .. Intrinsic Functions ..\n      INTRINSIC          ABS, MAX\n*     ..\n*     .. Statement Functions ..\n      REAL               CABS1\n*     ..\n*     .. Statement Function Definitions ..\n      CABS1( ZDUM ) = ABS( REAL( ZDUM ) ) + ABS( AIMAG( ZDUM ) )\n*     ..\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  info, __out__ = NumRu::Lapack.cla_gbrcond_c( trans, kl, ku, ab, afb, ipiv, c, capply, work, rwork, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 10)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 10)", argc);
  rblapack_trans = argv[0];
  rblapack_kl = argv[1];
  rblapack_ku = argv[2];
  rblapack_ab = argv[3];
  rblapack_afb = argv[4];
  rblapack_ipiv = argv[5];
  rblapack_c = argv[6];
  rblapack_capply = argv[7];
  rblapack_work = argv[8];
  rblapack_rwork = argv[9];
  if (rb_options != Qnil) {
  }

  if (!NA_IsNArray(rblapack_ipiv))
    rb_raise(rb_eArgError, "ipiv (6th argument) must be NArray");
  if (NA_RANK(rblapack_ipiv) != 1)
    rb_raise(rb_eArgError, "rank of ipiv (6th argument) must be %d", 1);
  n = NA_SHAPE0(rblapack_ipiv);
  if (NA_TYPE(rblapack_ipiv) != NA_LINT)
    rblapack_ipiv = na_change_type(rblapack_ipiv, NA_LINT);
  ipiv = NA_PTR_TYPE(rblapack_ipiv, integer*);
  if (!NA_IsNArray(rblapack_ab))
    rb_raise(rb_eArgError, "ab (4th argument) must be NArray");
  if (NA_RANK(rblapack_ab) != 2)
    rb_raise(rb_eArgError, "rank of ab (4th argument) must be %d", 2);
  if (NA_SHAPE1(rblapack_ab) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of ab must be the same as shape 0 of ipiv");
  ldab = NA_SHAPE0(rblapack_ab);
  if (NA_TYPE(rblapack_ab) != NA_SCOMPLEX)
    rblapack_ab = na_change_type(rblapack_ab, NA_SCOMPLEX);
  ab = NA_PTR_TYPE(rblapack_ab, complex*);
  kl = NUM2INT(rblapack_kl);
  if (!NA_IsNArray(rblapack_rwork))
    rb_raise(rb_eArgError, "rwork (10th argument) must be NArray");
  if (NA_RANK(rblapack_rwork) != 1)
    rb_raise(rb_eArgError, "rank of rwork (10th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_rwork) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of rwork must be the same as shape 0 of ipiv");
  if (NA_TYPE(rblapack_rwork) != NA_SFLOAT)
    rblapack_rwork = na_change_type(rblapack_rwork, NA_SFLOAT);
  rwork = NA_PTR_TYPE(rblapack_rwork, real*);
  if (!NA_IsNArray(rblapack_c))
    rb_raise(rb_eArgError, "c (7th argument) must be NArray");
  if (NA_RANK(rblapack_c) != 1)
    rb_raise(rb_eArgError, "rank of c (7th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_c) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of c must be the same as shape 0 of ipiv");
  if (NA_TYPE(rblapack_c) != NA_SFLOAT)
    rblapack_c = na_change_type(rblapack_c, NA_SFLOAT);
  c = NA_PTR_TYPE(rblapack_c, real*);
  ku = NUM2INT(rblapack_ku);
  if (!NA_IsNArray(rblapack_afb))
    rb_raise(rb_eArgError, "afb (5th argument) must be NArray");
  if (NA_RANK(rblapack_afb) != 2)
    rb_raise(rb_eArgError, "rank of afb (5th argument) must be %d", 2);
  if (NA_SHAPE1(rblapack_afb) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of afb must be the same as shape 0 of ipiv");
  ldafb = NA_SHAPE0(rblapack_afb);
  if (NA_TYPE(rblapack_afb) != NA_SCOMPLEX)
    rblapack_afb = na_change_type(rblapack_afb, NA_SCOMPLEX);
  afb = NA_PTR_TYPE(rblapack_afb, complex*);
  trans = StringValueCStr(rblapack_trans)[0];
  capply = (rblapack_capply == Qtrue);
  if (!NA_IsNArray(rblapack_work))
    rb_raise(rb_eArgError, "work (9th argument) must be NArray");
  if (NA_RANK(rblapack_work) != 1)
    rb_raise(rb_eArgError, "rank of work (9th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_work) != (2*n))
    rb_raise(rb_eRuntimeError, "shape 0 of work must be %d", 2*n);
  if (NA_TYPE(rblapack_work) != NA_SCOMPLEX)
    rblapack_work = na_change_type(rblapack_work, NA_SCOMPLEX);
  work = NA_PTR_TYPE(rblapack_work, complex*);

  __out__ = cla_gbrcond_c_(&trans, &n, &kl, &ku, ab, &ldab, afb, &ldafb, ipiv, c, &capply, &info, work, rwork);

  rblapack_info = INT2NUM(info);
  rblapack___out__ = rb_float_new((double)__out__);
  return rb_ary_new3(2, rblapack_info, rblapack___out__);
}

void
init_lapack_cla_gbrcond_c(VALUE mLapack){
  rb_define_module_function(mLapack, "cla_gbrcond_c", rblapack_cla_gbrcond_c, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}

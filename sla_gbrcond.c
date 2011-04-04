#include "rb_lapack.h"

extern real sla_gbrcond_(char *trans, integer *n, integer *kl, integer *ku, real *ab, integer *ldab, real *afb, integer *ldafb, integer *ipiv, integer *cmode, real *c, integer *info, real *work, integer *iwork);

static VALUE
rb_sla_gbrcond(int argc, VALUE *argv, VALUE self){
  VALUE rb_trans;
  char trans; 
  VALUE rb_kl;
  integer kl; 
  VALUE rb_ku;
  integer ku; 
  VALUE rb_ab;
  real *ab; 
  VALUE rb_afb;
  real *afb; 
  VALUE rb_ipiv;
  integer *ipiv; 
  VALUE rb_cmode;
  integer cmode; 
  VALUE rb_c;
  real *c; 
  VALUE rb_work;
  real *work; 
  VALUE rb_iwork;
  integer *iwork; 
  VALUE rb_info;
  integer info; 
  VALUE rb___out__;
  real __out__; 

  integer ldab;
  integer n;
  integer ldafb;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  info, __out__ = NumRu::Lapack.sla_gbrcond( trans, kl, ku, ab, afb, ipiv, cmode, c, work, iwork)\n    or\n  NumRu::Lapack.sla_gbrcond  # print help\n\n\nFORTRAN MANUAL\n      REAL FUNCTION SLA_GBRCOND( TRANS, N, KL, KU, AB, LDAB, AFB, LDAFB, IPIV, CMODE, C, INFO, WORK, IWORK )\n\n*  Purpose\n*  =======\n*\n*     SLA_GBRCOND Estimates the Skeel condition number of  op(A) * op2(C)\n*     where op2 is determined by CMODE as follows\n*     CMODE =  1    op2(C) = C\n*     CMODE =  0    op2(C) = I\n*     CMODE = -1    op2(C) = inv(C)\n*     The Skeel condition number  cond(A) = norminf( |inv(A)||A| )\n*     is computed by computing scaling factors R such that\n*     diag(R)*A*op2(C) is row equilibrated and computing the standard\n*     infinity-norm condition number.\n*\n\n*  Arguments\n*  ==========\n*\n*     TRANS   (input) CHARACTER*1\n*     Specifies the form of the system of equations:\n*       = 'N':  A * X = B     (No transpose)\n*       = 'T':  A**T * X = B  (Transpose)\n*       = 'C':  A**H * X = B  (Conjugate Transpose = Transpose)\n*\n*     N       (input) INTEGER\n*     The number of linear equations, i.e., the order of the\n*     matrix A.  N >= 0.\n*\n*     KL      (input) INTEGER\n*     The number of subdiagonals within the band of A.  KL >= 0.\n*\n*     KU      (input) INTEGER\n*     The number of superdiagonals within the band of A.  KU >= 0.\n*\n*     AB      (input) REAL array, dimension (LDAB,N)\n*     On entry, the matrix A in band storage, in rows 1 to KL+KU+1.\n*     The j-th column of A is stored in the j-th column of the\n*     array AB as follows:\n*     AB(KU+1+i-j,j) = A(i,j) for max(1,j-KU)<=i<=min(N,j+kl)\n*\n*     LDAB    (input) INTEGER\n*     The leading dimension of the array AB.  LDAB >= KL+KU+1.\n*\n*     AFB     (input) REAL array, dimension (LDAFB,N)\n*     Details of the LU factorization of the band matrix A, as\n*     computed by SGBTRF.  U is stored as an upper triangular\n*     band matrix with KL+KU superdiagonals in rows 1 to KL+KU+1,\n*     and the multipliers used during the factorization are stored\n*     in rows KL+KU+2 to 2*KL+KU+1.\n*\n*     LDAFB   (input) INTEGER\n*     The leading dimension of the array AFB.  LDAFB >= 2*KL+KU+1.\n*\n*     IPIV    (input) INTEGER array, dimension (N)\n*     The pivot indices from the factorization A = P*L*U\n*     as computed by SGBTRF; row i of the matrix was interchanged\n*     with row IPIV(i).\n*\n*     CMODE   (input) INTEGER\n*     Determines op2(C) in the formula op(A) * op2(C) as follows:\n*     CMODE =  1    op2(C) = C\n*     CMODE =  0    op2(C) = I\n*     CMODE = -1    op2(C) = inv(C)\n*\n*     C       (input) REAL array, dimension (N)\n*     The vector C in the formula op(A) * op2(C).\n*\n*     INFO    (output) INTEGER\n*       = 0:  Successful exit.\n*     i > 0:  The ith argument is invalid.\n*\n*     WORK    (input) REAL array, dimension (5*N).\n*     Workspace.\n*\n*     IWORK   (input) INTEGER array, dimension (N).\n*     Workspace.\n*\n\n*  =====================================================================\n*\n*     .. Local Scalars ..\n      LOGICAL            NOTRANS\n      INTEGER            KASE, I, J, KD, KE\n      REAL               AINVNM, TMP\n*     ..\n*     .. Local Arrays ..\n      INTEGER            ISAVE( 3 )\n*     ..\n*     .. External Functions ..\n      LOGICAL            LSAME\n      EXTERNAL           LSAME\n*     ..\n*     .. External Subroutines ..\n      EXTERNAL           SLACN2, SGBTRS, XERBLA\n*     ..\n*     .. Intrinsic Functions ..\n      INTRINSIC          ABS, MAX\n*     ..\n\n");
    return Qnil;
  }
  if (argc != 10)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 10)", argc);
  rb_trans = argv[0];
  rb_kl = argv[1];
  rb_ku = argv[2];
  rb_ab = argv[3];
  rb_afb = argv[4];
  rb_ipiv = argv[5];
  rb_cmode = argv[6];
  rb_c = argv[7];
  rb_work = argv[8];
  rb_iwork = argv[9];

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
  if (NA_TYPE(rb_ab) != NA_SFLOAT)
    rb_ab = na_change_type(rb_ab, NA_SFLOAT);
  ab = NA_PTR_TYPE(rb_ab, real*);
  kl = NUM2INT(rb_kl);
  if (!NA_IsNArray(rb_c))
    rb_raise(rb_eArgError, "c (8th argument) must be NArray");
  if (NA_RANK(rb_c) != 1)
    rb_raise(rb_eArgError, "rank of c (8th argument) must be %d", 1);
  if (NA_SHAPE0(rb_c) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of c must be the same as shape 0 of ipiv");
  if (NA_TYPE(rb_c) != NA_SFLOAT)
    rb_c = na_change_type(rb_c, NA_SFLOAT);
  c = NA_PTR_TYPE(rb_c, real*);
  if (!NA_IsNArray(rb_iwork))
    rb_raise(rb_eArgError, "iwork (10th argument) must be NArray");
  if (NA_RANK(rb_iwork) != 1)
    rb_raise(rb_eArgError, "rank of iwork (10th argument) must be %d", 1);
  if (NA_SHAPE0(rb_iwork) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of iwork must be the same as shape 0 of ipiv");
  if (NA_TYPE(rb_iwork) != NA_LINT)
    rb_iwork = na_change_type(rb_iwork, NA_LINT);
  iwork = NA_PTR_TYPE(rb_iwork, integer*);
  ku = NUM2INT(rb_ku);
  cmode = NUM2INT(rb_cmode);
  trans = StringValueCStr(rb_trans)[0];
  if (!NA_IsNArray(rb_afb))
    rb_raise(rb_eArgError, "afb (5th argument) must be NArray");
  if (NA_RANK(rb_afb) != 2)
    rb_raise(rb_eArgError, "rank of afb (5th argument) must be %d", 2);
  if (NA_SHAPE1(rb_afb) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of afb must be the same as shape 0 of ipiv");
  ldafb = NA_SHAPE0(rb_afb);
  if (NA_TYPE(rb_afb) != NA_SFLOAT)
    rb_afb = na_change_type(rb_afb, NA_SFLOAT);
  afb = NA_PTR_TYPE(rb_afb, real*);
  if (!NA_IsNArray(rb_work))
    rb_raise(rb_eArgError, "work (9th argument) must be NArray");
  if (NA_RANK(rb_work) != 1)
    rb_raise(rb_eArgError, "rank of work (9th argument) must be %d", 1);
  if (NA_SHAPE0(rb_work) != (5*n))
    rb_raise(rb_eRuntimeError, "shape 0 of work must be %d", 5*n);
  if (NA_TYPE(rb_work) != NA_SFLOAT)
    rb_work = na_change_type(rb_work, NA_SFLOAT);
  work = NA_PTR_TYPE(rb_work, real*);

  __out__ = sla_gbrcond_(&trans, &n, &kl, &ku, ab, &ldab, afb, &ldafb, ipiv, &cmode, c, &info, work, iwork);

  rb_info = INT2NUM(info);
  rb___out__ = rb_float_new((double)__out__);
  return rb_ary_new3(2, rb_info, rb___out__);
}

void
init_lapack_sla_gbrcond(VALUE mLapack){
  rb_define_module_function(mLapack, "sla_gbrcond", rb_sla_gbrcond, -1);
}

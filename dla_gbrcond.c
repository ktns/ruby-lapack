#include "rb_lapack.h"

extern doublereal dla_gbrcond_(char *trans, integer *n, integer *kl, integer *ku, doublereal *ab, integer *ldab, doublereal *afb, integer *ldafb, integer *ipiv, integer *cmode, doublereal *c, integer *info, doublereal *work, integer *iwork);

static VALUE sHelp, sUsage;

static VALUE
rblapack_dla_gbrcond(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_trans;
  char trans; 
  VALUE rblapack_kl;
  integer kl; 
  VALUE rblapack_ku;
  integer ku; 
  VALUE rblapack_ab;
  doublereal *ab; 
  VALUE rblapack_afb;
  doublereal *afb; 
  VALUE rblapack_ipiv;
  integer *ipiv; 
  VALUE rblapack_cmode;
  integer cmode; 
  VALUE rblapack_c;
  doublereal *c; 
  VALUE rblapack_work;
  doublereal *work; 
  VALUE rblapack_iwork;
  integer *iwork; 
  VALUE rblapack_info;
  integer info; 
  VALUE rblapack___out__;
  doublereal __out__; 

  integer ldab;
  integer n;
  integer ldafb;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  info, __out__ = NumRu::Lapack.dla_gbrcond( trans, kl, ku, ab, afb, ipiv, cmode, c, work, iwork, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      DOUBLE PRECISION FUNCTION DLA_GBRCOND( TRANS, N, KL, KU, AB, LDAB, AFB, LDAFB, IPIV, CMODE, C, INFO, WORK, IWORK )\n\n*  Purpose\n*  =======\n*\n*     DLA_GBRCOND Estimates the Skeel condition number of  op(A) * op2(C)\n*     where op2 is determined by CMODE as follows\n*     CMODE =  1    op2(C) = C\n*     CMODE =  0    op2(C) = I\n*     CMODE = -1    op2(C) = inv(C)\n*     The Skeel condition number  cond(A) = norminf( |inv(A)||A| )\n*     is computed by computing scaling factors R such that\n*     diag(R)*A*op2(C) is row equilibrated and computing the standard\n*     infinity-norm condition number.\n*\n\n*  Arguments\n*  =========\n*\n*     TRANS   (input) CHARACTER*1\n*     Specifies the form of the system of equations:\n*       = 'N':  A * X = B     (No transpose)\n*       = 'T':  A**T * X = B  (Transpose)\n*       = 'C':  A**H * X = B  (Conjugate Transpose = Transpose)\n*\n*     N       (input) INTEGER\n*     The number of linear equations, i.e., the order of the\n*     matrix A.  N >= 0.\n*\n*     KL      (input) INTEGER\n*     The number of subdiagonals within the band of A.  KL >= 0.\n*\n*     KU      (input) INTEGER\n*     The number of superdiagonals within the band of A.  KU >= 0.\n*\n*     AB      (input) DOUBLE PRECISION array, dimension (LDAB,N)\n*     On entry, the matrix A in band storage, in rows 1 to KL+KU+1.\n*     The j-th column of A is stored in the j-th column of the\n*     array AB as follows:\n*     AB(KU+1+i-j,j) = A(i,j) for max(1,j-KU)<=i<=min(N,j+kl)\n*\n*     LDAB    (input) INTEGER\n*     The leading dimension of the array AB.  LDAB >= KL+KU+1.\n*\n*     AFB     (input) DOUBLE PRECISION array, dimension (LDAFB,N)\n*     Details of the LU factorization of the band matrix A, as\n*     computed by DGBTRF.  U is stored as an upper triangular\n*     band matrix with KL+KU superdiagonals in rows 1 to KL+KU+1,\n*     and the multipliers used during the factorization are stored\n*     in rows KL+KU+2 to 2*KL+KU+1.\n*\n*     LDAFB   (input) INTEGER\n*     The leading dimension of the array AFB.  LDAFB >= 2*KL+KU+1.\n*\n*     IPIV    (input) INTEGER array, dimension (N)\n*     The pivot indices from the factorization A = P*L*U\n*     as computed by DGBTRF; row i of the matrix was interchanged\n*     with row IPIV(i).\n*\n*     CMODE   (input) INTEGER\n*     Determines op2(C) in the formula op(A) * op2(C) as follows:\n*     CMODE =  1    op2(C) = C\n*     CMODE =  0    op2(C) = I\n*     CMODE = -1    op2(C) = inv(C)\n*\n*     C       (input) DOUBLE PRECISION array, dimension (N)\n*     The vector C in the formula op(A) * op2(C).\n*\n*     INFO    (output) INTEGER\n*       = 0:  Successful exit.\n*     i > 0:  The ith argument is invalid.\n*\n*     WORK    (input) DOUBLE PRECISION array, dimension (5*N).\n*     Workspace.\n*\n*     IWORK   (input) INTEGER array, dimension (N).\n*     Workspace.\n*\n\n*  =====================================================================\n*\n*     .. Local Scalars ..\n      LOGICAL            NOTRANS\n      INTEGER            KASE, I, J, KD, KE\n      DOUBLE PRECISION   AINVNM, TMP\n*     ..\n*     .. Local Arrays ..\n      INTEGER            ISAVE( 3 )\n*     ..\n*     .. External Functions ..\n      LOGICAL            LSAME\n      EXTERNAL           LSAME\n*     ..\n*     .. External Subroutines ..\n      EXTERNAL           DLACN2, DGBTRS, XERBLA\n*     ..\n*     .. Intrinsic Functions ..\n      INTRINSIC          ABS, MAX\n*     ..\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  info, __out__ = NumRu::Lapack.dla_gbrcond( trans, kl, ku, ab, afb, ipiv, cmode, c, work, iwork, [:usage => usage, :help => help])\n");
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
  rblapack_cmode = argv[6];
  rblapack_c = argv[7];
  rblapack_work = argv[8];
  rblapack_iwork = argv[9];
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
  if (NA_TYPE(rblapack_ab) != NA_DFLOAT)
    rblapack_ab = na_change_type(rblapack_ab, NA_DFLOAT);
  ab = NA_PTR_TYPE(rblapack_ab, doublereal*);
  kl = NUM2INT(rblapack_kl);
  if (!NA_IsNArray(rblapack_c))
    rb_raise(rb_eArgError, "c (8th argument) must be NArray");
  if (NA_RANK(rblapack_c) != 1)
    rb_raise(rb_eArgError, "rank of c (8th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_c) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of c must be the same as shape 0 of ipiv");
  if (NA_TYPE(rblapack_c) != NA_DFLOAT)
    rblapack_c = na_change_type(rblapack_c, NA_DFLOAT);
  c = NA_PTR_TYPE(rblapack_c, doublereal*);
  if (!NA_IsNArray(rblapack_iwork))
    rb_raise(rb_eArgError, "iwork (10th argument) must be NArray");
  if (NA_RANK(rblapack_iwork) != 1)
    rb_raise(rb_eArgError, "rank of iwork (10th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_iwork) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of iwork must be the same as shape 0 of ipiv");
  if (NA_TYPE(rblapack_iwork) != NA_LINT)
    rblapack_iwork = na_change_type(rblapack_iwork, NA_LINT);
  iwork = NA_PTR_TYPE(rblapack_iwork, integer*);
  ku = NUM2INT(rblapack_ku);
  cmode = NUM2INT(rblapack_cmode);
  trans = StringValueCStr(rblapack_trans)[0];
  if (!NA_IsNArray(rblapack_afb))
    rb_raise(rb_eArgError, "afb (5th argument) must be NArray");
  if (NA_RANK(rblapack_afb) != 2)
    rb_raise(rb_eArgError, "rank of afb (5th argument) must be %d", 2);
  if (NA_SHAPE1(rblapack_afb) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of afb must be the same as shape 0 of ipiv");
  ldafb = NA_SHAPE0(rblapack_afb);
  if (NA_TYPE(rblapack_afb) != NA_DFLOAT)
    rblapack_afb = na_change_type(rblapack_afb, NA_DFLOAT);
  afb = NA_PTR_TYPE(rblapack_afb, doublereal*);
  if (!NA_IsNArray(rblapack_work))
    rb_raise(rb_eArgError, "work (9th argument) must be NArray");
  if (NA_RANK(rblapack_work) != 1)
    rb_raise(rb_eArgError, "rank of work (9th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_work) != (5*n))
    rb_raise(rb_eRuntimeError, "shape 0 of work must be %d", 5*n);
  if (NA_TYPE(rblapack_work) != NA_DFLOAT)
    rblapack_work = na_change_type(rblapack_work, NA_DFLOAT);
  work = NA_PTR_TYPE(rblapack_work, doublereal*);

  __out__ = dla_gbrcond_(&trans, &n, &kl, &ku, ab, &ldab, afb, &ldafb, ipiv, &cmode, c, &info, work, iwork);

  rblapack_info = INT2NUM(info);
  rblapack___out__ = rb_float_new((double)__out__);
  return rb_ary_new3(2, rblapack_info, rblapack___out__);
}

void
init_lapack_dla_gbrcond(VALUE mLapack){
  rb_define_module_function(mLapack, "dla_gbrcond", rblapack_dla_gbrcond, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}

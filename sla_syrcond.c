#include "rb_lapack.h"

extern VOID sla_syrcond_(real *__out__, char *uplo, integer *n, real *a, integer *lda, real *af, integer *ldaf, integer *ipiv, integer *cmode, real *c, integer *info, real *work, integer *iwork);
static VALUE
rb_sla_syrcond(int argc, VALUE *argv, VALUE self){
  VALUE rb_uplo;
  char uplo; 
  VALUE rb_a;
  real *a; 
  VALUE rb_af;
  real *af; 
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

  integer lda;
  integer n;
  integer ldaf;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  info, __out__ = NumRu::Lapack.sla_syrcond( uplo, a, af, ipiv, cmode, c, work, iwork)\n    or\n  NumRu::Lapack.sla_syrcond  # print help\n\n\nFORTRAN MANUAL\n      REAL FUNCTION SLA_SYRCOND( UPLO, N, A, LDA, AF, LDAF, IPIV, CMODE, C, INFO, WORK, IWORK )\n\n*  Purpose\n*  =======\n*\n*     SLA_SYRCOND estimates the Skeel condition number of  op(A) * op2(C)\n*     where op2 is determined by CMODE as follows\n*     CMODE =  1    op2(C) = C\n*     CMODE =  0    op2(C) = I\n*     CMODE = -1    op2(C) = inv(C)\n*     The Skeel condition number cond(A) = norminf( |inv(A)||A| )\n*     is computed by computing scaling factors R such that\n*     diag(R)*A*op2(C) is row equilibrated and computing the standard\n*     infinity-norm condition number.\n*\n\n*  Arguments\n*  ==========\n*\n*     UPLO    (input) CHARACTER*1\n*       = 'U':  Upper triangle of A is stored;\n*       = 'L':  Lower triangle of A is stored.\n*\n*     N       (input) INTEGER\n*     The number of linear equations, i.e., the order of the\n*     matrix A.  N >= 0.\n*\n*     A       (input) REAL array, dimension (LDA,N)\n*     On entry, the N-by-N matrix A.\n*\n*     LDA     (input) INTEGER\n*     The leading dimension of the array A.  LDA >= max(1,N).\n*\n*     AF      (input) REAL array, dimension (LDAF,N)\n*     The block diagonal matrix D and the multipliers used to\n*     obtain the factor U or L as computed by SSYTRF.\n*\n*     LDAF    (input) INTEGER\n*     The leading dimension of the array AF.  LDAF >= max(1,N).\n*\n*     IPIV    (input) INTEGER array, dimension (N)\n*     Details of the interchanges and the block structure of D\n*     as determined by SSYTRF.\n*\n*     CMODE   (input) INTEGER\n*     Determines op2(C) in the formula op(A) * op2(C) as follows:\n*     CMODE =  1    op2(C) = C\n*     CMODE =  0    op2(C) = I\n*     CMODE = -1    op2(C) = inv(C)\n*\n*     C       (input) REAL array, dimension (N)\n*     The vector C in the formula op(A) * op2(C).\n*\n*     INFO    (output) INTEGER\n*       = 0:  Successful exit.\n*     i > 0:  The ith argument is invalid.\n*\n*     WORK    (input) REAL array, dimension (3*N).\n*     Workspace.\n*\n*     IWORK   (input) INTEGER array, dimension (N).\n*     Workspace.\n*\n\n*  =====================================================================\n*\n*     .. Local Scalars ..\n      CHARACTER          NORMIN\n      INTEGER            KASE, I, J\n      REAL               AINVNM, SMLNUM, TMP\n      LOGICAL            UP\n*     ..\n*     .. Local Arrays ..\n      INTEGER            ISAVE( 3 )\n*     ..\n*     .. External Functions ..\n      LOGICAL            LSAME\n      INTEGER            ISAMAX\n      REAL               SLAMCH\n      EXTERNAL           LSAME, ISAMAX, SLAMCH\n*     ..\n*     .. External Subroutines ..\n      EXTERNAL           SLACN2, SLATRS, SRSCL, XERBLA, SSYTRS\n*     ..\n*     .. Intrinsic Functions ..\n      INTRINSIC          ABS, MAX\n*     ..\n\n");
    return Qnil;
  }
  if (argc != 8)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 8)", argc);
  rb_uplo = argv[0];
  rb_a = argv[1];
  rb_af = argv[2];
  rb_ipiv = argv[3];
  rb_cmode = argv[4];
  rb_c = argv[5];
  rb_work = argv[6];
  rb_iwork = argv[7];

  uplo = StringValueCStr(rb_uplo)[0];
  cmode = NUM2INT(rb_cmode);
  if (!NA_IsNArray(rb_a))
    rb_raise(rb_eArgError, "a (2th argument) must be NArray");
  if (NA_RANK(rb_a) != 2)
    rb_raise(rb_eArgError, "rank of a (2th argument) must be %d", 2);
  lda = NA_SHAPE0(rb_a);
  n = NA_SHAPE1(rb_a);
  if (NA_TYPE(rb_a) != NA_SFLOAT)
    rb_a = na_change_type(rb_a, NA_SFLOAT);
  a = NA_PTR_TYPE(rb_a, real*);
  if (!NA_IsNArray(rb_af))
    rb_raise(rb_eArgError, "af (3th argument) must be NArray");
  if (NA_RANK(rb_af) != 2)
    rb_raise(rb_eArgError, "rank of af (3th argument) must be %d", 2);
  ldaf = NA_SHAPE0(rb_af);
  if (NA_SHAPE1(rb_af) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of af must be the same as shape 1 of a");
  if (NA_TYPE(rb_af) != NA_SFLOAT)
    rb_af = na_change_type(rb_af, NA_SFLOAT);
  af = NA_PTR_TYPE(rb_af, real*);
  if (!NA_IsNArray(rb_ipiv))
    rb_raise(rb_eArgError, "ipiv (4th argument) must be NArray");
  if (NA_RANK(rb_ipiv) != 1)
    rb_raise(rb_eArgError, "rank of ipiv (4th argument) must be %d", 1);
  if (NA_SHAPE0(rb_ipiv) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of ipiv must be the same as shape 1 of a");
  if (NA_TYPE(rb_ipiv) != NA_LINT)
    rb_ipiv = na_change_type(rb_ipiv, NA_LINT);
  ipiv = NA_PTR_TYPE(rb_ipiv, integer*);
  if (!NA_IsNArray(rb_c))
    rb_raise(rb_eArgError, "c (6th argument) must be NArray");
  if (NA_RANK(rb_c) != 1)
    rb_raise(rb_eArgError, "rank of c (6th argument) must be %d", 1);
  if (NA_SHAPE0(rb_c) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of c must be the same as shape 1 of a");
  if (NA_TYPE(rb_c) != NA_SFLOAT)
    rb_c = na_change_type(rb_c, NA_SFLOAT);
  c = NA_PTR_TYPE(rb_c, real*);
  if (!NA_IsNArray(rb_work))
    rb_raise(rb_eArgError, "work (7th argument) must be NArray");
  if (NA_RANK(rb_work) != 1)
    rb_raise(rb_eArgError, "rank of work (7th argument) must be %d", 1);
  if (NA_SHAPE0(rb_work) != (3*n))
    rb_raise(rb_eRuntimeError, "shape 0 of work must be %d", 3*n);
  if (NA_TYPE(rb_work) != NA_SFLOAT)
    rb_work = na_change_type(rb_work, NA_SFLOAT);
  work = NA_PTR_TYPE(rb_work, real*);
  if (!NA_IsNArray(rb_iwork))
    rb_raise(rb_eArgError, "iwork (8th argument) must be NArray");
  if (NA_RANK(rb_iwork) != 1)
    rb_raise(rb_eArgError, "rank of iwork (8th argument) must be %d", 1);
  if (NA_SHAPE0(rb_iwork) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of iwork must be the same as shape 1 of a");
  if (NA_TYPE(rb_iwork) != NA_LINT)
    rb_iwork = na_change_type(rb_iwork, NA_LINT);
  iwork = NA_PTR_TYPE(rb_iwork, integer*);

  sla_syrcond_(&__out__, &uplo, &n, a, &lda, af, &ldaf, ipiv, &cmode, c, &info, work, iwork);

  rb_info = INT2NUM(info);
  rb___out__ = rb_float_new((double)__out__);
  return rb_ary_new3(2, rb_info, rb___out__);
}

void
init_lapack_sla_syrcond(VALUE mLapack){
  rb_define_module_function(mLapack, "sla_syrcond", rb_sla_syrcond, -1);
}

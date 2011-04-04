#include "rb_lapack.h"

extern real cla_hercond_c_(char *uplo, integer *n, complex *a, integer *lda, complex *af, integer *ldaf, integer *ipiv, real *c, logical *capply, integer *info, complex *work, real *rwork);

static VALUE
rb_cla_hercond_c(int argc, VALUE *argv, VALUE self){
  VALUE rb_uplo;
  char uplo; 
  VALUE rb_a;
  complex *a; 
  VALUE rb_af;
  complex *af; 
  VALUE rb_ipiv;
  integer *ipiv; 
  VALUE rb_c;
  real *c; 
  VALUE rb_capply;
  logical capply; 
  VALUE rb_work;
  complex *work; 
  VALUE rb_rwork;
  real *rwork; 
  VALUE rb_info;
  integer info; 
  VALUE rb___out__;
  real __out__; 

  integer lda;
  integer n;
  integer ldaf;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  info, __out__ = NumRu::Lapack.cla_hercond_c( uplo, a, af, ipiv, c, capply, work, rwork)\n    or\n  NumRu::Lapack.cla_hercond_c  # print help\n\n\nFORTRAN MANUAL\n      REAL FUNCTION CLA_HERCOND_C( UPLO, N, A, LDA, AF, LDAF, IPIV, C, CAPPLY, INFO, WORK, RWORK )\n\n*  Purpose\n*  =======\n*\n*     CLA_HERCOND_C computes the infinity norm condition number of\n*     op(A) * inv(diag(C)) where C is a REAL vector.\n*\n\n*  Arguments\n*  =========\n*\n*     UPLO    (input) CHARACTER*1\n*       = 'U':  Upper triangle of A is stored;\n*       = 'L':  Lower triangle of A is stored.\n*\n*     N       (input) INTEGER\n*     The number of linear equations, i.e., the order of the\n*     matrix A.  N >= 0.\n*\n*     A       (input) COMPLEX array, dimension (LDA,N)\n*     On entry, the N-by-N matrix A\n*\n*     LDA     (input) INTEGER\n*     The leading dimension of the array A.  LDA >= max(1,N).\n*\n*     AF      (input) COMPLEX array, dimension (LDAF,N)\n*     The block diagonal matrix D and the multipliers used to\n*     obtain the factor U or L as computed by CHETRF.\n*\n*     LDAF    (input) INTEGER\n*     The leading dimension of the array AF.  LDAF >= max(1,N).\n*\n*     IPIV    (input) INTEGER array, dimension (N)\n*     Details of the interchanges and the block structure of D\n*     as determined by CHETRF.\n*\n*     C       (input) REAL array, dimension (N)\n*     The vector C in the formula op(A) * inv(diag(C)).\n*\n*     CAPPLY  (input) LOGICAL\n*     If .TRUE. then access the vector C in the formula above.\n*\n*     INFO    (output) INTEGER\n*       = 0:  Successful exit.\n*     i > 0:  The ith argument is invalid.\n*\n*     WORK    (input) COMPLEX array, dimension (2*N).\n*     Workspace.\n*\n*     RWORK   (input) REAL array, dimension (N).\n*     Workspace.\n*\n\n*  =====================================================================\n*\n*     .. Local Scalars ..\n      INTEGER            KASE, I, J\n      REAL               AINVNM, ANORM, TMP\n      LOGICAL            UP\n      COMPLEX            ZDUM\n*     ..\n*     .. Local Arrays ..\n      INTEGER            ISAVE( 3 )\n*     ..\n*     .. External Functions ..\n      LOGICAL            LSAME\n      EXTERNAL           LSAME\n*     ..\n*     .. External Subroutines ..\n      EXTERNAL           CLACN2, CHETRS, XERBLA\n*     ..\n*     .. Intrinsic Functions ..\n      INTRINSIC          ABS, MAX\n*     ..\n*     .. Statement Functions ..\n      REAL               CABS1\n*     ..\n*     .. Statement Function Definitions ..\n      CABS1( ZDUM ) = ABS( REAL( ZDUM ) ) + ABS( AIMAG( ZDUM ) )\n*     ..\n\n");
    return Qnil;
  }
  if (argc != 8)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 8)", argc);
  rb_uplo = argv[0];
  rb_a = argv[1];
  rb_af = argv[2];
  rb_ipiv = argv[3];
  rb_c = argv[4];
  rb_capply = argv[5];
  rb_work = argv[6];
  rb_rwork = argv[7];

  if (!NA_IsNArray(rb_ipiv))
    rb_raise(rb_eArgError, "ipiv (4th argument) must be NArray");
  if (NA_RANK(rb_ipiv) != 1)
    rb_raise(rb_eArgError, "rank of ipiv (4th argument) must be %d", 1);
  n = NA_SHAPE0(rb_ipiv);
  if (NA_TYPE(rb_ipiv) != NA_LINT)
    rb_ipiv = na_change_type(rb_ipiv, NA_LINT);
  ipiv = NA_PTR_TYPE(rb_ipiv, integer*);
  if (!NA_IsNArray(rb_a))
    rb_raise(rb_eArgError, "a (2th argument) must be NArray");
  if (NA_RANK(rb_a) != 2)
    rb_raise(rb_eArgError, "rank of a (2th argument) must be %d", 2);
  if (NA_SHAPE1(rb_a) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of a must be the same as shape 0 of ipiv");
  lda = NA_SHAPE0(rb_a);
  if (NA_TYPE(rb_a) != NA_SCOMPLEX)
    rb_a = na_change_type(rb_a, NA_SCOMPLEX);
  a = NA_PTR_TYPE(rb_a, complex*);
  if (!NA_IsNArray(rb_rwork))
    rb_raise(rb_eArgError, "rwork (8th argument) must be NArray");
  if (NA_RANK(rb_rwork) != 1)
    rb_raise(rb_eArgError, "rank of rwork (8th argument) must be %d", 1);
  if (NA_SHAPE0(rb_rwork) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of rwork must be the same as shape 0 of ipiv");
  if (NA_TYPE(rb_rwork) != NA_SFLOAT)
    rb_rwork = na_change_type(rb_rwork, NA_SFLOAT);
  rwork = NA_PTR_TYPE(rb_rwork, real*);
  if (!NA_IsNArray(rb_c))
    rb_raise(rb_eArgError, "c (5th argument) must be NArray");
  if (NA_RANK(rb_c) != 1)
    rb_raise(rb_eArgError, "rank of c (5th argument) must be %d", 1);
  if (NA_SHAPE0(rb_c) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of c must be the same as shape 0 of ipiv");
  if (NA_TYPE(rb_c) != NA_SFLOAT)
    rb_c = na_change_type(rb_c, NA_SFLOAT);
  c = NA_PTR_TYPE(rb_c, real*);
  uplo = StringValueCStr(rb_uplo)[0];
  if (!NA_IsNArray(rb_af))
    rb_raise(rb_eArgError, "af (3th argument) must be NArray");
  if (NA_RANK(rb_af) != 2)
    rb_raise(rb_eArgError, "rank of af (3th argument) must be %d", 2);
  if (NA_SHAPE1(rb_af) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of af must be the same as shape 0 of ipiv");
  ldaf = NA_SHAPE0(rb_af);
  if (NA_TYPE(rb_af) != NA_SCOMPLEX)
    rb_af = na_change_type(rb_af, NA_SCOMPLEX);
  af = NA_PTR_TYPE(rb_af, complex*);
  capply = (rb_capply == Qtrue);
  if (!NA_IsNArray(rb_work))
    rb_raise(rb_eArgError, "work (7th argument) must be NArray");
  if (NA_RANK(rb_work) != 1)
    rb_raise(rb_eArgError, "rank of work (7th argument) must be %d", 1);
  if (NA_SHAPE0(rb_work) != (2*n))
    rb_raise(rb_eRuntimeError, "shape 0 of work must be %d", 2*n);
  if (NA_TYPE(rb_work) != NA_SCOMPLEX)
    rb_work = na_change_type(rb_work, NA_SCOMPLEX);
  work = NA_PTR_TYPE(rb_work, complex*);

  __out__ = cla_hercond_c_(&uplo, &n, a, &lda, af, &ldaf, ipiv, c, &capply, &info, work, rwork);

  rb_info = INT2NUM(info);
  rb___out__ = rb_float_new((double)__out__);
  return rb_ary_new3(2, rb_info, rb___out__);
}

void
init_lapack_cla_hercond_c(VALUE mLapack){
  rb_define_module_function(mLapack, "cla_hercond_c", rb_cla_hercond_c, -1);
}

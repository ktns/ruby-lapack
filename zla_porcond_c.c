#include "rb_lapack.h"

extern VOID zla_porcond_c_(doublereal *__out__, char *uplo, integer *n, doublecomplex *a, integer *lda, doublecomplex *af, integer *ldaf, doublereal *c, logical *capply, integer *info, doublecomplex *work, doublereal *rwork);
static VALUE
rb_zla_porcond_c(int argc, VALUE *argv, VALUE self){
  VALUE rb_uplo;
  char uplo; 
  VALUE rb_a;
  doublecomplex *a; 
  VALUE rb_af;
  doublecomplex *af; 
  VALUE rb_c;
  doublereal *c; 
  VALUE rb_capply;
  logical capply; 
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
    printf("%s\n", "USAGE:\n  info, __out__ = NumRu::Lapack.zla_porcond_c( uplo, a, af, c, capply, work, rwork)\n    or\n  NumRu::Lapack.zla_porcond_c  # print help\n\n\nFORTRAN MANUAL\n      DOUBLE PRECISION FUNCTION ZLA_PORCOND_C( UPLO, N, A, LDA, AF,  LDAF, C, CAPPLY, INFO, WORK, RWORK )\n\n*  Purpose\n*  =======\n*\n*     ZLA_PORCOND_C Computes the infinity norm condition number of\n*     op(A) * inv(diag(C)) where C is a DOUBLE PRECISION vector\n*\n\n*  Arguments\n*  =========\n*\n*     UPLO    (input) CHARACTER*1\n*       = 'U':  Upper triangle of A is stored;\n*       = 'L':  Lower triangle of A is stored.\n*\n*     N       (input) INTEGER\n*     The number of linear equations, i.e., the order of the\n*     matrix A.  N >= 0.\n*\n*     A       (input) COMPLEX*16 array, dimension (LDA,N)\n*     On entry, the N-by-N matrix A\n*\n*     LDA     (input) INTEGER\n*     The leading dimension of the array A.  LDA >= max(1,N).\n*\n*     AF      (input) COMPLEX*16 array, dimension (LDAF,N)\n*     The triangular factor U or L from the Cholesky factorization\n*     A = U**T*U or A = L*L**T, as computed by ZPOTRF.\n*\n*     LDAF    (input) INTEGER\n*     The leading dimension of the array AF.  LDAF >= max(1,N).\n*\n*     C       (input) DOUBLE PRECISION array, dimension (N)\n*     The vector C in the formula op(A) * inv(diag(C)).\n*\n*     CAPPLY  (input) LOGICAL\n*     If .TRUE. then access the vector C in the formula above.\n*\n*     INFO    (output) INTEGER\n*       = 0:  Successful exit.\n*     i > 0:  The ith argument is invalid.\n*\n*     WORK    (input) COMPLEX*16 array, dimension (2*N).\n*     Workspace.\n*\n*     RWORK   (input) DOUBLE PRECISION array, dimension (N).\n*     Workspace.\n*\n\n*  =====================================================================\n*\n*     .. Local Scalars ..\n      INTEGER            KASE\n      DOUBLE PRECISION   AINVNM, ANORM, TMP\n      INTEGER            I, J\n      LOGICAL            UP\n      COMPLEX*16         ZDUM\n*     ..\n*     .. Local Arrays ..\n      INTEGER            ISAVE( 3 )\n*     ..\n*     .. External Functions ..\n      LOGICAL            LSAME\n      EXTERNAL           LSAME\n*     ..\n*     .. External Subroutines ..\n      EXTERNAL           ZLACN2, ZPOTRS, XERBLA\n*     ..\n*     .. Intrinsic Functions ..\n      INTRINSIC          ABS, MAX, REAL, DIMAG\n*     ..\n*     .. Statement Functions ..\n      DOUBLE PRECISION CABS1\n*     ..\n*     .. Statement Function Definitions ..\n      CABS1( ZDUM ) = ABS( DBLE( ZDUM ) ) + ABS( DIMAG( ZDUM ) )\n*     ..\n\n");
    return Qnil;
  }
  if (argc != 7)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 7)", argc);
  rb_uplo = argv[0];
  rb_a = argv[1];
  rb_af = argv[2];
  rb_c = argv[3];
  rb_capply = argv[4];
  rb_work = argv[5];
  rb_rwork = argv[6];

  uplo = StringValueCStr(rb_uplo)[0];
  capply = (rb_capply == Qtrue);
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
  if (!NA_IsNArray(rb_c))
    rb_raise(rb_eArgError, "c (4th argument) must be NArray");
  if (NA_RANK(rb_c) != 1)
    rb_raise(rb_eArgError, "rank of c (4th argument) must be %d", 1);
  if (NA_SHAPE0(rb_c) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of c must be the same as shape 1 of a");
  if (NA_TYPE(rb_c) != NA_DFLOAT)
    rb_c = na_change_type(rb_c, NA_DFLOAT);
  c = NA_PTR_TYPE(rb_c, doublereal*);
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

  zla_porcond_c_(&__out__, &uplo, &n, a, &lda, af, &ldaf, c, &capply, &info, work, rwork);

  rb_info = INT2NUM(info);
  rb___out__ = rb_float_new((double)__out__);
  return rb_ary_new3(2, rb_info, rb___out__);
}

void
init_lapack_zla_porcond_c(VALUE mLapack){
  rb_define_module_function(mLapack, "zla_porcond_c", rb_zla_porcond_c, -1);
}

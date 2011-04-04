#include "rb_lapack.h"

extern doublereal dla_syrpvgrw_(char *uplo, integer *n, integer *info, doublereal *a, integer *lda, doublereal *af, integer *ldaf, integer *ipiv, doublereal *work);

static VALUE
rb_dla_syrpvgrw(int argc, VALUE *argv, VALUE self){
  VALUE rb_uplo;
  char uplo; 
  VALUE rb_info;
  integer info; 
  VALUE rb_a;
  doublereal *a; 
  VALUE rb_af;
  doublereal *af; 
  VALUE rb_ipiv;
  integer *ipiv; 
  VALUE rb_work;
  doublereal *work; 
  VALUE rb___out__;
  doublereal __out__; 

  integer lda;
  integer n;
  integer ldaf;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  __out__ = NumRu::Lapack.dla_syrpvgrw( uplo, info, a, af, ipiv, work)\n    or\n  NumRu::Lapack.dla_syrpvgrw  # print help\n\n\nFORTRAN MANUAL\n      DOUBLE PRECISION FUNCTION DLA_SYRPVGRW( UPLO, N, INFO, A, LDA, AF, LDAF, IPIV, WORK )\n\n*  Purpose\n*  =======\n* \n*  DLA_SYRPVGRW computes the reciprocal pivot growth factor\n*  norm(A)/norm(U). The \"max absolute element\" norm is used. If this is\n*  much less than 1, the stability of the LU factorization of the\n*  (equilibrated) matrix A could be poor. This also means that the\n*  solution X, estimated condition numbers, and error bounds could be\n*  unreliable.\n*\n\n*  Arguments\n*  =========\n*\n*     UPLO    (input) CHARACTER*1\n*       = 'U':  Upper triangle of A is stored;\n*       = 'L':  Lower triangle of A is stored.\n*\n*     N       (input) INTEGER\n*     The number of linear equations, i.e., the order of the\n*     matrix A.  N >= 0.\n*\n*     INFO    (input) INTEGER\n*     The value of INFO returned from DSYTRF, .i.e., the pivot in\n*     column INFO is exactly 0.\n*\n*     NCOLS   (input) INTEGER\n*     The number of columns of the matrix A. NCOLS >= 0.\n*\n*     A       (input) DOUBLE PRECISION array, dimension (LDA,N)\n*     On entry, the N-by-N matrix A.\n*\n*     LDA     (input) INTEGER\n*     The leading dimension of the array A.  LDA >= max(1,N).\n*\n*     AF      (input) DOUBLE PRECISION array, dimension (LDAF,N)\n*     The block diagonal matrix D and the multipliers used to\n*     obtain the factor U or L as computed by DSYTRF.\n*\n*     LDAF    (input) INTEGER\n*     The leading dimension of the array AF.  LDAF >= max(1,N).\n*\n*     IPIV    (input) INTEGER array, dimension (N)\n*     Details of the interchanges and the block structure of D\n*     as determined by DSYTRF.\n*\n*     WORK    (input) DOUBLE PRECISION array, dimension (2*N)\n*\n\n*  =====================================================================\n*\n*     .. Local Scalars ..\n      INTEGER            NCOLS, I, J, K, KP\n      DOUBLE PRECISION   AMAX, UMAX, RPVGRW, TMP\n      LOGICAL            UPPER\n*     ..\n*     .. Intrinsic Functions ..\n      INTRINSIC          ABS, MAX, MIN\n*     ..\n*     .. External Functions ..\n      EXTERNAL           LSAME, DLASET\n      LOGICAL            LSAME\n*     ..\n\n");
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

  if (!NA_IsNArray(rb_ipiv))
    rb_raise(rb_eArgError, "ipiv (5th argument) must be NArray");
  if (NA_RANK(rb_ipiv) != 1)
    rb_raise(rb_eArgError, "rank of ipiv (5th argument) must be %d", 1);
  n = NA_SHAPE0(rb_ipiv);
  if (NA_TYPE(rb_ipiv) != NA_LINT)
    rb_ipiv = na_change_type(rb_ipiv, NA_LINT);
  ipiv = NA_PTR_TYPE(rb_ipiv, integer*);
  if (!NA_IsNArray(rb_a))
    rb_raise(rb_eArgError, "a (3th argument) must be NArray");
  if (NA_RANK(rb_a) != 2)
    rb_raise(rb_eArgError, "rank of a (3th argument) must be %d", 2);
  if (NA_SHAPE1(rb_a) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of a must be the same as shape 0 of ipiv");
  lda = NA_SHAPE0(rb_a);
  if (NA_TYPE(rb_a) != NA_DFLOAT)
    rb_a = na_change_type(rb_a, NA_DFLOAT);
  a = NA_PTR_TYPE(rb_a, doublereal*);
  if (!NA_IsNArray(rb_af))
    rb_raise(rb_eArgError, "af (4th argument) must be NArray");
  if (NA_RANK(rb_af) != 2)
    rb_raise(rb_eArgError, "rank of af (4th argument) must be %d", 2);
  if (NA_SHAPE1(rb_af) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of af must be the same as shape 0 of ipiv");
  ldaf = NA_SHAPE0(rb_af);
  if (NA_TYPE(rb_af) != NA_DFLOAT)
    rb_af = na_change_type(rb_af, NA_DFLOAT);
  af = NA_PTR_TYPE(rb_af, doublereal*);
  uplo = StringValueCStr(rb_uplo)[0];
  info = NUM2INT(rb_info);
  if (!NA_IsNArray(rb_work))
    rb_raise(rb_eArgError, "work (6th argument) must be NArray");
  if (NA_RANK(rb_work) != 1)
    rb_raise(rb_eArgError, "rank of work (6th argument) must be %d", 1);
  if (NA_SHAPE0(rb_work) != (2*n))
    rb_raise(rb_eRuntimeError, "shape 0 of work must be %d", 2*n);
  if (NA_TYPE(rb_work) != NA_DFLOAT)
    rb_work = na_change_type(rb_work, NA_DFLOAT);
  work = NA_PTR_TYPE(rb_work, doublereal*);

  __out__ = dla_syrpvgrw_(&uplo, &n, &info, a, &lda, af, &ldaf, ipiv, work);

  rb___out__ = rb_float_new((double)__out__);
  return rb___out__;
}

void
init_lapack_dla_syrpvgrw(VALUE mLapack){
  rb_define_module_function(mLapack, "dla_syrpvgrw", rb_dla_syrpvgrw, -1);
}

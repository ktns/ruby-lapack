#include "rb_lapack.h"

extern doublereal zla_porpvgrw_(char *uplo, integer *ncols, doublecomplex *a, integer *lda, doublecomplex *af, integer *ldaf, doublecomplex *work);

static VALUE
rb_zla_porpvgrw(int argc, VALUE *argv, VALUE self){
  VALUE rb_uplo;
  char uplo; 
  VALUE rb_ncols;
  integer ncols; 
  VALUE rb_a;
  doublecomplex *a; 
  VALUE rb_af;
  doublecomplex *af; 
  VALUE rb_work;
  doublecomplex *work; 
  VALUE rb___out__;
  doublereal __out__; 

  integer lda;
  integer n;
  integer ldaf;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  __out__ = NumRu::Lapack.zla_porpvgrw( uplo, ncols, a, af, work)\n    or\n  NumRu::Lapack.zla_porpvgrw  # print help\n\n\nFORTRAN MANUAL\n      DOUBLE PRECISION FUNCTION ZLA_PORPVGRW( UPLO, NCOLS, A, LDA, AF,  LDAF, WORK )\n\n*  Purpose\n*  =======\n* \n*  ZLA_PORPVGRW computes the reciprocal pivot growth factor\n*  norm(A)/norm(U). The \"max absolute element\" norm is used. If this is\n*  much less than 1, the stability of the LU factorization of the\n*  (equilibrated) matrix A could be poor. This also means that the\n*  solution X, estimated condition numbers, and error bounds could be\n*  unreliable.\n*\n\n*  Arguments\n*  =========\n*\n*     UPLO    (input) CHARACTER*1\n*       = 'U':  Upper triangle of A is stored;\n*       = 'L':  Lower triangle of A is stored.\n*\n*     NCOLS   (input) INTEGER\n*     The number of columns of the matrix A. NCOLS >= 0.\n*\n*     A       (input) COMPLEX*16 array, dimension (LDA,N)\n*     On entry, the N-by-N matrix A.\n*\n*     LDA     (input) INTEGER\n*     The leading dimension of the array A.  LDA >= max(1,N).\n*\n*     AF      (input) COMPLEX*16 array, dimension (LDAF,N)\n*     The triangular factor U or L from the Cholesky factorization\n*     A = U**T*U or A = L*L**T, as computed by ZPOTRF.\n*\n*     LDAF    (input) INTEGER\n*     The leading dimension of the array AF.  LDAF >= max(1,N).\n*\n*     WORK    (input) COMPLEX*16 array, dimension (2*N)\n*\n\n*  =====================================================================\n*\n*     .. Local Scalars ..\n      INTEGER            I, J\n      DOUBLE PRECISION   AMAX, UMAX, RPVGRW\n      LOGICAL            UPPER\n      COMPLEX*16         ZDUM\n*     ..\n*     .. External Functions ..\n      EXTERNAL           LSAME, ZLASET\n      LOGICAL            LSAME\n*     ..\n*     .. Intrinsic Functions ..\n      INTRINSIC          ABS, MAX, MIN, REAL, DIMAG\n*     ..\n*     .. Statement Functions ..\n      DOUBLE PRECISION   CABS1\n*     ..\n*     .. Statement Function Definitions ..\n      CABS1( ZDUM ) = ABS( DBLE( ZDUM ) ) + ABS( DIMAG( ZDUM ) )\n*     ..\n\n");
    return Qnil;
  }
  if (argc != 5)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 5)", argc);
  rb_uplo = argv[0];
  rb_ncols = argv[1];
  rb_a = argv[2];
  rb_af = argv[3];
  rb_work = argv[4];

  if (!NA_IsNArray(rb_a))
    rb_raise(rb_eArgError, "a (3th argument) must be NArray");
  if (NA_RANK(rb_a) != 2)
    rb_raise(rb_eArgError, "rank of a (3th argument) must be %d", 2);
  n = NA_SHAPE1(rb_a);
  lda = NA_SHAPE0(rb_a);
  if (NA_TYPE(rb_a) != NA_DCOMPLEX)
    rb_a = na_change_type(rb_a, NA_DCOMPLEX);
  a = NA_PTR_TYPE(rb_a, doublecomplex*);
  ncols = NUM2INT(rb_ncols);
  uplo = StringValueCStr(rb_uplo)[0];
  if (!NA_IsNArray(rb_af))
    rb_raise(rb_eArgError, "af (4th argument) must be NArray");
  if (NA_RANK(rb_af) != 2)
    rb_raise(rb_eArgError, "rank of af (4th argument) must be %d", 2);
  if (NA_SHAPE1(rb_af) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of af must be the same as shape 1 of a");
  ldaf = NA_SHAPE0(rb_af);
  if (NA_TYPE(rb_af) != NA_DCOMPLEX)
    rb_af = na_change_type(rb_af, NA_DCOMPLEX);
  af = NA_PTR_TYPE(rb_af, doublecomplex*);
  if (!NA_IsNArray(rb_work))
    rb_raise(rb_eArgError, "work (5th argument) must be NArray");
  if (NA_RANK(rb_work) != 1)
    rb_raise(rb_eArgError, "rank of work (5th argument) must be %d", 1);
  if (NA_SHAPE0(rb_work) != (2*n))
    rb_raise(rb_eRuntimeError, "shape 0 of work must be %d", 2*n);
  if (NA_TYPE(rb_work) != NA_DCOMPLEX)
    rb_work = na_change_type(rb_work, NA_DCOMPLEX);
  work = NA_PTR_TYPE(rb_work, doublecomplex*);

  __out__ = zla_porpvgrw_(&uplo, &ncols, a, &lda, af, &ldaf, work);

  rb___out__ = rb_float_new((double)__out__);
  return rb___out__;
}

void
init_lapack_zla_porpvgrw(VALUE mLapack){
  rb_define_module_function(mLapack, "zla_porpvgrw", rb_zla_porpvgrw, -1);
}

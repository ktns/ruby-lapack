#include "rb_lapack.h"

extern real sla_porpvgrw_(char *uplo, integer *ncols, real *a, integer *lda, real *af, integer *ldaf, real *work);

static VALUE
rb_sla_porpvgrw(int argc, VALUE *argv, VALUE self){
  VALUE rb_uplo;
  char uplo; 
  VALUE rb_ncols;
  integer ncols; 
  VALUE rb_a;
  real *a; 
  VALUE rb_af;
  real *af; 
  VALUE rb_work;
  real *work; 
  VALUE rb___out__;
  real __out__; 

  integer lda;
  integer n;
  integer ldaf;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  __out__ = NumRu::Lapack.sla_porpvgrw( uplo, ncols, a, af, work)\n    or\n  NumRu::Lapack.sla_porpvgrw  # print help\n\n\nFORTRAN MANUAL\n      REAL FUNCTION SLA_PORPVGRW( UPLO, NCOLS, A, LDA, AF, LDAF, WORK )\n\n*  Purpose\n*  =======\n* \n*  SLA_PORPVGRW computes the reciprocal pivot growth factor\n*  norm(A)/norm(U). The \"max absolute element\" norm is used. If this is\n*  much less than 1, the stability of the LU factorization of the\n*  (equilibrated) matrix A could be poor. This also means that the\n*  solution X, estimated condition numbers, and error bounds could be\n*  unreliable.\n*\n\n*  Arguments\n*  =========\n*\n*     UPLO    (input) CHARACTER*1\n*       = 'U':  Upper triangle of A is stored;\n*       = 'L':  Lower triangle of A is stored.\n*\n*     NCOLS   (input) INTEGER\n*     The number of columns of the matrix A. NCOLS >= 0.\n*\n*     A       (input) REAL array, dimension (LDA,N)\n*     On entry, the N-by-N matrix A.\n*\n*     LDA     (input) INTEGER\n*     The leading dimension of the array A.  LDA >= max(1,N).\n*\n*     AF      (input) REAL array, dimension (LDAF,N)\n*     The triangular factor U or L from the Cholesky factorization\n*     A = U**T*U or A = L*L**T, as computed by SPOTRF.\n*\n*     LDAF    (input) INTEGER\n*     The leading dimension of the array AF.  LDAF >= max(1,N).\n*\n*     WORK    (input) REAL array, dimension (2*N)\n*\n\n*  =====================================================================\n*\n*     .. Local Scalars ..\n      INTEGER            I, J\n      REAL               AMAX, UMAX, RPVGRW\n      LOGICAL            UPPER\n*     ..\n*     .. Intrinsic Functions ..\n      INTRINSIC          ABS, MAX, MIN\n*     ..\n*     .. External Functions ..\n      EXTERNAL           LSAME, SLASET\n      LOGICAL            LSAME\n*     ..\n\n");
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
  if (NA_TYPE(rb_a) != NA_SFLOAT)
    rb_a = na_change_type(rb_a, NA_SFLOAT);
  a = NA_PTR_TYPE(rb_a, real*);
  ncols = NUM2INT(rb_ncols);
  uplo = StringValueCStr(rb_uplo)[0];
  if (!NA_IsNArray(rb_af))
    rb_raise(rb_eArgError, "af (4th argument) must be NArray");
  if (NA_RANK(rb_af) != 2)
    rb_raise(rb_eArgError, "rank of af (4th argument) must be %d", 2);
  if (NA_SHAPE1(rb_af) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of af must be the same as shape 1 of a");
  ldaf = NA_SHAPE0(rb_af);
  if (NA_TYPE(rb_af) != NA_SFLOAT)
    rb_af = na_change_type(rb_af, NA_SFLOAT);
  af = NA_PTR_TYPE(rb_af, real*);
  if (!NA_IsNArray(rb_work))
    rb_raise(rb_eArgError, "work (5th argument) must be NArray");
  if (NA_RANK(rb_work) != 1)
    rb_raise(rb_eArgError, "rank of work (5th argument) must be %d", 1);
  if (NA_SHAPE0(rb_work) != (2*n))
    rb_raise(rb_eRuntimeError, "shape 0 of work must be %d", 2*n);
  if (NA_TYPE(rb_work) != NA_SFLOAT)
    rb_work = na_change_type(rb_work, NA_SFLOAT);
  work = NA_PTR_TYPE(rb_work, real*);

  __out__ = sla_porpvgrw_(&uplo, &ncols, a, &lda, af, &ldaf, work);

  rb___out__ = rb_float_new((double)__out__);
  return rb___out__;
}

void
init_lapack_sla_porpvgrw(VALUE mLapack){
  rb_define_module_function(mLapack, "sla_porpvgrw", rb_sla_porpvgrw, -1);
}

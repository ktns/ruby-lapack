#include "rb_lapack.h"

extern real sla_rpvgrw_(integer *n, integer *ncols, real *a, integer *lda, real *af, integer *ldaf);

static VALUE
rb_sla_rpvgrw(int argc, VALUE *argv, VALUE self){
  VALUE rb_ncols;
  integer ncols; 
  VALUE rb_a;
  real *a; 
  VALUE rb_af;
  real *af; 
  VALUE rb___out__;
  real __out__; 

  integer lda;
  integer n;
  integer ldaf;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  __out__ = NumRu::Lapack.sla_rpvgrw( ncols, a, af)\n    or\n  NumRu::Lapack.sla_rpvgrw  # print help\n\n\nFORTRAN MANUAL\n      REAL FUNCTION SLA_RPVGRW( N, NCOLS, A, LDA, AF, LDAF )\n\n*  Purpose\n*  =======\n*\n*  SLA_RPVGRW computes the reciprocal pivot growth factor\n*  norm(A)/norm(U). The \"max absolute element\" norm is used. If this is\n*  much less than 1, the stability of the LU factorization of the\n*  (equilibrated) matrix A could be poor. This also means that the\n*  solution X, estimated condition numbers, and error bounds could be\n*  unreliable.\n*\n\n*  Arguments\n*  =========\n*\n*     N       (input) INTEGER\n*     The number of linear equations, i.e., the order of the\n*     matrix A.  N >= 0.\n*\n*     NCOLS   (input) INTEGER\n*     The number of columns of the matrix A. NCOLS >= 0.\n*\n*     A       (input) REAL array, dimension (LDA,N)\n*     On entry, the N-by-N matrix A.\n*\n*     LDA     (input) INTEGER\n*     The leading dimension of the array A.  LDA >= max(1,N).\n*\n*     AF      (input) REAL array, dimension (LDAF,N)\n*     The factors L and U from the factorization\n*     A = P*L*U as computed by SGETRF.\n*\n*     LDAF    (input) INTEGER\n*     The leading dimension of the array AF.  LDAF >= max(1,N).\n*\n\n*  =====================================================================\n*\n*     .. Local Scalars ..\n      INTEGER            I, J\n      REAL               AMAX, UMAX, RPVGRW\n*     ..\n*     .. Intrinsic Functions ..\n      INTRINSIC          ABS, MAX, MIN\n*     ..\n\n");
    return Qnil;
  }
  if (argc != 3)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 3)", argc);
  rb_ncols = argv[0];
  rb_a = argv[1];
  rb_af = argv[2];

  if (!NA_IsNArray(rb_a))
    rb_raise(rb_eArgError, "a (2th argument) must be NArray");
  if (NA_RANK(rb_a) != 2)
    rb_raise(rb_eArgError, "rank of a (2th argument) must be %d", 2);
  n = NA_SHAPE1(rb_a);
  lda = NA_SHAPE0(rb_a);
  if (NA_TYPE(rb_a) != NA_SFLOAT)
    rb_a = na_change_type(rb_a, NA_SFLOAT);
  a = NA_PTR_TYPE(rb_a, real*);
  ncols = NUM2INT(rb_ncols);
  if (!NA_IsNArray(rb_af))
    rb_raise(rb_eArgError, "af (3th argument) must be NArray");
  if (NA_RANK(rb_af) != 2)
    rb_raise(rb_eArgError, "rank of af (3th argument) must be %d", 2);
  if (NA_SHAPE1(rb_af) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of af must be the same as shape 1 of a");
  ldaf = NA_SHAPE0(rb_af);
  if (NA_TYPE(rb_af) != NA_SFLOAT)
    rb_af = na_change_type(rb_af, NA_SFLOAT);
  af = NA_PTR_TYPE(rb_af, real*);

  __out__ = sla_rpvgrw_(&n, &ncols, a, &lda, af, &ldaf);

  rb___out__ = rb_float_new((double)__out__);
  return rb___out__;
}

void
init_lapack_sla_rpvgrw(VALUE mLapack){
  rb_define_module_function(mLapack, "sla_rpvgrw", rb_sla_rpvgrw, -1);
}

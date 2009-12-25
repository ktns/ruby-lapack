#include "rb_lapack.h"

extern VOID dla_rpvgrw_(doublereal *__out__, integer *n, integer *ncols, doublereal *a, integer *lda, doublereal *af, integer *ldaf);
static VALUE
rb_dla_rpvgrw(int argc, VALUE *argv, VALUE self){
  VALUE rb_ncols;
  integer ncols; 
  VALUE rb_a;
  doublereal *a; 
  VALUE rb_af;
  doublereal *af; 
  VALUE rb___out__;
  doublereal __out__; 

  integer lda;
  integer n;
  integer ldaf;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  __out__ = NumRu::Lapack.dla_rpvgrw( ncols, a, af)\n    or\n  NumRu::Lapack.dla_rpvgrw  # print help\n\n\nFORTRAN MANUAL\n      DOUBLE PRECISION FUNCTION DLA_RPVGRW( N, NCOLS, A, LDA, AF, LDAF )\n\n*  Purpose\n*  =======\n* \n*  DLA_RPVGRW computes the reciprocal pivot growth factor\n*  norm(A)/norm(U). The \"max absolute element\" norm is used. If this is\n*  much less than 1, the stability of the LU factorization of the\n*  (equilibrated) matrix A could be poor. This also means that the\n*  solution X, estimated condition numbers, and error bounds could be\n*  unreliable.\n*\n\n*  Arguments\n*  =========\n*\n*     N       (input) INTEGER\n*     The number of linear equations, i.e., the order of the\n*     matrix A.  N >= 0.\n*\n*     NCOLS   (input) INTEGER\n*     The number of columns of the matrix A. NCOLS >= 0.\n*\n*     A       (input) DOUBLE PRECISION array, dimension (LDA,N)\n*     On entry, the N-by-N matrix A.\n*\n*     LDA     (input) INTEGER\n*     The leading dimension of the array A.  LDA >= max(1,N).\n*\n*     AF      (input) DOUBLE PRECISION array, dimension (LDAF,N)\n*     The factors L and U from the factorization\n*     A = P*L*U as computed by DGETRF.\n*\n*     LDAF    (input) INTEGER\n*     The leading dimension of the array AF.  LDAF >= max(1,N).\n*\n\n*  =====================================================================\n*\n*     .. Local Scalars ..\n      INTEGER            I, J\n      DOUBLE PRECISION   AMAX, UMAX, RPVGRW\n*     ..\n*     .. Intrinsic Functions ..\n      INTRINSIC          ABS, MAX, MIN\n*     ..\n\n");
    return Qnil;
  }
  if (argc != 3)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 3)", argc);
  rb_ncols = argv[0];
  rb_a = argv[1];
  rb_af = argv[2];

  ncols = NUM2INT(rb_ncols);
  if (!NA_IsNArray(rb_a))
    rb_raise(rb_eArgError, "a (2th argument) must be NArray");
  if (NA_RANK(rb_a) != 2)
    rb_raise(rb_eArgError, "rank of a (2th argument) must be %d", 2);
  lda = NA_SHAPE0(rb_a);
  n = NA_SHAPE1(rb_a);
  if (NA_TYPE(rb_a) != NA_DFLOAT)
    rb_a = na_change_type(rb_a, NA_DFLOAT);
  a = NA_PTR_TYPE(rb_a, doublereal*);
  if (!NA_IsNArray(rb_af))
    rb_raise(rb_eArgError, "af (3th argument) must be NArray");
  if (NA_RANK(rb_af) != 2)
    rb_raise(rb_eArgError, "rank of af (3th argument) must be %d", 2);
  ldaf = NA_SHAPE0(rb_af);
  if (NA_SHAPE1(rb_af) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of af must be the same as shape 1 of a");
  if (NA_TYPE(rb_af) != NA_DFLOAT)
    rb_af = na_change_type(rb_af, NA_DFLOAT);
  af = NA_PTR_TYPE(rb_af, doublereal*);

  dla_rpvgrw_(&__out__, &n, &ncols, a, &lda, af, &ldaf);

  rb___out__ = rb_float_new((double)__out__);
  return rb___out__;
}

void
init_lapack_dla_rpvgrw(VALUE mLapack){
  rb_define_module_function(mLapack, "dla_rpvgrw", rb_dla_rpvgrw, -1);
}

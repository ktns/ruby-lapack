#include "rb_lapack.h"

extern VOID cla_geamv_(integer *trans, integer *m, integer *n, real *alpha, complex *a, integer *lda, complex *x, integer *incx, real *beta, real *y, integer *incy);

static VALUE
rb_cla_geamv(int argc, VALUE *argv, VALUE self){
  VALUE rb_trans;
  integer trans; 
  VALUE rb_m;
  integer m; 
  VALUE rb_alpha;
  real alpha; 
  VALUE rb_a;
  complex *a; 
  VALUE rb_x;
  complex *x; 
  VALUE rb_incx;
  integer incx; 
  VALUE rb_beta;
  real beta; 
  VALUE rb_y;
  real *y; 
  VALUE rb_incy;
  integer incy; 
  VALUE rb_y_out__;
  real *y_out__;

  integer lda;
  integer n;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  y = NumRu::Lapack.cla_geamv( trans, m, alpha, a, x, incx, beta, y, incy)\n    or\n  NumRu::Lapack.cla_geamv  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE CLA_GEAMV ( TRANS, M, N, ALPHA, A, LDA, X, INCX, BETA, Y, INCY )\n\n*  Purpose\n*  =======\n*\n*  CLA_GEAMV  performs one of the matrix-vector operations\n*\n*          y := alpha*abs(A)*abs(x) + beta*abs(y),\n*     or   y := alpha*abs(A)'*abs(x) + beta*abs(y),\n*\n*  where alpha and beta are scalars, x and y are vectors and A is an\n*  m by n matrix.\n*\n*  This function is primarily used in calculating error bounds.\n*  To protect against underflow during evaluation, components in\n*  the resulting vector are perturbed away from zero by (N+1)\n*  times the underflow threshold.  To prevent unnecessarily large\n*  errors for block-structure embedded in general matrices,\n*  \"symbolically\" zero components are not perturbed.  A zero\n*  entry is considered \"symbolic\" if all multiplications involved\n*  in computing that entry have at least one zero multiplicand.\n*\n\n*  Arguments\n*  ==========\n*\n*  TRANS   (input) INTEGER\n*           On entry, TRANS specifies the operation to be performed as\n*           follows:\n*\n*             BLAS_NO_TRANS      y := alpha*abs(A)*abs(x) + beta*abs(y)\n*             BLAS_TRANS         y := alpha*abs(A')*abs(x) + beta*abs(y)\n*             BLAS_CONJ_TRANS    y := alpha*abs(A')*abs(x) + beta*abs(y)\n*\n*           Unchanged on exit.\n*\n*  M       (input) INTEGER\n*           On entry, M specifies the number of rows of the matrix A.\n*           M must be at least zero.\n*           Unchanged on exit.\n*\n*  N       (input) INTEGER\n*           On entry, N specifies the number of columns of the matrix A.\n*           N must be at least zero.\n*           Unchanged on exit.\n*\n*  ALPHA   (input) REAL\n*           On entry, ALPHA specifies the scalar alpha.\n*           Unchanged on exit.\n*\n*  A       (input) COMPLEX array, dimension (LDA,n)\n*           Before entry, the leading m by n part of the array A must\n*           contain the matrix of coefficients.\n*           Unchanged on exit.\n*\n*  LDA     (input) INTEGER\n*           On entry, LDA specifies the first dimension of A as declared\n*           in the calling (sub) program. LDA must be at least\n*           max( 1, m ).\n*           Unchanged on exit.\n*\n*  X       (input) COMPLEX array, dimension\n*           ( 1 + ( n - 1 )*abs( INCX ) ) when TRANS = 'N' or 'n'\n*           and at least\n*           ( 1 + ( m - 1 )*abs( INCX ) ) otherwise.\n*           Before entry, the incremented array X must contain the\n*           vector x.\n*           Unchanged on exit.\n*\n*  INCX    (input) INTEGER\n*           On entry, INCX specifies the increment for the elements of\n*           X. INCX must not be zero.\n*           Unchanged on exit.\n*\n*  BETA    (input) REAL\n*           On entry, BETA specifies the scalar beta. When BETA is\n*           supplied as zero then Y need not be set on input.\n*           Unchanged on exit.\n*\n*  Y       (input/output) REAL array, dimension\n*           ( 1 + ( m - 1 )*abs( INCY ) ) when TRANS = 'N' or 'n'\n*           and at least\n*           ( 1 + ( n - 1 )*abs( INCY ) ) otherwise.\n*           Before entry with BETA non-zero, the incremented array Y\n*           must contain the vector y. On exit, Y is overwritten by the\n*           updated vector y.\n*\n*  INCY    (input) INTEGER\n*           On entry, INCY specifies the increment for the elements of\n*           Y. INCY must not be zero.\n*           Unchanged on exit.\n*\n*\n*  Level 2 Blas routine.\n*\n\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 9)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 9)", argc);
  rb_trans = argv[0];
  rb_m = argv[1];
  rb_alpha = argv[2];
  rb_a = argv[3];
  rb_x = argv[4];
  rb_incx = argv[5];
  rb_beta = argv[6];
  rb_y = argv[7];
  rb_incy = argv[8];

  if (!NA_IsNArray(rb_a))
    rb_raise(rb_eArgError, "a (4th argument) must be NArray");
  if (NA_RANK(rb_a) != 2)
    rb_raise(rb_eArgError, "rank of a (4th argument) must be %d", 2);
  n = NA_SHAPE1(rb_a);
  lda = NA_SHAPE0(rb_a);
  if (NA_TYPE(rb_a) != NA_SCOMPLEX)
    rb_a = na_change_type(rb_a, NA_SCOMPLEX);
  a = NA_PTR_TYPE(rb_a, complex*);
  trans = NUM2INT(rb_trans);
  m = NUM2INT(rb_m);
  alpha = (real)NUM2DBL(rb_alpha);
  incy = NUM2INT(rb_incy);
  incx = NUM2INT(rb_incx);
  beta = (real)NUM2DBL(rb_beta);
  if (!NA_IsNArray(rb_y))
    rb_raise(rb_eArgError, "y (8th argument) must be NArray");
  if (NA_RANK(rb_y) != 1)
    rb_raise(rb_eArgError, "rank of y (8th argument) must be %d", 1);
  if (NA_SHAPE0(rb_y) != (((lsame_(&trans,"N")) || (lsame_(&trans,"n"))) ? 1 + ( m - 1 )*abs( incy ) : 1 + ( n - 1 )*abs( incy )))
    rb_raise(rb_eRuntimeError, "shape 0 of y must be %d", ((lsame_(&trans,"N")) || (lsame_(&trans,"n"))) ? 1 + ( m - 1 )*abs( incy ) : 1 + ( n - 1 )*abs( incy ));
  if (NA_TYPE(rb_y) != NA_SFLOAT)
    rb_y = na_change_type(rb_y, NA_SFLOAT);
  y = NA_PTR_TYPE(rb_y, real*);
  if (!NA_IsNArray(rb_x))
    rb_raise(rb_eArgError, "x (5th argument) must be NArray");
  if (NA_RANK(rb_x) != 1)
    rb_raise(rb_eArgError, "rank of x (5th argument) must be %d", 1);
  if (NA_SHAPE0(rb_x) != (((lsame_(&trans,"N")) || (lsame_(&trans,"n"))) ? 1 + ( n - 1 )*abs( incx ) : 1 + ( m - 1 )*abs( incx )))
    rb_raise(rb_eRuntimeError, "shape 0 of x must be %d", ((lsame_(&trans,"N")) || (lsame_(&trans,"n"))) ? 1 + ( n - 1 )*abs( incx ) : 1 + ( m - 1 )*abs( incx ));
  if (NA_TYPE(rb_x) != NA_SCOMPLEX)
    rb_x = na_change_type(rb_x, NA_SCOMPLEX);
  x = NA_PTR_TYPE(rb_x, complex*);
  {
    int shape[1];
    shape[0] = ((lsame_(&trans,"N")) || (lsame_(&trans,"n"))) ? 1 + ( m - 1 )*abs( incy ) : 1 + ( n - 1 )*abs( incy );
    rb_y_out__ = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  y_out__ = NA_PTR_TYPE(rb_y_out__, real*);
  MEMCPY(y_out__, y, real, NA_TOTAL(rb_y));
  rb_y = rb_y_out__;
  y = y_out__;

  cla_geamv_(&trans, &m, &n, &alpha, a, &lda, x, &incx, &beta, y, &incy);

  return rb_y;
}

void
init_lapack_cla_geamv(VALUE mLapack){
  rb_define_module_function(mLapack, "cla_geamv", rb_cla_geamv, -1);
}

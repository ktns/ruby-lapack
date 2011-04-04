#include "rb_lapack.h"

extern VOID dla_geamv_(integer *trans, integer *m, integer *n, doublereal *alpha, doublereal *a, integer *lda, doublereal *x, integer *incx, doublereal *beta, doublereal *y, integer *incy);

static VALUE
rb_dla_geamv(int argc, VALUE *argv, VALUE self){
  VALUE rb_trans;
  integer trans; 
  VALUE rb_m;
  integer m; 
  VALUE rb_alpha;
  doublereal alpha; 
  VALUE rb_a;
  doublereal *a; 
  VALUE rb_x;
  doublereal *x; 
  VALUE rb_incx;
  integer incx; 
  VALUE rb_beta;
  doublereal beta; 
  VALUE rb_y;
  doublereal *y; 
  VALUE rb_incy;
  integer incy; 
  VALUE rb_y_out__;
  doublereal *y_out__;

  integer lda;
  integer n;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  y = NumRu::Lapack.dla_geamv( trans, m, alpha, a, x, incx, beta, y, incy)\n    or\n  NumRu::Lapack.dla_geamv  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE DLA_GEAMV ( TRANS, M, N, ALPHA, A, LDA, X, INCX, BETA, Y, INCY )\n\n*  Purpose\n*  =======\n*\n*  DLA_GEAMV  performs one of the matrix-vector operations\n*\n*          y := alpha*abs(A)*abs(x) + beta*abs(y),\n*     or   y := alpha*abs(A)'*abs(x) + beta*abs(y),\n*\n*  where alpha and beta are scalars, x and y are vectors and A is an\n*  m by n matrix.\n*\n*  This function is primarily used in calculating error bounds.\n*  To protect against underflow during evaluation, components in\n*  the resulting vector are perturbed away from zero by (N+1)\n*  times the underflow threshold.  To prevent unnecessarily large\n*  errors for block-structure embedded in general matrices,\n*  \"symbolically\" zero components are not perturbed.  A zero\n*  entry is considered \"symbolic\" if all multiplications involved\n*  in computing that entry have at least one zero multiplicand.\n*\n\n*  Arguments\n*  ==========\n*\n*  TRANS   (input) INTEGER\n*           On entry, TRANS specifies the operation to be performed as\n*           follows:\n*\n*             BLAS_NO_TRANS      y := alpha*abs(A)*abs(x) + beta*abs(y)\n*             BLAS_TRANS         y := alpha*abs(A')*abs(x) + beta*abs(y)\n*             BLAS_CONJ_TRANS    y := alpha*abs(A')*abs(x) + beta*abs(y)\n*\n*           Unchanged on exit.\n*\n*  M       (input) INTEGER\n*           On entry, M specifies the number of rows of the matrix A.\n*           M must be at least zero.\n*           Unchanged on exit.\n*\n*  N       (input) INTEGER\n*           On entry, N specifies the number of columns of the matrix A.\n*           N must be at least zero.\n*           Unchanged on exit.\n*\n*  ALPHA  - DOUBLE PRECISION\n*           On entry, ALPHA specifies the scalar alpha.\n*           Unchanged on exit.\n*\n*  A      - DOUBLE PRECISION   array of DIMENSION ( LDA, n )\n*           Before entry, the leading m by n part of the array A must\n*           contain the matrix of coefficients.\n*           Unchanged on exit.\n*\n*  LDA     (input) INTEGER\n*           On entry, LDA specifies the first dimension of A as declared\n*           in the calling (sub) program. LDA must be at least\n*           max( 1, m ).\n*           Unchanged on exit.\n*\n*  X       (input) DOUBLE PRECISION array, dimension\n*           ( 1 + ( n - 1 )*abs( INCX ) ) when TRANS = 'N' or 'n'\n*           and at least\n*           ( 1 + ( m - 1 )*abs( INCX ) ) otherwise.\n*           Before entry, the incremented array X must contain the\n*           vector x.\n*           Unchanged on exit.\n*\n*  INCX    (input) INTEGER\n*           On entry, INCX specifies the increment for the elements of\n*           X. INCX must not be zero.\n*           Unchanged on exit.\n*\n*  BETA   - DOUBLE PRECISION\n*           On entry, BETA specifies the scalar beta. When BETA is\n*           supplied as zero then Y need not be set on input.\n*           Unchanged on exit.\n*\n*  Y      - DOUBLE PRECISION\n*           Array of DIMENSION at least\n*           ( 1 + ( m - 1 )*abs( INCY ) ) when TRANS = 'N' or 'n'\n*           and at least\n*           ( 1 + ( n - 1 )*abs( INCY ) ) otherwise.\n*           Before entry with BETA non-zero, the incremented array Y\n*           must contain the vector y. On exit, Y is overwritten by the\n*           updated vector y.\n*\n*  INCY    (input) INTEGER\n*           On entry, INCY specifies the increment for the elements of\n*           Y. INCY must not be zero.\n*           Unchanged on exit.\n*\n*  Level 2 Blas routine.\n*\n\n*  =====================================================================\n*\n\n");
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
  if (lda != (MAX(1, m)))
    rb_raise(rb_eRuntimeError, "shape 0 of a must be %d", MAX(1, m));
  if (NA_TYPE(rb_a) != NA_DFLOAT)
    rb_a = na_change_type(rb_a, NA_DFLOAT);
  a = NA_PTR_TYPE(rb_a, doublereal*);
  trans = NUM2INT(rb_trans);
  m = NUM2INT(rb_m);
  alpha = NUM2DBL(rb_alpha);
  beta = NUM2DBL(rb_beta);
  incy = NUM2INT(rb_incy);
  incx = NUM2INT(rb_incx);
  lda = MAX(1, m);
  if (!NA_IsNArray(rb_y))
    rb_raise(rb_eArgError, "y (8th argument) must be NArray");
  if (NA_RANK(rb_y) != 1)
    rb_raise(rb_eArgError, "rank of y (8th argument) must be %d", 1);
  if (NA_SHAPE0(rb_y) != (lsame_(&trans,"N") ? 1+(m-1)*abs(incy) : 1+(n-1)*abs(incy)))
    rb_raise(rb_eRuntimeError, "shape 0 of y must be %d", lsame_(&trans,"N") ? 1+(m-1)*abs(incy) : 1+(n-1)*abs(incy));
  if (NA_TYPE(rb_y) != NA_DFLOAT)
    rb_y = na_change_type(rb_y, NA_DFLOAT);
  y = NA_PTR_TYPE(rb_y, doublereal*);
  if (!NA_IsNArray(rb_x))
    rb_raise(rb_eArgError, "x (5th argument) must be NArray");
  if (NA_RANK(rb_x) != 1)
    rb_raise(rb_eArgError, "rank of x (5th argument) must be %d", 1);
  if (NA_SHAPE0(rb_x) != (((lsame_(&trans,"N")) || (lsame_(&trans,"n"))) ? 1 + ( n - 1 )*abs( incx ) : 1 + ( m - 1 )*abs( incx )))
    rb_raise(rb_eRuntimeError, "shape 0 of x must be %d", ((lsame_(&trans,"N")) || (lsame_(&trans,"n"))) ? 1 + ( n - 1 )*abs( incx ) : 1 + ( m - 1 )*abs( incx ));
  if (NA_TYPE(rb_x) != NA_DFLOAT)
    rb_x = na_change_type(rb_x, NA_DFLOAT);
  x = NA_PTR_TYPE(rb_x, doublereal*);
  {
    int shape[1];
    shape[0] = lsame_(&trans,"N") ? 1+(m-1)*abs(incy) : 1+(n-1)*abs(incy);
    rb_y_out__ = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  y_out__ = NA_PTR_TYPE(rb_y_out__, doublereal*);
  MEMCPY(y_out__, y, doublereal, NA_TOTAL(rb_y));
  rb_y = rb_y_out__;
  y = y_out__;

  dla_geamv_(&trans, &m, &n, &alpha, a, &lda, x, &incx, &beta, y, &incy);

  return rb_y;
}

void
init_lapack_dla_geamv(VALUE mLapack){
  rb_define_module_function(mLapack, "dla_geamv", rb_dla_geamv, -1);
}

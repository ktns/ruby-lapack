#include "rb_lapack.h"

extern VOID cla_gbamv_(integer *trans, integer *m, integer *n, integer *kl, integer *ku, real *alpha, real *ab, integer *ldab, real *x, integer *incx, real *beta, real *y, integer *incy);

static VALUE
rb_cla_gbamv(int argc, VALUE *argv, VALUE self){
  VALUE rb_trans;
  integer trans; 
  VALUE rb_m;
  integer m; 
  VALUE rb_kl;
  integer kl; 
  VALUE rb_ku;
  integer ku; 
  VALUE rb_alpha;
  real alpha; 
  VALUE rb_ab;
  real *ab; 
  VALUE rb_x;
  real *x; 
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

  integer ldab;
  integer n;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  y = NumRu::Lapack.cla_gbamv( trans, m, kl, ku, alpha, ab, x, incx, beta, y, incy)\n    or\n  NumRu::Lapack.cla_gbamv  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE CLA_GBAMV( TRANS, M, N, KL, KU, ALPHA, AB, LDAB, X, INCX, BETA, Y, INCY )\n\n*  Purpose\n*  =======\n*\n*  SLA_GBAMV  performs one of the matrix-vector operations\n*\n*          y := alpha*abs(A)*abs(x) + beta*abs(y),\n*     or   y := alpha*abs(A)'*abs(x) + beta*abs(y),\n*\n*  where alpha and beta are scalars, x and y are vectors and A is an\n*  m by n matrix.\n*\n*  This function is primarily used in calculating error bounds.\n*  To protect against underflow during evaluation, components in\n*  the resulting vector are perturbed away from zero by (N+1)\n*  times the underflow threshold.  To prevent unnecessarily large\n*  errors for block-structure embedded in general matrices,\n*  \"symbolically\" zero components are not perturbed.  A zero\n*  entry is considered \"symbolic\" if all multiplications involved\n*  in computing that entry have at least one zero multiplicand.\n*\n\n*  Arguments\n*  ==========\n*\n*  TRANS   (input) INTEGER\n*           On entry, TRANS specifies the operation to be performed as\n*           follows:\n*\n*             BLAS_NO_TRANS      y := alpha*abs(A)*abs(x) + beta*abs(y)\n*             BLAS_TRANS         y := alpha*abs(A')*abs(x) + beta*abs(y)\n*             BLAS_CONJ_TRANS    y := alpha*abs(A')*abs(x) + beta*abs(y)\n*\n*           Unchanged on exit.\n*\n*  M      (input) INTEGER\n*           On entry, M specifies the number of rows of the matrix A.\n*           M must be at least zero.\n*           Unchanged on exit.\n*\n*  N      (input) INTEGER\n*           On entry, N specifies the number of columns of the matrix A.\n*           N must be at least zero.\n*           Unchanged on exit.\n*\n*  KL     (input) INTEGER\n*           The number of subdiagonals within the band of A.  KL >= 0.\n*\n*  KU     (input) INTEGER\n*           The number of superdiagonals within the band of A.  KU >= 0.\n*\n*  ALPHA  (input) REAL\n*           On entry, ALPHA specifies the scalar alpha.\n*           Unchanged on exit.\n*\n*  A      (input) REAL array, dimension (LDA,n)\n*           Before entry, the leading m by n part of the array A must\n*           contain the matrix of coefficients.\n*           Unchanged on exit.\n*\n*  LDA    (input) INTEGER\n*           On entry, LDA specifies the first dimension of A as declared\n*           in the calling (sub) program. LDA must be at least\n*           max( 1, m ).\n*           Unchanged on exit.\n*\n*  X      (input) REAL array, dimension at least\n*           ( 1 + ( n - 1 )*abs( INCX ) ) when TRANS = 'N' or 'n'\n*           and at least\n*           ( 1 + ( m - 1 )*abs( INCX ) ) otherwise.\n*           Before entry, the incremented array X must contain the\n*           vector x.\n*           Unchanged on exit.\n*\n*  INCX   (input) INTEGER\n*           On entry, INCX specifies the increment for the elements of\n*           X. INCX must not be zero.\n*           Unchanged on exit.\n*\n*  BETA   (input) REAL\n*           On entry, BETA specifies the scalar beta. When BETA is\n*           supplied as zero then Y need not be set on input.\n*           Unchanged on exit.\n*\n*  Y      (input/output) REAL array, dimension at least\n*           ( 1 + ( m - 1 )*abs( INCY ) ) when TRANS = 'N' or 'n'\n*           and at least\n*           ( 1 + ( n - 1 )*abs( INCY ) ) otherwise.\n*           Before entry with BETA non-zero, the incremented array Y\n*           must contain the vector y. On exit, Y is overwritten by the\n*           updated vector y.\n*\n*  INCY   (input) INTEGER\n*           On entry, INCY specifies the increment for the elements of\n*           Y. INCY must not be zero.\n*           Unchanged on exit.\n*\n*\n*  Level 2 Blas routine.\n*\n\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 11)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 11)", argc);
  rb_trans = argv[0];
  rb_m = argv[1];
  rb_kl = argv[2];
  rb_ku = argv[3];
  rb_alpha = argv[4];
  rb_ab = argv[5];
  rb_x = argv[6];
  rb_incx = argv[7];
  rb_beta = argv[8];
  rb_y = argv[9];
  rb_incy = argv[10];

  if (!NA_IsNArray(rb_ab))
    rb_raise(rb_eArgError, "ab (6th argument) must be NArray");
  if (NA_RANK(rb_ab) != 2)
    rb_raise(rb_eArgError, "rank of ab (6th argument) must be %d", 2);
  n = NA_SHAPE1(rb_ab);
  ldab = NA_SHAPE0(rb_ab);
  if (ldab != (MAX(1,m)))
    rb_raise(rb_eRuntimeError, "shape 0 of ab must be %d", MAX(1,m));
  if (NA_TYPE(rb_ab) != NA_SFLOAT)
    rb_ab = na_change_type(rb_ab, NA_SFLOAT);
  ab = NA_PTR_TYPE(rb_ab, real*);
  kl = NUM2INT(rb_kl);
  trans = NUM2INT(rb_trans);
  m = NUM2INT(rb_m);
  ku = NUM2INT(rb_ku);
  beta = (real)NUM2DBL(rb_beta);
  incy = NUM2INT(rb_incy);
  alpha = (real)NUM2DBL(rb_alpha);
  incx = NUM2INT(rb_incx);
  ldab = MAX(1,m);
  if (!NA_IsNArray(rb_y))
    rb_raise(rb_eArgError, "y (10th argument) must be NArray");
  if (NA_RANK(rb_y) != 1)
    rb_raise(rb_eArgError, "rank of y (10th argument) must be %d", 1);
  if (NA_SHAPE0(rb_y) != (((lsame_(&trans,"N")) || (lsame_(&trans,"n"))) ? 1 + ( m - 1 )*abs( incy ) : 1 + ( n - 1 )*abs( incy )))
    rb_raise(rb_eRuntimeError, "shape 0 of y must be %d", ((lsame_(&trans,"N")) || (lsame_(&trans,"n"))) ? 1 + ( m - 1 )*abs( incy ) : 1 + ( n - 1 )*abs( incy ));
  if (NA_TYPE(rb_y) != NA_SFLOAT)
    rb_y = na_change_type(rb_y, NA_SFLOAT);
  y = NA_PTR_TYPE(rb_y, real*);
  if (!NA_IsNArray(rb_x))
    rb_raise(rb_eArgError, "x (7th argument) must be NArray");
  if (NA_RANK(rb_x) != 1)
    rb_raise(rb_eArgError, "rank of x (7th argument) must be %d", 1);
  if (NA_SHAPE0(rb_x) != (((lsame_(&trans,"N")) || (lsame_(&trans,"n"))) ? 1 + ( n - 1 )*abs( incx ) : 1 + ( m - 1 )*abs( incx )))
    rb_raise(rb_eRuntimeError, "shape 0 of x must be %d", ((lsame_(&trans,"N")) || (lsame_(&trans,"n"))) ? 1 + ( n - 1 )*abs( incx ) : 1 + ( m - 1 )*abs( incx ));
  if (NA_TYPE(rb_x) != NA_SFLOAT)
    rb_x = na_change_type(rb_x, NA_SFLOAT);
  x = NA_PTR_TYPE(rb_x, real*);
  {
    int shape[1];
    shape[0] = ((lsame_(&trans,"N")) || (lsame_(&trans,"n"))) ? 1 + ( m - 1 )*abs( incy ) : 1 + ( n - 1 )*abs( incy );
    rb_y_out__ = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  y_out__ = NA_PTR_TYPE(rb_y_out__, real*);
  MEMCPY(y_out__, y, real, NA_TOTAL(rb_y));
  rb_y = rb_y_out__;
  y = y_out__;

  cla_gbamv_(&trans, &m, &n, &kl, &ku, &alpha, ab, &ldab, x, &incx, &beta, y, &incy);

  return rb_y;
}

void
init_lapack_cla_gbamv(VALUE mLapack){
  rb_define_module_function(mLapack, "cla_gbamv", rb_cla_gbamv, -1);
}

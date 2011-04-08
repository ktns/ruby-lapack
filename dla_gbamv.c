#include "rb_lapack.h"

extern VOID dla_gbamv_(integer *trans, integer *m, integer *n, integer *kl, integer *ku, doublereal *alpha, doublereal *ab, integer *ldab, doublereal *x, integer *incx, doublereal *beta, doublereal *y, integer *incy);

static VALUE sHelp, sUsage;

static VALUE
rblapack_dla_gbamv(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_trans;
  integer trans; 
  VALUE rblapack_m;
  integer m; 
  VALUE rblapack_n;
  integer n; 
  VALUE rblapack_kl;
  integer kl; 
  VALUE rblapack_ku;
  integer ku; 
  VALUE rblapack_alpha;
  doublereal alpha; 
  VALUE rblapack_ab;
  doublereal *ab; 
  VALUE rblapack_x;
  doublereal *x; 
  VALUE rblapack_incx;
  integer incx; 
  VALUE rblapack_beta;
  doublereal beta; 
  VALUE rblapack_y;
  doublereal *y; 
  VALUE rblapack_incy;
  integer incy; 
  VALUE rblapack_y_out__;
  doublereal *y_out__;

  integer ldab;
  integer lda;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  y = NumRu::Lapack.dla_gbamv( trans, m, n, kl, ku, alpha, ab, x, incx, beta, y, incy, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE DLA_GBAMV( TRANS, M, N, KL, KU, ALPHA, AB, LDAB, X, INCX, BETA, Y, INCY )\n\n*  Purpose\n*  =======\n*\n*  DLA_GBAMV  performs one of the matrix-vector operations\n*\n*          y := alpha*abs(A)*abs(x) + beta*abs(y),\n*     or   y := alpha*abs(A)'*abs(x) + beta*abs(y),\n*\n*  where alpha and beta are scalars, x and y are vectors and A is an\n*  m by n matrix.\n*\n*  This function is primarily used in calculating error bounds.\n*  To protect against underflow during evaluation, components in\n*  the resulting vector are perturbed away from zero by (N+1)\n*  times the underflow threshold.  To prevent unnecessarily large\n*  errors for block-structure embedded in general matrices,\n*  \"symbolically\" zero components are not perturbed.  A zero\n*  entry is considered \"symbolic\" if all multiplications involved\n*  in computing that entry have at least one zero multiplicand.\n*\n\n*  Arguments\n*  ==========\n*\n*  TRANS   (input) INTEGER\n*           On entry, TRANS specifies the operation to be performed as\n*           follows:\n*\n*             BLAS_NO_TRANS      y := alpha*abs(A)*abs(x) + beta*abs(y)\n*             BLAS_TRANS         y := alpha*abs(A')*abs(x) + beta*abs(y)\n*             BLAS_CONJ_TRANS    y := alpha*abs(A')*abs(x) + beta*abs(y)\n*\n*           Unchanged on exit.\n*\n*  M       (input) INTEGER\n*           On entry, M specifies the number of rows of the matrix A.\n*           M must be at least zero.\n*           Unchanged on exit.\n*\n*  N       (input) INTEGER\n*           On entry, N specifies the number of columns of the matrix A.\n*           N must be at least zero.\n*           Unchanged on exit.\n*\n*  KL      (input) INTEGER\n*           The number of subdiagonals within the band of A.  KL >= 0.\n*\n*  KU      (input) INTEGER\n*           The number of superdiagonals within the band of A.  KU >= 0.\n*\n*  ALPHA  - DOUBLE PRECISION\n*           On entry, ALPHA specifies the scalar alpha.\n*           Unchanged on exit.\n*\n*  A      - DOUBLE PRECISION   array of DIMENSION ( LDA, n )\n*           Before entry, the leading m by n part of the array A must\n*           contain the matrix of coefficients.\n*           Unchanged on exit.\n*\n*  LDA     (input) INTEGER\n*           On entry, LDA specifies the first dimension of A as declared\n*           in the calling (sub) program. LDA must be at least\n*           max( 1, m ).\n*           Unchanged on exit.\n*\n*  X       (input) DOUBLE PRECISION array, dimension\n*           ( 1 + ( n - 1 )*abs( INCX ) ) when TRANS = 'N' or 'n'\n*           and at least\n*           ( 1 + ( m - 1 )*abs( INCX ) ) otherwise.\n*           Before entry, the incremented array X must contain the\n*           vector x.\n*           Unchanged on exit.\n*\n*  INCX    (input) INTEGER\n*           On entry, INCX specifies the increment for the elements of\n*           X. INCX must not be zero.\n*           Unchanged on exit.\n*\n*  BETA   - DOUBLE PRECISION\n*           On entry, BETA specifies the scalar beta. When BETA is\n*           supplied as zero then Y need not be set on input.\n*           Unchanged on exit.\n*\n*  Y       (input/output) DOUBLE PRECISION  array, dimension\n*           ( 1 + ( m - 1 )*abs( INCY ) ) when TRANS = 'N' or 'n'\n*           and at least\n*           ( 1 + ( n - 1 )*abs( INCY ) ) otherwise.\n*           Before entry with BETA non-zero, the incremented array Y\n*           must contain the vector y. On exit, Y is overwritten by the\n*           updated vector y.\n*\n*  INCY    (input) INTEGER\n*           On entry, INCY specifies the increment for the elements of\n*           Y. INCY must not be zero.\n*           Unchanged on exit.\n*\n*\n*  Level 2 Blas routine.\n*\n\n*  =====================================================================\n\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  y = NumRu::Lapack.dla_gbamv( trans, m, n, kl, ku, alpha, ab, x, incx, beta, y, incy, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 12)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 12)", argc);
  rblapack_trans = argv[0];
  rblapack_m = argv[1];
  rblapack_n = argv[2];
  rblapack_kl = argv[3];
  rblapack_ku = argv[4];
  rblapack_alpha = argv[5];
  rblapack_ab = argv[6];
  rblapack_x = argv[7];
  rblapack_incx = argv[8];
  rblapack_beta = argv[9];
  rblapack_y = argv[10];
  rblapack_incy = argv[11];
  if (rb_options != Qnil) {
  }

  if (!NA_IsNArray(rblapack_ab))
    rb_raise(rb_eArgError, "ab (7th argument) must be NArray");
  if (NA_RANK(rblapack_ab) != 1)
    rb_raise(rb_eArgError, "rank of ab (7th argument) must be %d", 1);
  ldab = NA_SHAPE0(rblapack_ab);
  if (NA_TYPE(rblapack_ab) != NA_DFLOAT)
    rblapack_ab = na_change_type(rblapack_ab, NA_DFLOAT);
  ab = NA_PTR_TYPE(rblapack_ab, doublereal*);
  kl = NUM2INT(rblapack_kl);
  trans = NUM2INT(rblapack_trans);
  m = NUM2INT(rblapack_m);
  ku = NUM2INT(rblapack_ku);
  n = NUM2INT(rblapack_n);
  alpha = NUM2DBL(rblapack_alpha);
  incx = NUM2INT(rblapack_incx);
  incy = NUM2INT(rblapack_incy);
  beta = NUM2DBL(rblapack_beta);
  lda = max( 1, m );
  if (!NA_IsNArray(rblapack_y))
    rb_raise(rb_eArgError, "y (11th argument) must be NArray");
  if (NA_RANK(rblapack_y) != 1)
    rb_raise(rb_eArgError, "rank of y (11th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_y) != (((lsame_(&trans,"N")) || (lsame_(&trans,"n"))) ? 1 + ( m - 1 )*abs( incy ) : 1 + ( n - 1 )*abs( incy )))
    rb_raise(rb_eRuntimeError, "shape 0 of y must be %d", ((lsame_(&trans,"N")) || (lsame_(&trans,"n"))) ? 1 + ( m - 1 )*abs( incy ) : 1 + ( n - 1 )*abs( incy ));
  if (NA_TYPE(rblapack_y) != NA_DFLOAT)
    rblapack_y = na_change_type(rblapack_y, NA_DFLOAT);
  y = NA_PTR_TYPE(rblapack_y, doublereal*);
  if (!NA_IsNArray(rblapack_x))
    rb_raise(rb_eArgError, "x (8th argument) must be NArray");
  if (NA_RANK(rblapack_x) != 1)
    rb_raise(rb_eArgError, "rank of x (8th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_x) != (((lsame_(&trans,"N")) || (lsame_(&trans,"n"))) ? 1 + ( n - 1 )*abs( incx ) : 1 + ( m - 1 )*abs( incx )))
    rb_raise(rb_eRuntimeError, "shape 0 of x must be %d", ((lsame_(&trans,"N")) || (lsame_(&trans,"n"))) ? 1 + ( n - 1 )*abs( incx ) : 1 + ( m - 1 )*abs( incx ));
  if (NA_TYPE(rblapack_x) != NA_DFLOAT)
    rblapack_x = na_change_type(rblapack_x, NA_DFLOAT);
  x = NA_PTR_TYPE(rblapack_x, doublereal*);
  {
    int shape[1];
    shape[0] = ((lsame_(&trans,"N")) || (lsame_(&trans,"n"))) ? 1 + ( m - 1 )*abs( incy ) : 1 + ( n - 1 )*abs( incy );
    rblapack_y_out__ = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  y_out__ = NA_PTR_TYPE(rblapack_y_out__, doublereal*);
  MEMCPY(y_out__, y, doublereal, NA_TOTAL(rblapack_y));
  rblapack_y = rblapack_y_out__;
  y = y_out__;

  dla_gbamv_(&trans, &m, &n, &kl, &ku, &alpha, ab, &ldab, x, &incx, &beta, y, &incy);

  return rblapack_y;
}

void
init_lapack_dla_gbamv(VALUE mLapack){
  rb_define_module_function(mLapack, "dla_gbamv", rblapack_dla_gbamv, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}

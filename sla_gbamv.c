#include "rb_lapack.h"

extern VOID sla_gbamv_(integer *trans, integer *m, integer *n, integer *kl, integer *ku, real *alpha, real *ab, integer *ldab, real *x, integer *incx, real *beta, real *y, integer *incy);

static VALUE
rb_sla_gbamv(int argc, VALUE *argv, VALUE self){
  VALUE rb_trans;
  integer trans; 
  VALUE rb_m;
  integer m; 
  VALUE rb_n;
  integer n; 
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
  integer lda;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  y = NumRu::Lapack.sla_gbamv( trans, m, n, kl, ku, alpha, ab, x, incx, beta, y, incy)\n    or\n  NumRu::Lapack.sla_gbamv  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 12)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 12)", argc);
  rb_trans = argv[0];
  rb_m = argv[1];
  rb_n = argv[2];
  rb_kl = argv[3];
  rb_ku = argv[4];
  rb_alpha = argv[5];
  rb_ab = argv[6];
  rb_x = argv[7];
  rb_incx = argv[8];
  rb_beta = argv[9];
  rb_y = argv[10];
  rb_incy = argv[11];

  if (!NA_IsNArray(rb_ab))
    rb_raise(rb_eArgError, "ab (7th argument) must be NArray");
  if (NA_RANK(rb_ab) != 1)
    rb_raise(rb_eArgError, "rank of ab (7th argument) must be %d", 1);
  ldab = NA_SHAPE0(rb_ab);
  if (NA_TYPE(rb_ab) != NA_SFLOAT)
    rb_ab = na_change_type(rb_ab, NA_SFLOAT);
  ab = NA_PTR_TYPE(rb_ab, real*);
  kl = NUM2INT(rb_kl);
  trans = NUM2INT(rb_trans);
  m = NUM2INT(rb_m);
  ku = NUM2INT(rb_ku);
  n = NUM2INT(rb_n);
  alpha = (real)NUM2DBL(rb_alpha);
  incx = NUM2INT(rb_incx);
  incy = NUM2INT(rb_incy);
  beta = (real)NUM2DBL(rb_beta);
  lda = MAX(1,m);
  if (!NA_IsNArray(rb_y))
    rb_raise(rb_eArgError, "y (11th argument) must be NArray");
  if (NA_RANK(rb_y) != 1)
    rb_raise(rb_eArgError, "rank of y (11th argument) must be %d", 1);
  if (NA_SHAPE0(rb_y) != (((lsame_(&trans,"N")) || (lsame_(&trans,"n"))) ? 1 + ( m - 1 )*abs( incy ) : 1 + ( n - 1 )*abs( incy )))
    rb_raise(rb_eRuntimeError, "shape 0 of y must be %d", ((lsame_(&trans,"N")) || (lsame_(&trans,"n"))) ? 1 + ( m - 1 )*abs( incy ) : 1 + ( n - 1 )*abs( incy ));
  if (NA_TYPE(rb_y) != NA_SFLOAT)
    rb_y = na_change_type(rb_y, NA_SFLOAT);
  y = NA_PTR_TYPE(rb_y, real*);
  if (!NA_IsNArray(rb_x))
    rb_raise(rb_eArgError, "x (8th argument) must be NArray");
  if (NA_RANK(rb_x) != 1)
    rb_raise(rb_eArgError, "rank of x (8th argument) must be %d", 1);
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

  sla_gbamv_(&trans, &m, &n, &kl, &ku, &alpha, ab, &ldab, x, &incx, &beta, y, &incy);

  return rb_y;
}

void
init_lapack_sla_gbamv(VALUE mLapack){
  rb_define_module_function(mLapack, "sla_gbamv", rb_sla_gbamv, -1);
}

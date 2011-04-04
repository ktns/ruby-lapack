#include "rb_lapack.h"

extern VOID zla_gbamv_(integer *trans, integer *m, integer *n, integer *kl, integer *ku, doublereal *alpha, doublereal *ab, integer *ldab, doublereal *x, integer *incx, doublereal *beta, doublereal *y, integer *incy);

static VALUE
rb_zla_gbamv(int argc, VALUE *argv, VALUE self){
  VALUE rb_trans;
  integer trans; 
  VALUE rb_m;
  integer m; 
  VALUE rb_kl;
  integer kl; 
  VALUE rb_ku;
  integer ku; 
  VALUE rb_alpha;
  doublereal alpha; 
  VALUE rb_ab;
  doublereal *ab; 
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

  integer ldab;
  integer n;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  y = NumRu::Lapack.zla_gbamv( trans, m, kl, ku, alpha, ab, x, incx, beta, y, incy)\n    or\n  NumRu::Lapack.zla_gbamv  # print help\n\n\nFORTRAN MANUAL\n\n");
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
  if (ldab != (MAX(1, m)))
    rb_raise(rb_eRuntimeError, "shape 0 of ab must be %d", MAX(1, m));
  if (NA_TYPE(rb_ab) != NA_DFLOAT)
    rb_ab = na_change_type(rb_ab, NA_DFLOAT);
  ab = NA_PTR_TYPE(rb_ab, doublereal*);
  kl = NUM2INT(rb_kl);
  trans = NUM2INT(rb_trans);
  m = NUM2INT(rb_m);
  ku = NUM2INT(rb_ku);
  beta = NUM2DBL(rb_beta);
  incy = NUM2INT(rb_incy);
  alpha = NUM2DBL(rb_alpha);
  incx = NUM2INT(rb_incx);
  ldab = MAX(1, m);
  if (!NA_IsNArray(rb_y))
    rb_raise(rb_eArgError, "y (10th argument) must be NArray");
  if (NA_RANK(rb_y) != 1)
    rb_raise(rb_eArgError, "rank of y (10th argument) must be %d", 1);
  if (NA_SHAPE0(rb_y) != (((lsame_(&trans,"N")) || (lsame_(&trans,"n"))) ? 1 + ( m - 1 )*abs( incy ) : 1 + ( n - 1 )*abs( incy )))
    rb_raise(rb_eRuntimeError, "shape 0 of y must be %d", ((lsame_(&trans,"N")) || (lsame_(&trans,"n"))) ? 1 + ( m - 1 )*abs( incy ) : 1 + ( n - 1 )*abs( incy ));
  if (NA_TYPE(rb_y) != NA_DFLOAT)
    rb_y = na_change_type(rb_y, NA_DFLOAT);
  y = NA_PTR_TYPE(rb_y, doublereal*);
  if (!NA_IsNArray(rb_x))
    rb_raise(rb_eArgError, "x (7th argument) must be NArray");
  if (NA_RANK(rb_x) != 1)
    rb_raise(rb_eArgError, "rank of x (7th argument) must be %d", 1);
  if (NA_SHAPE0(rb_x) != (((lsame_(&trans,"N")) || (lsame_(&trans,"n"))) ? 1 + ( n - 1 )*abs( incx ) : 1 + ( m - 1 )*abs( incx )))
    rb_raise(rb_eRuntimeError, "shape 0 of x must be %d", ((lsame_(&trans,"N")) || (lsame_(&trans,"n"))) ? 1 + ( n - 1 )*abs( incx ) : 1 + ( m - 1 )*abs( incx ));
  if (NA_TYPE(rb_x) != NA_DFLOAT)
    rb_x = na_change_type(rb_x, NA_DFLOAT);
  x = NA_PTR_TYPE(rb_x, doublereal*);
  {
    int shape[1];
    shape[0] = ((lsame_(&trans,"N")) || (lsame_(&trans,"n"))) ? 1 + ( m - 1 )*abs( incy ) : 1 + ( n - 1 )*abs( incy );
    rb_y_out__ = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  y_out__ = NA_PTR_TYPE(rb_y_out__, doublereal*);
  MEMCPY(y_out__, y, doublereal, NA_TOTAL(rb_y));
  rb_y = rb_y_out__;
  y = y_out__;

  zla_gbamv_(&trans, &m, &n, &kl, &ku, &alpha, ab, &ldab, x, &incx, &beta, y, &incy);

  return rb_y;
}

void
init_lapack_zla_gbamv(VALUE mLapack){
  rb_define_module_function(mLapack, "zla_gbamv", rb_zla_gbamv, -1);
}
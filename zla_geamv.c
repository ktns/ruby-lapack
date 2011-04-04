#include "rb_lapack.h"

extern VOID zla_geamv_(integer *trans, integer *m, integer *n, doublereal *alpha, doublereal *a, integer *lda, doublereal *x, integer *incx, doublereal *beta, doublereal *y, integer *incy);

static VALUE
rb_zla_geamv(int argc, VALUE *argv, VALUE self){
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
    printf("%s\n", "USAGE:\n  y = NumRu::Lapack.zla_geamv( trans, m, alpha, a, x, incx, beta, y, incy)\n    or\n  NumRu::Lapack.zla_geamv  # print help\n\n\nFORTRAN MANUAL\n\n");
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
  if (NA_SHAPE0(rb_y) != (((lsame_(&trans,"N")) || (lsame_(&trans,"n"))) ? 1 + ( m - 1 )*abs( incy ) : 1 + ( n - 1 )*abs( incy )))
    rb_raise(rb_eRuntimeError, "shape 0 of y must be %d", ((lsame_(&trans,"N")) || (lsame_(&trans,"n"))) ? 1 + ( m - 1 )*abs( incy ) : 1 + ( n - 1 )*abs( incy ));
  if (NA_TYPE(rb_y) != NA_DFLOAT)
    rb_y = na_change_type(rb_y, NA_DFLOAT);
  y = NA_PTR_TYPE(rb_y, doublereal*);
  if (!NA_IsNArray(rb_x))
    rb_raise(rb_eArgError, "x (5th argument) must be NArray");
  if (NA_RANK(rb_x) != 1)
    rb_raise(rb_eArgError, "rank of x (5th argument) must be %d", 1);
  if (NA_SHAPE0(rb_x) != (lsame_(&trans,"N") ? 1+(n-1)*abs(incx) : 1+(m-1)*abs(incx)))
    rb_raise(rb_eRuntimeError, "shape 0 of x must be %d", lsame_(&trans,"N") ? 1+(n-1)*abs(incx) : 1+(m-1)*abs(incx));
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

  zla_geamv_(&trans, &m, &n, &alpha, a, &lda, x, &incx, &beta, y, &incy);

  return rb_y;
}

void
init_lapack_zla_geamv(VALUE mLapack){
  rb_define_module_function(mLapack, "zla_geamv", rb_zla_geamv, -1);
}

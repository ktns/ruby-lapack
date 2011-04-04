#include "rb_lapack.h"

extern VOID zlaein_(logical *rightv, logical *noinit, integer *n, doublecomplex *h, integer *ldh, doublecomplex *w, doublecomplex *v, doublecomplex *b, integer *ldb, doublereal *rwork, doublereal *eps3, doublereal *smlnum, integer *info);

static VALUE
rb_zlaein(int argc, VALUE *argv, VALUE self){
  VALUE rb_rightv;
  logical rightv; 
  VALUE rb_noinit;
  logical noinit; 
  VALUE rb_h;
  doublecomplex *h; 
  VALUE rb_w;
  doublecomplex w; 
  VALUE rb_v;
  doublecomplex *v; 
  VALUE rb_eps3;
  doublereal eps3; 
  VALUE rb_smlnum;
  doublereal smlnum; 
  VALUE rb_info;
  integer info; 
  VALUE rb_v_out__;
  doublecomplex *v_out__;
  doublecomplex *b;
  doublereal *rwork;

  integer ldh;
  integer n;
  integer ldb;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  info, v = NumRu::Lapack.zlaein( rightv, noinit, h, w, v, eps3, smlnum)\n    or\n  NumRu::Lapack.zlaein  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 7)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 7)", argc);
  rb_rightv = argv[0];
  rb_noinit = argv[1];
  rb_h = argv[2];
  rb_w = argv[3];
  rb_v = argv[4];
  rb_eps3 = argv[5];
  rb_smlnum = argv[6];

  smlnum = NUM2DBL(rb_smlnum);
  if (!NA_IsNArray(rb_v))
    rb_raise(rb_eArgError, "v (5th argument) must be NArray");
  if (NA_RANK(rb_v) != 1)
    rb_raise(rb_eArgError, "rank of v (5th argument) must be %d", 1);
  n = NA_SHAPE0(rb_v);
  if (NA_TYPE(rb_v) != NA_DCOMPLEX)
    rb_v = na_change_type(rb_v, NA_DCOMPLEX);
  v = NA_PTR_TYPE(rb_v, doublecomplex*);
  w.r = NUM2DBL(rb_funcall(rb_w, rb_intern("real"), 0));
  w.i = NUM2DBL(rb_funcall(rb_w, rb_intern("imag"), 0));
  eps3 = NUM2DBL(rb_eps3);
  noinit = (rb_noinit == Qtrue);
  if (!NA_IsNArray(rb_h))
    rb_raise(rb_eArgError, "h (3th argument) must be NArray");
  if (NA_RANK(rb_h) != 2)
    rb_raise(rb_eArgError, "rank of h (3th argument) must be %d", 2);
  if (NA_SHAPE1(rb_h) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of h must be the same as shape 0 of v");
  ldh = NA_SHAPE0(rb_h);
  if (NA_TYPE(rb_h) != NA_DCOMPLEX)
    rb_h = na_change_type(rb_h, NA_DCOMPLEX);
  h = NA_PTR_TYPE(rb_h, doublecomplex*);
  rightv = (rb_rightv == Qtrue);
  ldb = MAX(1,n);
  {
    int shape[1];
    shape[0] = n;
    rb_v_out__ = na_make_object(NA_DCOMPLEX, 1, shape, cNArray);
  }
  v_out__ = NA_PTR_TYPE(rb_v_out__, doublecomplex*);
  MEMCPY(v_out__, v, doublecomplex, NA_TOTAL(rb_v));
  rb_v = rb_v_out__;
  v = v_out__;
  b = ALLOC_N(doublecomplex, (ldb)*(n));
  rwork = ALLOC_N(doublereal, (n));

  zlaein_(&rightv, &noinit, &n, h, &ldh, &w, v, b, &ldb, rwork, &eps3, &smlnum, &info);

  free(b);
  free(rwork);
  rb_info = INT2NUM(info);
  return rb_ary_new3(2, rb_info, rb_v);
}

void
init_lapack_zlaein(VALUE mLapack){
  rb_define_module_function(mLapack, "zlaein", rb_zlaein, -1);
}

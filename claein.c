#include "rb_lapack.h"

extern VOID claein_(logical *rightv, logical *noinit, integer *n, complex *h, integer *ldh, complex *w, complex *v, complex *b, integer *ldb, real *rwork, real *eps3, real *smlnum, integer *info);

static VALUE
rb_claein(int argc, VALUE *argv, VALUE self){
  VALUE rb_rightv;
  logical rightv; 
  VALUE rb_noinit;
  logical noinit; 
  VALUE rb_h;
  complex *h; 
  VALUE rb_w;
  complex w; 
  VALUE rb_v;
  complex *v; 
  VALUE rb_eps3;
  real eps3; 
  VALUE rb_smlnum;
  real smlnum; 
  VALUE rb_info;
  integer info; 
  VALUE rb_v_out__;
  complex *v_out__;
  complex *b;
  real *rwork;

  integer ldh;
  integer n;
  integer ldb;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  info, v = NumRu::Lapack.claein( rightv, noinit, h, w, v, eps3, smlnum)\n    or\n  NumRu::Lapack.claein  # print help\n\n\nFORTRAN MANUAL\n\n");
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

  smlnum = (real)NUM2DBL(rb_smlnum);
  if (!NA_IsNArray(rb_v))
    rb_raise(rb_eArgError, "v (5th argument) must be NArray");
  if (NA_RANK(rb_v) != 1)
    rb_raise(rb_eArgError, "rank of v (5th argument) must be %d", 1);
  n = NA_SHAPE0(rb_v);
  if (NA_TYPE(rb_v) != NA_SCOMPLEX)
    rb_v = na_change_type(rb_v, NA_SCOMPLEX);
  v = NA_PTR_TYPE(rb_v, complex*);
  w.r = (real)NUM2DBL(rb_funcall(rb_w, rb_intern("real"), 0));
  w.i = (real)NUM2DBL(rb_funcall(rb_w, rb_intern("imag"), 0));
  eps3 = (real)NUM2DBL(rb_eps3);
  noinit = (rb_noinit == Qtrue);
  if (!NA_IsNArray(rb_h))
    rb_raise(rb_eArgError, "h (3th argument) must be NArray");
  if (NA_RANK(rb_h) != 2)
    rb_raise(rb_eArgError, "rank of h (3th argument) must be %d", 2);
  if (NA_SHAPE1(rb_h) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of h must be the same as shape 0 of v");
  ldh = NA_SHAPE0(rb_h);
  if (NA_TYPE(rb_h) != NA_SCOMPLEX)
    rb_h = na_change_type(rb_h, NA_SCOMPLEX);
  h = NA_PTR_TYPE(rb_h, complex*);
  rightv = (rb_rightv == Qtrue);
  ldb = MAX(1,n);
  {
    int shape[1];
    shape[0] = n;
    rb_v_out__ = na_make_object(NA_SCOMPLEX, 1, shape, cNArray);
  }
  v_out__ = NA_PTR_TYPE(rb_v_out__, complex*);
  MEMCPY(v_out__, v, complex, NA_TOTAL(rb_v));
  rb_v = rb_v_out__;
  v = v_out__;
  b = ALLOC_N(complex, (ldb)*(n));
  rwork = ALLOC_N(real, (n));

  claein_(&rightv, &noinit, &n, h, &ldh, &w, v, b, &ldb, rwork, &eps3, &smlnum, &info);

  free(b);
  free(rwork);
  rb_info = INT2NUM(info);
  return rb_ary_new3(2, rb_info, rb_v);
}

void
init_lapack_claein(VALUE mLapack){
  rb_define_module_function(mLapack, "claein", rb_claein, -1);
}

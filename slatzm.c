#include "rb_lapack.h"

extern VOID slatzm_(char *side, integer *m, integer *n, real *v, integer *incv, real *tau, real *c1, real *c2, integer *ldc, real *work);

static VALUE
rb_slatzm(int argc, VALUE *argv, VALUE self){
  VALUE rb_side;
  char side; 
  VALUE rb_m;
  integer m; 
  VALUE rb_n;
  integer n; 
  VALUE rb_v;
  real *v; 
  VALUE rb_incv;
  integer incv; 
  VALUE rb_tau;
  real tau; 
  VALUE rb_c1;
  real *c1; 
  VALUE rb_c2;
  real *c2; 
  VALUE rb_c1_out__;
  real *c1_out__;
  VALUE rb_c2_out__;
  real *c2_out__;
  real *work;

  integer ldc;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  c1, c2 = NumRu::Lapack.slatzm( side, m, n, v, incv, tau, c1, c2)\n    or\n  NumRu::Lapack.slatzm  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 8)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 8)", argc);
  rb_side = argv[0];
  rb_m = argv[1];
  rb_n = argv[2];
  rb_v = argv[3];
  rb_incv = argv[4];
  rb_tau = argv[5];
  rb_c1 = argv[6];
  rb_c2 = argv[7];

  tau = (real)NUM2DBL(rb_tau);
  side = StringValueCStr(rb_side)[0];
  m = NUM2INT(rb_m);
  incv = NUM2INT(rb_incv);
  n = NUM2INT(rb_n);
  if (!NA_IsNArray(rb_c2))
    rb_raise(rb_eArgError, "c2 (8th argument) must be NArray");
  if (NA_RANK(rb_c2) != 2)
    rb_raise(rb_eArgError, "rank of c2 (8th argument) must be %d", 2);
  if (NA_SHAPE1(rb_c2) != (lsame_(&side,"L") ? n : lsame_(&side,"R") ? n-1 : 0))
    rb_raise(rb_eRuntimeError, "shape 1 of c2 must be %d", lsame_(&side,"L") ? n : lsame_(&side,"R") ? n-1 : 0);
  ldc = NA_SHAPE0(rb_c2);
  if (NA_TYPE(rb_c2) != NA_SFLOAT)
    rb_c2 = na_change_type(rb_c2, NA_SFLOAT);
  c2 = NA_PTR_TYPE(rb_c2, real*);
  if (!NA_IsNArray(rb_c1))
    rb_raise(rb_eArgError, "c1 (7th argument) must be NArray");
  if (NA_RANK(rb_c1) != 2)
    rb_raise(rb_eArgError, "rank of c1 (7th argument) must be %d", 2);
  if (NA_SHAPE1(rb_c1) != (lsame_(&side,"L") ? n : lsame_(&side,"R") ? 1 : 0))
    rb_raise(rb_eRuntimeError, "shape 1 of c1 must be %d", lsame_(&side,"L") ? n : lsame_(&side,"R") ? 1 : 0);
  if (NA_SHAPE0(rb_c1) != (lsame_(&side,"L") ? ldc : lsame_(&side,"R") ? m : 0))
    rb_raise(rb_eRuntimeError, "shape 0 of c1 must be %d", lsame_(&side,"L") ? ldc : lsame_(&side,"R") ? m : 0);
  if (NA_TYPE(rb_c1) != NA_SFLOAT)
    rb_c1 = na_change_type(rb_c1, NA_SFLOAT);
  c1 = NA_PTR_TYPE(rb_c1, real*);
  if (!NA_IsNArray(rb_v))
    rb_raise(rb_eArgError, "v (4th argument) must be NArray");
  if (NA_RANK(rb_v) != 1)
    rb_raise(rb_eArgError, "rank of v (4th argument) must be %d", 1);
  if (NA_SHAPE0(rb_v) != (1 + (m-1)*abs(incv)))
    rb_raise(rb_eRuntimeError, "shape 0 of v must be %d", 1 + (m-1)*abs(incv));
  if (NA_TYPE(rb_v) != NA_SFLOAT)
    rb_v = na_change_type(rb_v, NA_SFLOAT);
  v = NA_PTR_TYPE(rb_v, real*);
  {
    int shape[2];
    shape[0] = lsame_(&side,"L") ? ldc : lsame_(&side,"R") ? m : 0;
    shape[1] = lsame_(&side,"L") ? n : lsame_(&side,"R") ? 1 : 0;
    rb_c1_out__ = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  c1_out__ = NA_PTR_TYPE(rb_c1_out__, real*);
  MEMCPY(c1_out__, c1, real, NA_TOTAL(rb_c1));
  rb_c1 = rb_c1_out__;
  c1 = c1_out__;
  {
    int shape[2];
    shape[0] = ldc;
    shape[1] = lsame_(&side,"L") ? n : lsame_(&side,"R") ? n-1 : 0;
    rb_c2_out__ = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  c2_out__ = NA_PTR_TYPE(rb_c2_out__, real*);
  MEMCPY(c2_out__, c2, real, NA_TOTAL(rb_c2));
  rb_c2 = rb_c2_out__;
  c2 = c2_out__;
  work = ALLOC_N(real, (lsame_(&side,"L") ? n : lsame_(&side,"R") ? m : 0));

  slatzm_(&side, &m, &n, v, &incv, &tau, c1, c2, &ldc, work);

  free(work);
  return rb_ary_new3(2, rb_c1, rb_c2);
}

void
init_lapack_slatzm(VALUE mLapack){
  rb_define_module_function(mLapack, "slatzm", rb_slatzm, -1);
}

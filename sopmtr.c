#include "rb_lapack.h"

extern VOID sopmtr_(char *side, char *uplo, char *trans, integer *m, integer *n, real *ap, real *tau, real *c, integer *ldc, real *work, integer *info);

static VALUE
rb_sopmtr(int argc, VALUE *argv, VALUE self){
  VALUE rb_side;
  char side; 
  VALUE rb_uplo;
  char uplo; 
  VALUE rb_trans;
  char trans; 
  VALUE rb_m;
  integer m; 
  VALUE rb_ap;
  real *ap; 
  VALUE rb_tau;
  real *tau; 
  VALUE rb_c;
  real *c; 
  VALUE rb_info;
  integer info; 
  VALUE rb_c_out__;
  real *c_out__;
  real *work;

  integer ldc;
  integer n;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  info, c = NumRu::Lapack.sopmtr( side, uplo, trans, m, ap, tau, c)\n    or\n  NumRu::Lapack.sopmtr  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 7)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 7)", argc);
  rb_side = argv[0];
  rb_uplo = argv[1];
  rb_trans = argv[2];
  rb_m = argv[3];
  rb_ap = argv[4];
  rb_tau = argv[5];
  rb_c = argv[6];

  side = StringValueCStr(rb_side)[0];
  m = NUM2INT(rb_m);
  if (!NA_IsNArray(rb_c))
    rb_raise(rb_eArgError, "c (7th argument) must be NArray");
  if (NA_RANK(rb_c) != 2)
    rb_raise(rb_eArgError, "rank of c (7th argument) must be %d", 2);
  n = NA_SHAPE1(rb_c);
  ldc = NA_SHAPE0(rb_c);
  if (NA_TYPE(rb_c) != NA_SFLOAT)
    rb_c = na_change_type(rb_c, NA_SFLOAT);
  c = NA_PTR_TYPE(rb_c, real*);
  trans = StringValueCStr(rb_trans)[0];
  uplo = StringValueCStr(rb_uplo)[0];
  if (!NA_IsNArray(rb_tau))
    rb_raise(rb_eArgError, "tau (6th argument) must be NArray");
  if (NA_RANK(rb_tau) != 1)
    rb_raise(rb_eArgError, "rank of tau (6th argument) must be %d", 1);
  if (NA_SHAPE0(rb_tau) != (m-1))
    rb_raise(rb_eRuntimeError, "shape 0 of tau must be %d", m-1);
  if (NA_TYPE(rb_tau) != NA_SFLOAT)
    rb_tau = na_change_type(rb_tau, NA_SFLOAT);
  tau = NA_PTR_TYPE(rb_tau, real*);
  if (!NA_IsNArray(rb_ap))
    rb_raise(rb_eArgError, "ap (5th argument) must be NArray");
  if (NA_RANK(rb_ap) != 1)
    rb_raise(rb_eArgError, "rank of ap (5th argument) must be %d", 1);
  if (NA_SHAPE0(rb_ap) != (m*(m+1)/2))
    rb_raise(rb_eRuntimeError, "shape 0 of ap must be %d", m*(m+1)/2);
  if (NA_TYPE(rb_ap) != NA_SFLOAT)
    rb_ap = na_change_type(rb_ap, NA_SFLOAT);
  ap = NA_PTR_TYPE(rb_ap, real*);
  {
    int shape[2];
    shape[0] = ldc;
    shape[1] = n;
    rb_c_out__ = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  c_out__ = NA_PTR_TYPE(rb_c_out__, real*);
  MEMCPY(c_out__, c, real, NA_TOTAL(rb_c));
  rb_c = rb_c_out__;
  c = c_out__;
  work = ALLOC_N(real, (lsame_(&side,"L") ? n : lsame_(&side,"R") ? m : 0));

  sopmtr_(&side, &uplo, &trans, &m, &n, ap, tau, c, &ldc, work, &info);

  free(work);
  rb_info = INT2NUM(info);
  return rb_ary_new3(2, rb_info, rb_c);
}

void
init_lapack_sopmtr(VALUE mLapack){
  rb_define_module_function(mLapack, "sopmtr", rb_sopmtr, -1);
}

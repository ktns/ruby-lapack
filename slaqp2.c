#include "rb_lapack.h"

extern VOID slaqp2_(integer *m, integer *n, integer *offset, real *a, integer *lda, integer *jpvt, real *tau, real *vn1, real *vn2, real *work);

static VALUE
rb_slaqp2(int argc, VALUE *argv, VALUE self){
  VALUE rb_m;
  integer m; 
  VALUE rb_offset;
  integer offset; 
  VALUE rb_a;
  real *a; 
  VALUE rb_jpvt;
  integer *jpvt; 
  VALUE rb_vn1;
  real *vn1; 
  VALUE rb_vn2;
  real *vn2; 
  VALUE rb_tau;
  real *tau; 
  VALUE rb_a_out__;
  real *a_out__;
  VALUE rb_jpvt_out__;
  integer *jpvt_out__;
  VALUE rb_vn1_out__;
  real *vn1_out__;
  VALUE rb_vn2_out__;
  real *vn2_out__;
  real *work;

  integer lda;
  integer n;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  tau, a, jpvt, vn1, vn2 = NumRu::Lapack.slaqp2( m, offset, a, jpvt, vn1, vn2)\n    or\n  NumRu::Lapack.slaqp2  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 6)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 6)", argc);
  rb_m = argv[0];
  rb_offset = argv[1];
  rb_a = argv[2];
  rb_jpvt = argv[3];
  rb_vn1 = argv[4];
  rb_vn2 = argv[5];

  if (!NA_IsNArray(rb_a))
    rb_raise(rb_eArgError, "a (3th argument) must be NArray");
  if (NA_RANK(rb_a) != 2)
    rb_raise(rb_eArgError, "rank of a (3th argument) must be %d", 2);
  n = NA_SHAPE1(rb_a);
  lda = NA_SHAPE0(rb_a);
  if (NA_TYPE(rb_a) != NA_SFLOAT)
    rb_a = na_change_type(rb_a, NA_SFLOAT);
  a = NA_PTR_TYPE(rb_a, real*);
  m = NUM2INT(rb_m);
  if (!NA_IsNArray(rb_vn1))
    rb_raise(rb_eArgError, "vn1 (5th argument) must be NArray");
  if (NA_RANK(rb_vn1) != 1)
    rb_raise(rb_eArgError, "rank of vn1 (5th argument) must be %d", 1);
  if (NA_SHAPE0(rb_vn1) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of vn1 must be the same as shape 1 of a");
  if (NA_TYPE(rb_vn1) != NA_SFLOAT)
    rb_vn1 = na_change_type(rb_vn1, NA_SFLOAT);
  vn1 = NA_PTR_TYPE(rb_vn1, real*);
  if (!NA_IsNArray(rb_vn2))
    rb_raise(rb_eArgError, "vn2 (6th argument) must be NArray");
  if (NA_RANK(rb_vn2) != 1)
    rb_raise(rb_eArgError, "rank of vn2 (6th argument) must be %d", 1);
  if (NA_SHAPE0(rb_vn2) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of vn2 must be the same as shape 1 of a");
  if (NA_TYPE(rb_vn2) != NA_SFLOAT)
    rb_vn2 = na_change_type(rb_vn2, NA_SFLOAT);
  vn2 = NA_PTR_TYPE(rb_vn2, real*);
  if (!NA_IsNArray(rb_jpvt))
    rb_raise(rb_eArgError, "jpvt (4th argument) must be NArray");
  if (NA_RANK(rb_jpvt) != 1)
    rb_raise(rb_eArgError, "rank of jpvt (4th argument) must be %d", 1);
  if (NA_SHAPE0(rb_jpvt) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of jpvt must be the same as shape 1 of a");
  if (NA_TYPE(rb_jpvt) != NA_LINT)
    rb_jpvt = na_change_type(rb_jpvt, NA_LINT);
  jpvt = NA_PTR_TYPE(rb_jpvt, integer*);
  offset = NUM2INT(rb_offset);
  {
    int shape[1];
    shape[0] = MIN(m,n);
    rb_tau = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  tau = NA_PTR_TYPE(rb_tau, real*);
  {
    int shape[2];
    shape[0] = lda;
    shape[1] = n;
    rb_a_out__ = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  a_out__ = NA_PTR_TYPE(rb_a_out__, real*);
  MEMCPY(a_out__, a, real, NA_TOTAL(rb_a));
  rb_a = rb_a_out__;
  a = a_out__;
  {
    int shape[1];
    shape[0] = n;
    rb_jpvt_out__ = na_make_object(NA_LINT, 1, shape, cNArray);
  }
  jpvt_out__ = NA_PTR_TYPE(rb_jpvt_out__, integer*);
  MEMCPY(jpvt_out__, jpvt, integer, NA_TOTAL(rb_jpvt));
  rb_jpvt = rb_jpvt_out__;
  jpvt = jpvt_out__;
  {
    int shape[1];
    shape[0] = n;
    rb_vn1_out__ = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  vn1_out__ = NA_PTR_TYPE(rb_vn1_out__, real*);
  MEMCPY(vn1_out__, vn1, real, NA_TOTAL(rb_vn1));
  rb_vn1 = rb_vn1_out__;
  vn1 = vn1_out__;
  {
    int shape[1];
    shape[0] = n;
    rb_vn2_out__ = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  vn2_out__ = NA_PTR_TYPE(rb_vn2_out__, real*);
  MEMCPY(vn2_out__, vn2, real, NA_TOTAL(rb_vn2));
  rb_vn2 = rb_vn2_out__;
  vn2 = vn2_out__;
  work = ALLOC_N(real, (n));

  slaqp2_(&m, &n, &offset, a, &lda, jpvt, tau, vn1, vn2, work);

  free(work);
  return rb_ary_new3(5, rb_tau, rb_a, rb_jpvt, rb_vn1, rb_vn2);
}

void
init_lapack_slaqp2(VALUE mLapack){
  rb_define_module_function(mLapack, "slaqp2", rb_slaqp2, -1);
}

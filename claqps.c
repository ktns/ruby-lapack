#include "rb_lapack.h"

extern VOID claqps_(integer *m, integer *n, integer *offset, integer *nb, integer *kb, complex *a, integer *lda, integer *jpvt, complex *tau, real *vn1, real *vn2, complex *auxv, complex *f, integer *ldf);

static VALUE
rb_claqps(int argc, VALUE *argv, VALUE self){
  VALUE rb_m;
  integer m; 
  VALUE rb_offset;
  integer offset; 
  VALUE rb_a;
  complex *a; 
  VALUE rb_jpvt;
  integer *jpvt; 
  VALUE rb_vn1;
  real *vn1; 
  VALUE rb_vn2;
  real *vn2; 
  VALUE rb_auxv;
  complex *auxv; 
  VALUE rb_f;
  complex *f; 
  VALUE rb_kb;
  integer kb; 
  VALUE rb_tau;
  complex *tau; 
  VALUE rb_a_out__;
  complex *a_out__;
  VALUE rb_jpvt_out__;
  integer *jpvt_out__;
  VALUE rb_vn1_out__;
  real *vn1_out__;
  VALUE rb_vn2_out__;
  real *vn2_out__;
  VALUE rb_auxv_out__;
  complex *auxv_out__;
  VALUE rb_f_out__;
  complex *f_out__;

  integer lda;
  integer n;
  integer nb;
  integer ldf;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  kb, tau, a, jpvt, vn1, vn2, auxv, f = NumRu::Lapack.claqps( m, offset, a, jpvt, vn1, vn2, auxv, f)\n    or\n  NumRu::Lapack.claqps  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 8)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 8)", argc);
  rb_m = argv[0];
  rb_offset = argv[1];
  rb_a = argv[2];
  rb_jpvt = argv[3];
  rb_vn1 = argv[4];
  rb_vn2 = argv[5];
  rb_auxv = argv[6];
  rb_f = argv[7];

  if (!NA_IsNArray(rb_auxv))
    rb_raise(rb_eArgError, "auxv (7th argument) must be NArray");
  if (NA_RANK(rb_auxv) != 1)
    rb_raise(rb_eArgError, "rank of auxv (7th argument) must be %d", 1);
  nb = NA_SHAPE0(rb_auxv);
  if (NA_TYPE(rb_auxv) != NA_SCOMPLEX)
    rb_auxv = na_change_type(rb_auxv, NA_SCOMPLEX);
  auxv = NA_PTR_TYPE(rb_auxv, complex*);
  offset = NUM2INT(rb_offset);
  if (!NA_IsNArray(rb_a))
    rb_raise(rb_eArgError, "a (3th argument) must be NArray");
  if (NA_RANK(rb_a) != 2)
    rb_raise(rb_eArgError, "rank of a (3th argument) must be %d", 2);
  n = NA_SHAPE1(rb_a);
  lda = NA_SHAPE0(rb_a);
  if (NA_TYPE(rb_a) != NA_SCOMPLEX)
    rb_a = na_change_type(rb_a, NA_SCOMPLEX);
  a = NA_PTR_TYPE(rb_a, complex*);
  m = NUM2INT(rb_m);
  if (!NA_IsNArray(rb_jpvt))
    rb_raise(rb_eArgError, "jpvt (4th argument) must be NArray");
  if (NA_RANK(rb_jpvt) != 1)
    rb_raise(rb_eArgError, "rank of jpvt (4th argument) must be %d", 1);
  if (NA_SHAPE0(rb_jpvt) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of jpvt must be the same as shape 1 of a");
  if (NA_TYPE(rb_jpvt) != NA_LINT)
    rb_jpvt = na_change_type(rb_jpvt, NA_LINT);
  jpvt = NA_PTR_TYPE(rb_jpvt, integer*);
  if (!NA_IsNArray(rb_vn2))
    rb_raise(rb_eArgError, "vn2 (6th argument) must be NArray");
  if (NA_RANK(rb_vn2) != 1)
    rb_raise(rb_eArgError, "rank of vn2 (6th argument) must be %d", 1);
  if (NA_SHAPE0(rb_vn2) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of vn2 must be the same as shape 1 of a");
  if (NA_TYPE(rb_vn2) != NA_SFLOAT)
    rb_vn2 = na_change_type(rb_vn2, NA_SFLOAT);
  vn2 = NA_PTR_TYPE(rb_vn2, real*);
  if (!NA_IsNArray(rb_f))
    rb_raise(rb_eArgError, "f (8th argument) must be NArray");
  if (NA_RANK(rb_f) != 2)
    rb_raise(rb_eArgError, "rank of f (8th argument) must be %d", 2);
  if (NA_SHAPE1(rb_f) != nb)
    rb_raise(rb_eRuntimeError, "shape 1 of f must be the same as shape 0 of auxv");
  ldf = NA_SHAPE0(rb_f);
  if (NA_TYPE(rb_f) != NA_SCOMPLEX)
    rb_f = na_change_type(rb_f, NA_SCOMPLEX);
  f = NA_PTR_TYPE(rb_f, complex*);
  if (!NA_IsNArray(rb_vn1))
    rb_raise(rb_eArgError, "vn1 (5th argument) must be NArray");
  if (NA_RANK(rb_vn1) != 1)
    rb_raise(rb_eArgError, "rank of vn1 (5th argument) must be %d", 1);
  if (NA_SHAPE0(rb_vn1) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of vn1 must be the same as shape 1 of a");
  if (NA_TYPE(rb_vn1) != NA_SFLOAT)
    rb_vn1 = na_change_type(rb_vn1, NA_SFLOAT);
  vn1 = NA_PTR_TYPE(rb_vn1, real*);
  kb = nb;
  {
    int shape[1];
    shape[0] = kb;
    rb_tau = na_make_object(NA_SCOMPLEX, 1, shape, cNArray);
  }
  tau = NA_PTR_TYPE(rb_tau, complex*);
  {
    int shape[2];
    shape[0] = lda;
    shape[1] = n;
    rb_a_out__ = na_make_object(NA_SCOMPLEX, 2, shape, cNArray);
  }
  a_out__ = NA_PTR_TYPE(rb_a_out__, complex*);
  MEMCPY(a_out__, a, complex, NA_TOTAL(rb_a));
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
  {
    int shape[1];
    shape[0] = nb;
    rb_auxv_out__ = na_make_object(NA_SCOMPLEX, 1, shape, cNArray);
  }
  auxv_out__ = NA_PTR_TYPE(rb_auxv_out__, complex*);
  MEMCPY(auxv_out__, auxv, complex, NA_TOTAL(rb_auxv));
  rb_auxv = rb_auxv_out__;
  auxv = auxv_out__;
  {
    int shape[2];
    shape[0] = ldf;
    shape[1] = nb;
    rb_f_out__ = na_make_object(NA_SCOMPLEX, 2, shape, cNArray);
  }
  f_out__ = NA_PTR_TYPE(rb_f_out__, complex*);
  MEMCPY(f_out__, f, complex, NA_TOTAL(rb_f));
  rb_f = rb_f_out__;
  f = f_out__;

  claqps_(&m, &n, &offset, &nb, &kb, a, &lda, jpvt, tau, vn1, vn2, auxv, f, &ldf);

  rb_kb = INT2NUM(kb);
  return rb_ary_new3(8, rb_kb, rb_tau, rb_a, rb_jpvt, rb_vn1, rb_vn2, rb_auxv, rb_f);
}

void
init_lapack_claqps(VALUE mLapack){
  rb_define_module_function(mLapack, "claqps", rb_claqps, -1);
}

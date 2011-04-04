#include "rb_lapack.h"

extern VOID sorbdb_(char *trans, char *signs, integer *m, integer *p, integer *q, real *x11, integer *ldx11, real *x12, integer *ldx12, real *x21, integer *ldx21, real *x22, integer *ldx22, real *theta, real *phi, real *taup1, real *taup2, real *tauq1, real *tauq2, real *work, integer *lwork, integer *info);

static VALUE
rb_sorbdb(int argc, VALUE *argv, VALUE self){
  VALUE rb_trans;
  char trans; 
  VALUE rb_signs;
  char signs; 
  VALUE rb_m;
  integer m; 
  VALUE rb_x11;
  real *x11; 
  VALUE rb_x12;
  real *x12; 
  VALUE rb_x21;
  real *x21; 
  VALUE rb_x22;
  real *x22; 
  VALUE rb_lwork;
  integer lwork; 
  VALUE rb_theta;
  real *theta; 
  VALUE rb_phi;
  real *phi; 
  VALUE rb_taup1;
  real *taup1; 
  VALUE rb_taup2;
  real *taup2; 
  VALUE rb_tauq1;
  real *tauq1; 
  VALUE rb_tauq2;
  real *tauq2; 
  VALUE rb_info;
  integer info; 
  VALUE rb_x11_out__;
  real *x11_out__;
  VALUE rb_x12_out__;
  real *x12_out__;
  VALUE rb_x21_out__;
  real *x21_out__;
  VALUE rb_x22_out__;
  real *x22_out__;
  real *work;

  integer ldx11;
  integer q;
  integer ldx12;
  integer ldx21;
  integer ldx22;
  integer p;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  theta, phi, taup1, taup2, tauq1, tauq2, info, x11, x12, x21, x22 = NumRu::Lapack.sorbdb( trans, signs, m, x11, x12, x21, x22, lwork)\n    or\n  NumRu::Lapack.sorbdb  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 8)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 8)", argc);
  rb_trans = argv[0];
  rb_signs = argv[1];
  rb_m = argv[2];
  rb_x11 = argv[3];
  rb_x12 = argv[4];
  rb_x21 = argv[5];
  rb_x22 = argv[6];
  rb_lwork = argv[7];

  trans = StringValueCStr(rb_trans)[0];
  lwork = NUM2INT(rb_lwork);
  signs = StringValueCStr(rb_signs)[0];
  if (!NA_IsNArray(rb_x21))
    rb_raise(rb_eArgError, "x21 (6th argument) must be NArray");
  if (NA_RANK(rb_x21) != 2)
    rb_raise(rb_eArgError, "rank of x21 (6th argument) must be %d", 2);
  q = NA_SHAPE1(rb_x21);
  ldx21 = NA_SHAPE0(rb_x21);
  if (ldx21 != (p))
    rb_raise(rb_eRuntimeError, "shape 0 of x21 must be %d", p);
  p = ldx21;
  if (NA_TYPE(rb_x21) != NA_SFLOAT)
    rb_x21 = na_change_type(rb_x21, NA_SFLOAT);
  x21 = NA_PTR_TYPE(rb_x21, real*);
  if (!NA_IsNArray(rb_x11))
    rb_raise(rb_eArgError, "x11 (4th argument) must be NArray");
  if (NA_RANK(rb_x11) != 2)
    rb_raise(rb_eArgError, "rank of x11 (4th argument) must be %d", 2);
  if (NA_SHAPE1(rb_x11) != q)
    rb_raise(rb_eRuntimeError, "shape 1 of x11 must be the same as shape 1 of x21");
  ldx11 = NA_SHAPE0(rb_x11);
  if (ldx11 != (p))
    rb_raise(rb_eRuntimeError, "shape 0 of x11 must be %d", p);
  p = ldx11;
  if (NA_TYPE(rb_x11) != NA_SFLOAT)
    rb_x11 = na_change_type(rb_x11, NA_SFLOAT);
  x11 = NA_PTR_TYPE(rb_x11, real*);
  m = NUM2INT(rb_m);
  if (!NA_IsNArray(rb_x22))
    rb_raise(rb_eArgError, "x22 (7th argument) must be NArray");
  if (NA_RANK(rb_x22) != 2)
    rb_raise(rb_eArgError, "rank of x22 (7th argument) must be %d", 2);
  if (NA_SHAPE1(rb_x22) != (m-q))
    rb_raise(rb_eRuntimeError, "shape 1 of x22 must be %d", m-q);
  ldx22 = NA_SHAPE0(rb_x22);
  if (ldx22 != (p))
    rb_raise(rb_eRuntimeError, "shape 0 of x22 must be %d", p);
  p = ldx22;
  if (NA_TYPE(rb_x22) != NA_SFLOAT)
    rb_x22 = na_change_type(rb_x22, NA_SFLOAT);
  x22 = NA_PTR_TYPE(rb_x22, real*);
  if (!NA_IsNArray(rb_x12))
    rb_raise(rb_eArgError, "x12 (5th argument) must be NArray");
  if (NA_RANK(rb_x12) != 2)
    rb_raise(rb_eArgError, "rank of x12 (5th argument) must be %d", 2);
  if (NA_SHAPE1(rb_x12) != (m-q))
    rb_raise(rb_eRuntimeError, "shape 1 of x12 must be %d", m-q);
  ldx12 = NA_SHAPE0(rb_x12);
  if (ldx12 != (p))
    rb_raise(rb_eRuntimeError, "shape 0 of x12 must be %d", p);
  p = ldx12;
  if (NA_TYPE(rb_x12) != NA_SFLOAT)
    rb_x12 = na_change_type(rb_x12, NA_SFLOAT);
  x12 = NA_PTR_TYPE(rb_x12, real*);
  ldx12 = p;
  ldx22 = p;
  ldx21 = p;
  ldx11 = p;
  {
    int shape[1];
    shape[0] = q;
    rb_theta = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  theta = NA_PTR_TYPE(rb_theta, real*);
  {
    int shape[1];
    shape[0] = q-1;
    rb_phi = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  phi = NA_PTR_TYPE(rb_phi, real*);
  {
    int shape[1];
    shape[0] = p;
    rb_taup1 = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  taup1 = NA_PTR_TYPE(rb_taup1, real*);
  {
    int shape[1];
    shape[0] = m-p;
    rb_taup2 = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  taup2 = NA_PTR_TYPE(rb_taup2, real*);
  {
    int shape[1];
    shape[0] = q;
    rb_tauq1 = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  tauq1 = NA_PTR_TYPE(rb_tauq1, real*);
  {
    int shape[1];
    shape[0] = m-q;
    rb_tauq2 = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  tauq2 = NA_PTR_TYPE(rb_tauq2, real*);
  {
    int shape[2];
    shape[0] = ldx11;
    shape[1] = q;
    rb_x11_out__ = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  x11_out__ = NA_PTR_TYPE(rb_x11_out__, real*);
  MEMCPY(x11_out__, x11, real, NA_TOTAL(rb_x11));
  rb_x11 = rb_x11_out__;
  x11 = x11_out__;
  {
    int shape[2];
    shape[0] = ldx12;
    shape[1] = m-q;
    rb_x12_out__ = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  x12_out__ = NA_PTR_TYPE(rb_x12_out__, real*);
  MEMCPY(x12_out__, x12, real, NA_TOTAL(rb_x12));
  rb_x12 = rb_x12_out__;
  x12 = x12_out__;
  {
    int shape[2];
    shape[0] = ldx21;
    shape[1] = q;
    rb_x21_out__ = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  x21_out__ = NA_PTR_TYPE(rb_x21_out__, real*);
  MEMCPY(x21_out__, x21, real, NA_TOTAL(rb_x21));
  rb_x21 = rb_x21_out__;
  x21 = x21_out__;
  {
    int shape[2];
    shape[0] = ldx22;
    shape[1] = m-q;
    rb_x22_out__ = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  x22_out__ = NA_PTR_TYPE(rb_x22_out__, real*);
  MEMCPY(x22_out__, x22, real, NA_TOTAL(rb_x22));
  rb_x22 = rb_x22_out__;
  x22 = x22_out__;
  work = ALLOC_N(real, (MAX(1,lwork)));

  sorbdb_(&trans, &signs, &m, &p, &q, x11, &ldx11, x12, &ldx12, x21, &ldx21, x22, &ldx22, theta, phi, taup1, taup2, tauq1, tauq2, work, &lwork, &info);

  free(work);
  rb_info = INT2NUM(info);
  return rb_ary_new3(11, rb_theta, rb_phi, rb_taup1, rb_taup2, rb_tauq1, rb_tauq2, rb_info, rb_x11, rb_x12, rb_x21, rb_x22);
}

void
init_lapack_sorbdb(VALUE mLapack){
  rb_define_module_function(mLapack, "sorbdb", rb_sorbdb, -1);
}
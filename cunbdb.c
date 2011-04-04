#include "rb_lapack.h"

extern VOID cunbdb_(char *trans, char *signs, integer *m, integer *p, integer *q, complex *x11, integer *ldx11, complex *x12, integer *ldx12, complex *x21, integer *ldx21, complex *x22, integer *ldx22, real *theta, real *phi, complex *taup1, complex *taup2, complex *tauq1, complex *tauq2, complex *work, integer *lwork, integer *info);

static VALUE
rb_cunbdb(int argc, VALUE *argv, VALUE self){
  VALUE rb_trans;
  char trans; 
  VALUE rb_signs;
  char signs; 
  VALUE rb_m;
  integer m; 
  VALUE rb_x11;
  complex *x11; 
  VALUE rb_x12;
  complex *x12; 
  VALUE rb_x21;
  complex *x21; 
  VALUE rb_x22;
  complex *x22; 
  VALUE rb_lwork;
  integer lwork; 
  VALUE rb_theta;
  real *theta; 
  VALUE rb_phi;
  real *phi; 
  VALUE rb_taup1;
  complex *taup1; 
  VALUE rb_taup2;
  complex *taup2; 
  VALUE rb_tauq1;
  complex *tauq1; 
  VALUE rb_tauq2;
  complex *tauq2; 
  VALUE rb_info;
  integer info; 
  VALUE rb_x11_out__;
  complex *x11_out__;
  VALUE rb_x12_out__;
  complex *x12_out__;
  VALUE rb_x21_out__;
  complex *x21_out__;
  VALUE rb_x22_out__;
  complex *x22_out__;
  complex *work;

  integer ldx11;
  integer q;
  integer ldx12;
  integer ldx21;
  integer ldx22;
  integer p;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  theta, phi, taup1, taup2, tauq1, tauq2, info, x11, x12, x21, x22 = NumRu::Lapack.cunbdb( trans, signs, m, x11, x12, x21, x22, lwork)\n    or\n  NumRu::Lapack.cunbdb  # print help\n\n\nFORTRAN MANUAL\n\n");
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
  if (!NA_IsNArray(rb_x21))
    rb_raise(rb_eArgError, "x21 (6th argument) must be NArray");
  if (NA_RANK(rb_x21) != 2)
    rb_raise(rb_eArgError, "rank of x21 (6th argument) must be %d", 2);
  q = NA_SHAPE1(rb_x21);
  ldx21 = NA_SHAPE0(rb_x21);
  if (ldx21 != (p))
    rb_raise(rb_eRuntimeError, "shape 0 of x21 must be %d", p);
  p = ldx21;
  if (NA_TYPE(rb_x21) != NA_SCOMPLEX)
    rb_x21 = na_change_type(rb_x21, NA_SCOMPLEX);
  x21 = NA_PTR_TYPE(rb_x21, complex*);
  signs = StringValueCStr(rb_signs)[0];
  if (!NA_IsNArray(rb_x11))
    rb_raise(rb_eArgError, "x11 (4th argument) must be NArray");
  if (NA_RANK(rb_x11) != 2)
    rb_raise(rb_eArgError, "rank of x11 (4th argument) must be %d", 2);
  if (NA_SHAPE1(rb_x11) != q)
    rb_raise(rb_eRuntimeError, "shape 1 of x11 must be the same as shape 1 of x21");
  ldx11 = NA_SHAPE0(rb_x11);
  if (NA_TYPE(rb_x11) != NA_SCOMPLEX)
    rb_x11 = na_change_type(rb_x11, NA_SCOMPLEX);
  x11 = NA_PTR_TYPE(rb_x11, complex*);
  m = NUM2INT(rb_m);
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
  if (NA_TYPE(rb_x12) != NA_SCOMPLEX)
    rb_x12 = na_change_type(rb_x12, NA_SCOMPLEX);
  x12 = NA_PTR_TYPE(rb_x12, complex*);
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
  if (NA_TYPE(rb_x22) != NA_SCOMPLEX)
    rb_x22 = na_change_type(rb_x22, NA_SCOMPLEX);
  x22 = NA_PTR_TYPE(rb_x22, complex*);
  ldx12 = p;
  ldx21 = p;
  ldx22 = p;
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
    rb_taup1 = na_make_object(NA_SCOMPLEX, 1, shape, cNArray);
  }
  taup1 = NA_PTR_TYPE(rb_taup1, complex*);
  {
    int shape[1];
    shape[0] = m-p;
    rb_taup2 = na_make_object(NA_SCOMPLEX, 1, shape, cNArray);
  }
  taup2 = NA_PTR_TYPE(rb_taup2, complex*);
  {
    int shape[1];
    shape[0] = q;
    rb_tauq1 = na_make_object(NA_SCOMPLEX, 1, shape, cNArray);
  }
  tauq1 = NA_PTR_TYPE(rb_tauq1, complex*);
  {
    int shape[1];
    shape[0] = m-q;
    rb_tauq2 = na_make_object(NA_SCOMPLEX, 1, shape, cNArray);
  }
  tauq2 = NA_PTR_TYPE(rb_tauq2, complex*);
  {
    int shape[2];
    shape[0] = ldx11;
    shape[1] = q;
    rb_x11_out__ = na_make_object(NA_SCOMPLEX, 2, shape, cNArray);
  }
  x11_out__ = NA_PTR_TYPE(rb_x11_out__, complex*);
  MEMCPY(x11_out__, x11, complex, NA_TOTAL(rb_x11));
  rb_x11 = rb_x11_out__;
  x11 = x11_out__;
  {
    int shape[2];
    shape[0] = ldx12;
    shape[1] = m-q;
    rb_x12_out__ = na_make_object(NA_SCOMPLEX, 2, shape, cNArray);
  }
  x12_out__ = NA_PTR_TYPE(rb_x12_out__, complex*);
  MEMCPY(x12_out__, x12, complex, NA_TOTAL(rb_x12));
  rb_x12 = rb_x12_out__;
  x12 = x12_out__;
  {
    int shape[2];
    shape[0] = ldx21;
    shape[1] = q;
    rb_x21_out__ = na_make_object(NA_SCOMPLEX, 2, shape, cNArray);
  }
  x21_out__ = NA_PTR_TYPE(rb_x21_out__, complex*);
  MEMCPY(x21_out__, x21, complex, NA_TOTAL(rb_x21));
  rb_x21 = rb_x21_out__;
  x21 = x21_out__;
  {
    int shape[2];
    shape[0] = ldx22;
    shape[1] = m-q;
    rb_x22_out__ = na_make_object(NA_SCOMPLEX, 2, shape, cNArray);
  }
  x22_out__ = NA_PTR_TYPE(rb_x22_out__, complex*);
  MEMCPY(x22_out__, x22, complex, NA_TOTAL(rb_x22));
  rb_x22 = rb_x22_out__;
  x22 = x22_out__;
  work = ALLOC_N(complex, (MAX(1,lwork)));

  cunbdb_(&trans, &signs, &m, &p, &q, x11, &ldx11, x12, &ldx12, x21, &ldx21, x22, &ldx22, theta, phi, taup1, taup2, tauq1, tauq2, work, &lwork, &info);

  free(work);
  rb_info = INT2NUM(info);
  return rb_ary_new3(11, rb_theta, rb_phi, rb_taup1, rb_taup2, rb_tauq1, rb_tauq2, rb_info, rb_x11, rb_x12, rb_x21, rb_x22);
}

void
init_lapack_cunbdb(VALUE mLapack){
  rb_define_module_function(mLapack, "cunbdb", rb_cunbdb, -1);
}

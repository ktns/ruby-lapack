#include "rb_lapack.h"

extern VOID sbbcsd_(char *jobu1, char *jobu2, char *jobv1t, char *jobv2t, char *trans, integer *m, integer *p, integer *q, real *theta, real *phi, real *u1, integer *ldu1, real *u2, integer *ldu2, real *v1t, integer *ldv1t, real *v2t, integer *ldv2t, real *b11d, real *b11e, real *b12d, real *b12e, real *b21d, real *b21e, real *b22d, real *b22e, real *work, integer *lwork, integer *info);

static VALUE
rb_sbbcsd(int argc, VALUE *argv, VALUE self){
  VALUE rb_jobu1;
  char jobu1; 
  VALUE rb_jobu2;
  char jobu2; 
  VALUE rb_jobv1t;
  char jobv1t; 
  VALUE rb_jobv2t;
  char jobv2t; 
  VALUE rb_trans;
  char trans; 
  VALUE rb_m;
  integer m; 
  VALUE rb_theta;
  real *theta; 
  VALUE rb_phi;
  real *phi; 
  VALUE rb_u1;
  real *u1; 
  VALUE rb_u2;
  real *u2; 
  VALUE rb_v1t;
  real *v1t; 
  VALUE rb_v2t;
  real *v2t; 
  VALUE rb_lwork;
  integer lwork; 
  VALUE rb_b11d;
  real *b11d; 
  VALUE rb_b11e;
  real *b11e; 
  VALUE rb_b12d;
  real *b12d; 
  VALUE rb_b12e;
  real *b12e; 
  VALUE rb_b21d;
  real *b21d; 
  VALUE rb_b21e;
  real *b21e; 
  VALUE rb_b22d;
  real *b22d; 
  VALUE rb_b22e;
  real *b22e; 
  VALUE rb_info;
  integer info; 
  VALUE rb_theta_out__;
  real *theta_out__;
  VALUE rb_u1_out__;
  real *u1_out__;
  VALUE rb_u2_out__;
  real *u2_out__;
  VALUE rb_v1t_out__;
  real *v1t_out__;
  VALUE rb_v2t_out__;
  real *v2t_out__;
  real *work;

  integer q;
  integer ldu1;
  integer p;
  integer ldu2;
  integer ldv1t;
  integer ldv2t;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  b11d, b11e, b12d, b12e, b21d, b21e, b22d, b22e, info, theta, u1, u2, v1t, v2t = NumRu::Lapack.sbbcsd( jobu1, jobu2, jobv1t, jobv2t, trans, m, theta, phi, u1, u2, v1t, v2t, lwork)\n    or\n  NumRu::Lapack.sbbcsd  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 13)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 13)", argc);
  rb_jobu1 = argv[0];
  rb_jobu2 = argv[1];
  rb_jobv1t = argv[2];
  rb_jobv2t = argv[3];
  rb_trans = argv[4];
  rb_m = argv[5];
  rb_theta = argv[6];
  rb_phi = argv[7];
  rb_u1 = argv[8];
  rb_u2 = argv[9];
  rb_v1t = argv[10];
  rb_v2t = argv[11];
  rb_lwork = argv[12];

  if (!NA_IsNArray(rb_theta))
    rb_raise(rb_eArgError, "theta (7th argument) must be NArray");
  if (NA_RANK(rb_theta) != 1)
    rb_raise(rb_eArgError, "rank of theta (7th argument) must be %d", 1);
  q = NA_SHAPE0(rb_theta);
  if (NA_TYPE(rb_theta) != NA_SFLOAT)
    rb_theta = na_change_type(rb_theta, NA_SFLOAT);
  theta = NA_PTR_TYPE(rb_theta, real*);
  jobu1 = StringValueCStr(rb_jobu1)[0];
  trans = StringValueCStr(rb_trans)[0];
  m = NUM2INT(rb_m);
  jobu2 = StringValueCStr(rb_jobu2)[0];
  if (!NA_IsNArray(rb_v1t))
    rb_raise(rb_eArgError, "v1t (11th argument) must be NArray");
  if (NA_RANK(rb_v1t) != 2)
    rb_raise(rb_eArgError, "rank of v1t (11th argument) must be %d", 2);
  if (NA_SHAPE1(rb_v1t) != q)
    rb_raise(rb_eRuntimeError, "shape 1 of v1t must be the same as shape 0 of theta");
  ldv1t = NA_SHAPE0(rb_v1t);
  if (NA_TYPE(rb_v1t) != NA_SFLOAT)
    rb_v1t = na_change_type(rb_v1t, NA_SFLOAT);
  v1t = NA_PTR_TYPE(rb_v1t, real*);
  jobv2t = StringValueCStr(rb_jobv2t)[0];
  if (!NA_IsNArray(rb_u1))
    rb_raise(rb_eArgError, "u1 (9th argument) must be NArray");
  if (NA_RANK(rb_u1) != 2)
    rb_raise(rb_eArgError, "rank of u1 (9th argument) must be %d", 2);
  p = NA_SHAPE1(rb_u1);
  ldu1 = NA_SHAPE0(rb_u1);
  if (NA_TYPE(rb_u1) != NA_SFLOAT)
    rb_u1 = na_change_type(rb_u1, NA_SFLOAT);
  u1 = NA_PTR_TYPE(rb_u1, real*);
  jobv1t = StringValueCStr(rb_jobv1t)[0];
  lwork = NUM2INT(rb_lwork);
  if (!NA_IsNArray(rb_u2))
    rb_raise(rb_eArgError, "u2 (10th argument) must be NArray");
  if (NA_RANK(rb_u2) != 2)
    rb_raise(rb_eArgError, "rank of u2 (10th argument) must be %d", 2);
  if (NA_SHAPE1(rb_u2) != (m-p))
    rb_raise(rb_eRuntimeError, "shape 1 of u2 must be %d", m-p);
  ldu2 = NA_SHAPE0(rb_u2);
  if (NA_TYPE(rb_u2) != NA_SFLOAT)
    rb_u2 = na_change_type(rb_u2, NA_SFLOAT);
  u2 = NA_PTR_TYPE(rb_u2, real*);
  if (!NA_IsNArray(rb_phi))
    rb_raise(rb_eArgError, "phi (8th argument) must be NArray");
  if (NA_RANK(rb_phi) != 1)
    rb_raise(rb_eArgError, "rank of phi (8th argument) must be %d", 1);
  if (NA_SHAPE0(rb_phi) != (q-1))
    rb_raise(rb_eRuntimeError, "shape 0 of phi must be %d", q-1);
  if (NA_TYPE(rb_phi) != NA_SFLOAT)
    rb_phi = na_change_type(rb_phi, NA_SFLOAT);
  phi = NA_PTR_TYPE(rb_phi, real*);
  if (!NA_IsNArray(rb_v2t))
    rb_raise(rb_eArgError, "v2t (12th argument) must be NArray");
  if (NA_RANK(rb_v2t) != 2)
    rb_raise(rb_eArgError, "rank of v2t (12th argument) must be %d", 2);
  if (NA_SHAPE1(rb_v2t) != (m-q))
    rb_raise(rb_eRuntimeError, "shape 1 of v2t must be %d", m-q);
  ldv2t = NA_SHAPE0(rb_v2t);
  if (NA_TYPE(rb_v2t) != NA_SFLOAT)
    rb_v2t = na_change_type(rb_v2t, NA_SFLOAT);
  v2t = NA_PTR_TYPE(rb_v2t, real*);
  {
    int shape[1];
    shape[0] = q;
    rb_b11d = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  b11d = NA_PTR_TYPE(rb_b11d, real*);
  {
    int shape[1];
    shape[0] = q-1;
    rb_b11e = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  b11e = NA_PTR_TYPE(rb_b11e, real*);
  {
    int shape[1];
    shape[0] = q;
    rb_b12d = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  b12d = NA_PTR_TYPE(rb_b12d, real*);
  {
    int shape[1];
    shape[0] = q-1;
    rb_b12e = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  b12e = NA_PTR_TYPE(rb_b12e, real*);
  {
    int shape[1];
    shape[0] = q;
    rb_b21d = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  b21d = NA_PTR_TYPE(rb_b21d, real*);
  {
    int shape[1];
    shape[0] = q-1;
    rb_b21e = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  b21e = NA_PTR_TYPE(rb_b21e, real*);
  {
    int shape[1];
    shape[0] = q;
    rb_b22d = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  b22d = NA_PTR_TYPE(rb_b22d, real*);
  {
    int shape[1];
    shape[0] = q-1;
    rb_b22e = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  b22e = NA_PTR_TYPE(rb_b22e, real*);
  {
    int shape[1];
    shape[0] = q;
    rb_theta_out__ = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  theta_out__ = NA_PTR_TYPE(rb_theta_out__, real*);
  MEMCPY(theta_out__, theta, real, NA_TOTAL(rb_theta));
  rb_theta = rb_theta_out__;
  theta = theta_out__;
  {
    int shape[2];
    shape[0] = ldu1;
    shape[1] = p;
    rb_u1_out__ = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  u1_out__ = NA_PTR_TYPE(rb_u1_out__, real*);
  MEMCPY(u1_out__, u1, real, NA_TOTAL(rb_u1));
  rb_u1 = rb_u1_out__;
  u1 = u1_out__;
  {
    int shape[2];
    shape[0] = ldu2;
    shape[1] = m-p;
    rb_u2_out__ = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  u2_out__ = NA_PTR_TYPE(rb_u2_out__, real*);
  MEMCPY(u2_out__, u2, real, NA_TOTAL(rb_u2));
  rb_u2 = rb_u2_out__;
  u2 = u2_out__;
  {
    int shape[2];
    shape[0] = ldv1t;
    shape[1] = q;
    rb_v1t_out__ = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  v1t_out__ = NA_PTR_TYPE(rb_v1t_out__, real*);
  MEMCPY(v1t_out__, v1t, real, NA_TOTAL(rb_v1t));
  rb_v1t = rb_v1t_out__;
  v1t = v1t_out__;
  {
    int shape[2];
    shape[0] = ldv2t;
    shape[1] = m-q;
    rb_v2t_out__ = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  v2t_out__ = NA_PTR_TYPE(rb_v2t_out__, real*);
  MEMCPY(v2t_out__, v2t, real, NA_TOTAL(rb_v2t));
  rb_v2t = rb_v2t_out__;
  v2t = v2t_out__;
  work = ALLOC_N(real, (MAX(1,lwork)));

  sbbcsd_(&jobu1, &jobu2, &jobv1t, &jobv2t, &trans, &m, &p, &q, theta, phi, u1, &ldu1, u2, &ldu2, v1t, &ldv1t, v2t, &ldv2t, b11d, b11e, b12d, b12e, b21d, b21e, b22d, b22e, work, &lwork, &info);

  free(work);
  rb_info = INT2NUM(info);
  return rb_ary_new3(14, rb_b11d, rb_b11e, rb_b12d, rb_b12e, rb_b21d, rb_b21e, rb_b22d, rb_b22e, rb_info, rb_theta, rb_u1, rb_u2, rb_v1t, rb_v2t);
}

void
init_lapack_sbbcsd(VALUE mLapack){
  rb_define_module_function(mLapack, "sbbcsd", rb_sbbcsd, -1);
}

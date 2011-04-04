#include "rb_lapack.h"

extern VOID slasd1_(integer *nl, integer *nr, integer *sqre, real *d, real *alpha, real *beta, real *u, integer *ldu, real *vt, integer *ldvt, integer *idxq, integer *iwork, real *work, integer *info);

static VALUE
rb_slasd1(int argc, VALUE *argv, VALUE self){
  VALUE rb_nl;
  integer nl; 
  VALUE rb_nr;
  integer nr; 
  VALUE rb_sqre;
  integer sqre; 
  VALUE rb_d;
  real *d; 
  VALUE rb_alpha;
  real alpha; 
  VALUE rb_beta;
  real beta; 
  VALUE rb_u;
  real *u; 
  VALUE rb_vt;
  real *vt; 
  VALUE rb_idxq;
  integer *idxq; 
  VALUE rb_info;
  integer info; 
  VALUE rb_d_out__;
  real *d_out__;
  VALUE rb_u_out__;
  real *u_out__;
  VALUE rb_vt_out__;
  real *vt_out__;
  integer *iwork;
  real *work;

  integer ldu;
  integer n;
  integer ldvt;
  integer m;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  idxq, info, d, alpha, beta, u, vt = NumRu::Lapack.slasd1( nl, nr, sqre, d, alpha, beta, u, vt)\n    or\n  NumRu::Lapack.slasd1  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 8)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 8)", argc);
  rb_nl = argv[0];
  rb_nr = argv[1];
  rb_sqre = argv[2];
  rb_d = argv[3];
  rb_alpha = argv[4];
  rb_beta = argv[5];
  rb_u = argv[6];
  rb_vt = argv[7];

  alpha = (real)NUM2DBL(rb_alpha);
  if (!NA_IsNArray(rb_u))
    rb_raise(rb_eArgError, "u (7th argument) must be NArray");
  if (NA_RANK(rb_u) != 2)
    rb_raise(rb_eArgError, "rank of u (7th argument) must be %d", 2);
  n = NA_SHAPE1(rb_u);
  ldu = NA_SHAPE0(rb_u);
  if (NA_TYPE(rb_u) != NA_SFLOAT)
    rb_u = na_change_type(rb_u, NA_SFLOAT);
  u = NA_PTR_TYPE(rb_u, real*);
  nr = NUM2INT(rb_nr);
  beta = (real)NUM2DBL(rb_beta);
  nl = NUM2INT(rb_nl);
  if (!NA_IsNArray(rb_vt))
    rb_raise(rb_eArgError, "vt (8th argument) must be NArray");
  if (NA_RANK(rb_vt) != 2)
    rb_raise(rb_eArgError, "rank of vt (8th argument) must be %d", 2);
  m = NA_SHAPE1(rb_vt);
  if (m != (n + sqre))
    rb_raise(rb_eRuntimeError, "shape 1 of vt must be %d", n + sqre);
  ldvt = NA_SHAPE0(rb_vt);
  if (NA_TYPE(rb_vt) != NA_SFLOAT)
    rb_vt = na_change_type(rb_vt, NA_SFLOAT);
  vt = NA_PTR_TYPE(rb_vt, real*);
  sqre = NUM2INT(rb_sqre);
  if (!NA_IsNArray(rb_d))
    rb_raise(rb_eArgError, "d (4th argument) must be NArray");
  if (NA_RANK(rb_d) != 1)
    rb_raise(rb_eArgError, "rank of d (4th argument) must be %d", 1);
  if (NA_SHAPE0(rb_d) != (nl+nr+1))
    rb_raise(rb_eRuntimeError, "shape 0 of d must be %d", nl+nr+1);
  if (NA_TYPE(rb_d) != NA_SFLOAT)
    rb_d = na_change_type(rb_d, NA_SFLOAT);
  d = NA_PTR_TYPE(rb_d, real*);
  m = n + sqre;
  {
    int shape[1];
    shape[0] = n;
    rb_idxq = na_make_object(NA_LINT, 1, shape, cNArray);
  }
  idxq = NA_PTR_TYPE(rb_idxq, integer*);
  {
    int shape[1];
    shape[0] = nl+nr+1;
    rb_d_out__ = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  d_out__ = NA_PTR_TYPE(rb_d_out__, real*);
  MEMCPY(d_out__, d, real, NA_TOTAL(rb_d));
  rb_d = rb_d_out__;
  d = d_out__;
  {
    int shape[2];
    shape[0] = ldu;
    shape[1] = n;
    rb_u_out__ = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  u_out__ = NA_PTR_TYPE(rb_u_out__, real*);
  MEMCPY(u_out__, u, real, NA_TOTAL(rb_u));
  rb_u = rb_u_out__;
  u = u_out__;
  {
    int shape[2];
    shape[0] = ldvt;
    shape[1] = m;
    rb_vt_out__ = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  vt_out__ = NA_PTR_TYPE(rb_vt_out__, real*);
  MEMCPY(vt_out__, vt, real, NA_TOTAL(rb_vt));
  rb_vt = rb_vt_out__;
  vt = vt_out__;
  iwork = ALLOC_N(integer, (4*n));
  work = ALLOC_N(real, (3*pow(m,2)+2*m));

  slasd1_(&nl, &nr, &sqre, d, &alpha, &beta, u, &ldu, vt, &ldvt, idxq, iwork, work, &info);

  free(iwork);
  free(work);
  rb_info = INT2NUM(info);
  rb_alpha = rb_float_new((double)alpha);
  rb_beta = rb_float_new((double)beta);
  return rb_ary_new3(7, rb_idxq, rb_info, rb_d, rb_alpha, rb_beta, rb_u, rb_vt);
}

void
init_lapack_slasd1(VALUE mLapack){
  rb_define_module_function(mLapack, "slasd1", rb_slasd1, -1);
}

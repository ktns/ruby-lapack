#include "rb_lapack.h"

extern VOID dlasd1_(integer *nl, integer *nr, integer *sqre, doublereal *d, doublereal *alpha, doublereal *beta, doublereal *u, integer *ldu, doublereal *vt, integer *ldvt, integer *idxq, integer *iwork, doublereal *work, integer *info);

static VALUE
rb_dlasd1(int argc, VALUE *argv, VALUE self){
  VALUE rb_nl;
  integer nl; 
  VALUE rb_nr;
  integer nr; 
  VALUE rb_sqre;
  integer sqre; 
  VALUE rb_d;
  doublereal *d; 
  VALUE rb_alpha;
  doublereal alpha; 
  VALUE rb_beta;
  doublereal beta; 
  VALUE rb_u;
  doublereal *u; 
  VALUE rb_vt;
  doublereal *vt; 
  VALUE rb_idxq;
  integer *idxq; 
  VALUE rb_info;
  integer info; 
  VALUE rb_d_out__;
  doublereal *d_out__;
  VALUE rb_u_out__;
  doublereal *u_out__;
  VALUE rb_vt_out__;
  doublereal *vt_out__;
  integer *iwork;
  doublereal *work;

  integer n;
  integer ldu;
  integer ldvt;
  integer m;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  idxq, info, d, alpha, beta, u, vt = NumRu::Lapack.dlasd1( nl, nr, sqre, d, alpha, beta, u, vt)\n    or\n  NumRu::Lapack.dlasd1  # print help\n\n\nFORTRAN MANUAL\n\n");
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

  alpha = NUM2DBL(rb_alpha);
  if (!NA_IsNArray(rb_u))
    rb_raise(rb_eArgError, "u (7th argument) must be NArray");
  if (NA_RANK(rb_u) != 2)
    rb_raise(rb_eArgError, "rank of u (7th argument) must be %d", 2);
  n = NA_SHAPE1(rb_u);
  if (n != (nl+nr+1))
    rb_raise(rb_eRuntimeError, "shape 1 of u must be %d", nl+nr+1);
  ldu = NA_SHAPE0(rb_u);
  if (NA_TYPE(rb_u) != NA_DFLOAT)
    rb_u = na_change_type(rb_u, NA_DFLOAT);
  u = NA_PTR_TYPE(rb_u, doublereal*);
  if (!NA_IsNArray(rb_d))
    rb_raise(rb_eArgError, "d (4th argument) must be NArray");
  if (NA_RANK(rb_d) != 1)
    rb_raise(rb_eArgError, "rank of d (4th argument) must be %d", 1);
  if (NA_SHAPE0(rb_d) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of d must be the same as shape 1 of u");
  if (NA_TYPE(rb_d) != NA_DFLOAT)
    rb_d = na_change_type(rb_d, NA_DFLOAT);
  d = NA_PTR_TYPE(rb_d, doublereal*);
  nr = NUM2INT(rb_nr);
  beta = NUM2DBL(rb_beta);
  nl = NUM2INT(rb_nl);
  if (!NA_IsNArray(rb_vt))
    rb_raise(rb_eArgError, "vt (8th argument) must be NArray");
  if (NA_RANK(rb_vt) != 2)
    rb_raise(rb_eArgError, "rank of vt (8th argument) must be %d", 2);
  m = NA_SHAPE1(rb_vt);
  if (m != (n + sqre))
    rb_raise(rb_eRuntimeError, "shape 1 of vt must be %d", n + sqre);
  ldvt = NA_SHAPE0(rb_vt);
  if (NA_TYPE(rb_vt) != NA_DFLOAT)
    rb_vt = na_change_type(rb_vt, NA_DFLOAT);
  vt = NA_PTR_TYPE(rb_vt, doublereal*);
  sqre = NUM2INT(rb_sqre);
  n = nl+nr+1;
  m = n + sqre;
  {
    int shape[1];
    shape[0] = n;
    rb_idxq = na_make_object(NA_LINT, 1, shape, cNArray);
  }
  idxq = NA_PTR_TYPE(rb_idxq, integer*);
  {
    int shape[1];
    shape[0] = n;
    rb_d_out__ = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  d_out__ = NA_PTR_TYPE(rb_d_out__, doublereal*);
  MEMCPY(d_out__, d, doublereal, NA_TOTAL(rb_d));
  rb_d = rb_d_out__;
  d = d_out__;
  {
    int shape[2];
    shape[0] = ldu;
    shape[1] = n;
    rb_u_out__ = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  u_out__ = NA_PTR_TYPE(rb_u_out__, doublereal*);
  MEMCPY(u_out__, u, doublereal, NA_TOTAL(rb_u));
  rb_u = rb_u_out__;
  u = u_out__;
  {
    int shape[2];
    shape[0] = ldvt;
    shape[1] = m;
    rb_vt_out__ = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  vt_out__ = NA_PTR_TYPE(rb_vt_out__, doublereal*);
  MEMCPY(vt_out__, vt, doublereal, NA_TOTAL(rb_vt));
  rb_vt = rb_vt_out__;
  vt = vt_out__;
  iwork = ALLOC_N(integer, (4 * n));
  work = ALLOC_N(doublereal, (3*pow(m,2) + 2*m));

  dlasd1_(&nl, &nr, &sqre, d, &alpha, &beta, u, &ldu, vt, &ldvt, idxq, iwork, work, &info);

  free(iwork);
  free(work);
  rb_info = INT2NUM(info);
  rb_alpha = rb_float_new((double)alpha);
  rb_beta = rb_float_new((double)beta);
  return rb_ary_new3(7, rb_idxq, rb_info, rb_d, rb_alpha, rb_beta, rb_u, rb_vt);
}

void
init_lapack_dlasd1(VALUE mLapack){
  rb_define_module_function(mLapack, "dlasd1", rb_dlasd1, -1);
}

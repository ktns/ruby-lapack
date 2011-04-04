#include "rb_lapack.h"

extern VOID slasdq_(char *uplo, integer *sqre, integer *n, integer *ncvt, integer *nru, integer *ncc, real *d, real *e, real *vt, integer *ldvt, real *u, integer *ldu, real *c, integer *ldc, real *work, integer *info);

static VALUE
rb_slasdq(int argc, VALUE *argv, VALUE self){
  VALUE rb_uplo;
  char uplo; 
  VALUE rb_sqre;
  integer sqre; 
  VALUE rb_nru;
  integer nru; 
  VALUE rb_d;
  real *d; 
  VALUE rb_e;
  real *e; 
  VALUE rb_vt;
  real *vt; 
  VALUE rb_u;
  real *u; 
  VALUE rb_c;
  real *c; 
  VALUE rb_info;
  integer info; 
  VALUE rb_d_out__;
  real *d_out__;
  VALUE rb_e_out__;
  real *e_out__;
  VALUE rb_vt_out__;
  real *vt_out__;
  VALUE rb_u_out__;
  real *u_out__;
  VALUE rb_c_out__;
  real *c_out__;
  real *work;

  integer n;
  integer ldvt;
  integer ncvt;
  integer ldu;
  integer ldc;
  integer ncc;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  info, d, e, vt, u, c = NumRu::Lapack.slasdq( uplo, sqre, nru, d, e, vt, u, c)\n    or\n  NumRu::Lapack.slasdq  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 8)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 8)", argc);
  rb_uplo = argv[0];
  rb_sqre = argv[1];
  rb_nru = argv[2];
  rb_d = argv[3];
  rb_e = argv[4];
  rb_vt = argv[5];
  rb_u = argv[6];
  rb_c = argv[7];

  if (!NA_IsNArray(rb_c))
    rb_raise(rb_eArgError, "c (8th argument) must be NArray");
  if (NA_RANK(rb_c) != 2)
    rb_raise(rb_eArgError, "rank of c (8th argument) must be %d", 2);
  ncc = NA_SHAPE1(rb_c);
  ldc = NA_SHAPE0(rb_c);
  if (NA_TYPE(rb_c) != NA_SFLOAT)
    rb_c = na_change_type(rb_c, NA_SFLOAT);
  c = NA_PTR_TYPE(rb_c, real*);
  if (!NA_IsNArray(rb_d))
    rb_raise(rb_eArgError, "d (4th argument) must be NArray");
  if (NA_RANK(rb_d) != 1)
    rb_raise(rb_eArgError, "rank of d (4th argument) must be %d", 1);
  n = NA_SHAPE0(rb_d);
  if (NA_TYPE(rb_d) != NA_SFLOAT)
    rb_d = na_change_type(rb_d, NA_SFLOAT);
  d = NA_PTR_TYPE(rb_d, real*);
  uplo = StringValueCStr(rb_uplo)[0];
  nru = NUM2INT(rb_nru);
  if (!NA_IsNArray(rb_u))
    rb_raise(rb_eArgError, "u (7th argument) must be NArray");
  if (NA_RANK(rb_u) != 2)
    rb_raise(rb_eArgError, "rank of u (7th argument) must be %d", 2);
  if (NA_SHAPE1(rb_u) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of u must be the same as shape 0 of d");
  ldu = NA_SHAPE0(rb_u);
  if (NA_TYPE(rb_u) != NA_SFLOAT)
    rb_u = na_change_type(rb_u, NA_SFLOAT);
  u = NA_PTR_TYPE(rb_u, real*);
  if (!NA_IsNArray(rb_vt))
    rb_raise(rb_eArgError, "vt (6th argument) must be NArray");
  if (NA_RANK(rb_vt) != 2)
    rb_raise(rb_eArgError, "rank of vt (6th argument) must be %d", 2);
  ncvt = NA_SHAPE1(rb_vt);
  ldvt = NA_SHAPE0(rb_vt);
  if (NA_TYPE(rb_vt) != NA_SFLOAT)
    rb_vt = na_change_type(rb_vt, NA_SFLOAT);
  vt = NA_PTR_TYPE(rb_vt, real*);
  sqre = NUM2INT(rb_sqre);
  if (!NA_IsNArray(rb_e))
    rb_raise(rb_eArgError, "e (5th argument) must be NArray");
  if (NA_RANK(rb_e) != 1)
    rb_raise(rb_eArgError, "rank of e (5th argument) must be %d", 1);
  if (NA_SHAPE0(rb_e) != (sqre==0 ? n-1 : sqre==1 ? n : 0))
    rb_raise(rb_eRuntimeError, "shape 0 of e must be %d", sqre==0 ? n-1 : sqre==1 ? n : 0);
  if (NA_TYPE(rb_e) != NA_SFLOAT)
    rb_e = na_change_type(rb_e, NA_SFLOAT);
  e = NA_PTR_TYPE(rb_e, real*);
  {
    int shape[1];
    shape[0] = n;
    rb_d_out__ = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  d_out__ = NA_PTR_TYPE(rb_d_out__, real*);
  MEMCPY(d_out__, d, real, NA_TOTAL(rb_d));
  rb_d = rb_d_out__;
  d = d_out__;
  {
    int shape[1];
    shape[0] = sqre==0 ? n-1 : sqre==1 ? n : 0;
    rb_e_out__ = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  e_out__ = NA_PTR_TYPE(rb_e_out__, real*);
  MEMCPY(e_out__, e, real, NA_TOTAL(rb_e));
  rb_e = rb_e_out__;
  e = e_out__;
  {
    int shape[2];
    shape[0] = ldvt;
    shape[1] = ncvt;
    rb_vt_out__ = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  vt_out__ = NA_PTR_TYPE(rb_vt_out__, real*);
  MEMCPY(vt_out__, vt, real, NA_TOTAL(rb_vt));
  rb_vt = rb_vt_out__;
  vt = vt_out__;
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
    shape[0] = ldc;
    shape[1] = ncc;
    rb_c_out__ = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  c_out__ = NA_PTR_TYPE(rb_c_out__, real*);
  MEMCPY(c_out__, c, real, NA_TOTAL(rb_c));
  rb_c = rb_c_out__;
  c = c_out__;
  work = ALLOC_N(real, (4*n));

  slasdq_(&uplo, &sqre, &n, &ncvt, &nru, &ncc, d, e, vt, &ldvt, u, &ldu, c, &ldc, work, &info);

  free(work);
  rb_info = INT2NUM(info);
  return rb_ary_new3(6, rb_info, rb_d, rb_e, rb_vt, rb_u, rb_c);
}

void
init_lapack_slasdq(VALUE mLapack){
  rb_define_module_function(mLapack, "slasdq", rb_slasdq, -1);
}

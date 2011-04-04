#include "rb_lapack.h"

extern VOID dbdsqr_(char *uplo, integer *n, integer *ncvt, integer *nru, integer *ncc, doublereal *d, doublereal *e, doublereal *vt, integer *ldvt, doublereal *u, integer *ldu, doublereal *c, integer *ldc, doublereal *work, integer *info);

static VALUE
rb_dbdsqr(int argc, VALUE *argv, VALUE self){
  VALUE rb_uplo;
  char uplo; 
  VALUE rb_nru;
  integer nru; 
  VALUE rb_d;
  doublereal *d; 
  VALUE rb_e;
  doublereal *e; 
  VALUE rb_vt;
  doublereal *vt; 
  VALUE rb_u;
  doublereal *u; 
  VALUE rb_c;
  doublereal *c; 
  VALUE rb_info;
  integer info; 
  VALUE rb_d_out__;
  doublereal *d_out__;
  VALUE rb_e_out__;
  doublereal *e_out__;
  VALUE rb_vt_out__;
  doublereal *vt_out__;
  VALUE rb_u_out__;
  doublereal *u_out__;
  VALUE rb_c_out__;
  doublereal *c_out__;
  doublereal *work;

  integer n;
  integer ldvt;
  integer ncvt;
  integer ldu;
  integer ldc;
  integer ncc;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  info, d, e, vt, u, c = NumRu::Lapack.dbdsqr( uplo, nru, d, e, vt, u, c)\n    or\n  NumRu::Lapack.dbdsqr  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 7)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 7)", argc);
  rb_uplo = argv[0];
  rb_nru = argv[1];
  rb_d = argv[2];
  rb_e = argv[3];
  rb_vt = argv[4];
  rb_u = argv[5];
  rb_c = argv[6];

  if (!NA_IsNArray(rb_c))
    rb_raise(rb_eArgError, "c (7th argument) must be NArray");
  if (NA_RANK(rb_c) != 2)
    rb_raise(rb_eArgError, "rank of c (7th argument) must be %d", 2);
  ncc = NA_SHAPE1(rb_c);
  ldc = NA_SHAPE0(rb_c);
  if (NA_TYPE(rb_c) != NA_DFLOAT)
    rb_c = na_change_type(rb_c, NA_DFLOAT);
  c = NA_PTR_TYPE(rb_c, doublereal*);
  if (!NA_IsNArray(rb_d))
    rb_raise(rb_eArgError, "d (3th argument) must be NArray");
  if (NA_RANK(rb_d) != 1)
    rb_raise(rb_eArgError, "rank of d (3th argument) must be %d", 1);
  n = NA_SHAPE0(rb_d);
  if (NA_TYPE(rb_d) != NA_DFLOAT)
    rb_d = na_change_type(rb_d, NA_DFLOAT);
  d = NA_PTR_TYPE(rb_d, doublereal*);
  uplo = StringValueCStr(rb_uplo)[0];
  if (!NA_IsNArray(rb_u))
    rb_raise(rb_eArgError, "u (6th argument) must be NArray");
  if (NA_RANK(rb_u) != 2)
    rb_raise(rb_eArgError, "rank of u (6th argument) must be %d", 2);
  if (NA_SHAPE1(rb_u) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of u must be the same as shape 0 of d");
  ldu = NA_SHAPE0(rb_u);
  if (NA_TYPE(rb_u) != NA_DFLOAT)
    rb_u = na_change_type(rb_u, NA_DFLOAT);
  u = NA_PTR_TYPE(rb_u, doublereal*);
  if (!NA_IsNArray(rb_vt))
    rb_raise(rb_eArgError, "vt (5th argument) must be NArray");
  if (NA_RANK(rb_vt) != 2)
    rb_raise(rb_eArgError, "rank of vt (5th argument) must be %d", 2);
  ncvt = NA_SHAPE1(rb_vt);
  ldvt = NA_SHAPE0(rb_vt);
  if (NA_TYPE(rb_vt) != NA_DFLOAT)
    rb_vt = na_change_type(rb_vt, NA_DFLOAT);
  vt = NA_PTR_TYPE(rb_vt, doublereal*);
  nru = NUM2INT(rb_nru);
  if (!NA_IsNArray(rb_e))
    rb_raise(rb_eArgError, "e (4th argument) must be NArray");
  if (NA_RANK(rb_e) != 1)
    rb_raise(rb_eArgError, "rank of e (4th argument) must be %d", 1);
  if (NA_SHAPE0(rb_e) != (n-1))
    rb_raise(rb_eRuntimeError, "shape 0 of e must be %d", n-1);
  if (NA_TYPE(rb_e) != NA_DFLOAT)
    rb_e = na_change_type(rb_e, NA_DFLOAT);
  e = NA_PTR_TYPE(rb_e, doublereal*);
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
    int shape[1];
    shape[0] = n-1;
    rb_e_out__ = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  e_out__ = NA_PTR_TYPE(rb_e_out__, doublereal*);
  MEMCPY(e_out__, e, doublereal, NA_TOTAL(rb_e));
  rb_e = rb_e_out__;
  e = e_out__;
  {
    int shape[2];
    shape[0] = ldvt;
    shape[1] = ncvt;
    rb_vt_out__ = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  vt_out__ = NA_PTR_TYPE(rb_vt_out__, doublereal*);
  MEMCPY(vt_out__, vt, doublereal, NA_TOTAL(rb_vt));
  rb_vt = rb_vt_out__;
  vt = vt_out__;
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
    shape[0] = ldc;
    shape[1] = ncc;
    rb_c_out__ = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  c_out__ = NA_PTR_TYPE(rb_c_out__, doublereal*);
  MEMCPY(c_out__, c, doublereal, NA_TOTAL(rb_c));
  rb_c = rb_c_out__;
  c = c_out__;
  work = ALLOC_N(doublereal, (4*n));

  dbdsqr_(&uplo, &n, &ncvt, &nru, &ncc, d, e, vt, &ldvt, u, &ldu, c, &ldc, work, &info);

  free(work);
  rb_info = INT2NUM(info);
  return rb_ary_new3(6, rb_info, rb_d, rb_e, rb_vt, rb_u, rb_c);
}

void
init_lapack_dbdsqr(VALUE mLapack){
  rb_define_module_function(mLapack, "dbdsqr", rb_dbdsqr, -1);
}

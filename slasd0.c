#include "rb_lapack.h"

extern VOID slasd0_(integer *n, integer *sqre, real *d, real *e, real *u, integer *ldu, real *vt, integer *ldvt, integer *smlsiz, integer *iwork, real *work, integer *info);

static VALUE
rb_slasd0(int argc, VALUE *argv, VALUE self){
  VALUE rb_sqre;
  integer sqre; 
  VALUE rb_d;
  real *d; 
  VALUE rb_e;
  real *e; 
  VALUE rb_smlsiz;
  integer smlsiz; 
  VALUE rb_u;
  real *u; 
  VALUE rb_vt;
  real *vt; 
  VALUE rb_info;
  integer info; 
  VALUE rb_d_out__;
  real *d_out__;
  integer *iwork;
  real *work;

  integer n;
  integer ldu;
  integer ldvt;
  integer m;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  u, vt, info, d = NumRu::Lapack.slasd0( sqre, d, e, smlsiz)\n    or\n  NumRu::Lapack.slasd0  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 4)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 4)", argc);
  rb_sqre = argv[0];
  rb_d = argv[1];
  rb_e = argv[2];
  rb_smlsiz = argv[3];

  smlsiz = NUM2INT(rb_smlsiz);
  sqre = NUM2INT(rb_sqre);
  if (!NA_IsNArray(rb_d))
    rb_raise(rb_eArgError, "d (2th argument) must be NArray");
  if (NA_RANK(rb_d) != 1)
    rb_raise(rb_eArgError, "rank of d (2th argument) must be %d", 1);
  n = NA_SHAPE0(rb_d);
  if (NA_TYPE(rb_d) != NA_SFLOAT)
    rb_d = na_change_type(rb_d, NA_SFLOAT);
  d = NA_PTR_TYPE(rb_d, real*);
  ldu = n;
  m = sqre == 0 ? n : sqre == 1 ? n+1 : 0;
  ldvt = m;
  if (!NA_IsNArray(rb_e))
    rb_raise(rb_eArgError, "e (3th argument) must be NArray");
  if (NA_RANK(rb_e) != 1)
    rb_raise(rb_eArgError, "rank of e (3th argument) must be %d", 1);
  if (NA_SHAPE0(rb_e) != (m-1))
    rb_raise(rb_eRuntimeError, "shape 0 of e must be %d", m-1);
  if (NA_TYPE(rb_e) != NA_SFLOAT)
    rb_e = na_change_type(rb_e, NA_SFLOAT);
  e = NA_PTR_TYPE(rb_e, real*);
  {
    int shape[2];
    shape[0] = ldu;
    shape[1] = n;
    rb_u = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  u = NA_PTR_TYPE(rb_u, real*);
  {
    int shape[2];
    shape[0] = ldvt;
    shape[1] = m;
    rb_vt = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  vt = NA_PTR_TYPE(rb_vt, real*);
  {
    int shape[1];
    shape[0] = n;
    rb_d_out__ = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  d_out__ = NA_PTR_TYPE(rb_d_out__, real*);
  MEMCPY(d_out__, d, real, NA_TOTAL(rb_d));
  rb_d = rb_d_out__;
  d = d_out__;
  iwork = ALLOC_N(integer, (8*n));
  work = ALLOC_N(real, (3*pow(m,2)+2*m));

  slasd0_(&n, &sqre, d, e, u, &ldu, vt, &ldvt, &smlsiz, iwork, work, &info);

  free(iwork);
  free(work);
  rb_info = INT2NUM(info);
  return rb_ary_new3(4, rb_u, rb_vt, rb_info, rb_d);
}

void
init_lapack_slasd0(VALUE mLapack){
  rb_define_module_function(mLapack, "slasd0", rb_slasd0, -1);
}

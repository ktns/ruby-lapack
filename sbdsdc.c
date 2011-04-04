#include "rb_lapack.h"

extern VOID sbdsdc_(char *uplo, char *compq, integer *n, real *d, real *e, real *u, integer *ldu, real *vt, integer *ldvt, real *q, integer *iq, real *work, integer *iwork, integer *info);

static VALUE
rb_sbdsdc(int argc, VALUE *argv, VALUE self){
  VALUE rb_uplo;
  char uplo; 
  VALUE rb_compq;
  char compq; 
  VALUE rb_d;
  real *d; 
  VALUE rb_e;
  real *e; 
  VALUE rb_ldu;
  integer ldu; 
  VALUE rb_ldvt;
  integer ldvt; 
  VALUE rb_u;
  real *u; 
  VALUE rb_vt;
  real *vt; 
  VALUE rb_q;
  real *q; 
  VALUE rb_iq;
  integer *iq; 
  VALUE rb_info;
  integer info; 
  VALUE rb_d_out__;
  real *d_out__;
  VALUE rb_e_out__;
  real *e_out__;
  real *work;
  integer *iwork;

  integer n;
  integer c__9;
  integer c__0;
  integer ldq;
  integer ldiq;
  integer lwork;
  integer smlsiz;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  u, vt, q, iq, info, d, e = NumRu::Lapack.sbdsdc( uplo, compq, d, e, ldu, ldvt)\n    or\n  NumRu::Lapack.sbdsdc  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 6)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 6)", argc);
  rb_uplo = argv[0];
  rb_compq = argv[1];
  rb_d = argv[2];
  rb_e = argv[3];
  rb_ldu = argv[4];
  rb_ldvt = argv[5];

  c__9 = 9;
  c__0 = 0;
  uplo = StringValueCStr(rb_uplo)[0];
  if (!NA_IsNArray(rb_d))
    rb_raise(rb_eArgError, "d (3th argument) must be NArray");
  if (NA_RANK(rb_d) != 1)
    rb_raise(rb_eArgError, "rank of d (3th argument) must be %d", 1);
  n = NA_SHAPE0(rb_d);
  if (NA_TYPE(rb_d) != NA_SFLOAT)
    rb_d = na_change_type(rb_d, NA_SFLOAT);
  d = NA_PTR_TYPE(rb_d, real*);
  compq = StringValueCStr(rb_compq)[0];
  ldvt = lsame_(&compq,"I") ? MAX(1,n) : 0;
  ldiq = lsame_(&compq,"P") ? n*(3+3*(int)(log(((double)n)/(smlsiz+1))/log(2.0))) : 0;
  lwork = lsame_(&compq,"N") ? 4*n : lsame_(&compq,"P") ? 6*n : lsame_(&compq,"I") ? 3*n*n+4*n : 0;
  if (!NA_IsNArray(rb_e))
    rb_raise(rb_eArgError, "e (4th argument) must be NArray");
  if (NA_RANK(rb_e) != 1)
    rb_raise(rb_eArgError, "rank of e (4th argument) must be %d", 1);
  if (NA_SHAPE0(rb_e) != (n-1))
    rb_raise(rb_eRuntimeError, "shape 0 of e must be %d", n-1);
  if (NA_TYPE(rb_e) != NA_SFLOAT)
    rb_e = na_change_type(rb_e, NA_SFLOAT);
  e = NA_PTR_TYPE(rb_e, real*);
  ldu = lsame_(&compq,"I") ? MAX(1,n) : 0;
  smlsiz = ilaenv_(&c__9, "SBDSDC", " ", &c__0, &c__0, &c__0, &c__0);
  ldq = lsame_(&compq,"P") ? n*(11+2*smlsiz+8*(int)(log(((double)n)/(smlsiz+1))/log(2.0))) : 0;
  {
    int shape[2];
    shape[0] = lsame_(&compq,"I") ? ldu : 0;
    shape[1] = lsame_(&compq,"I") ? n : 0;
    rb_u = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  u = NA_PTR_TYPE(rb_u, real*);
  {
    int shape[2];
    shape[0] = lsame_(&compq,"I") ? ldvt : 0;
    shape[1] = lsame_(&compq,"I") ? n : 0;
    rb_vt = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  vt = NA_PTR_TYPE(rb_vt, real*);
  {
    int shape[1];
    shape[0] = lsame_(&compq,"I") ? ldq : 0;
    rb_q = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  q = NA_PTR_TYPE(rb_q, real*);
  {
    int shape[1];
    shape[0] = lsame_(&compq,"I") ? ldiq : 0;
    rb_iq = na_make_object(NA_LINT, 1, shape, cNArray);
  }
  iq = NA_PTR_TYPE(rb_iq, integer*);
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
    shape[0] = n-1;
    rb_e_out__ = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  e_out__ = NA_PTR_TYPE(rb_e_out__, real*);
  MEMCPY(e_out__, e, real, NA_TOTAL(rb_e));
  rb_e = rb_e_out__;
  e = e_out__;
  work = ALLOC_N(real, (MAX(1,lwork)));
  iwork = ALLOC_N(integer, (8*n));

  sbdsdc_(&uplo, &compq, &n, d, e, u, &ldu, vt, &ldvt, q, iq, work, iwork, &info);

  free(work);
  free(iwork);
  rb_info = INT2NUM(info);
  return rb_ary_new3(7, rb_u, rb_vt, rb_q, rb_iq, rb_info, rb_d, rb_e);
}

void
init_lapack_sbdsdc(VALUE mLapack){
  rb_define_module_function(mLapack, "sbdsdc", rb_sbdsdc, -1);
}

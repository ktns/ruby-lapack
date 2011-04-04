#include "rb_lapack.h"

extern VOID zuncsd_(char *jobu1, char *jobu2, char *jobv1t, char *jobv2t, char *trans, char *signs, integer *m, integer *p, integer *q, doublecomplex *x11, integer *ldx11, doublecomplex *x12, integer *ldx12, doublecomplex *x21, integer *ldx21, doublecomplex *x22, integer *ldx22, doublereal *theta, doublecomplex *u1, integer *ldu1, doublecomplex *u2, integer *ldu2, doublecomplex *v1t, integer *ldv1t, doublecomplex *v2t, integer *ldv2t, doublecomplex *work, integer *lwork, doublereal *rwork, integer *lrwork, integer *iwork, integer *info);

static VALUE
rb_zuncsd(int argc, VALUE *argv, VALUE self){
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
  VALUE rb_signs;
  char signs; 
  VALUE rb_m;
  integer m; 
  VALUE rb_x11;
  doublecomplex *x11; 
  VALUE rb_ldx11;
  integer ldx11; 
  VALUE rb_x12;
  doublecomplex *x12; 
  VALUE rb_ldx12;
  integer ldx12; 
  VALUE rb_x21;
  doublecomplex *x21; 
  VALUE rb_ldx21;
  integer ldx21; 
  VALUE rb_x22;
  doublecomplex *x22; 
  VALUE rb_ldx22;
  integer ldx22; 
  VALUE rb_ldu1;
  integer ldu1; 
  VALUE rb_ldu2;
  integer ldu2; 
  VALUE rb_ldv1t;
  integer ldv1t; 
  VALUE rb_ldv2t;
  integer ldv2t; 
  VALUE rb_lwork;
  integer lwork; 
  VALUE rb_lrwork;
  integer lrwork; 
  VALUE rb_theta;
  doublereal *theta; 
  VALUE rb_u1;
  doublecomplex *u1; 
  VALUE rb_u2;
  doublecomplex *u2; 
  VALUE rb_v1t;
  doublecomplex *v1t; 
  VALUE rb_v2t;
  doublecomplex *v2t; 
  VALUE rb_info;
  integer info; 
  doublecomplex *work;
  doublereal *rwork;
  integer *iwork;

  integer p;
  integer q;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  theta, u1, u2, v1t, v2t, info = NumRu::Lapack.zuncsd( jobu1, jobu2, jobv1t, jobv2t, trans, signs, m, x11, ldx11, x12, ldx12, x21, ldx21, x22, ldx22, ldu1, ldu2, ldv1t, ldv2t, lwork, lrwork)\n    or\n  NumRu::Lapack.zuncsd  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 21)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 21)", argc);
  rb_jobu1 = argv[0];
  rb_jobu2 = argv[1];
  rb_jobv1t = argv[2];
  rb_jobv2t = argv[3];
  rb_trans = argv[4];
  rb_signs = argv[5];
  rb_m = argv[6];
  rb_x11 = argv[7];
  rb_ldx11 = argv[8];
  rb_x12 = argv[9];
  rb_ldx12 = argv[10];
  rb_x21 = argv[11];
  rb_ldx21 = argv[12];
  rb_x22 = argv[13];
  rb_ldx22 = argv[14];
  rb_ldu1 = argv[15];
  rb_ldu2 = argv[16];
  rb_ldv1t = argv[17];
  rb_ldv2t = argv[18];
  rb_lwork = argv[19];
  rb_lrwork = argv[20];

  trans = StringValueCStr(rb_trans)[0];
  jobv1t = StringValueCStr(rb_jobv1t)[0];
  jobv2t = StringValueCStr(rb_jobv2t)[0];
  lwork = NUM2INT(rb_lwork);
  if (!NA_IsNArray(rb_x21))
    rb_raise(rb_eArgError, "x21 (12th argument) must be NArray");
  if (NA_RANK(rb_x21) != 2)
    rb_raise(rb_eArgError, "rank of x21 (12th argument) must be %d", 2);
  q = NA_SHAPE1(rb_x21);
  p = NA_SHAPE0(rb_x21);
  if (NA_TYPE(rb_x21) != NA_DCOMPLEX)
    rb_x21 = na_change_type(rb_x21, NA_DCOMPLEX);
  x21 = NA_PTR_TYPE(rb_x21, doublecomplex*);
  signs = StringValueCStr(rb_signs)[0];
  jobu1 = StringValueCStr(rb_jobu1)[0];
  lrwork = NUM2INT(rb_lrwork);
  jobu2 = StringValueCStr(rb_jobu2)[0];
  if (!NA_IsNArray(rb_x11))
    rb_raise(rb_eArgError, "x11 (8th argument) must be NArray");
  if (NA_RANK(rb_x11) != 2)
    rb_raise(rb_eArgError, "rank of x11 (8th argument) must be %d", 2);
  if (NA_SHAPE1(rb_x11) != q)
    rb_raise(rb_eRuntimeError, "shape 1 of x11 must be the same as shape 1 of x21");
  if (NA_SHAPE0(rb_x11) != p)
    rb_raise(rb_eRuntimeError, "shape 0 of x11 must be the same as shape 0 of x21");
  if (NA_TYPE(rb_x11) != NA_DCOMPLEX)
    rb_x11 = na_change_type(rb_x11, NA_DCOMPLEX);
  x11 = NA_PTR_TYPE(rb_x11, doublecomplex*);
  m = NUM2INT(rb_m);
  if (!NA_IsNArray(rb_x22))
    rb_raise(rb_eArgError, "x22 (14th argument) must be NArray");
  if (NA_RANK(rb_x22) != 2)
    rb_raise(rb_eArgError, "rank of x22 (14th argument) must be %d", 2);
  if (NA_SHAPE1(rb_x22) != (m-q))
    rb_raise(rb_eRuntimeError, "shape 1 of x22 must be %d", m-q);
  if (NA_SHAPE0(rb_x22) != p)
    rb_raise(rb_eRuntimeError, "shape 0 of x22 must be the same as shape 0 of x21");
  if (NA_TYPE(rb_x22) != NA_DCOMPLEX)
    rb_x22 = na_change_type(rb_x22, NA_DCOMPLEX);
  x22 = NA_PTR_TYPE(rb_x22, doublecomplex*);
  ldu1 = lsame_(&jobu1,"Y") ? MAX(1,p) : 0;
  ldu2 = lsame_(&jobu2,"Y") ? MAX(1,m-p) : 0;
  if (!NA_IsNArray(rb_x12))
    rb_raise(rb_eArgError, "x12 (10th argument) must be NArray");
  if (NA_RANK(rb_x12) != 2)
    rb_raise(rb_eArgError, "rank of x12 (10th argument) must be %d", 2);
  if (NA_SHAPE1(rb_x12) != (m-q))
    rb_raise(rb_eRuntimeError, "shape 1 of x12 must be %d", m-q);
  if (NA_SHAPE0(rb_x12) != p)
    rb_raise(rb_eRuntimeError, "shape 0 of x12 must be the same as shape 0 of x21");
  if (NA_TYPE(rb_x12) != NA_DCOMPLEX)
    rb_x12 = na_change_type(rb_x12, NA_DCOMPLEX);
  x12 = NA_PTR_TYPE(rb_x12, doublecomplex*);
  ldx12 = p;
  ldv1t = lsame_(&jobv1t,"Y") ? MAX(1,q) : 0;
  ldx11 = p;
  ldx22 = p;
  ldx21 = p;
  ldv2t = lsame_(&jobv2t,"Y") ? MAX(1,m-q) : 0;
  {
    int shape[1];
    shape[0] = MIN(MIN(MIN(p,m-p),q),m-q);
    rb_theta = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  theta = NA_PTR_TYPE(rb_theta, doublereal*);
  {
    int shape[1];
    shape[0] = p;
    rb_u1 = na_make_object(NA_DCOMPLEX, 1, shape, cNArray);
  }
  u1 = NA_PTR_TYPE(rb_u1, doublecomplex*);
  {
    int shape[1];
    shape[0] = m-p;
    rb_u2 = na_make_object(NA_DCOMPLEX, 1, shape, cNArray);
  }
  u2 = NA_PTR_TYPE(rb_u2, doublecomplex*);
  {
    int shape[1];
    shape[0] = q;
    rb_v1t = na_make_object(NA_DCOMPLEX, 1, shape, cNArray);
  }
  v1t = NA_PTR_TYPE(rb_v1t, doublecomplex*);
  {
    int shape[1];
    shape[0] = m-q;
    rb_v2t = na_make_object(NA_DCOMPLEX, 1, shape, cNArray);
  }
  v2t = NA_PTR_TYPE(rb_v2t, doublecomplex*);
  work = ALLOC_N(doublecomplex, (MAX(1,lwork)));
  rwork = ALLOC_N(doublereal, (MAX(1,lrwork)));
  iwork = ALLOC_N(integer, (m-q));

  zuncsd_(&jobu1, &jobu2, &jobv1t, &jobv2t, &trans, &signs, &m, &p, &q, x11, &ldx11, x12, &ldx12, x21, &ldx21, x22, &ldx22, theta, u1, &ldu1, u2, &ldu2, v1t, &ldv1t, v2t, &ldv2t, work, &lwork, rwork, &lrwork, iwork, &info);

  free(work);
  free(rwork);
  free(iwork);
  rb_info = INT2NUM(info);
  return rb_ary_new3(6, rb_theta, rb_u1, rb_u2, rb_v1t, rb_v2t, rb_info);
}

void
init_lapack_zuncsd(VALUE mLapack){
  rb_define_module_function(mLapack, "zuncsd", rb_zuncsd, -1);
}

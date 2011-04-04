#include "rb_lapack.h"

extern VOID zunmbr_(char *vect, char *side, char *trans, integer *m, integer *n, integer *k, doublecomplex *a, integer *lda, doublecomplex *tau, doublecomplex *c, integer *ldc, doublecomplex *work, integer *lwork, integer *info);

static VALUE
rb_zunmbr(int argc, VALUE *argv, VALUE self){
  VALUE rb_vect;
  char vect; 
  VALUE rb_side;
  char side; 
  VALUE rb_trans;
  char trans; 
  VALUE rb_m;
  integer m; 
  VALUE rb_k;
  integer k; 
  VALUE rb_a;
  doublecomplex *a; 
  VALUE rb_tau;
  doublecomplex *tau; 
  VALUE rb_c;
  doublecomplex *c; 
  VALUE rb_lwork;
  integer lwork; 
  VALUE rb_work;
  doublecomplex *work; 
  VALUE rb_info;
  integer info; 
  VALUE rb_c_out__;
  doublecomplex *c_out__;

  integer lda;
  integer ldc;
  integer n;
  integer nq;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  work, info, c = NumRu::Lapack.zunmbr( vect, side, trans, m, k, a, tau, c, lwork)\n    or\n  NumRu::Lapack.zunmbr  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 9)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 9)", argc);
  rb_vect = argv[0];
  rb_side = argv[1];
  rb_trans = argv[2];
  rb_m = argv[3];
  rb_k = argv[4];
  rb_a = argv[5];
  rb_tau = argv[6];
  rb_c = argv[7];
  rb_lwork = argv[8];

  k = NUM2INT(rb_k);
  lwork = NUM2INT(rb_lwork);
  side = StringValueCStr(rb_side)[0];
  m = NUM2INT(rb_m);
  vect = StringValueCStr(rb_vect)[0];
  if (!NA_IsNArray(rb_c))
    rb_raise(rb_eArgError, "c (8th argument) must be NArray");
  if (NA_RANK(rb_c) != 2)
    rb_raise(rb_eArgError, "rank of c (8th argument) must be %d", 2);
  n = NA_SHAPE1(rb_c);
  ldc = NA_SHAPE0(rb_c);
  if (NA_TYPE(rb_c) != NA_DCOMPLEX)
    rb_c = na_change_type(rb_c, NA_DCOMPLEX);
  c = NA_PTR_TYPE(rb_c, doublecomplex*);
  trans = StringValueCStr(rb_trans)[0];
  nq = lsame_(&side,"L") ? m : lsame_(&side,"R") ? n : 0;
  if (!NA_IsNArray(rb_tau))
    rb_raise(rb_eArgError, "tau (7th argument) must be NArray");
  if (NA_RANK(rb_tau) != 1)
    rb_raise(rb_eArgError, "rank of tau (7th argument) must be %d", 1);
  if (NA_SHAPE0(rb_tau) != (MIN(nq,k)))
    rb_raise(rb_eRuntimeError, "shape 0 of tau must be %d", MIN(nq,k));
  if (NA_TYPE(rb_tau) != NA_DCOMPLEX)
    rb_tau = na_change_type(rb_tau, NA_DCOMPLEX);
  tau = NA_PTR_TYPE(rb_tau, doublecomplex*);
  if (!NA_IsNArray(rb_a))
    rb_raise(rb_eArgError, "a (6th argument) must be NArray");
  if (NA_RANK(rb_a) != 2)
    rb_raise(rb_eArgError, "rank of a (6th argument) must be %d", 2);
  if (NA_SHAPE1(rb_a) != (MIN(nq,k)))
    rb_raise(rb_eRuntimeError, "shape 1 of a must be %d", MIN(nq,k));
  lda = NA_SHAPE0(rb_a);
  if (NA_TYPE(rb_a) != NA_DCOMPLEX)
    rb_a = na_change_type(rb_a, NA_DCOMPLEX);
  a = NA_PTR_TYPE(rb_a, doublecomplex*);
  {
    int shape[1];
    shape[0] = MAX(1,lwork);
    rb_work = na_make_object(NA_DCOMPLEX, 1, shape, cNArray);
  }
  work = NA_PTR_TYPE(rb_work, doublecomplex*);
  {
    int shape[2];
    shape[0] = ldc;
    shape[1] = n;
    rb_c_out__ = na_make_object(NA_DCOMPLEX, 2, shape, cNArray);
  }
  c_out__ = NA_PTR_TYPE(rb_c_out__, doublecomplex*);
  MEMCPY(c_out__, c, doublecomplex, NA_TOTAL(rb_c));
  rb_c = rb_c_out__;
  c = c_out__;

  zunmbr_(&vect, &side, &trans, &m, &n, &k, a, &lda, tau, c, &ldc, work, &lwork, &info);

  rb_info = INT2NUM(info);
  return rb_ary_new3(3, rb_work, rb_info, rb_c);
}

void
init_lapack_zunmbr(VALUE mLapack){
  rb_define_module_function(mLapack, "zunmbr", rb_zunmbr, -1);
}

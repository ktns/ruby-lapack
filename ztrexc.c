#include "rb_lapack.h"

extern VOID ztrexc_(char *compq, integer *n, doublecomplex *t, integer *ldt, doublecomplex *q, integer *ldq, integer *ifst, integer *ilst, integer *info);

static VALUE
rb_ztrexc(int argc, VALUE *argv, VALUE self){
  VALUE rb_compq;
  char compq; 
  VALUE rb_t;
  doublecomplex *t; 
  VALUE rb_q;
  doublecomplex *q; 
  VALUE rb_ifst;
  integer ifst; 
  VALUE rb_ilst;
  integer ilst; 
  VALUE rb_info;
  integer info; 
  VALUE rb_t_out__;
  doublecomplex *t_out__;
  VALUE rb_q_out__;
  doublecomplex *q_out__;

  integer ldt;
  integer n;
  integer ldq;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  info, t, q = NumRu::Lapack.ztrexc( compq, t, q, ifst, ilst)\n    or\n  NumRu::Lapack.ztrexc  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 5)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 5)", argc);
  rb_compq = argv[0];
  rb_t = argv[1];
  rb_q = argv[2];
  rb_ifst = argv[3];
  rb_ilst = argv[4];

  compq = StringValueCStr(rb_compq)[0];
  ilst = NUM2INT(rb_ilst);
  if (!NA_IsNArray(rb_q))
    rb_raise(rb_eArgError, "q (3th argument) must be NArray");
  if (NA_RANK(rb_q) != 2)
    rb_raise(rb_eArgError, "rank of q (3th argument) must be %d", 2);
  n = NA_SHAPE1(rb_q);
  ldq = NA_SHAPE0(rb_q);
  if (NA_TYPE(rb_q) != NA_DCOMPLEX)
    rb_q = na_change_type(rb_q, NA_DCOMPLEX);
  q = NA_PTR_TYPE(rb_q, doublecomplex*);
  if (!NA_IsNArray(rb_t))
    rb_raise(rb_eArgError, "t (2th argument) must be NArray");
  if (NA_RANK(rb_t) != 2)
    rb_raise(rb_eArgError, "rank of t (2th argument) must be %d", 2);
  if (NA_SHAPE1(rb_t) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of t must be the same as shape 1 of q");
  ldt = NA_SHAPE0(rb_t);
  if (NA_TYPE(rb_t) != NA_DCOMPLEX)
    rb_t = na_change_type(rb_t, NA_DCOMPLEX);
  t = NA_PTR_TYPE(rb_t, doublecomplex*);
  ifst = NUM2INT(rb_ifst);
  {
    int shape[2];
    shape[0] = ldt;
    shape[1] = n;
    rb_t_out__ = na_make_object(NA_DCOMPLEX, 2, shape, cNArray);
  }
  t_out__ = NA_PTR_TYPE(rb_t_out__, doublecomplex*);
  MEMCPY(t_out__, t, doublecomplex, NA_TOTAL(rb_t));
  rb_t = rb_t_out__;
  t = t_out__;
  {
    int shape[2];
    shape[0] = ldq;
    shape[1] = n;
    rb_q_out__ = na_make_object(NA_DCOMPLEX, 2, shape, cNArray);
  }
  q_out__ = NA_PTR_TYPE(rb_q_out__, doublecomplex*);
  MEMCPY(q_out__, q, doublecomplex, NA_TOTAL(rb_q));
  rb_q = rb_q_out__;
  q = q_out__;

  ztrexc_(&compq, &n, t, &ldt, q, &ldq, &ifst, &ilst, &info);

  rb_info = INT2NUM(info);
  return rb_ary_new3(3, rb_info, rb_t, rb_q);
}

void
init_lapack_ztrexc(VALUE mLapack){
  rb_define_module_function(mLapack, "ztrexc", rb_ztrexc, -1);
}

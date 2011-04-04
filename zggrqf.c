#include "rb_lapack.h"

extern VOID zggrqf_(integer *m, integer *p, integer *n, doublecomplex *a, integer *lda, doublecomplex *taua, doublecomplex *b, integer *ldb, doublecomplex *taub, doublecomplex *work, integer *lwork, integer *info);

static VALUE
rb_zggrqf(int argc, VALUE *argv, VALUE self){
  VALUE rb_m;
  integer m; 
  VALUE rb_p;
  integer p; 
  VALUE rb_a;
  doublecomplex *a; 
  VALUE rb_b;
  doublecomplex *b; 
  VALUE rb_lwork;
  integer lwork; 
  VALUE rb_taua;
  doublecomplex *taua; 
  VALUE rb_taub;
  doublecomplex *taub; 
  VALUE rb_work;
  doublecomplex *work; 
  VALUE rb_info;
  integer info; 
  VALUE rb_a_out__;
  doublecomplex *a_out__;
  VALUE rb_b_out__;
  doublecomplex *b_out__;

  integer lda;
  integer n;
  integer ldb;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  taua, taub, work, info, a, b = NumRu::Lapack.zggrqf( m, p, a, b, lwork)\n    or\n  NumRu::Lapack.zggrqf  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 5)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 5)", argc);
  rb_m = argv[0];
  rb_p = argv[1];
  rb_a = argv[2];
  rb_b = argv[3];
  rb_lwork = argv[4];

  if (!NA_IsNArray(rb_a))
    rb_raise(rb_eArgError, "a (3th argument) must be NArray");
  if (NA_RANK(rb_a) != 2)
    rb_raise(rb_eArgError, "rank of a (3th argument) must be %d", 2);
  n = NA_SHAPE1(rb_a);
  lda = NA_SHAPE0(rb_a);
  if (NA_TYPE(rb_a) != NA_DCOMPLEX)
    rb_a = na_change_type(rb_a, NA_DCOMPLEX);
  a = NA_PTR_TYPE(rb_a, doublecomplex*);
  if (!NA_IsNArray(rb_b))
    rb_raise(rb_eArgError, "b (4th argument) must be NArray");
  if (NA_RANK(rb_b) != 2)
    rb_raise(rb_eArgError, "rank of b (4th argument) must be %d", 2);
  if (NA_SHAPE1(rb_b) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of b must be the same as shape 1 of a");
  ldb = NA_SHAPE0(rb_b);
  if (NA_TYPE(rb_b) != NA_DCOMPLEX)
    rb_b = na_change_type(rb_b, NA_DCOMPLEX);
  b = NA_PTR_TYPE(rb_b, doublecomplex*);
  m = NUM2INT(rb_m);
  lwork = NUM2INT(rb_lwork);
  p = NUM2INT(rb_p);
  {
    int shape[1];
    shape[0] = MIN(m,n);
    rb_taua = na_make_object(NA_DCOMPLEX, 1, shape, cNArray);
  }
  taua = NA_PTR_TYPE(rb_taua, doublecomplex*);
  {
    int shape[1];
    shape[0] = MIN(p,n);
    rb_taub = na_make_object(NA_DCOMPLEX, 1, shape, cNArray);
  }
  taub = NA_PTR_TYPE(rb_taub, doublecomplex*);
  {
    int shape[1];
    shape[0] = MAX(1,lwork);
    rb_work = na_make_object(NA_DCOMPLEX, 1, shape, cNArray);
  }
  work = NA_PTR_TYPE(rb_work, doublecomplex*);
  {
    int shape[2];
    shape[0] = lda;
    shape[1] = n;
    rb_a_out__ = na_make_object(NA_DCOMPLEX, 2, shape, cNArray);
  }
  a_out__ = NA_PTR_TYPE(rb_a_out__, doublecomplex*);
  MEMCPY(a_out__, a, doublecomplex, NA_TOTAL(rb_a));
  rb_a = rb_a_out__;
  a = a_out__;
  {
    int shape[2];
    shape[0] = ldb;
    shape[1] = n;
    rb_b_out__ = na_make_object(NA_DCOMPLEX, 2, shape, cNArray);
  }
  b_out__ = NA_PTR_TYPE(rb_b_out__, doublecomplex*);
  MEMCPY(b_out__, b, doublecomplex, NA_TOTAL(rb_b));
  rb_b = rb_b_out__;
  b = b_out__;

  zggrqf_(&m, &p, &n, a, &lda, taua, b, &ldb, taub, work, &lwork, &info);

  rb_info = INT2NUM(info);
  return rb_ary_new3(6, rb_taua, rb_taub, rb_work, rb_info, rb_a, rb_b);
}

void
init_lapack_zggrqf(VALUE mLapack){
  rb_define_module_function(mLapack, "zggrqf", rb_zggrqf, -1);
}

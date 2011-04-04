#include "rb_lapack.h"

extern VOID cggqrf_(integer *n, integer *m, integer *p, complex *a, integer *lda, complex *taua, complex *b, integer *ldb, complex *taub, complex *work, integer *lwork, integer *info);

static VALUE
rb_cggqrf(int argc, VALUE *argv, VALUE self){
  VALUE rb_n;
  integer n; 
  VALUE rb_a;
  complex *a; 
  VALUE rb_b;
  complex *b; 
  VALUE rb_lwork;
  integer lwork; 
  VALUE rb_taua;
  complex *taua; 
  VALUE rb_taub;
  complex *taub; 
  VALUE rb_work;
  complex *work; 
  VALUE rb_info;
  integer info; 
  VALUE rb_a_out__;
  complex *a_out__;
  VALUE rb_b_out__;
  complex *b_out__;

  integer lda;
  integer m;
  integer ldb;
  integer p;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  taua, taub, work, info, a, b = NumRu::Lapack.cggqrf( n, a, b, lwork)\n    or\n  NumRu::Lapack.cggqrf  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 4)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 4)", argc);
  rb_n = argv[0];
  rb_a = argv[1];
  rb_b = argv[2];
  rb_lwork = argv[3];

  if (!NA_IsNArray(rb_a))
    rb_raise(rb_eArgError, "a (2th argument) must be NArray");
  if (NA_RANK(rb_a) != 2)
    rb_raise(rb_eArgError, "rank of a (2th argument) must be %d", 2);
  m = NA_SHAPE1(rb_a);
  lda = NA_SHAPE0(rb_a);
  if (NA_TYPE(rb_a) != NA_SCOMPLEX)
    rb_a = na_change_type(rb_a, NA_SCOMPLEX);
  a = NA_PTR_TYPE(rb_a, complex*);
  if (!NA_IsNArray(rb_b))
    rb_raise(rb_eArgError, "b (3th argument) must be NArray");
  if (NA_RANK(rb_b) != 2)
    rb_raise(rb_eArgError, "rank of b (3th argument) must be %d", 2);
  p = NA_SHAPE1(rb_b);
  ldb = NA_SHAPE0(rb_b);
  if (NA_TYPE(rb_b) != NA_SCOMPLEX)
    rb_b = na_change_type(rb_b, NA_SCOMPLEX);
  b = NA_PTR_TYPE(rb_b, complex*);
  n = NUM2INT(rb_n);
  lwork = NUM2INT(rb_lwork);
  {
    int shape[1];
    shape[0] = MIN(n,m);
    rb_taua = na_make_object(NA_SCOMPLEX, 1, shape, cNArray);
  }
  taua = NA_PTR_TYPE(rb_taua, complex*);
  {
    int shape[1];
    shape[0] = MIN(n,p);
    rb_taub = na_make_object(NA_SCOMPLEX, 1, shape, cNArray);
  }
  taub = NA_PTR_TYPE(rb_taub, complex*);
  {
    int shape[1];
    shape[0] = MAX(1,lwork);
    rb_work = na_make_object(NA_SCOMPLEX, 1, shape, cNArray);
  }
  work = NA_PTR_TYPE(rb_work, complex*);
  {
    int shape[2];
    shape[0] = lda;
    shape[1] = m;
    rb_a_out__ = na_make_object(NA_SCOMPLEX, 2, shape, cNArray);
  }
  a_out__ = NA_PTR_TYPE(rb_a_out__, complex*);
  MEMCPY(a_out__, a, complex, NA_TOTAL(rb_a));
  rb_a = rb_a_out__;
  a = a_out__;
  {
    int shape[2];
    shape[0] = ldb;
    shape[1] = p;
    rb_b_out__ = na_make_object(NA_SCOMPLEX, 2, shape, cNArray);
  }
  b_out__ = NA_PTR_TYPE(rb_b_out__, complex*);
  MEMCPY(b_out__, b, complex, NA_TOTAL(rb_b));
  rb_b = rb_b_out__;
  b = b_out__;

  cggqrf_(&n, &m, &p, a, &lda, taua, b, &ldb, taub, work, &lwork, &info);

  rb_info = INT2NUM(info);
  return rb_ary_new3(6, rb_taua, rb_taub, rb_work, rb_info, rb_a, rb_b);
}

void
init_lapack_cggqrf(VALUE mLapack){
  rb_define_module_function(mLapack, "cggqrf", rb_cggqrf, -1);
}

#include "rb_lapack.h"

extern VOID cunmtr_(char *side, char *uplo, char *trans, integer *m, integer *n, complex *a, integer *lda, complex *tau, complex *c, integer *ldc, complex *work, integer *lwork, integer *info);

static VALUE
rb_cunmtr(int argc, VALUE *argv, VALUE self){
  VALUE rb_side;
  char side; 
  VALUE rb_uplo;
  char uplo; 
  VALUE rb_trans;
  char trans; 
  VALUE rb_a;
  complex *a; 
  VALUE rb_tau;
  complex *tau; 
  VALUE rb_c;
  complex *c; 
  VALUE rb_lwork;
  integer lwork; 
  VALUE rb_work;
  complex *work; 
  VALUE rb_info;
  integer info; 
  VALUE rb_c_out__;
  complex *c_out__;

  integer lda;
  integer m;
  integer ldc;
  integer n;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  work, info, c = NumRu::Lapack.cunmtr( side, uplo, trans, a, tau, c, lwork)\n    or\n  NumRu::Lapack.cunmtr  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 7)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 7)", argc);
  rb_side = argv[0];
  rb_uplo = argv[1];
  rb_trans = argv[2];
  rb_a = argv[3];
  rb_tau = argv[4];
  rb_c = argv[5];
  rb_lwork = argv[6];

  if (!NA_IsNArray(rb_a))
    rb_raise(rb_eArgError, "a (4th argument) must be NArray");
  if (NA_RANK(rb_a) != 2)
    rb_raise(rb_eArgError, "rank of a (4th argument) must be %d", 2);
  m = NA_SHAPE1(rb_a);
  lda = NA_SHAPE0(rb_a);
  if (NA_TYPE(rb_a) != NA_SCOMPLEX)
    rb_a = na_change_type(rb_a, NA_SCOMPLEX);
  a = NA_PTR_TYPE(rb_a, complex*);
  side = StringValueCStr(rb_side)[0];
  if (!NA_IsNArray(rb_c))
    rb_raise(rb_eArgError, "c (6th argument) must be NArray");
  if (NA_RANK(rb_c) != 2)
    rb_raise(rb_eArgError, "rank of c (6th argument) must be %d", 2);
  n = NA_SHAPE1(rb_c);
  ldc = NA_SHAPE0(rb_c);
  if (NA_TYPE(rb_c) != NA_SCOMPLEX)
    rb_c = na_change_type(rb_c, NA_SCOMPLEX);
  c = NA_PTR_TYPE(rb_c, complex*);
  trans = StringValueCStr(rb_trans)[0];
  uplo = StringValueCStr(rb_uplo)[0];
  lwork = NUM2INT(rb_lwork);
  if (!NA_IsNArray(rb_tau))
    rb_raise(rb_eArgError, "tau (5th argument) must be NArray");
  if (NA_RANK(rb_tau) != 1)
    rb_raise(rb_eArgError, "rank of tau (5th argument) must be %d", 1);
  if (NA_SHAPE0(rb_tau) != (m-1))
    rb_raise(rb_eRuntimeError, "shape 0 of tau must be %d", m-1);
  if (NA_TYPE(rb_tau) != NA_SCOMPLEX)
    rb_tau = na_change_type(rb_tau, NA_SCOMPLEX);
  tau = NA_PTR_TYPE(rb_tau, complex*);
  {
    int shape[1];
    shape[0] = MAX(1,lwork);
    rb_work = na_make_object(NA_SCOMPLEX, 1, shape, cNArray);
  }
  work = NA_PTR_TYPE(rb_work, complex*);
  {
    int shape[2];
    shape[0] = ldc;
    shape[1] = n;
    rb_c_out__ = na_make_object(NA_SCOMPLEX, 2, shape, cNArray);
  }
  c_out__ = NA_PTR_TYPE(rb_c_out__, complex*);
  MEMCPY(c_out__, c, complex, NA_TOTAL(rb_c));
  rb_c = rb_c_out__;
  c = c_out__;

  cunmtr_(&side, &uplo, &trans, &m, &n, a, &lda, tau, c, &ldc, work, &lwork, &info);

  rb_info = INT2NUM(info);
  return rb_ary_new3(3, rb_work, rb_info, rb_c);
}

void
init_lapack_cunmtr(VALUE mLapack){
  rb_define_module_function(mLapack, "cunmtr", rb_cunmtr, -1);
}

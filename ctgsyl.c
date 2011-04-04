#include "rb_lapack.h"

extern VOID ctgsyl_(char *trans, integer *ijob, integer *m, integer *n, complex *a, integer *lda, complex *b, integer *ldb, complex *c, integer *ldc, complex *d, integer *ldd, complex *e, integer *lde, complex *f, integer *ldf, real *scale, real *dif, complex *work, integer *lwork, integer *iwork, integer *info);

static VALUE
rb_ctgsyl(int argc, VALUE *argv, VALUE self){
  VALUE rb_trans;
  char trans; 
  VALUE rb_ijob;
  integer ijob; 
  VALUE rb_a;
  complex *a; 
  VALUE rb_b;
  complex *b; 
  VALUE rb_c;
  complex *c; 
  VALUE rb_d;
  complex *d; 
  VALUE rb_e;
  complex *e; 
  VALUE rb_f;
  complex *f; 
  VALUE rb_lwork;
  integer lwork; 
  VALUE rb_scale;
  real scale; 
  VALUE rb_dif;
  real dif; 
  VALUE rb_work;
  complex *work; 
  VALUE rb_info;
  integer info; 
  VALUE rb_c_out__;
  complex *c_out__;
  VALUE rb_f_out__;
  complex *f_out__;
  integer *iwork;

  integer lda;
  integer m;
  integer ldb;
  integer n;
  integer ldc;
  integer ldd;
  integer lde;
  integer ldf;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  scale, dif, work, info, c, f = NumRu::Lapack.ctgsyl( trans, ijob, a, b, c, d, e, f, lwork)\n    or\n  NumRu::Lapack.ctgsyl  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 9)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 9)", argc);
  rb_trans = argv[0];
  rb_ijob = argv[1];
  rb_a = argv[2];
  rb_b = argv[3];
  rb_c = argv[4];
  rb_d = argv[5];
  rb_e = argv[6];
  rb_f = argv[7];
  rb_lwork = argv[8];

  ijob = NUM2INT(rb_ijob);
  if (!NA_IsNArray(rb_a))
    rb_raise(rb_eArgError, "a (3th argument) must be NArray");
  if (NA_RANK(rb_a) != 2)
    rb_raise(rb_eArgError, "rank of a (3th argument) must be %d", 2);
  m = NA_SHAPE1(rb_a);
  lda = NA_SHAPE0(rb_a);
  if (NA_TYPE(rb_a) != NA_SCOMPLEX)
    rb_a = na_change_type(rb_a, NA_SCOMPLEX);
  a = NA_PTR_TYPE(rb_a, complex*);
  if (!NA_IsNArray(rb_b))
    rb_raise(rb_eArgError, "b (4th argument) must be NArray");
  if (NA_RANK(rb_b) != 2)
    rb_raise(rb_eArgError, "rank of b (4th argument) must be %d", 2);
  n = NA_SHAPE1(rb_b);
  ldb = NA_SHAPE0(rb_b);
  if (NA_TYPE(rb_b) != NA_SCOMPLEX)
    rb_b = na_change_type(rb_b, NA_SCOMPLEX);
  b = NA_PTR_TYPE(rb_b, complex*);
  if (!NA_IsNArray(rb_c))
    rb_raise(rb_eArgError, "c (5th argument) must be NArray");
  if (NA_RANK(rb_c) != 2)
    rb_raise(rb_eArgError, "rank of c (5th argument) must be %d", 2);
  if (NA_SHAPE1(rb_c) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of c must be the same as shape 1 of b");
  ldc = NA_SHAPE0(rb_c);
  if (NA_TYPE(rb_c) != NA_SCOMPLEX)
    rb_c = na_change_type(rb_c, NA_SCOMPLEX);
  c = NA_PTR_TYPE(rb_c, complex*);
  if (!NA_IsNArray(rb_d))
    rb_raise(rb_eArgError, "d (6th argument) must be NArray");
  if (NA_RANK(rb_d) != 2)
    rb_raise(rb_eArgError, "rank of d (6th argument) must be %d", 2);
  if (NA_SHAPE1(rb_d) != m)
    rb_raise(rb_eRuntimeError, "shape 1 of d must be the same as shape 1 of a");
  ldd = NA_SHAPE0(rb_d);
  if (NA_TYPE(rb_d) != NA_SCOMPLEX)
    rb_d = na_change_type(rb_d, NA_SCOMPLEX);
  d = NA_PTR_TYPE(rb_d, complex*);
  lwork = NUM2INT(rb_lwork);
  if (!NA_IsNArray(rb_e))
    rb_raise(rb_eArgError, "e (7th argument) must be NArray");
  if (NA_RANK(rb_e) != 2)
    rb_raise(rb_eArgError, "rank of e (7th argument) must be %d", 2);
  if (NA_SHAPE1(rb_e) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of e must be the same as shape 1 of b");
  lde = NA_SHAPE0(rb_e);
  if (NA_TYPE(rb_e) != NA_SCOMPLEX)
    rb_e = na_change_type(rb_e, NA_SCOMPLEX);
  e = NA_PTR_TYPE(rb_e, complex*);
  if (!NA_IsNArray(rb_f))
    rb_raise(rb_eArgError, "f (8th argument) must be NArray");
  if (NA_RANK(rb_f) != 2)
    rb_raise(rb_eArgError, "rank of f (8th argument) must be %d", 2);
  if (NA_SHAPE1(rb_f) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of f must be the same as shape 1 of b");
  ldf = NA_SHAPE0(rb_f);
  if (NA_TYPE(rb_f) != NA_SCOMPLEX)
    rb_f = na_change_type(rb_f, NA_SCOMPLEX);
  f = NA_PTR_TYPE(rb_f, complex*);
  trans = StringValueCStr(rb_trans)[0];
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
  {
    int shape[2];
    shape[0] = ldf;
    shape[1] = n;
    rb_f_out__ = na_make_object(NA_SCOMPLEX, 2, shape, cNArray);
  }
  f_out__ = NA_PTR_TYPE(rb_f_out__, complex*);
  MEMCPY(f_out__, f, complex, NA_TOTAL(rb_f));
  rb_f = rb_f_out__;
  f = f_out__;
  iwork = ALLOC_N(integer, (m+n+2));

  ctgsyl_(&trans, &ijob, &m, &n, a, &lda, b, &ldb, c, &ldc, d, &ldd, e, &lde, f, &ldf, &scale, &dif, work, &lwork, iwork, &info);

  free(iwork);
  rb_scale = rb_float_new((double)scale);
  rb_dif = rb_float_new((double)dif);
  rb_info = INT2NUM(info);
  return rb_ary_new3(6, rb_scale, rb_dif, rb_work, rb_info, rb_c, rb_f);
}

void
init_lapack_ctgsyl(VALUE mLapack){
  rb_define_module_function(mLapack, "ctgsyl", rb_ctgsyl, -1);
}

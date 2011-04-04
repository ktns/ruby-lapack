#include "rb_lapack.h"

extern VOID dtgsy2_(char *trans, integer *ijob, integer *m, integer *n, doublereal *a, integer *lda, doublereal *b, integer *ldb, doublereal *c, integer *ldc, doublereal *d, integer *ldd, doublereal *e, integer *lde, doublereal *f, integer *ldf, doublereal *scale, doublereal *rdsum, doublereal *rdscal, integer *iwork, integer *pq, integer *info);

static VALUE
rb_dtgsy2(int argc, VALUE *argv, VALUE self){
  VALUE rb_trans;
  char trans; 
  VALUE rb_ijob;
  integer ijob; 
  VALUE rb_a;
  doublereal *a; 
  VALUE rb_b;
  doublereal *b; 
  VALUE rb_c;
  doublereal *c; 
  VALUE rb_d;
  doublereal *d; 
  VALUE rb_e;
  doublereal *e; 
  VALUE rb_f;
  doublereal *f; 
  VALUE rb_rdsum;
  doublereal rdsum; 
  VALUE rb_rdscal;
  doublereal rdscal; 
  VALUE rb_scale;
  doublereal scale; 
  VALUE rb_pq;
  integer pq; 
  VALUE rb_info;
  integer info; 
  VALUE rb_c_out__;
  doublereal *c_out__;
  VALUE rb_f_out__;
  doublereal *f_out__;
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
    printf("%s\n", "USAGE:\n  scale, pq, info, c, f, rdsum, rdscal = NumRu::Lapack.dtgsy2( trans, ijob, a, b, c, d, e, f, rdsum, rdscal)\n    or\n  NumRu::Lapack.dtgsy2  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 10)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 10)", argc);
  rb_trans = argv[0];
  rb_ijob = argv[1];
  rb_a = argv[2];
  rb_b = argv[3];
  rb_c = argv[4];
  rb_d = argv[5];
  rb_e = argv[6];
  rb_f = argv[7];
  rb_rdsum = argv[8];
  rb_rdscal = argv[9];

  rdscal = NUM2DBL(rb_rdscal);
  ijob = NUM2INT(rb_ijob);
  if (!NA_IsNArray(rb_a))
    rb_raise(rb_eArgError, "a (3th argument) must be NArray");
  if (NA_RANK(rb_a) != 2)
    rb_raise(rb_eArgError, "rank of a (3th argument) must be %d", 2);
  m = NA_SHAPE1(rb_a);
  lda = NA_SHAPE0(rb_a);
  if (NA_TYPE(rb_a) != NA_DFLOAT)
    rb_a = na_change_type(rb_a, NA_DFLOAT);
  a = NA_PTR_TYPE(rb_a, doublereal*);
  if (!NA_IsNArray(rb_b))
    rb_raise(rb_eArgError, "b (4th argument) must be NArray");
  if (NA_RANK(rb_b) != 2)
    rb_raise(rb_eArgError, "rank of b (4th argument) must be %d", 2);
  n = NA_SHAPE1(rb_b);
  ldb = NA_SHAPE0(rb_b);
  if (NA_TYPE(rb_b) != NA_DFLOAT)
    rb_b = na_change_type(rb_b, NA_DFLOAT);
  b = NA_PTR_TYPE(rb_b, doublereal*);
  if (!NA_IsNArray(rb_c))
    rb_raise(rb_eArgError, "c (5th argument) must be NArray");
  if (NA_RANK(rb_c) != 2)
    rb_raise(rb_eArgError, "rank of c (5th argument) must be %d", 2);
  if (NA_SHAPE1(rb_c) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of c must be the same as shape 1 of b");
  ldc = NA_SHAPE0(rb_c);
  if (NA_TYPE(rb_c) != NA_DFLOAT)
    rb_c = na_change_type(rb_c, NA_DFLOAT);
  c = NA_PTR_TYPE(rb_c, doublereal*);
  if (!NA_IsNArray(rb_d))
    rb_raise(rb_eArgError, "d (6th argument) must be NArray");
  if (NA_RANK(rb_d) != 2)
    rb_raise(rb_eArgError, "rank of d (6th argument) must be %d", 2);
  if (NA_SHAPE1(rb_d) != m)
    rb_raise(rb_eRuntimeError, "shape 1 of d must be the same as shape 1 of a");
  ldd = NA_SHAPE0(rb_d);
  if (NA_TYPE(rb_d) != NA_DFLOAT)
    rb_d = na_change_type(rb_d, NA_DFLOAT);
  d = NA_PTR_TYPE(rb_d, doublereal*);
  rdsum = NUM2DBL(rb_rdsum);
  if (!NA_IsNArray(rb_e))
    rb_raise(rb_eArgError, "e (7th argument) must be NArray");
  if (NA_RANK(rb_e) != 2)
    rb_raise(rb_eArgError, "rank of e (7th argument) must be %d", 2);
  if (NA_SHAPE1(rb_e) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of e must be the same as shape 1 of b");
  lde = NA_SHAPE0(rb_e);
  if (NA_TYPE(rb_e) != NA_DFLOAT)
    rb_e = na_change_type(rb_e, NA_DFLOAT);
  e = NA_PTR_TYPE(rb_e, doublereal*);
  if (!NA_IsNArray(rb_f))
    rb_raise(rb_eArgError, "f (8th argument) must be NArray");
  if (NA_RANK(rb_f) != 2)
    rb_raise(rb_eArgError, "rank of f (8th argument) must be %d", 2);
  if (NA_SHAPE1(rb_f) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of f must be the same as shape 1 of b");
  ldf = NA_SHAPE0(rb_f);
  if (NA_TYPE(rb_f) != NA_DFLOAT)
    rb_f = na_change_type(rb_f, NA_DFLOAT);
  f = NA_PTR_TYPE(rb_f, doublereal*);
  trans = StringValueCStr(rb_trans)[0];
  {
    int shape[2];
    shape[0] = ldc;
    shape[1] = n;
    rb_c_out__ = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  c_out__ = NA_PTR_TYPE(rb_c_out__, doublereal*);
  MEMCPY(c_out__, c, doublereal, NA_TOTAL(rb_c));
  rb_c = rb_c_out__;
  c = c_out__;
  {
    int shape[2];
    shape[0] = ldf;
    shape[1] = n;
    rb_f_out__ = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  f_out__ = NA_PTR_TYPE(rb_f_out__, doublereal*);
  MEMCPY(f_out__, f, doublereal, NA_TOTAL(rb_f));
  rb_f = rb_f_out__;
  f = f_out__;
  iwork = ALLOC_N(integer, (m+n+2));

  dtgsy2_(&trans, &ijob, &m, &n, a, &lda, b, &ldb, c, &ldc, d, &ldd, e, &lde, f, &ldf, &scale, &rdsum, &rdscal, iwork, &pq, &info);

  free(iwork);
  rb_scale = rb_float_new((double)scale);
  rb_pq = INT2NUM(pq);
  rb_info = INT2NUM(info);
  rb_rdsum = rb_float_new((double)rdsum);
  rb_rdscal = rb_float_new((double)rdscal);
  return rb_ary_new3(7, rb_scale, rb_pq, rb_info, rb_c, rb_f, rb_rdsum, rb_rdscal);
}

void
init_lapack_dtgsy2(VALUE mLapack){
  rb_define_module_function(mLapack, "dtgsy2", rb_dtgsy2, -1);
}

#include "rb_lapack.h"

extern VOID cgebd2_(integer *m, integer *n, complex *a, integer *lda, real *d, real *e, complex *tauq, complex *taup, complex *work, integer *info);

static VALUE
rb_cgebd2(int argc, VALUE *argv, VALUE self){
  VALUE rb_m;
  integer m; 
  VALUE rb_a;
  complex *a; 
  VALUE rb_d;
  real *d; 
  VALUE rb_e;
  real *e; 
  VALUE rb_tauq;
  complex *tauq; 
  VALUE rb_taup;
  complex *taup; 
  VALUE rb_info;
  integer info; 
  VALUE rb_a_out__;
  complex *a_out__;
  complex *work;

  integer lda;
  integer n;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  d, e, tauq, taup, info, a = NumRu::Lapack.cgebd2( m, a)\n    or\n  NumRu::Lapack.cgebd2  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 2)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 2)", argc);
  rb_m = argv[0];
  rb_a = argv[1];

  if (!NA_IsNArray(rb_a))
    rb_raise(rb_eArgError, "a (2th argument) must be NArray");
  if (NA_RANK(rb_a) != 2)
    rb_raise(rb_eArgError, "rank of a (2th argument) must be %d", 2);
  n = NA_SHAPE1(rb_a);
  lda = NA_SHAPE0(rb_a);
  if (NA_TYPE(rb_a) != NA_SCOMPLEX)
    rb_a = na_change_type(rb_a, NA_SCOMPLEX);
  a = NA_PTR_TYPE(rb_a, complex*);
  m = NUM2INT(rb_m);
  {
    int shape[1];
    shape[0] = MIN(m,n);
    rb_d = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  d = NA_PTR_TYPE(rb_d, real*);
  {
    int shape[1];
    shape[0] = MIN(m,n)-1;
    rb_e = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  e = NA_PTR_TYPE(rb_e, real*);
  {
    int shape[1];
    shape[0] = MIN(m,n);
    rb_tauq = na_make_object(NA_SCOMPLEX, 1, shape, cNArray);
  }
  tauq = NA_PTR_TYPE(rb_tauq, complex*);
  {
    int shape[1];
    shape[0] = MIN(m,n);
    rb_taup = na_make_object(NA_SCOMPLEX, 1, shape, cNArray);
  }
  taup = NA_PTR_TYPE(rb_taup, complex*);
  {
    int shape[2];
    shape[0] = lda;
    shape[1] = n;
    rb_a_out__ = na_make_object(NA_SCOMPLEX, 2, shape, cNArray);
  }
  a_out__ = NA_PTR_TYPE(rb_a_out__, complex*);
  MEMCPY(a_out__, a, complex, NA_TOTAL(rb_a));
  rb_a = rb_a_out__;
  a = a_out__;
  work = ALLOC_N(complex, (MAX(m,n)));

  cgebd2_(&m, &n, a, &lda, d, e, tauq, taup, work, &info);

  free(work);
  rb_info = INT2NUM(info);
  return rb_ary_new3(6, rb_d, rb_e, rb_tauq, rb_taup, rb_info, rb_a);
}

void
init_lapack_cgebd2(VALUE mLapack){
  rb_define_module_function(mLapack, "cgebd2", rb_cgebd2, -1);
}

#include "rb_lapack.h"

extern VOID ssbtrd_(char *vect, char *uplo, integer *n, integer *kd, real *ab, integer *ldab, real *d, real *e, real *q, integer *ldq, real *work, integer *info);

static VALUE
rb_ssbtrd(int argc, VALUE *argv, VALUE self){
  VALUE rb_vect;
  char vect; 
  VALUE rb_uplo;
  char uplo; 
  VALUE rb_kd;
  integer kd; 
  VALUE rb_ab;
  real *ab; 
  VALUE rb_q;
  real *q; 
  VALUE rb_d;
  real *d; 
  VALUE rb_e;
  real *e; 
  VALUE rb_info;
  integer info; 
  VALUE rb_ab_out__;
  real *ab_out__;
  VALUE rb_q_out__;
  real *q_out__;
  real *work;

  integer ldab;
  integer n;
  integer ldq;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  d, e, info, ab, q = NumRu::Lapack.ssbtrd( vect, uplo, kd, ab, q)\n    or\n  NumRu::Lapack.ssbtrd  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 5)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 5)", argc);
  rb_vect = argv[0];
  rb_uplo = argv[1];
  rb_kd = argv[2];
  rb_ab = argv[3];
  rb_q = argv[4];

  if (!NA_IsNArray(rb_ab))
    rb_raise(rb_eArgError, "ab (4th argument) must be NArray");
  if (NA_RANK(rb_ab) != 2)
    rb_raise(rb_eArgError, "rank of ab (4th argument) must be %d", 2);
  n = NA_SHAPE1(rb_ab);
  ldab = NA_SHAPE0(rb_ab);
  if (NA_TYPE(rb_ab) != NA_SFLOAT)
    rb_ab = na_change_type(rb_ab, NA_SFLOAT);
  ab = NA_PTR_TYPE(rb_ab, real*);
  vect = StringValueCStr(rb_vect)[0];
  kd = NUM2INT(rb_kd);
  if (!NA_IsNArray(rb_q))
    rb_raise(rb_eArgError, "q (5th argument) must be NArray");
  if (NA_RANK(rb_q) != 2)
    rb_raise(rb_eArgError, "rank of q (5th argument) must be %d", 2);
  if (NA_SHAPE1(rb_q) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of q must be the same as shape 1 of ab");
  ldq = NA_SHAPE0(rb_q);
  if (NA_TYPE(rb_q) != NA_SFLOAT)
    rb_q = na_change_type(rb_q, NA_SFLOAT);
  q = NA_PTR_TYPE(rb_q, real*);
  uplo = StringValueCStr(rb_uplo)[0];
  {
    int shape[1];
    shape[0] = n;
    rb_d = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  d = NA_PTR_TYPE(rb_d, real*);
  {
    int shape[1];
    shape[0] = n-1;
    rb_e = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  e = NA_PTR_TYPE(rb_e, real*);
  {
    int shape[2];
    shape[0] = ldab;
    shape[1] = n;
    rb_ab_out__ = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  ab_out__ = NA_PTR_TYPE(rb_ab_out__, real*);
  MEMCPY(ab_out__, ab, real, NA_TOTAL(rb_ab));
  rb_ab = rb_ab_out__;
  ab = ab_out__;
  {
    int shape[2];
    shape[0] = ldq;
    shape[1] = n;
    rb_q_out__ = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  q_out__ = NA_PTR_TYPE(rb_q_out__, real*);
  MEMCPY(q_out__, q, real, NA_TOTAL(rb_q));
  rb_q = rb_q_out__;
  q = q_out__;
  work = ALLOC_N(real, (n));

  ssbtrd_(&vect, &uplo, &n, &kd, ab, &ldab, d, e, q, &ldq, work, &info);

  free(work);
  rb_info = INT2NUM(info);
  return rb_ary_new3(5, rb_d, rb_e, rb_info, rb_ab, rb_q);
}

void
init_lapack_ssbtrd(VALUE mLapack){
  rb_define_module_function(mLapack, "ssbtrd", rb_ssbtrd, -1);
}

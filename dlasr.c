#include "rb_lapack.h"

extern VOID dlasr_(char *side, char *pivot, char *direct, integer *m, integer *n, doublereal *c, doublereal *s, doublereal *a, integer *lda);

static VALUE
rb_dlasr(int argc, VALUE *argv, VALUE self){
  VALUE rb_side;
  char side; 
  VALUE rb_pivot;
  char pivot; 
  VALUE rb_direct;
  char direct; 
  VALUE rb_m;
  integer m; 
  VALUE rb_c;
  doublereal *c; 
  VALUE rb_s;
  doublereal *s; 
  VALUE rb_a;
  doublereal *a; 
  VALUE rb_a_out__;
  doublereal *a_out__;

  integer lda;
  integer n;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  a = NumRu::Lapack.dlasr( side, pivot, direct, m, c, s, a)\n    or\n  NumRu::Lapack.dlasr  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 7)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 7)", argc);
  rb_side = argv[0];
  rb_pivot = argv[1];
  rb_direct = argv[2];
  rb_m = argv[3];
  rb_c = argv[4];
  rb_s = argv[5];
  rb_a = argv[6];

  direct = StringValueCStr(rb_direct)[0];
  if (!NA_IsNArray(rb_a))
    rb_raise(rb_eArgError, "a (7th argument) must be NArray");
  if (NA_RANK(rb_a) != 2)
    rb_raise(rb_eArgError, "rank of a (7th argument) must be %d", 2);
  n = NA_SHAPE1(rb_a);
  lda = NA_SHAPE0(rb_a);
  if (NA_TYPE(rb_a) != NA_DFLOAT)
    rb_a = na_change_type(rb_a, NA_DFLOAT);
  a = NA_PTR_TYPE(rb_a, doublereal*);
  side = StringValueCStr(rb_side)[0];
  pivot = StringValueCStr(rb_pivot)[0];
  m = NUM2INT(rb_m);
  if (!NA_IsNArray(rb_c))
    rb_raise(rb_eArgError, "c (5th argument) must be NArray");
  if (NA_RANK(rb_c) != 1)
    rb_raise(rb_eArgError, "rank of c (5th argument) must be %d", 1);
  if (NA_SHAPE0(rb_c) != (m-1))
    rb_raise(rb_eRuntimeError, "shape 0 of c must be %d", m-1);
  if (NA_TYPE(rb_c) != NA_DFLOAT)
    rb_c = na_change_type(rb_c, NA_DFLOAT);
  c = NA_PTR_TYPE(rb_c, doublereal*);
  if (!NA_IsNArray(rb_s))
    rb_raise(rb_eArgError, "s (6th argument) must be NArray");
  if (NA_RANK(rb_s) != 1)
    rb_raise(rb_eArgError, "rank of s (6th argument) must be %d", 1);
  if (NA_SHAPE0(rb_s) != (m-1))
    rb_raise(rb_eRuntimeError, "shape 0 of s must be %d", m-1);
  if (NA_TYPE(rb_s) != NA_DFLOAT)
    rb_s = na_change_type(rb_s, NA_DFLOAT);
  s = NA_PTR_TYPE(rb_s, doublereal*);
  {
    int shape[2];
    shape[0] = lda;
    shape[1] = n;
    rb_a_out__ = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  a_out__ = NA_PTR_TYPE(rb_a_out__, doublereal*);
  MEMCPY(a_out__, a, doublereal, NA_TOTAL(rb_a));
  rb_a = rb_a_out__;
  a = a_out__;

  dlasr_(&side, &pivot, &direct, &m, &n, c, s, a, &lda);

  return rb_a;
}

void
init_lapack_dlasr(VALUE mLapack){
  rb_define_module_function(mLapack, "dlasr", rb_dlasr, -1);
}

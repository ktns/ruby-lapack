#include "rb_lapack.h"

extern VOID dlatrd_(char *uplo, integer *n, integer *nb, doublereal *a, integer *lda, doublereal *e, doublereal *tau, doublereal *w, integer *ldw);

static VALUE
rb_dlatrd(int argc, VALUE *argv, VALUE self){
  VALUE rb_uplo;
  char uplo; 
  VALUE rb_nb;
  integer nb; 
  VALUE rb_a;
  doublereal *a; 
  VALUE rb_e;
  doublereal *e; 
  VALUE rb_tau;
  doublereal *tau; 
  VALUE rb_w;
  doublereal *w; 
  VALUE rb_a_out__;
  doublereal *a_out__;

  integer lda;
  integer n;
  integer ldw;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  e, tau, w, a = NumRu::Lapack.dlatrd( uplo, nb, a)\n    or\n  NumRu::Lapack.dlatrd  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 3)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 3)", argc);
  rb_uplo = argv[0];
  rb_nb = argv[1];
  rb_a = argv[2];

  if (!NA_IsNArray(rb_a))
    rb_raise(rb_eArgError, "a (3th argument) must be NArray");
  if (NA_RANK(rb_a) != 2)
    rb_raise(rb_eArgError, "rank of a (3th argument) must be %d", 2);
  n = NA_SHAPE1(rb_a);
  lda = NA_SHAPE0(rb_a);
  if (NA_TYPE(rb_a) != NA_DFLOAT)
    rb_a = na_change_type(rb_a, NA_DFLOAT);
  a = NA_PTR_TYPE(rb_a, doublereal*);
  uplo = StringValueCStr(rb_uplo)[0];
  nb = NUM2INT(rb_nb);
  ldw = MAX(1,n);
  {
    int shape[1];
    shape[0] = n-1;
    rb_e = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  e = NA_PTR_TYPE(rb_e, doublereal*);
  {
    int shape[1];
    shape[0] = n-1;
    rb_tau = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  tau = NA_PTR_TYPE(rb_tau, doublereal*);
  {
    int shape[2];
    shape[0] = ldw;
    shape[1] = MAX(n,nb);
    rb_w = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  w = NA_PTR_TYPE(rb_w, doublereal*);
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

  dlatrd_(&uplo, &n, &nb, a, &lda, e, tau, w, &ldw);

  return rb_ary_new3(4, rb_e, rb_tau, rb_w, rb_a);
}

void
init_lapack_dlatrd(VALUE mLapack){
  rb_define_module_function(mLapack, "dlatrd", rb_dlatrd, -1);
}

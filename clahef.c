#include "rb_lapack.h"

extern VOID clahef_(char *uplo, integer *n, integer *nb, integer *kb, complex *a, integer *lda, integer *ipiv, complex *w, integer *ldw, integer *info);

static VALUE
rb_clahef(int argc, VALUE *argv, VALUE self){
  VALUE rb_uplo;
  char uplo; 
  VALUE rb_nb;
  integer nb; 
  VALUE rb_a;
  complex *a; 
  VALUE rb_kb;
  integer kb; 
  VALUE rb_ipiv;
  integer *ipiv; 
  VALUE rb_info;
  integer info; 
  VALUE rb_a_out__;
  complex *a_out__;
  complex *w;

  integer lda;
  integer n;
  integer ldw;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  kb, ipiv, info, a = NumRu::Lapack.clahef( uplo, nb, a)\n    or\n  NumRu::Lapack.clahef  # print help\n\n\nFORTRAN MANUAL\n\n");
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
  if (NA_TYPE(rb_a) != NA_SCOMPLEX)
    rb_a = na_change_type(rb_a, NA_SCOMPLEX);
  a = NA_PTR_TYPE(rb_a, complex*);
  uplo = StringValueCStr(rb_uplo)[0];
  nb = NUM2INT(rb_nb);
  ldw = MAX(1,n);
  {
    int shape[1];
    shape[0] = n;
    rb_ipiv = na_make_object(NA_LINT, 1, shape, cNArray);
  }
  ipiv = NA_PTR_TYPE(rb_ipiv, integer*);
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
  w = ALLOC_N(complex, (ldw)*(MAX(n,nb)));

  clahef_(&uplo, &n, &nb, &kb, a, &lda, ipiv, w, &ldw, &info);

  free(w);
  rb_kb = INT2NUM(kb);
  rb_info = INT2NUM(info);
  return rb_ary_new3(4, rb_kb, rb_ipiv, rb_info, rb_a);
}

void
init_lapack_clahef(VALUE mLapack){
  rb_define_module_function(mLapack, "clahef", rb_clahef, -1);
}

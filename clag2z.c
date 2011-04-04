#include "rb_lapack.h"

extern VOID clag2z_(integer *m, integer *n, complex *sa, integer *ldsa, doublecomplex *a, integer *lda, integer *info);

static VALUE
rb_clag2z(int argc, VALUE *argv, VALUE self){
  VALUE rb_m;
  integer m; 
  VALUE rb_sa;
  complex *sa; 
  VALUE rb_a;
  doublecomplex *a; 
  VALUE rb_info;
  integer info; 

  integer ldsa;
  integer n;
  integer lda;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  a, info = NumRu::Lapack.clag2z( m, sa)\n    or\n  NumRu::Lapack.clag2z  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 2)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 2)", argc);
  rb_m = argv[0];
  rb_sa = argv[1];

  m = NUM2INT(rb_m);
  if (!NA_IsNArray(rb_sa))
    rb_raise(rb_eArgError, "sa (2th argument) must be NArray");
  if (NA_RANK(rb_sa) != 2)
    rb_raise(rb_eArgError, "rank of sa (2th argument) must be %d", 2);
  n = NA_SHAPE1(rb_sa);
  ldsa = NA_SHAPE0(rb_sa);
  if (NA_TYPE(rb_sa) != NA_SCOMPLEX)
    rb_sa = na_change_type(rb_sa, NA_SCOMPLEX);
  sa = NA_PTR_TYPE(rb_sa, complex*);
  lda = MAX(1,m);
  {
    int shape[2];
    shape[0] = lda;
    shape[1] = n;
    rb_a = na_make_object(NA_DCOMPLEX, 2, shape, cNArray);
  }
  a = NA_PTR_TYPE(rb_a, doublecomplex*);

  clag2z_(&m, &n, sa, &ldsa, a, &lda, &info);

  rb_info = INT2NUM(info);
  return rb_ary_new3(2, rb_a, rb_info);
}

void
init_lapack_clag2z(VALUE mLapack){
  rb_define_module_function(mLapack, "clag2z", rb_clag2z, -1);
}

#include "rb_lapack.h"

extern VOID sopgtr_(char *uplo, integer *n, real *ap, real *tau, real *q, integer *ldq, real *work, integer *info);

static VALUE
rb_sopgtr(int argc, VALUE *argv, VALUE self){
  VALUE rb_uplo;
  char uplo; 
  VALUE rb_ap;
  real *ap; 
  VALUE rb_tau;
  real *tau; 
  VALUE rb_q;
  real *q; 
  VALUE rb_info;
  integer info; 
  real *work;

  integer ldap;
  integer ldtau;
  integer ldq;
  integer n;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  q, info = NumRu::Lapack.sopgtr( uplo, ap, tau)\n    or\n  NumRu::Lapack.sopgtr  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 3)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 3)", argc);
  rb_uplo = argv[0];
  rb_ap = argv[1];
  rb_tau = argv[2];

  uplo = StringValueCStr(rb_uplo)[0];
  if (!NA_IsNArray(rb_tau))
    rb_raise(rb_eArgError, "tau (3th argument) must be NArray");
  if (NA_RANK(rb_tau) != 1)
    rb_raise(rb_eArgError, "rank of tau (3th argument) must be %d", 1);
  ldtau = NA_SHAPE0(rb_tau);
  if (NA_TYPE(rb_tau) != NA_SFLOAT)
    rb_tau = na_change_type(rb_tau, NA_SFLOAT);
  tau = NA_PTR_TYPE(rb_tau, real*);
  if (!NA_IsNArray(rb_ap))
    rb_raise(rb_eArgError, "ap (2th argument) must be NArray");
  if (NA_RANK(rb_ap) != 1)
    rb_raise(rb_eArgError, "rank of ap (2th argument) must be %d", 1);
  ldap = NA_SHAPE0(rb_ap);
  if (NA_TYPE(rb_ap) != NA_SFLOAT)
    rb_ap = na_change_type(rb_ap, NA_SFLOAT);
  ap = NA_PTR_TYPE(rb_ap, real*);
  n = ldtau+1;
  ldq = MAX(1,n);
  {
    int shape[2];
    shape[0] = ldq;
    shape[1] = n;
    rb_q = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  q = NA_PTR_TYPE(rb_q, real*);
  work = ALLOC_N(real, (n-1));

  sopgtr_(&uplo, &n, ap, tau, q, &ldq, work, &info);

  free(work);
  rb_info = INT2NUM(info);
  return rb_ary_new3(2, rb_q, rb_info);
}

void
init_lapack_sopgtr(VALUE mLapack){
  rb_define_module_function(mLapack, "sopgtr", rb_sopgtr, -1);
}

#include "rb_lapack.h"

extern doublereal dla_gercond_(char *trans, integer *n, doublereal *a, integer *lda, doublereal *af, integer *ldaf, integer *ipiv, integer *cmode, doublereal *c, integer *info, doublereal *work, integer *iwork);

static VALUE
rb_dla_gercond(int argc, VALUE *argv, VALUE self){
  VALUE rb_trans;
  char trans; 
  VALUE rb_a;
  doublereal *a; 
  VALUE rb_af;
  doublereal *af; 
  VALUE rb_ipiv;
  integer *ipiv; 
  VALUE rb_cmode;
  integer cmode; 
  VALUE rb_c;
  doublereal *c; 
  VALUE rb_work;
  doublereal *work; 
  VALUE rb_iwork;
  integer *iwork; 
  VALUE rb_info;
  integer info; 
  VALUE rb___out__;
  doublereal __out__; 

  integer lda;
  integer n;
  integer ldaf;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  info, __out__ = NumRu::Lapack.dla_gercond( trans, a, af, ipiv, cmode, c, work, iwork)\n    or\n  NumRu::Lapack.dla_gercond  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 8)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 8)", argc);
  rb_trans = argv[0];
  rb_a = argv[1];
  rb_af = argv[2];
  rb_ipiv = argv[3];
  rb_cmode = argv[4];
  rb_c = argv[5];
  rb_work = argv[6];
  rb_iwork = argv[7];

  if (!NA_IsNArray(rb_ipiv))
    rb_raise(rb_eArgError, "ipiv (4th argument) must be NArray");
  if (NA_RANK(rb_ipiv) != 1)
    rb_raise(rb_eArgError, "rank of ipiv (4th argument) must be %d", 1);
  n = NA_SHAPE0(rb_ipiv);
  if (NA_TYPE(rb_ipiv) != NA_LINT)
    rb_ipiv = na_change_type(rb_ipiv, NA_LINT);
  ipiv = NA_PTR_TYPE(rb_ipiv, integer*);
  if (!NA_IsNArray(rb_a))
    rb_raise(rb_eArgError, "a (2th argument) must be NArray");
  if (NA_RANK(rb_a) != 2)
    rb_raise(rb_eArgError, "rank of a (2th argument) must be %d", 2);
  if (NA_SHAPE1(rb_a) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of a must be the same as shape 0 of ipiv");
  lda = NA_SHAPE0(rb_a);
  if (NA_TYPE(rb_a) != NA_DFLOAT)
    rb_a = na_change_type(rb_a, NA_DFLOAT);
  a = NA_PTR_TYPE(rb_a, doublereal*);
  if (!NA_IsNArray(rb_c))
    rb_raise(rb_eArgError, "c (6th argument) must be NArray");
  if (NA_RANK(rb_c) != 1)
    rb_raise(rb_eArgError, "rank of c (6th argument) must be %d", 1);
  if (NA_SHAPE0(rb_c) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of c must be the same as shape 0 of ipiv");
  if (NA_TYPE(rb_c) != NA_DFLOAT)
    rb_c = na_change_type(rb_c, NA_DFLOAT);
  c = NA_PTR_TYPE(rb_c, doublereal*);
  if (!NA_IsNArray(rb_af))
    rb_raise(rb_eArgError, "af (3th argument) must be NArray");
  if (NA_RANK(rb_af) != 2)
    rb_raise(rb_eArgError, "rank of af (3th argument) must be %d", 2);
  if (NA_SHAPE1(rb_af) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of af must be the same as shape 0 of ipiv");
  ldaf = NA_SHAPE0(rb_af);
  if (NA_TYPE(rb_af) != NA_DFLOAT)
    rb_af = na_change_type(rb_af, NA_DFLOAT);
  af = NA_PTR_TYPE(rb_af, doublereal*);
  trans = StringValueCStr(rb_trans)[0];
  cmode = NUM2INT(rb_cmode);
  if (!NA_IsNArray(rb_iwork))
    rb_raise(rb_eArgError, "iwork (8th argument) must be NArray");
  if (NA_RANK(rb_iwork) != 1)
    rb_raise(rb_eArgError, "rank of iwork (8th argument) must be %d", 1);
  if (NA_SHAPE0(rb_iwork) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of iwork must be the same as shape 0 of ipiv");
  if (NA_TYPE(rb_iwork) != NA_LINT)
    rb_iwork = na_change_type(rb_iwork, NA_LINT);
  iwork = NA_PTR_TYPE(rb_iwork, integer*);
  if (!NA_IsNArray(rb_work))
    rb_raise(rb_eArgError, "work (7th argument) must be NArray");
  if (NA_RANK(rb_work) != 1)
    rb_raise(rb_eArgError, "rank of work (7th argument) must be %d", 1);
  if (NA_SHAPE0(rb_work) != (3*n))
    rb_raise(rb_eRuntimeError, "shape 0 of work must be %d", 3*n);
  if (NA_TYPE(rb_work) != NA_DFLOAT)
    rb_work = na_change_type(rb_work, NA_DFLOAT);
  work = NA_PTR_TYPE(rb_work, doublereal*);

  __out__ = dla_gercond_(&trans, &n, a, &lda, af, &ldaf, ipiv, &cmode, c, &info, work, iwork);

  rb_info = INT2NUM(info);
  rb___out__ = rb_float_new((double)__out__);
  return rb_ary_new3(2, rb_info, rb___out__);
}

void
init_lapack_dla_gercond(VALUE mLapack){
  rb_define_module_function(mLapack, "dla_gercond", rb_dla_gercond, -1);
}

#include "rb_lapack.h"

extern VOID zlatps_(char *uplo, char *trans, char *diag, char *normin, integer *n, doublecomplex *ap, doublecomplex *x, doublereal *scale, doublereal *cnorm, integer *info);

static VALUE
rb_zlatps(int argc, VALUE *argv, VALUE self){
  VALUE rb_uplo;
  char uplo; 
  VALUE rb_trans;
  char trans; 
  VALUE rb_diag;
  char diag; 
  VALUE rb_normin;
  char normin; 
  VALUE rb_ap;
  doublecomplex *ap; 
  VALUE rb_x;
  doublecomplex *x; 
  VALUE rb_cnorm;
  doublereal *cnorm; 
  VALUE rb_scale;
  doublereal scale; 
  VALUE rb_info;
  integer info; 
  VALUE rb_x_out__;
  doublecomplex *x_out__;
  VALUE rb_cnorm_out__;
  doublereal *cnorm_out__;

  integer n;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  scale, info, x, cnorm = NumRu::Lapack.zlatps( uplo, trans, diag, normin, ap, x, cnorm)\n    or\n  NumRu::Lapack.zlatps  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 7)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 7)", argc);
  rb_uplo = argv[0];
  rb_trans = argv[1];
  rb_diag = argv[2];
  rb_normin = argv[3];
  rb_ap = argv[4];
  rb_x = argv[5];
  rb_cnorm = argv[6];

  if (!NA_IsNArray(rb_cnorm))
    rb_raise(rb_eArgError, "cnorm (7th argument) must be NArray");
  if (NA_RANK(rb_cnorm) != 1)
    rb_raise(rb_eArgError, "rank of cnorm (7th argument) must be %d", 1);
  n = NA_SHAPE0(rb_cnorm);
  if (NA_TYPE(rb_cnorm) != NA_DFLOAT)
    rb_cnorm = na_change_type(rb_cnorm, NA_DFLOAT);
  cnorm = NA_PTR_TYPE(rb_cnorm, doublereal*);
  if (!NA_IsNArray(rb_x))
    rb_raise(rb_eArgError, "x (6th argument) must be NArray");
  if (NA_RANK(rb_x) != 1)
    rb_raise(rb_eArgError, "rank of x (6th argument) must be %d", 1);
  if (NA_SHAPE0(rb_x) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of x must be the same as shape 0 of cnorm");
  if (NA_TYPE(rb_x) != NA_DCOMPLEX)
    rb_x = na_change_type(rb_x, NA_DCOMPLEX);
  x = NA_PTR_TYPE(rb_x, doublecomplex*);
  diag = StringValueCStr(rb_diag)[0];
  trans = StringValueCStr(rb_trans)[0];
  normin = StringValueCStr(rb_normin)[0];
  uplo = StringValueCStr(rb_uplo)[0];
  if (!NA_IsNArray(rb_ap))
    rb_raise(rb_eArgError, "ap (5th argument) must be NArray");
  if (NA_RANK(rb_ap) != 1)
    rb_raise(rb_eArgError, "rank of ap (5th argument) must be %d", 1);
  if (NA_SHAPE0(rb_ap) != (n*(n+1)/2))
    rb_raise(rb_eRuntimeError, "shape 0 of ap must be %d", n*(n+1)/2);
  if (NA_TYPE(rb_ap) != NA_DCOMPLEX)
    rb_ap = na_change_type(rb_ap, NA_DCOMPLEX);
  ap = NA_PTR_TYPE(rb_ap, doublecomplex*);
  {
    int shape[1];
    shape[0] = n;
    rb_x_out__ = na_make_object(NA_DCOMPLEX, 1, shape, cNArray);
  }
  x_out__ = NA_PTR_TYPE(rb_x_out__, doublecomplex*);
  MEMCPY(x_out__, x, doublecomplex, NA_TOTAL(rb_x));
  rb_x = rb_x_out__;
  x = x_out__;
  {
    int shape[1];
    shape[0] = n;
    rb_cnorm_out__ = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  cnorm_out__ = NA_PTR_TYPE(rb_cnorm_out__, doublereal*);
  MEMCPY(cnorm_out__, cnorm, doublereal, NA_TOTAL(rb_cnorm));
  rb_cnorm = rb_cnorm_out__;
  cnorm = cnorm_out__;

  zlatps_(&uplo, &trans, &diag, &normin, &n, ap, x, &scale, cnorm, &info);

  rb_scale = rb_float_new((double)scale);
  rb_info = INT2NUM(info);
  return rb_ary_new3(4, rb_scale, rb_info, rb_x, rb_cnorm);
}

void
init_lapack_zlatps(VALUE mLapack){
  rb_define_module_function(mLapack, "zlatps", rb_zlatps, -1);
}

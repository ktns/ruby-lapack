#include "rb_lapack.h"

extern VOID chbevd_(char *jobz, char *uplo, integer *n, integer *kd, complex *ab, integer *ldab, real *w, complex *z, integer *ldz, complex *work, integer *lwork, real *rwork, integer *lrwork, integer *iwork, integer *liwork, integer *info);

static VALUE
rb_chbevd(int argc, VALUE *argv, VALUE self){
  VALUE rb_jobz;
  char jobz; 
  VALUE rb_uplo;
  char uplo; 
  VALUE rb_kd;
  integer kd; 
  VALUE rb_ab;
  complex *ab; 
  VALUE rb_lwork;
  integer lwork; 
  VALUE rb_lrwork;
  integer lrwork; 
  VALUE rb_liwork;
  integer liwork; 
  VALUE rb_w;
  real *w; 
  VALUE rb_z;
  complex *z; 
  VALUE rb_work;
  complex *work; 
  VALUE rb_rwork;
  real *rwork; 
  VALUE rb_iwork;
  integer *iwork; 
  VALUE rb_info;
  integer info; 
  VALUE rb_ab_out__;
  complex *ab_out__;

  integer ldab;
  integer n;
  integer ldz;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  w, z, work, rwork, iwork, info, ab = NumRu::Lapack.chbevd( jobz, uplo, kd, ab, lwork, lrwork, liwork)\n    or\n  NumRu::Lapack.chbevd  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 7)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 7)", argc);
  rb_jobz = argv[0];
  rb_uplo = argv[1];
  rb_kd = argv[2];
  rb_ab = argv[3];
  rb_lwork = argv[4];
  rb_lrwork = argv[5];
  rb_liwork = argv[6];

  if (!NA_IsNArray(rb_ab))
    rb_raise(rb_eArgError, "ab (4th argument) must be NArray");
  if (NA_RANK(rb_ab) != 2)
    rb_raise(rb_eArgError, "rank of ab (4th argument) must be %d", 2);
  n = NA_SHAPE1(rb_ab);
  ldab = NA_SHAPE0(rb_ab);
  if (NA_TYPE(rb_ab) != NA_SCOMPLEX)
    rb_ab = na_change_type(rb_ab, NA_SCOMPLEX);
  ab = NA_PTR_TYPE(rb_ab, complex*);
  uplo = StringValueCStr(rb_uplo)[0];
  liwork = NUM2INT(rb_liwork);
  lrwork = NUM2INT(rb_lrwork);
  jobz = StringValueCStr(rb_jobz)[0];
  lwork = NUM2INT(rb_lwork);
  kd = NUM2INT(rb_kd);
  ldz = lsame_(&jobz,"V") ? MAX(1,n) : 1;
  {
    int shape[1];
    shape[0] = n;
    rb_w = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  w = NA_PTR_TYPE(rb_w, real*);
  {
    int shape[2];
    shape[0] = ldz;
    shape[1] = n;
    rb_z = na_make_object(NA_SCOMPLEX, 2, shape, cNArray);
  }
  z = NA_PTR_TYPE(rb_z, complex*);
  {
    int shape[1];
    shape[0] = MAX(1,lwork);
    rb_work = na_make_object(NA_SCOMPLEX, 1, shape, cNArray);
  }
  work = NA_PTR_TYPE(rb_work, complex*);
  {
    int shape[1];
    shape[0] = MAX(1,lrwork);
    rb_rwork = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  rwork = NA_PTR_TYPE(rb_rwork, real*);
  {
    int shape[1];
    shape[0] = MAX(1,liwork);
    rb_iwork = na_make_object(NA_LINT, 1, shape, cNArray);
  }
  iwork = NA_PTR_TYPE(rb_iwork, integer*);
  {
    int shape[2];
    shape[0] = ldab;
    shape[1] = n;
    rb_ab_out__ = na_make_object(NA_SCOMPLEX, 2, shape, cNArray);
  }
  ab_out__ = NA_PTR_TYPE(rb_ab_out__, complex*);
  MEMCPY(ab_out__, ab, complex, NA_TOTAL(rb_ab));
  rb_ab = rb_ab_out__;
  ab = ab_out__;

  chbevd_(&jobz, &uplo, &n, &kd, ab, &ldab, w, z, &ldz, work, &lwork, rwork, &lrwork, iwork, &liwork, &info);

  rb_info = INT2NUM(info);
  return rb_ary_new3(7, rb_w, rb_z, rb_work, rb_rwork, rb_iwork, rb_info, rb_ab);
}

void
init_lapack_chbevd(VALUE mLapack){
  rb_define_module_function(mLapack, "chbevd", rb_chbevd, -1);
}

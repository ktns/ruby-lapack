#include "rb_lapack.h"

extern VOID chbevx_(char *jobz, char *range, char *uplo, integer *n, integer *kd, complex *ab, integer *ldab, complex *q, integer *ldq, real *vl, real *vu, integer *il, integer *iu, real *abstol, integer *m, real *w, complex *z, integer *ldz, complex *work, real *rwork, integer *iwork, integer *ifail, integer *info);

static VALUE
rb_chbevx(int argc, VALUE *argv, VALUE self){
  VALUE rb_jobz;
  char jobz; 
  VALUE rb_range;
  char range; 
  VALUE rb_uplo;
  char uplo; 
  VALUE rb_kd;
  integer kd; 
  VALUE rb_ab;
  complex *ab; 
  VALUE rb_vl;
  real vl; 
  VALUE rb_vu;
  real vu; 
  VALUE rb_il;
  integer il; 
  VALUE rb_iu;
  integer iu; 
  VALUE rb_abstol;
  real abstol; 
  VALUE rb_q;
  complex *q; 
  VALUE rb_m;
  integer m; 
  VALUE rb_w;
  real *w; 
  VALUE rb_z;
  complex *z; 
  VALUE rb_ifail;
  integer *ifail; 
  VALUE rb_info;
  integer info; 
  VALUE rb_ab_out__;
  complex *ab_out__;
  complex *work;
  real *rwork;
  integer *iwork;

  integer ldab;
  integer n;
  integer ldq;
  integer ldz;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  q, m, w, z, ifail, info, ab = NumRu::Lapack.chbevx( jobz, range, uplo, kd, ab, vl, vu, il, iu, abstol)\n    or\n  NumRu::Lapack.chbevx  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 10)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 10)", argc);
  rb_jobz = argv[0];
  rb_range = argv[1];
  rb_uplo = argv[2];
  rb_kd = argv[3];
  rb_ab = argv[4];
  rb_vl = argv[5];
  rb_vu = argv[6];
  rb_il = argv[7];
  rb_iu = argv[8];
  rb_abstol = argv[9];

  abstol = (real)NUM2DBL(rb_abstol);
  vl = (real)NUM2DBL(rb_vl);
  iu = NUM2INT(rb_iu);
  if (!NA_IsNArray(rb_ab))
    rb_raise(rb_eArgError, "ab (5th argument) must be NArray");
  if (NA_RANK(rb_ab) != 2)
    rb_raise(rb_eArgError, "rank of ab (5th argument) must be %d", 2);
  n = NA_SHAPE1(rb_ab);
  ldab = NA_SHAPE0(rb_ab);
  if (NA_TYPE(rb_ab) != NA_SCOMPLEX)
    rb_ab = na_change_type(rb_ab, NA_SCOMPLEX);
  ab = NA_PTR_TYPE(rb_ab, complex*);
  jobz = StringValueCStr(rb_jobz)[0];
  vu = (real)NUM2DBL(rb_vu);
  range = StringValueCStr(rb_range)[0];
  il = NUM2INT(rb_il);
  kd = NUM2INT(rb_kd);
  uplo = StringValueCStr(rb_uplo)[0];
  m = lsame_(&range,"A") ? n : lsame_(&range,"I") ? iu-il+1 : 0;
  ldz = lsame_(&jobz,"V") ? MAX(1,n) : 1;
  ldq = lsame_(&jobz,"V") ? MAX(1,n) : 0;
  {
    int shape[2];
    shape[0] = ldq;
    shape[1] = n;
    rb_q = na_make_object(NA_SCOMPLEX, 2, shape, cNArray);
  }
  q = NA_PTR_TYPE(rb_q, complex*);
  {
    int shape[1];
    shape[0] = n;
    rb_w = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  w = NA_PTR_TYPE(rb_w, real*);
  {
    int shape[2];
    shape[0] = ldz;
    shape[1] = MAX(1,m);
    rb_z = na_make_object(NA_SCOMPLEX, 2, shape, cNArray);
  }
  z = NA_PTR_TYPE(rb_z, complex*);
  {
    int shape[1];
    shape[0] = n;
    rb_ifail = na_make_object(NA_LINT, 1, shape, cNArray);
  }
  ifail = NA_PTR_TYPE(rb_ifail, integer*);
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
  work = ALLOC_N(complex, (n));
  rwork = ALLOC_N(real, (7*n));
  iwork = ALLOC_N(integer, (5*n));

  chbevx_(&jobz, &range, &uplo, &n, &kd, ab, &ldab, q, &ldq, &vl, &vu, &il, &iu, &abstol, &m, w, z, &ldz, work, rwork, iwork, ifail, &info);

  free(work);
  free(rwork);
  free(iwork);
  rb_m = INT2NUM(m);
  rb_info = INT2NUM(info);
  return rb_ary_new3(7, rb_q, rb_m, rb_w, rb_z, rb_ifail, rb_info, rb_ab);
}

void
init_lapack_chbevx(VALUE mLapack){
  rb_define_module_function(mLapack, "chbevx", rb_chbevx, -1);
}

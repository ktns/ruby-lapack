#include "rb_lapack.h"

extern VOID cgesdd_(char *jobz, integer *m, integer *n, complex *a, integer *lda, real *s, complex *u, integer *ldu, complex *vt, integer *ldvt, complex *work, integer *lwork, real *rwork, integer *iwork, integer *info);

static VALUE
rb_cgesdd(int argc, VALUE *argv, VALUE self){
  VALUE rb_jobz;
  char jobz; 
  VALUE rb_m;
  integer m; 
  VALUE rb_a;
  complex *a; 
  VALUE rb_lwork;
  integer lwork; 
  VALUE rb_s;
  real *s; 
  VALUE rb_u;
  complex *u; 
  VALUE rb_vt;
  complex *vt; 
  VALUE rb_work;
  complex *work; 
  VALUE rb_info;
  integer info; 
  VALUE rb_a_out__;
  complex *a_out__;
  real *rwork;
  integer *iwork;

  integer lda;
  integer n;
  integer ldu;
  integer ucol;
  integer ldvt;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  s, u, vt, work, info, a = NumRu::Lapack.cgesdd( jobz, m, a, lwork)\n    or\n  NumRu::Lapack.cgesdd  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 4)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 4)", argc);
  rb_jobz = argv[0];
  rb_m = argv[1];
  rb_a = argv[2];
  rb_lwork = argv[3];

  if (!NA_IsNArray(rb_a))
    rb_raise(rb_eArgError, "a (3th argument) must be NArray");
  if (NA_RANK(rb_a) != 2)
    rb_raise(rb_eArgError, "rank of a (3th argument) must be %d", 2);
  n = NA_SHAPE1(rb_a);
  lda = NA_SHAPE0(rb_a);
  if (NA_TYPE(rb_a) != NA_SCOMPLEX)
    rb_a = na_change_type(rb_a, NA_SCOMPLEX);
  a = NA_PTR_TYPE(rb_a, complex*);
  lwork = NUM2INT(rb_lwork);
  m = NUM2INT(rb_m);
  jobz = StringValueCStr(rb_jobz)[0];
  ucol = ((lsame_(&jobz,"A")) || (((lsame_(&jobz,"O")) && (m < n)))) ? m : lsame_(&jobz,"S") ? MIN(m,n) : 0;
  ldu = ((lsame_(&jobz,"S")) || ((('a') || (((lsame_(&jobz,"O")) && (m < n)))))) ? m : 1;
  ldvt = ((lsame_(&jobz,"A")) || (((lsame_(&jobz,"O")) && (m == n)))) ? n : lsame_(&jobz,"S") ? MIN(m,n) : 1;
  {
    int shape[1];
    shape[0] = MIN(m,n);
    rb_s = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  s = NA_PTR_TYPE(rb_s, real*);
  {
    int shape[2];
    shape[0] = ldu;
    shape[1] = ucol;
    rb_u = na_make_object(NA_SCOMPLEX, 2, shape, cNArray);
  }
  u = NA_PTR_TYPE(rb_u, complex*);
  {
    int shape[2];
    shape[0] = ldvt;
    shape[1] = n;
    rb_vt = na_make_object(NA_SCOMPLEX, 2, shape, cNArray);
  }
  vt = NA_PTR_TYPE(rb_vt, complex*);
  {
    int shape[1];
    shape[0] = MAX(1,lwork);
    rb_work = na_make_object(NA_SCOMPLEX, 1, shape, cNArray);
  }
  work = NA_PTR_TYPE(rb_work, complex*);
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
  rwork = ALLOC_N(real, (MAX(1, lsame_(&jobz,"N") ? 5*MIN(m,n) : MIN(m,n)*MAX(5*MIN(m,n)+7,2*MAX(m,n)+2*MIN(m,n)+1))));
  iwork = ALLOC_N(integer, (8*MIN(m,n)));

  cgesdd_(&jobz, &m, &n, a, &lda, s, u, &ldu, vt, &ldvt, work, &lwork, rwork, iwork, &info);

  free(rwork);
  free(iwork);
  rb_info = INT2NUM(info);
  return rb_ary_new3(6, rb_s, rb_u, rb_vt, rb_work, rb_info, rb_a);
}

void
init_lapack_cgesdd(VALUE mLapack){
  rb_define_module_function(mLapack, "cgesdd", rb_cgesdd, -1);
}

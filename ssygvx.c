#include "rb_lapack.h"

extern VOID ssygvx_(integer *itype, char *jobz, char *range, char *uplo, integer *n, real *a, integer *lda, real *b, integer *ldb, real *vl, real *vu, integer *il, integer *iu, real *abstol, integer *m, real *w, real *z, integer *ldz, real *work, integer *lwork, integer *iwork, integer *ifail, integer *info);

static VALUE
rb_ssygvx(int argc, VALUE *argv, VALUE self){
  VALUE rb_itype;
  integer itype; 
  VALUE rb_jobz;
  char jobz; 
  VALUE rb_range;
  char range; 
  VALUE rb_uplo;
  char uplo; 
  VALUE rb_a;
  real *a; 
  VALUE rb_b;
  real *b; 
  VALUE rb_ldb;
  integer ldb; 
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
  VALUE rb_ldz;
  integer ldz; 
  VALUE rb_lwork;
  integer lwork; 
  VALUE rb_m;
  integer m; 
  VALUE rb_w;
  real *w; 
  VALUE rb_z;
  real *z; 
  VALUE rb_work;
  real *work; 
  VALUE rb_ifail;
  integer *ifail; 
  VALUE rb_info;
  integer info; 
  VALUE rb_a_out__;
  real *a_out__;
  VALUE rb_b_out__;
  real *b_out__;
  integer *iwork;

  integer lda;
  integer n;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  m, w, z, work, ifail, info, a, b = NumRu::Lapack.ssygvx( itype, jobz, range, uplo, a, b, ldb, vl, vu, il, iu, abstol, ldz, lwork)\n    or\n  NumRu::Lapack.ssygvx  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 14)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 14)", argc);
  rb_itype = argv[0];
  rb_jobz = argv[1];
  rb_range = argv[2];
  rb_uplo = argv[3];
  rb_a = argv[4];
  rb_b = argv[5];
  rb_ldb = argv[6];
  rb_vl = argv[7];
  rb_vu = argv[8];
  rb_il = argv[9];
  rb_iu = argv[10];
  rb_abstol = argv[11];
  rb_ldz = argv[12];
  rb_lwork = argv[13];

  abstol = (real)NUM2DBL(rb_abstol);
  vl = (real)NUM2DBL(rb_vl);
  iu = NUM2INT(rb_iu);
  if (!NA_IsNArray(rb_a))
    rb_raise(rb_eArgError, "a (5th argument) must be NArray");
  if (NA_RANK(rb_a) != 2)
    rb_raise(rb_eArgError, "rank of a (5th argument) must be %d", 2);
  n = NA_SHAPE1(rb_a);
  lda = NA_SHAPE0(rb_a);
  if (NA_TYPE(rb_a) != NA_SFLOAT)
    rb_a = na_change_type(rb_a, NA_SFLOAT);
  a = NA_PTR_TYPE(rb_a, real*);
  uplo = StringValueCStr(rb_uplo)[0];
  ldb = NUM2INT(rb_ldb);
  if (!NA_IsNArray(rb_b))
    rb_raise(rb_eArgError, "b (6th argument) must be NArray");
  if (NA_RANK(rb_b) != 2)
    rb_raise(rb_eArgError, "rank of b (6th argument) must be %d", 2);
  if (NA_SHAPE1(rb_b) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of b must be the same as shape 1 of a");
  if (NA_SHAPE0(rb_b) != lda)
    rb_raise(rb_eRuntimeError, "shape 0 of b must be the same as shape 0 of a");
  if (NA_TYPE(rb_b) != NA_SFLOAT)
    rb_b = na_change_type(rb_b, NA_SFLOAT);
  b = NA_PTR_TYPE(rb_b, real*);
  jobz = StringValueCStr(rb_jobz)[0];
  itype = NUM2INT(rb_itype);
  lwork = NUM2INT(rb_lwork);
  range = StringValueCStr(rb_range)[0];
  vu = (real)NUM2DBL(rb_vu);
  il = NUM2INT(rb_il);
  ldz = lsame_(&jobz,"V") ? MAX(1,n) : 1;
  m = lsame_(&range,"A") ? n : lsame_(&range,"I") ? iu-il+1 : 0;
  {
    int shape[1];
    shape[0] = n;
    rb_w = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  w = NA_PTR_TYPE(rb_w, real*);
  {
    int shape[2];
    shape[0] = lsame_(&jobz,"N") ? 0 : ldz;
    shape[1] = lsame_(&jobz,"N") ? 0 : MAX(1,m);
    rb_z = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  z = NA_PTR_TYPE(rb_z, real*);
  {
    int shape[1];
    shape[0] = MAX(1,lwork);
    rb_work = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  work = NA_PTR_TYPE(rb_work, real*);
  {
    int shape[1];
    shape[0] = n;
    rb_ifail = na_make_object(NA_LINT, 1, shape, cNArray);
  }
  ifail = NA_PTR_TYPE(rb_ifail, integer*);
  {
    int shape[2];
    shape[0] = lda;
    shape[1] = n;
    rb_a_out__ = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  a_out__ = NA_PTR_TYPE(rb_a_out__, real*);
  MEMCPY(a_out__, a, real, NA_TOTAL(rb_a));
  rb_a = rb_a_out__;
  a = a_out__;
  {
    int shape[2];
    shape[0] = lda;
    shape[1] = n;
    rb_b_out__ = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  b_out__ = NA_PTR_TYPE(rb_b_out__, real*);
  MEMCPY(b_out__, b, real, NA_TOTAL(rb_b));
  rb_b = rb_b_out__;
  b = b_out__;
  iwork = ALLOC_N(integer, (5*n));

  ssygvx_(&itype, &jobz, &range, &uplo, &n, a, &lda, b, &ldb, &vl, &vu, &il, &iu, &abstol, &m, w, z, &ldz, work, &lwork, iwork, ifail, &info);

  free(iwork);
  rb_m = INT2NUM(m);
  rb_info = INT2NUM(info);
  return rb_ary_new3(8, rb_m, rb_w, rb_z, rb_work, rb_ifail, rb_info, rb_a, rb_b);
}

void
init_lapack_ssygvx(VALUE mLapack){
  rb_define_module_function(mLapack, "ssygvx", rb_ssygvx, -1);
}

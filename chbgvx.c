#include "rb_lapack.h"

extern VOID chbgvx_(char *jobz, char *range, char *uplo, integer *n, integer *ka, integer *kb, complex *ab, integer *ldab, complex *bb, integer *ldbb, complex *q, integer *ldq, real *vl, real *vu, integer *il, integer *iu, real *abstol, integer *m, real *w, complex *z, integer *ldz, complex *work, real *rwork, integer *iwork, integer *ifail, integer *info);

static VALUE
rb_chbgvx(int argc, VALUE *argv, VALUE self){
  VALUE rb_jobz;
  char jobz; 
  VALUE rb_range;
  char range; 
  VALUE rb_uplo;
  char uplo; 
  VALUE rb_ka;
  integer ka; 
  VALUE rb_kb;
  integer kb; 
  VALUE rb_ab;
  complex *ab; 
  VALUE rb_bb;
  complex *bb; 
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
  VALUE rb_bb_out__;
  complex *bb_out__;
  complex *work;
  real *rwork;
  integer *iwork;

  integer ldab;
  integer n;
  integer ldbb;
  integer ldq;
  integer ldz;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  q, m, w, z, ifail, info, ab, bb = NumRu::Lapack.chbgvx( jobz, range, uplo, ka, kb, ab, bb, vl, vu, il, iu, abstol)\n    or\n  NumRu::Lapack.chbgvx  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 12)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 12)", argc);
  rb_jobz = argv[0];
  rb_range = argv[1];
  rb_uplo = argv[2];
  rb_ka = argv[3];
  rb_kb = argv[4];
  rb_ab = argv[5];
  rb_bb = argv[6];
  rb_vl = argv[7];
  rb_vu = argv[8];
  rb_il = argv[9];
  rb_iu = argv[10];
  rb_abstol = argv[11];

  abstol = (real)NUM2DBL(rb_abstol);
  vl = (real)NUM2DBL(rb_vl);
  if (!NA_IsNArray(rb_bb))
    rb_raise(rb_eArgError, "bb (7th argument) must be NArray");
  if (NA_RANK(rb_bb) != 2)
    rb_raise(rb_eArgError, "rank of bb (7th argument) must be %d", 2);
  n = NA_SHAPE1(rb_bb);
  ldbb = NA_SHAPE0(rb_bb);
  if (NA_TYPE(rb_bb) != NA_SCOMPLEX)
    rb_bb = na_change_type(rb_bb, NA_SCOMPLEX);
  bb = NA_PTR_TYPE(rb_bb, complex*);
  ka = NUM2INT(rb_ka);
  iu = NUM2INT(rb_iu);
  if (!NA_IsNArray(rb_ab))
    rb_raise(rb_eArgError, "ab (6th argument) must be NArray");
  if (NA_RANK(rb_ab) != 2)
    rb_raise(rb_eArgError, "rank of ab (6th argument) must be %d", 2);
  if (NA_SHAPE1(rb_ab) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of ab must be the same as shape 1 of bb");
  ldab = NA_SHAPE0(rb_ab);
  if (NA_TYPE(rb_ab) != NA_SCOMPLEX)
    rb_ab = na_change_type(rb_ab, NA_SCOMPLEX);
  ab = NA_PTR_TYPE(rb_ab, complex*);
  kb = NUM2INT(rb_kb);
  vu = (real)NUM2DBL(rb_vu);
  jobz = StringValueCStr(rb_jobz)[0];
  il = NUM2INT(rb_il);
  range = StringValueCStr(rb_range)[0];
  uplo = StringValueCStr(rb_uplo)[0];
  ldq = 1 ? jobz = 'n' : max(1,n) ? jobz = 'v' : 0;
  ldz = lsame_(&jobz,"V") ? n : 1;
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
    shape[1] = n;
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
  {
    int shape[2];
    shape[0] = ldbb;
    shape[1] = n;
    rb_bb_out__ = na_make_object(NA_SCOMPLEX, 2, shape, cNArray);
  }
  bb_out__ = NA_PTR_TYPE(rb_bb_out__, complex*);
  MEMCPY(bb_out__, bb, complex, NA_TOTAL(rb_bb));
  rb_bb = rb_bb_out__;
  bb = bb_out__;
  work = ALLOC_N(complex, (n));
  rwork = ALLOC_N(real, (7*n));
  iwork = ALLOC_N(integer, (5*n));

  chbgvx_(&jobz, &range, &uplo, &n, &ka, &kb, ab, &ldab, bb, &ldbb, q, &ldq, &vl, &vu, &il, &iu, &abstol, &m, w, z, &ldz, work, rwork, iwork, ifail, &info);

  free(work);
  free(rwork);
  free(iwork);
  rb_m = INT2NUM(m);
  rb_info = INT2NUM(info);
  return rb_ary_new3(8, rb_q, rb_m, rb_w, rb_z, rb_ifail, rb_info, rb_ab, rb_bb);
}

void
init_lapack_chbgvx(VALUE mLapack){
  rb_define_module_function(mLapack, "chbgvx", rb_chbgvx, -1);
}

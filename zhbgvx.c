#include "rb_lapack.h"

extern VOID zhbgvx_(char *jobz, char *range, char *uplo, integer *n, integer *ka, integer *kb, doublecomplex *ab, integer *ldab, doublecomplex *bb, integer *ldbb, doublecomplex *q, integer *ldq, doublereal *vl, doublereal *vu, integer *il, integer *iu, doublereal *abstol, integer *m, doublereal *w, doublecomplex *z, integer *ldz, doublecomplex *work, doublereal *rwork, integer *iwork, integer *ifail, integer *info);

static VALUE
rb_zhbgvx(int argc, VALUE *argv, VALUE self){
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
  doublecomplex *ab; 
  VALUE rb_bb;
  doublecomplex *bb; 
  VALUE rb_vl;
  doublereal vl; 
  VALUE rb_vu;
  doublereal vu; 
  VALUE rb_il;
  integer il; 
  VALUE rb_iu;
  integer iu; 
  VALUE rb_abstol;
  doublereal abstol; 
  VALUE rb_q;
  doublecomplex *q; 
  VALUE rb_m;
  integer m; 
  VALUE rb_w;
  doublereal *w; 
  VALUE rb_z;
  doublecomplex *z; 
  VALUE rb_ifail;
  integer *ifail; 
  VALUE rb_info;
  integer info; 
  VALUE rb_ab_out__;
  doublecomplex *ab_out__;
  VALUE rb_bb_out__;
  doublecomplex *bb_out__;
  doublecomplex *work;
  doublereal *rwork;
  integer *iwork;

  integer ldab;
  integer n;
  integer ldbb;
  integer ldq;
  integer ldz;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  q, m, w, z, ifail, info, ab, bb = NumRu::Lapack.zhbgvx( jobz, range, uplo, ka, kb, ab, bb, vl, vu, il, iu, abstol)\n    or\n  NumRu::Lapack.zhbgvx  # print help\n\n\nFORTRAN MANUAL\n\n");
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

  abstol = NUM2DBL(rb_abstol);
  vl = NUM2DBL(rb_vl);
  if (!NA_IsNArray(rb_bb))
    rb_raise(rb_eArgError, "bb (7th argument) must be NArray");
  if (NA_RANK(rb_bb) != 2)
    rb_raise(rb_eArgError, "rank of bb (7th argument) must be %d", 2);
  n = NA_SHAPE1(rb_bb);
  ldbb = NA_SHAPE0(rb_bb);
  if (NA_TYPE(rb_bb) != NA_DCOMPLEX)
    rb_bb = na_change_type(rb_bb, NA_DCOMPLEX);
  bb = NA_PTR_TYPE(rb_bb, doublecomplex*);
  ka = NUM2INT(rb_ka);
  iu = NUM2INT(rb_iu);
  if (!NA_IsNArray(rb_ab))
    rb_raise(rb_eArgError, "ab (6th argument) must be NArray");
  if (NA_RANK(rb_ab) != 2)
    rb_raise(rb_eArgError, "rank of ab (6th argument) must be %d", 2);
  if (NA_SHAPE1(rb_ab) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of ab must be the same as shape 1 of bb");
  ldab = NA_SHAPE0(rb_ab);
  if (NA_TYPE(rb_ab) != NA_DCOMPLEX)
    rb_ab = na_change_type(rb_ab, NA_DCOMPLEX);
  ab = NA_PTR_TYPE(rb_ab, doublecomplex*);
  kb = NUM2INT(rb_kb);
  vu = NUM2DBL(rb_vu);
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
    rb_q = na_make_object(NA_DCOMPLEX, 2, shape, cNArray);
  }
  q = NA_PTR_TYPE(rb_q, doublecomplex*);
  {
    int shape[1];
    shape[0] = n;
    rb_w = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  w = NA_PTR_TYPE(rb_w, doublereal*);
  {
    int shape[2];
    shape[0] = ldz;
    shape[1] = n;
    rb_z = na_make_object(NA_DCOMPLEX, 2, shape, cNArray);
  }
  z = NA_PTR_TYPE(rb_z, doublecomplex*);
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
    rb_ab_out__ = na_make_object(NA_DCOMPLEX, 2, shape, cNArray);
  }
  ab_out__ = NA_PTR_TYPE(rb_ab_out__, doublecomplex*);
  MEMCPY(ab_out__, ab, doublecomplex, NA_TOTAL(rb_ab));
  rb_ab = rb_ab_out__;
  ab = ab_out__;
  {
    int shape[2];
    shape[0] = ldbb;
    shape[1] = n;
    rb_bb_out__ = na_make_object(NA_DCOMPLEX, 2, shape, cNArray);
  }
  bb_out__ = NA_PTR_TYPE(rb_bb_out__, doublecomplex*);
  MEMCPY(bb_out__, bb, doublecomplex, NA_TOTAL(rb_bb));
  rb_bb = rb_bb_out__;
  bb = bb_out__;
  work = ALLOC_N(doublecomplex, (n));
  rwork = ALLOC_N(doublereal, (7*n));
  iwork = ALLOC_N(integer, (5*n));

  zhbgvx_(&jobz, &range, &uplo, &n, &ka, &kb, ab, &ldab, bb, &ldbb, q, &ldq, &vl, &vu, &il, &iu, &abstol, &m, w, z, &ldz, work, rwork, iwork, ifail, &info);

  free(work);
  free(rwork);
  free(iwork);
  rb_m = INT2NUM(m);
  rb_info = INT2NUM(info);
  return rb_ary_new3(8, rb_q, rb_m, rb_w, rb_z, rb_ifail, rb_info, rb_ab, rb_bb);
}

void
init_lapack_zhbgvx(VALUE mLapack){
  rb_define_module_function(mLapack, "zhbgvx", rb_zhbgvx, -1);
}

#include "rb_lapack.h"

extern VOID dsbgv_(char *jobz, char *uplo, integer *n, integer *ka, integer *kb, doublereal *ab, integer *ldab, doublereal *bb, integer *ldbb, doublereal *w, doublereal *z, integer *ldz, doublereal *work, integer *info);

static VALUE
rb_dsbgv(int argc, VALUE *argv, VALUE self){
  VALUE rb_jobz;
  char jobz; 
  VALUE rb_uplo;
  char uplo; 
  VALUE rb_ka;
  integer ka; 
  VALUE rb_kb;
  integer kb; 
  VALUE rb_ab;
  doublereal *ab; 
  VALUE rb_bb;
  doublereal *bb; 
  VALUE rb_w;
  doublereal *w; 
  VALUE rb_z;
  doublereal *z; 
  VALUE rb_info;
  integer info; 
  VALUE rb_ab_out__;
  doublereal *ab_out__;
  VALUE rb_bb_out__;
  doublereal *bb_out__;
  doublereal *work;

  integer ldab;
  integer n;
  integer ldbb;
  integer ldz;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  w, z, info, ab, bb = NumRu::Lapack.dsbgv( jobz, uplo, ka, kb, ab, bb)\n    or\n  NumRu::Lapack.dsbgv  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 6)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 6)", argc);
  rb_jobz = argv[0];
  rb_uplo = argv[1];
  rb_ka = argv[2];
  rb_kb = argv[3];
  rb_ab = argv[4];
  rb_bb = argv[5];

  if (!NA_IsNArray(rb_bb))
    rb_raise(rb_eArgError, "bb (6th argument) must be NArray");
  if (NA_RANK(rb_bb) != 2)
    rb_raise(rb_eArgError, "rank of bb (6th argument) must be %d", 2);
  n = NA_SHAPE1(rb_bb);
  ldbb = NA_SHAPE0(rb_bb);
  if (NA_TYPE(rb_bb) != NA_DFLOAT)
    rb_bb = na_change_type(rb_bb, NA_DFLOAT);
  bb = NA_PTR_TYPE(rb_bb, doublereal*);
  ka = NUM2INT(rb_ka);
  if (!NA_IsNArray(rb_ab))
    rb_raise(rb_eArgError, "ab (5th argument) must be NArray");
  if (NA_RANK(rb_ab) != 2)
    rb_raise(rb_eArgError, "rank of ab (5th argument) must be %d", 2);
  if (NA_SHAPE1(rb_ab) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of ab must be the same as shape 1 of bb");
  ldab = NA_SHAPE0(rb_ab);
  if (NA_TYPE(rb_ab) != NA_DFLOAT)
    rb_ab = na_change_type(rb_ab, NA_DFLOAT);
  ab = NA_PTR_TYPE(rb_ab, doublereal*);
  jobz = StringValueCStr(rb_jobz)[0];
  uplo = StringValueCStr(rb_uplo)[0];
  kb = NUM2INT(rb_kb);
  ldz = lsame_(&jobz,"V") ? n : 1;
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
    rb_z = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  z = NA_PTR_TYPE(rb_z, doublereal*);
  {
    int shape[2];
    shape[0] = ldab;
    shape[1] = n;
    rb_ab_out__ = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  ab_out__ = NA_PTR_TYPE(rb_ab_out__, doublereal*);
  MEMCPY(ab_out__, ab, doublereal, NA_TOTAL(rb_ab));
  rb_ab = rb_ab_out__;
  ab = ab_out__;
  {
    int shape[2];
    shape[0] = ldbb;
    shape[1] = n;
    rb_bb_out__ = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  bb_out__ = NA_PTR_TYPE(rb_bb_out__, doublereal*);
  MEMCPY(bb_out__, bb, doublereal, NA_TOTAL(rb_bb));
  rb_bb = rb_bb_out__;
  bb = bb_out__;
  work = ALLOC_N(doublereal, (3*n));

  dsbgv_(&jobz, &uplo, &n, &ka, &kb, ab, &ldab, bb, &ldbb, w, z, &ldz, work, &info);

  free(work);
  rb_info = INT2NUM(info);
  return rb_ary_new3(5, rb_w, rb_z, rb_info, rb_ab, rb_bb);
}

void
init_lapack_dsbgv(VALUE mLapack){
  rb_define_module_function(mLapack, "dsbgv", rb_dsbgv, -1);
}

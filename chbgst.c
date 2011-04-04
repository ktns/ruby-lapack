#include "rb_lapack.h"

extern VOID chbgst_(char *vect, char *uplo, integer *n, integer *ka, integer *kb, complex *ab, integer *ldab, complex *bb, integer *ldbb, complex *x, integer *ldx, complex *work, real *rwork, integer *info);

static VALUE
rb_chbgst(int argc, VALUE *argv, VALUE self){
  VALUE rb_vect;
  char vect; 
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
  VALUE rb_x;
  complex *x; 
  VALUE rb_info;
  integer info; 
  VALUE rb_ab_out__;
  complex *ab_out__;
  complex *work;
  real *rwork;

  integer ldab;
  integer n;
  integer ldbb;
  integer ldx;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  x, info, ab = NumRu::Lapack.chbgst( vect, uplo, ka, kb, ab, bb)\n    or\n  NumRu::Lapack.chbgst  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 6)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 6)", argc);
  rb_vect = argv[0];
  rb_uplo = argv[1];
  rb_ka = argv[2];
  rb_kb = argv[3];
  rb_ab = argv[4];
  rb_bb = argv[5];

  uplo = StringValueCStr(rb_uplo)[0];
  if (!NA_IsNArray(rb_bb))
    rb_raise(rb_eArgError, "bb (6th argument) must be NArray");
  if (NA_RANK(rb_bb) != 2)
    rb_raise(rb_eArgError, "rank of bb (6th argument) must be %d", 2);
  n = NA_SHAPE1(rb_bb);
  ldbb = NA_SHAPE0(rb_bb);
  if (NA_TYPE(rb_bb) != NA_SCOMPLEX)
    rb_bb = na_change_type(rb_bb, NA_SCOMPLEX);
  bb = NA_PTR_TYPE(rb_bb, complex*);
  ka = NUM2INT(rb_ka);
  vect = StringValueCStr(rb_vect)[0];
  kb = NUM2INT(rb_kb);
  if (!NA_IsNArray(rb_ab))
    rb_raise(rb_eArgError, "ab (5th argument) must be NArray");
  if (NA_RANK(rb_ab) != 2)
    rb_raise(rb_eArgError, "rank of ab (5th argument) must be %d", 2);
  if (NA_SHAPE1(rb_ab) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of ab must be the same as shape 1 of bb");
  ldab = NA_SHAPE0(rb_ab);
  if (NA_TYPE(rb_ab) != NA_SCOMPLEX)
    rb_ab = na_change_type(rb_ab, NA_SCOMPLEX);
  ab = NA_PTR_TYPE(rb_ab, complex*);
  ldx = lsame_(&vect,"V") ? MAX(1,n) : 1;
  {
    int shape[2];
    shape[0] = ldx;
    shape[1] = n;
    rb_x = na_make_object(NA_SCOMPLEX, 2, shape, cNArray);
  }
  x = NA_PTR_TYPE(rb_x, complex*);
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
  rwork = ALLOC_N(real, (n));

  chbgst_(&vect, &uplo, &n, &ka, &kb, ab, &ldab, bb, &ldbb, x, &ldx, work, rwork, &info);

  free(work);
  free(rwork);
  rb_info = INT2NUM(info);
  return rb_ary_new3(3, rb_x, rb_info, rb_ab);
}

void
init_lapack_chbgst(VALUE mLapack){
  rb_define_module_function(mLapack, "chbgst", rb_chbgst, -1);
}

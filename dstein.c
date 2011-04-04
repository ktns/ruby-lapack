#include "rb_lapack.h"

extern VOID dstein_(integer *n, doublereal *d, doublereal *e, integer *m, doublereal *w, integer *iblock, integer *isplit, doublereal *z, integer *ldz, doublereal *work, integer *iwork, integer *ifail, integer *info);

static VALUE
rb_dstein(int argc, VALUE *argv, VALUE self){
  VALUE rb_d;
  doublereal *d; 
  VALUE rb_e;
  doublereal *e; 
  VALUE rb_w;
  doublereal *w; 
  VALUE rb_iblock;
  integer *iblock; 
  VALUE rb_isplit;
  integer *isplit; 
  VALUE rb_z;
  doublereal *z; 
  VALUE rb_ifail;
  integer *ifail; 
  VALUE rb_info;
  integer info; 
  doublereal *work;
  integer *iwork;

  integer n;
  integer ldz;
  integer m;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  z, ifail, info = NumRu::Lapack.dstein( d, e, w, iblock, isplit)\n    or\n  NumRu::Lapack.dstein  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 5)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 5)", argc);
  rb_d = argv[0];
  rb_e = argv[1];
  rb_w = argv[2];
  rb_iblock = argv[3];
  rb_isplit = argv[4];

  if (!NA_IsNArray(rb_w))
    rb_raise(rb_eArgError, "w (3th argument) must be NArray");
  if (NA_RANK(rb_w) != 1)
    rb_raise(rb_eArgError, "rank of w (3th argument) must be %d", 1);
  n = NA_SHAPE0(rb_w);
  if (NA_TYPE(rb_w) != NA_DFLOAT)
    rb_w = na_change_type(rb_w, NA_DFLOAT);
  w = NA_PTR_TYPE(rb_w, doublereal*);
  if (!NA_IsNArray(rb_iblock))
    rb_raise(rb_eArgError, "iblock (4th argument) must be NArray");
  if (NA_RANK(rb_iblock) != 1)
    rb_raise(rb_eArgError, "rank of iblock (4th argument) must be %d", 1);
  if (NA_SHAPE0(rb_iblock) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of iblock must be the same as shape 0 of w");
  if (NA_TYPE(rb_iblock) != NA_LINT)
    rb_iblock = na_change_type(rb_iblock, NA_LINT);
  iblock = NA_PTR_TYPE(rb_iblock, integer*);
  if (!NA_IsNArray(rb_isplit))
    rb_raise(rb_eArgError, "isplit (5th argument) must be NArray");
  if (NA_RANK(rb_isplit) != 1)
    rb_raise(rb_eArgError, "rank of isplit (5th argument) must be %d", 1);
  if (NA_SHAPE0(rb_isplit) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of isplit must be the same as shape 0 of w");
  if (NA_TYPE(rb_isplit) != NA_LINT)
    rb_isplit = na_change_type(rb_isplit, NA_LINT);
  isplit = NA_PTR_TYPE(rb_isplit, integer*);
  if (!NA_IsNArray(rb_d))
    rb_raise(rb_eArgError, "d (1th argument) must be NArray");
  if (NA_RANK(rb_d) != 1)
    rb_raise(rb_eArgError, "rank of d (1th argument) must be %d", 1);
  if (NA_SHAPE0(rb_d) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of d must be the same as shape 0 of w");
  if (NA_TYPE(rb_d) != NA_DFLOAT)
    rb_d = na_change_type(rb_d, NA_DFLOAT);
  d = NA_PTR_TYPE(rb_d, doublereal*);
  if (!NA_IsNArray(rb_e))
    rb_raise(rb_eArgError, "e (2th argument) must be NArray");
  if (NA_RANK(rb_e) != 1)
    rb_raise(rb_eArgError, "rank of e (2th argument) must be %d", 1);
  if (NA_SHAPE0(rb_e) != (n-1))
    rb_raise(rb_eRuntimeError, "shape 0 of e must be %d", n-1);
  if (NA_TYPE(rb_e) != NA_DFLOAT)
    rb_e = na_change_type(rb_e, NA_DFLOAT);
  e = NA_PTR_TYPE(rb_e, doublereal*);
  m = n;
  ldz = MAX(1,n);
  {
    int shape[2];
    shape[0] = ldz;
    shape[1] = m;
    rb_z = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  z = NA_PTR_TYPE(rb_z, doublereal*);
  {
    int shape[1];
    shape[0] = m;
    rb_ifail = na_make_object(NA_LINT, 1, shape, cNArray);
  }
  ifail = NA_PTR_TYPE(rb_ifail, integer*);
  work = ALLOC_N(doublereal, (5*n));
  iwork = ALLOC_N(integer, (n));

  dstein_(&n, d, e, &m, w, iblock, isplit, z, &ldz, work, iwork, ifail, &info);

  free(work);
  free(iwork);
  rb_info = INT2NUM(info);
  return rb_ary_new3(3, rb_z, rb_ifail, rb_info);
}

void
init_lapack_dstein(VALUE mLapack){
  rb_define_module_function(mLapack, "dstein", rb_dstein, -1);
}

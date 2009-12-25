#include "rb_lapack.h"

static VALUE
rb_ilaslc(int argc, VALUE *argv, VALUE self){
  VALUE rb_m;
  integer m; 
  VALUE rb_a;
  real *a; 
  VALUE rb___out__;
  integer __out__; 

  integer lda;
  integer n;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  __out__ = NumRu::Lapack.ilaslc( m, a)\n    or\n  NumRu::Lapack.ilaslc  # print help\n\n\nFORTRAN MANUAL\n      INTEGER FUNCTION ILASLC(M, N, A, LDA)\n\n*  Purpose\n*  =======\n*\n*  ILASLC scans A for its last non-zero column.\n*\n\n*  Arguments\n*  =========\n*\n*  M       (input) INTEGER\n*          The number of rows of the matrix A.\n*\n*  N       (input) INTEGER\n*          The number of columns of the matrix A.\n*\n*  A       (input) REAL array, dimension (LDA,N)\n*          The m by n matrix A.\n*\n*  LDA     (input) INTEGER\n*          The leading dimension of the array A. LDA >= max(1,M).\n*\n\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 2)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 2)", argc);
  rb_m = argv[0];
  rb_a = argv[1];

  m = NUM2INT(rb_m);
  if (!NA_IsNArray(rb_a))
    rb_raise(rb_eArgError, "a (2th argument) must be NArray");
  if (NA_RANK(rb_a) != 2)
    rb_raise(rb_eArgError, "rank of a (2th argument) must be %d", 2);
  lda = NA_SHAPE0(rb_a);
  n = NA_SHAPE1(rb_a);
  if (NA_TYPE(rb_a) != NA_SFLOAT)
    rb_a = na_change_type(rb_a, NA_SFLOAT);
  a = NA_PTR_TYPE(rb_a, real*);

  __out__ = ilaslc_(&m, &n, a, &lda);

  rb___out__ = INT2NUM(__out__);
  return rb___out__;
}

void
init_lapack_ilaslc(VALUE mLapack){
  rb_define_module_function(mLapack, "ilaslc", rb_ilaslc, -1);
}

#include "rb_lapack.h"

static VALUE
rb_cpttrf(int argc, VALUE *argv, VALUE self){
  VALUE rb_d;
  real *d; 
  VALUE rb_e;
  complex *e; 
  VALUE rb_info;
  integer info; 
  VALUE rb_d_out__;
  real *d_out__;
  VALUE rb_e_out__;
  complex *e_out__;

  integer n;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  info, d, e = NumRu::Lapack.cpttrf( d, e)\n    or\n  NumRu::Lapack.cpttrf  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE CPTTRF( N, D, E, INFO )\n\n*  Purpose\n*  =======\n*\n*  CPTTRF computes the L*D*L' factorization of a complex Hermitian\n*  positive definite tridiagonal matrix A.  The factorization may also\n*  be regarded as having the form A = U'*D*U.\n*\n\n*  Arguments\n*  =========\n*\n*  N       (input) INTEGER\n*          The order of the matrix A.  N >= 0.\n*\n*  D       (input/output) REAL array, dimension (N)\n*          On entry, the n diagonal elements of the tridiagonal matrix\n*          A.  On exit, the n diagonal elements of the diagonal matrix\n*          D from the L*D*L' factorization of A.\n*\n*  E       (input/output) COMPLEX array, dimension (N-1)\n*          On entry, the (n-1) subdiagonal elements of the tridiagonal\n*          matrix A.  On exit, the (n-1) subdiagonal elements of the\n*          unit bidiagonal factor L from the L*D*L' factorization of A.\n*          E can also be regarded as the superdiagonal of the unit\n*          bidiagonal factor U from the U'*D*U factorization of A.\n*\n*  INFO    (output) INTEGER\n*          = 0: successful exit\n*          < 0: if INFO = -k, the k-th argument had an illegal value\n*          > 0: if INFO = k, the leading minor of order k is not\n*               positive definite; if k < N, the factorization could not\n*               be completed, while if k = N, the factorization was\n*               completed, but D(N) <= 0.\n*\n\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 2)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 2)", argc);
  rb_d = argv[0];
  rb_e = argv[1];

  if (!NA_IsNArray(rb_d))
    rb_raise(rb_eArgError, "d (1th argument) must be NArray");
  if (NA_RANK(rb_d) != 1)
    rb_raise(rb_eArgError, "rank of d (1th argument) must be %d", 1);
  n = NA_SHAPE0(rb_d);
  if (NA_TYPE(rb_d) != NA_SFLOAT)
    rb_d = na_change_type(rb_d, NA_SFLOAT);
  d = NA_PTR_TYPE(rb_d, real*);
  if (!NA_IsNArray(rb_e))
    rb_raise(rb_eArgError, "e (2th argument) must be NArray");
  if (NA_RANK(rb_e) != 1)
    rb_raise(rb_eArgError, "rank of e (2th argument) must be %d", 1);
  if (NA_SHAPE0(rb_e) != (n-1))
    rb_raise(rb_eRuntimeError, "shape 0 of e must be %d", n-1);
  if (NA_TYPE(rb_e) != NA_SCOMPLEX)
    rb_e = na_change_type(rb_e, NA_SCOMPLEX);
  e = NA_PTR_TYPE(rb_e, complex*);
  {
    int shape[1];
    shape[0] = n;
    rb_d_out__ = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  d_out__ = NA_PTR_TYPE(rb_d_out__, real*);
  MEMCPY(d_out__, d, real, NA_TOTAL(rb_d));
  rb_d = rb_d_out__;
  d = d_out__;
  {
    int shape[1];
    shape[0] = n-1;
    rb_e_out__ = na_make_object(NA_SCOMPLEX, 1, shape, cNArray);
  }
  e_out__ = NA_PTR_TYPE(rb_e_out__, complex*);
  MEMCPY(e_out__, e, complex, NA_TOTAL(rb_e));
  rb_e = rb_e_out__;
  e = e_out__;

  cpttrf_(&n, d, e, &info);

  rb_info = INT2NUM(info);
  return rb_ary_new3(3, rb_info, rb_d, rb_e);
}

void
init_lapack_cpttrf(VALUE mLapack){
  rb_define_module_function(mLapack, "cpttrf", rb_cpttrf, -1);
}

#include "rb_lapack.h"

static VALUE
rb_zlascl2(int argc, VALUE *argv, VALUE self){
  VALUE rb_d;
  doublereal *d; 
  VALUE rb_x;
  doublecomplex *x; 
  VALUE rb_x_out__;
  doublecomplex *x_out__;

  integer m;
  integer ldx;
  integer n;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  x = NumRu::Lapack.zlascl2( d, x)\n    or\n  NumRu::Lapack.zlascl2  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE ZLASCL2 ( M, N, D, X, LDX )\n\n*  Purpose\n*  =======\n*\n*  ZLASCL2 performs a diagonal scaling on a vector:\n*    x <-- D * x\n*  where the DOUBLE PRECISION diagonal matrix D is stored as a vector.\n*\n*  Eventually to be replaced by BLAS_zge_diag_scale in the new BLAS\n*  standard.\n*\n\n*  Arguments\n*  =========\n*\n*     M       (input) INTEGER\n*     The number of rows of D and X. M >= 0.\n*\n*     N       (input) INTEGER\n*     The number of columns of D and X. N >= 0.\n*\n*     D       (input) DOUBLE PRECISION array, length M\n*     Diagonal matrix D, stored as a vector of length M.\n*\n*     X       (input/output) COMPLEX*16 array, dimension (LDX,N)\n*     On entry, the vector X to be scaled by D.\n*     On exit, the scaled vector.\n*\n*     LDX     (input) INTEGER\n*     The leading dimension of the vector X. LDX >= 0.\n*\n\n*  =====================================================================\n*\n*     .. Local Scalars ..\n      INTEGER            I, J\n*     ..\n\n");
    return Qnil;
  }
  if (argc != 2)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 2)", argc);
  rb_d = argv[0];
  rb_x = argv[1];

  if (!NA_IsNArray(rb_d))
    rb_raise(rb_eArgError, "d (1th argument) must be NArray");
  if (NA_RANK(rb_d) != 1)
    rb_raise(rb_eArgError, "rank of d (1th argument) must be %d", 1);
  m = NA_SHAPE0(rb_d);
  if (NA_TYPE(rb_d) != NA_DFLOAT)
    rb_d = na_change_type(rb_d, NA_DFLOAT);
  d = NA_PTR_TYPE(rb_d, doublereal*);
  if (!NA_IsNArray(rb_x))
    rb_raise(rb_eArgError, "x (2th argument) must be NArray");
  if (NA_RANK(rb_x) != 2)
    rb_raise(rb_eArgError, "rank of x (2th argument) must be %d", 2);
  ldx = NA_SHAPE0(rb_x);
  n = NA_SHAPE1(rb_x);
  if (NA_TYPE(rb_x) != NA_DCOMPLEX)
    rb_x = na_change_type(rb_x, NA_DCOMPLEX);
  x = NA_PTR_TYPE(rb_x, doublecomplex*);
  {
    int shape[2];
    shape[0] = ldx;
    shape[1] = n;
    rb_x_out__ = na_make_object(NA_DCOMPLEX, 2, shape, cNArray);
  }
  x_out__ = NA_PTR_TYPE(rb_x_out__, doublecomplex*);
  MEMCPY(x_out__, x, doublecomplex, NA_TOTAL(rb_x));
  rb_x = rb_x_out__;
  x = x_out__;

  zlascl2_(&m, &n, d, x, &ldx);

  return rb_x;
}

void
init_lapack_zlascl2(VALUE mLapack){
  rb_define_module_function(mLapack, "zlascl2", rb_zlascl2, -1);
}

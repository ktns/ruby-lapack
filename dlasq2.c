#include "rb_lapack.h"

extern VOID dlasq2_(integer *n, doublereal *z, integer *info);

static VALUE
rb_dlasq2(int argc, VALUE *argv, VALUE self){
  VALUE rb_n;
  integer n; 
  VALUE rb_z;
  doublereal *z; 
  VALUE rb_info;
  integer info; 
  VALUE rb_z_out__;
  doublereal *z_out__;


  if (argc == 0) {
    printf("%s\n", "USAGE:\n  info, z = NumRu::Lapack.dlasq2( n, z)\n    or\n  NumRu::Lapack.dlasq2  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE DLASQ2( N, Z, INFO )\n\n*  Purpose\n*  =======\n*\n*  DLASQ2 computes all the eigenvalues of the symmetric positive \n*  definite tridiagonal matrix associated with the qd array Z to high\n*  relative accuracy are computed to high relative accuracy, in the\n*  absence of denormalization, underflow and overflow.\n*\n*  To see the relation of Z to the tridiagonal matrix, let L be a\n*  unit lower bidiagonal matrix with subdiagonals Z(2,4,6,,..) and\n*  let U be an upper bidiagonal matrix with 1's above and diagonal\n*  Z(1,3,5,,..). The tridiagonal is L*U or, if you prefer, the\n*  symmetric tridiagonal to which it is similar.\n*\n*  Note : DLASQ2 defines a logical variable, IEEE, which is true\n*  on machines which follow ieee-754 floating-point standard in their\n*  handling of infinities and NaNs, and false otherwise. This variable\n*  is passed to DLASQ3.\n*\n\n*  Arguments\n*  =========\n*\n*  N     (input) INTEGER\n*        The number of rows and columns in the matrix. N >= 0.\n*\n*  Z     (input/output) DOUBLE PRECISION array, dimension ( 4*N )\n*        On entry Z holds the qd array. On exit, entries 1 to N hold\n*        the eigenvalues in decreasing order, Z( 2*N+1 ) holds the\n*        trace, and Z( 2*N+2 ) holds the sum of the eigenvalues. If\n*        N > 2, then Z( 2*N+3 ) holds the iteration count, Z( 2*N+4 )\n*        holds NDIVS/NIN^2, and Z( 2*N+5 ) holds the percentage of\n*        shifts that failed.\n*\n*  INFO  (output) INTEGER\n*        = 0: successful exit\n*        < 0: if the i-th argument is a scalar and had an illegal\n*             value, then INFO = -i, if the i-th argument is an\n*             array and the j-entry had an illegal value, then\n*             INFO = -(i*100+j)\n*        > 0: the algorithm failed\n*              = 1, a split was marked by a positive value in E\n*              = 2, current block of Z not diagonalized after 30*N\n*                   iterations (in inner while loop)\n*              = 3, termination criterion of outer while loop not met \n*                   (program created more than N unreduced blocks)\n*\n\n*  Further Details\n*  ===============\n*  Local Variables: I0:N0 defines a current unreduced segment of Z.\n*  The shifts are accumulated in SIGMA. Iteration count is in ITER.\n*  Ping-pong is controlled by PP (alternates between 0 and 1).\n*\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 2)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 2)", argc);
  rb_n = argv[0];
  rb_z = argv[1];

  n = NUM2INT(rb_n);
  if (!NA_IsNArray(rb_z))
    rb_raise(rb_eArgError, "z (2th argument) must be NArray");
  if (NA_RANK(rb_z) != 1)
    rb_raise(rb_eArgError, "rank of z (2th argument) must be %d", 1);
  if (NA_SHAPE0(rb_z) != (4*n))
    rb_raise(rb_eRuntimeError, "shape 0 of z must be %d", 4*n);
  if (NA_TYPE(rb_z) != NA_DFLOAT)
    rb_z = na_change_type(rb_z, NA_DFLOAT);
  z = NA_PTR_TYPE(rb_z, doublereal*);
  {
    int shape[1];
    shape[0] = 4*n;
    rb_z_out__ = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  z_out__ = NA_PTR_TYPE(rb_z_out__, doublereal*);
  MEMCPY(z_out__, z, doublereal, NA_TOTAL(rb_z));
  rb_z = rb_z_out__;
  z = z_out__;

  dlasq2_(&n, z, &info);

  rb_info = INT2NUM(info);
  return rb_ary_new3(2, rb_info, rb_z);
}

void
init_lapack_dlasq2(VALUE mLapack){
  rb_define_module_function(mLapack, "dlasq2", rb_dlasq2, -1);
}

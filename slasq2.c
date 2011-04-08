#include "rb_lapack.h"

extern VOID slasq2_(integer *n, real *z, integer *info);

static VALUE sHelp, sUsage;

static VALUE
rblapack_slasq2(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_n;
  integer n; 
  VALUE rblapack_z;
  real *z; 
  VALUE rblapack_info;
  integer info; 
  VALUE rblapack_z_out__;
  real *z_out__;


  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  info, z = NumRu::Lapack.slasq2( n, z, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE SLASQ2( N, Z, INFO )\n\n*  Purpose\n*  =======\n*\n*  SLASQ2 computes all the eigenvalues of the symmetric positive \n*  definite tridiagonal matrix associated with the qd array Z to high\n*  relative accuracy are computed to high relative accuracy, in the\n*  absence of denormalization, underflow and overflow.\n*\n*  To see the relation of Z to the tridiagonal matrix, let L be a\n*  unit lower bidiagonal matrix with subdiagonals Z(2,4,6,,..) and\n*  let U be an upper bidiagonal matrix with 1's above and diagonal\n*  Z(1,3,5,,..). The tridiagonal is L*U or, if you prefer, the\n*  symmetric tridiagonal to which it is similar.\n*\n*  Note : SLASQ2 defines a logical variable, IEEE, which is true\n*  on machines which follow ieee-754 floating-point standard in their\n*  handling of infinities and NaNs, and false otherwise. This variable\n*  is passed to SLASQ3.\n*\n\n*  Arguments\n*  =========\n*\n*  N     (input) INTEGER\n*        The number of rows and columns in the matrix. N >= 0.\n*\n*  Z     (input/output) REAL array, dimension ( 4*N )\n*        On entry Z holds the qd array. On exit, entries 1 to N hold\n*        the eigenvalues in decreasing order, Z( 2*N+1 ) holds the\n*        trace, and Z( 2*N+2 ) holds the sum of the eigenvalues. If\n*        N > 2, then Z( 2*N+3 ) holds the iteration count, Z( 2*N+4 )\n*        holds NDIVS/NIN^2, and Z( 2*N+5 ) holds the percentage of\n*        shifts that failed.\n*\n*  INFO  (output) INTEGER\n*        = 0: successful exit\n*        < 0: if the i-th argument is a scalar and had an illegal\n*             value, then INFO = -i, if the i-th argument is an\n*             array and the j-entry had an illegal value, then\n*             INFO = -(i*100+j)\n*        > 0: the algorithm failed\n*              = 1, a split was marked by a positive value in E\n*              = 2, current block of Z not diagonalized after 30*N\n*                   iterations (in inner while loop)\n*              = 3, termination criterion of outer while loop not met \n*                   (program created more than N unreduced blocks)\n*\n\n*  Further Details\n*  ===============\n*  Local Variables: I0:N0 defines a current unreduced segment of Z.\n*  The shifts are accumulated in SIGMA. Iteration count is in ITER.\n*  Ping-pong is controlled by PP (alternates between 0 and 1).\n*\n*  =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  info, z = NumRu::Lapack.slasq2( n, z, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 2)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 2)", argc);
  rblapack_n = argv[0];
  rblapack_z = argv[1];
  if (rb_options != Qnil) {
  }

  n = NUM2INT(rblapack_n);
  if (!NA_IsNArray(rblapack_z))
    rb_raise(rb_eArgError, "z (2th argument) must be NArray");
  if (NA_RANK(rblapack_z) != 1)
    rb_raise(rb_eArgError, "rank of z (2th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_z) != (4*n))
    rb_raise(rb_eRuntimeError, "shape 0 of z must be %d", 4*n);
  if (NA_TYPE(rblapack_z) != NA_SFLOAT)
    rblapack_z = na_change_type(rblapack_z, NA_SFLOAT);
  z = NA_PTR_TYPE(rblapack_z, real*);
  {
    int shape[1];
    shape[0] = 4*n;
    rblapack_z_out__ = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  z_out__ = NA_PTR_TYPE(rblapack_z_out__, real*);
  MEMCPY(z_out__, z, real, NA_TOTAL(rblapack_z));
  rblapack_z = rblapack_z_out__;
  z = z_out__;

  slasq2_(&n, z, &info);

  rblapack_info = INT2NUM(info);
  return rb_ary_new3(2, rblapack_info, rblapack_z);
}

void
init_lapack_slasq2(VALUE mLapack){
  rb_define_module_function(mLapack, "slasq2", rblapack_slasq2, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}

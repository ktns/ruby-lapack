#include "rb_lapack.h"

extern VOID dlarscl2_(integer *m, integer *n, doublereal *d, doublereal *x, integer *ldx);

static VALUE sHelp, sUsage;

static VALUE
rblapack_dlarscl2(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_d;
  doublereal *d; 
  VALUE rblapack_x;
  doublereal *x; 
  VALUE rblapack_x_out__;
  doublereal *x_out__;

  integer m;
  integer ldx;
  integer n;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  x = NumRu::Lapack.dlarscl2( d, x, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE DLARSCL2 ( M, N, D, X, LDX )\n\n*  Purpose\n*  =======\n*\n*  DLARSCL2 performs a reciprocal diagonal scaling on an vector:\n*    x <-- inv(D) * x\n*  where the diagonal matrix D is stored as a vector.\n*\n*  Eventually to be replaced by BLAS_dge_diag_scale in the new BLAS\n*  standard.\n*\n\n*  Arguments\n*  =========\n*\n*     M       (input) INTEGER\n*     The number of rows of D and X. M >= 0.\n*\n*     N       (input) INTEGER\n*     The number of columns of D and X. N >= 0.\n*\n*     D       (input) DOUBLE PRECISION array, dimension (M)\n*     Diagonal matrix D, stored as a vector of length M.\n*\n*     X       (input/output) DOUBLE PRECISION array, dimension (LDX,N)\n*     On entry, the vector X to be scaled by D.\n*     On exit, the scaled vector.\n*\n*     LDX     (input) INTEGER\n*     The leading dimension of the vector X. LDX >= 0.\n*\n\n*  =====================================================================\n*\n*     .. Local Scalars ..\n      INTEGER            I, J\n*     ..\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  x = NumRu::Lapack.dlarscl2( d, x, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 2)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 2)", argc);
  rblapack_d = argv[0];
  rblapack_x = argv[1];
  if (rb_options != Qnil) {
  }

  if (!NA_IsNArray(rblapack_x))
    rb_raise(rb_eArgError, "x (2th argument) must be NArray");
  if (NA_RANK(rblapack_x) != 2)
    rb_raise(rb_eArgError, "rank of x (2th argument) must be %d", 2);
  n = NA_SHAPE1(rblapack_x);
  ldx = NA_SHAPE0(rblapack_x);
  if (NA_TYPE(rblapack_x) != NA_DFLOAT)
    rblapack_x = na_change_type(rblapack_x, NA_DFLOAT);
  x = NA_PTR_TYPE(rblapack_x, doublereal*);
  if (!NA_IsNArray(rblapack_d))
    rb_raise(rb_eArgError, "d (1th argument) must be NArray");
  if (NA_RANK(rblapack_d) != 1)
    rb_raise(rb_eArgError, "rank of d (1th argument) must be %d", 1);
  m = NA_SHAPE0(rblapack_d);
  if (NA_TYPE(rblapack_d) != NA_DFLOAT)
    rblapack_d = na_change_type(rblapack_d, NA_DFLOAT);
  d = NA_PTR_TYPE(rblapack_d, doublereal*);
  {
    int shape[2];
    shape[0] = ldx;
    shape[1] = n;
    rblapack_x_out__ = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  x_out__ = NA_PTR_TYPE(rblapack_x_out__, doublereal*);
  MEMCPY(x_out__, x, doublereal, NA_TOTAL(rblapack_x));
  rblapack_x = rblapack_x_out__;
  x = x_out__;

  dlarscl2_(&m, &n, d, x, &ldx);

  return rblapack_x;
}

void
init_lapack_dlarscl2(VALUE mLapack){
  rb_define_module_function(mLapack, "dlarscl2", rblapack_dlarscl2, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}

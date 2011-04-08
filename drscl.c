#include "rb_lapack.h"

extern VOID drscl_(integer *n, doublereal *sa, doublereal *sx, integer *incx);

static VALUE sHelp, sUsage;

static VALUE
rblapack_drscl(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_n;
  integer n; 
  VALUE rblapack_sa;
  doublereal sa; 
  VALUE rblapack_sx;
  doublereal *sx; 
  VALUE rblapack_incx;
  integer incx; 
  VALUE rblapack_sx_out__;
  doublereal *sx_out__;


  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  sx = NumRu::Lapack.drscl( n, sa, sx, incx, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE DRSCL( N, SA, SX, INCX )\n\n*  Purpose\n*  =======\n*\n*  DRSCL multiplies an n-element real vector x by the real scalar 1/a.\n*  This is done without overflow or underflow as long as\n*  the final result x/a does not overflow or underflow.\n*\n\n*  Arguments\n*  =========\n*\n*  N       (input) INTEGER\n*          The number of components of the vector x.\n*\n*  SA      (input) DOUBLE PRECISION\n*          The scalar a which is used to divide each component of x.\n*          SA must be >= 0, or the subroutine will divide by zero.\n*\n*  SX      (input/output) DOUBLE PRECISION array, dimension\n*                         (1+(N-1)*abs(INCX))\n*          The n-element vector x.\n*\n*  INCX    (input) INTEGER\n*          The increment between successive values of the vector SX.\n*          > 0:  SX(1) = X(1) and SX(1+(i-1)*INCX) = x(i),     1< i<= n\n*\n\n* =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  sx = NumRu::Lapack.drscl( n, sa, sx, incx, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 4)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 4)", argc);
  rblapack_n = argv[0];
  rblapack_sa = argv[1];
  rblapack_sx = argv[2];
  rblapack_incx = argv[3];
  if (rb_options != Qnil) {
  }

  sa = NUM2DBL(rblapack_sa);
  n = NUM2INT(rblapack_n);
  incx = NUM2INT(rblapack_incx);
  if (!NA_IsNArray(rblapack_sx))
    rb_raise(rb_eArgError, "sx (3th argument) must be NArray");
  if (NA_RANK(rblapack_sx) != 1)
    rb_raise(rb_eArgError, "rank of sx (3th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_sx) != (1+(n-1)*abs(incx)))
    rb_raise(rb_eRuntimeError, "shape 0 of sx must be %d", 1+(n-1)*abs(incx));
  if (NA_TYPE(rblapack_sx) != NA_DFLOAT)
    rblapack_sx = na_change_type(rblapack_sx, NA_DFLOAT);
  sx = NA_PTR_TYPE(rblapack_sx, doublereal*);
  {
    int shape[1];
    shape[0] = 1+(n-1)*abs(incx);
    rblapack_sx_out__ = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  sx_out__ = NA_PTR_TYPE(rblapack_sx_out__, doublereal*);
  MEMCPY(sx_out__, sx, doublereal, NA_TOTAL(rblapack_sx));
  rblapack_sx = rblapack_sx_out__;
  sx = sx_out__;

  drscl_(&n, &sa, sx, &incx);

  return rblapack_sx;
}

void
init_lapack_drscl(VALUE mLapack){
  rb_define_module_function(mLapack, "drscl", rblapack_drscl, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}

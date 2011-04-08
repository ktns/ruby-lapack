#include "rb_lapack.h"

extern VOID slar2v_(integer *n, real *x, real *y, real *z, integer *incx, real *c, real *s, integer *incc);

static VALUE sHelp, sUsage;

static VALUE
rblapack_slar2v(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_n;
  integer n; 
  VALUE rblapack_x;
  real *x; 
  VALUE rblapack_y;
  real *y; 
  VALUE rblapack_z;
  real *z; 
  VALUE rblapack_incx;
  integer incx; 
  VALUE rblapack_c;
  real *c; 
  VALUE rblapack_s;
  real *s; 
  VALUE rblapack_incc;
  integer incc; 
  VALUE rblapack_x_out__;
  real *x_out__;
  VALUE rblapack_y_out__;
  real *y_out__;
  VALUE rblapack_z_out__;
  real *z_out__;


  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  x, y, z = NumRu::Lapack.slar2v( n, x, y, z, incx, c, s, incc, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE SLAR2V( N, X, Y, Z, INCX, C, S, INCC )\n\n*  Purpose\n*  =======\n*\n*  SLAR2V applies a vector of real plane rotations from both sides to\n*  a sequence of 2-by-2 real symmetric matrices, defined by the elements\n*  of the vectors x, y and z. For i = 1,2,...,n\n*\n*     ( x(i)  z(i) ) := (  c(i)  s(i) ) ( x(i)  z(i) ) ( c(i) -s(i) )\n*     ( z(i)  y(i) )    ( -s(i)  c(i) ) ( z(i)  y(i) ) ( s(i)  c(i) )\n*\n\n*  Arguments\n*  =========\n*\n*  N       (input) INTEGER\n*          The number of plane rotations to be applied.\n*\n*  X       (input/output) REAL array,\n*                         dimension (1+(N-1)*INCX)\n*          The vector x.\n*\n*  Y       (input/output) REAL array,\n*                         dimension (1+(N-1)*INCX)\n*          The vector y.\n*\n*  Z       (input/output) REAL array,\n*                         dimension (1+(N-1)*INCX)\n*          The vector z.\n*\n*  INCX    (input) INTEGER\n*          The increment between elements of X, Y and Z. INCX > 0.\n*\n*  C       (input) REAL array, dimension (1+(N-1)*INCC)\n*          The cosines of the plane rotations.\n*\n*  S       (input) REAL array, dimension (1+(N-1)*INCC)\n*          The sines of the plane rotations.\n*\n*  INCC    (input) INTEGER\n*          The increment between elements of C and S. INCC > 0.\n*\n\n*  =====================================================================\n*\n*     .. Local Scalars ..\n      INTEGER            I, IC, IX\n      REAL               CI, SI, T1, T2, T3, T4, T5, T6, XI, YI, ZI\n*     ..\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  x, y, z = NumRu::Lapack.slar2v( n, x, y, z, incx, c, s, incc, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 8)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 8)", argc);
  rblapack_n = argv[0];
  rblapack_x = argv[1];
  rblapack_y = argv[2];
  rblapack_z = argv[3];
  rblapack_incx = argv[4];
  rblapack_c = argv[5];
  rblapack_s = argv[6];
  rblapack_incc = argv[7];
  if (rb_options != Qnil) {
  }

  n = NUM2INT(rblapack_n);
  incc = NUM2INT(rblapack_incc);
  incx = NUM2INT(rblapack_incx);
  if (!NA_IsNArray(rblapack_z))
    rb_raise(rb_eArgError, "z (4th argument) must be NArray");
  if (NA_RANK(rblapack_z) != 1)
    rb_raise(rb_eArgError, "rank of z (4th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_z) != (1+(n-1)*incx))
    rb_raise(rb_eRuntimeError, "shape 0 of z must be %d", 1+(n-1)*incx);
  if (NA_TYPE(rblapack_z) != NA_SFLOAT)
    rblapack_z = na_change_type(rblapack_z, NA_SFLOAT);
  z = NA_PTR_TYPE(rblapack_z, real*);
  if (!NA_IsNArray(rblapack_x))
    rb_raise(rb_eArgError, "x (2th argument) must be NArray");
  if (NA_RANK(rblapack_x) != 1)
    rb_raise(rb_eArgError, "rank of x (2th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_x) != (1+(n-1)*incx))
    rb_raise(rb_eRuntimeError, "shape 0 of x must be %d", 1+(n-1)*incx);
  if (NA_TYPE(rblapack_x) != NA_SFLOAT)
    rblapack_x = na_change_type(rblapack_x, NA_SFLOAT);
  x = NA_PTR_TYPE(rblapack_x, real*);
  if (!NA_IsNArray(rblapack_y))
    rb_raise(rb_eArgError, "y (3th argument) must be NArray");
  if (NA_RANK(rblapack_y) != 1)
    rb_raise(rb_eArgError, "rank of y (3th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_y) != (1+(n-1)*incx))
    rb_raise(rb_eRuntimeError, "shape 0 of y must be %d", 1+(n-1)*incx);
  if (NA_TYPE(rblapack_y) != NA_SFLOAT)
    rblapack_y = na_change_type(rblapack_y, NA_SFLOAT);
  y = NA_PTR_TYPE(rblapack_y, real*);
  if (!NA_IsNArray(rblapack_c))
    rb_raise(rb_eArgError, "c (6th argument) must be NArray");
  if (NA_RANK(rblapack_c) != 1)
    rb_raise(rb_eArgError, "rank of c (6th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_c) != (1+(n-1)*incc))
    rb_raise(rb_eRuntimeError, "shape 0 of c must be %d", 1+(n-1)*incc);
  if (NA_TYPE(rblapack_c) != NA_SFLOAT)
    rblapack_c = na_change_type(rblapack_c, NA_SFLOAT);
  c = NA_PTR_TYPE(rblapack_c, real*);
  if (!NA_IsNArray(rblapack_s))
    rb_raise(rb_eArgError, "s (7th argument) must be NArray");
  if (NA_RANK(rblapack_s) != 1)
    rb_raise(rb_eArgError, "rank of s (7th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_s) != (1+(n-1)*incc))
    rb_raise(rb_eRuntimeError, "shape 0 of s must be %d", 1+(n-1)*incc);
  if (NA_TYPE(rblapack_s) != NA_SFLOAT)
    rblapack_s = na_change_type(rblapack_s, NA_SFLOAT);
  s = NA_PTR_TYPE(rblapack_s, real*);
  {
    int shape[1];
    shape[0] = 1+(n-1)*incx;
    rblapack_x_out__ = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  x_out__ = NA_PTR_TYPE(rblapack_x_out__, real*);
  MEMCPY(x_out__, x, real, NA_TOTAL(rblapack_x));
  rblapack_x = rblapack_x_out__;
  x = x_out__;
  {
    int shape[1];
    shape[0] = 1+(n-1)*incx;
    rblapack_y_out__ = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  y_out__ = NA_PTR_TYPE(rblapack_y_out__, real*);
  MEMCPY(y_out__, y, real, NA_TOTAL(rblapack_y));
  rblapack_y = rblapack_y_out__;
  y = y_out__;
  {
    int shape[1];
    shape[0] = 1+(n-1)*incx;
    rblapack_z_out__ = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  z_out__ = NA_PTR_TYPE(rblapack_z_out__, real*);
  MEMCPY(z_out__, z, real, NA_TOTAL(rblapack_z));
  rblapack_z = rblapack_z_out__;
  z = z_out__;

  slar2v_(&n, x, y, z, &incx, c, s, &incc);

  return rb_ary_new3(3, rblapack_x, rblapack_y, rblapack_z);
}

void
init_lapack_slar2v(VALUE mLapack){
  rb_define_module_function(mLapack, "slar2v", rblapack_slar2v, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}

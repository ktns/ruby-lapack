#include "rb_lapack.h"

extern VOID clar2v_(integer *n, complex *x, complex *y, complex *z, integer *incx, real *c, complex *s, integer *incc);

static VALUE sHelp, sUsage;

static VALUE
rblapack_clar2v(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_n;
  integer n; 
  VALUE rblapack_x;
  complex *x; 
  VALUE rblapack_y;
  complex *y; 
  VALUE rblapack_z;
  complex *z; 
  VALUE rblapack_incx;
  integer incx; 
  VALUE rblapack_c;
  real *c; 
  VALUE rblapack_s;
  complex *s; 
  VALUE rblapack_incc;
  integer incc; 
  VALUE rblapack_x_out__;
  complex *x_out__;
  VALUE rblapack_y_out__;
  complex *y_out__;
  VALUE rblapack_z_out__;
  complex *z_out__;


  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  x, y, z = NumRu::Lapack.clar2v( n, x, y, z, incx, c, s, incc, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE CLAR2V( N, X, Y, Z, INCX, C, S, INCC )\n\n*  Purpose\n*  =======\n*\n*  CLAR2V applies a vector of complex plane rotations with real cosines\n*  from both sides to a sequence of 2-by-2 complex Hermitian matrices,\n*  defined by the elements of the vectors x, y and z. For i = 1,2,...,n\n*\n*     (       x(i)  z(i) ) :=\n*     ( conjg(z(i)) y(i) )\n*\n*       (  c(i) conjg(s(i)) ) (       x(i)  z(i) ) ( c(i) -conjg(s(i)) )\n*       ( -s(i)       c(i)  ) ( conjg(z(i)) y(i) ) ( s(i)        c(i)  )\n*\n\n*  Arguments\n*  =========\n*\n*  N       (input) INTEGER\n*          The number of plane rotations to be applied.\n*\n*  X       (input/output) COMPLEX array, dimension (1+(N-1)*INCX)\n*          The vector x; the elements of x are assumed to be real.\n*\n*  Y       (input/output) COMPLEX array, dimension (1+(N-1)*INCX)\n*          The vector y; the elements of y are assumed to be real.\n*\n*  Z       (input/output) COMPLEX array, dimension (1+(N-1)*INCX)\n*          The vector z.\n*\n*  INCX    (input) INTEGER\n*          The increment between elements of X, Y and Z. INCX > 0.\n*\n*  C       (input) REAL array, dimension (1+(N-1)*INCC)\n*          The cosines of the plane rotations.\n*\n*  S       (input) COMPLEX array, dimension (1+(N-1)*INCC)\n*          The sines of the plane rotations.\n*\n*  INCC    (input) INTEGER\n*          The increment between elements of C and S. INCC > 0.\n*\n\n*  =====================================================================\n*\n*     .. Local Scalars ..\n      INTEGER            I, IC, IX\n      REAL               CI, SII, SIR, T1I, T1R, T5, T6, XI, YI, ZII,\n     $                   ZIR\n      COMPLEX            SI, T2, T3, T4, ZI\n*     ..\n*     .. Intrinsic Functions ..\n      INTRINSIC          AIMAG, CMPLX, CONJG, REAL\n*     ..\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  x, y, z = NumRu::Lapack.clar2v( n, x, y, z, incx, c, s, incc, [:usage => usage, :help => help])\n");
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
  if (NA_TYPE(rblapack_z) != NA_SCOMPLEX)
    rblapack_z = na_change_type(rblapack_z, NA_SCOMPLEX);
  z = NA_PTR_TYPE(rblapack_z, complex*);
  if (!NA_IsNArray(rblapack_x))
    rb_raise(rb_eArgError, "x (2th argument) must be NArray");
  if (NA_RANK(rblapack_x) != 1)
    rb_raise(rb_eArgError, "rank of x (2th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_x) != (1+(n-1)*incx))
    rb_raise(rb_eRuntimeError, "shape 0 of x must be %d", 1+(n-1)*incx);
  if (NA_TYPE(rblapack_x) != NA_SCOMPLEX)
    rblapack_x = na_change_type(rblapack_x, NA_SCOMPLEX);
  x = NA_PTR_TYPE(rblapack_x, complex*);
  if (!NA_IsNArray(rblapack_y))
    rb_raise(rb_eArgError, "y (3th argument) must be NArray");
  if (NA_RANK(rblapack_y) != 1)
    rb_raise(rb_eArgError, "rank of y (3th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_y) != (1+(n-1)*incx))
    rb_raise(rb_eRuntimeError, "shape 0 of y must be %d", 1+(n-1)*incx);
  if (NA_TYPE(rblapack_y) != NA_SCOMPLEX)
    rblapack_y = na_change_type(rblapack_y, NA_SCOMPLEX);
  y = NA_PTR_TYPE(rblapack_y, complex*);
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
  if (NA_TYPE(rblapack_s) != NA_SCOMPLEX)
    rblapack_s = na_change_type(rblapack_s, NA_SCOMPLEX);
  s = NA_PTR_TYPE(rblapack_s, complex*);
  {
    int shape[1];
    shape[0] = 1+(n-1)*incx;
    rblapack_x_out__ = na_make_object(NA_SCOMPLEX, 1, shape, cNArray);
  }
  x_out__ = NA_PTR_TYPE(rblapack_x_out__, complex*);
  MEMCPY(x_out__, x, complex, NA_TOTAL(rblapack_x));
  rblapack_x = rblapack_x_out__;
  x = x_out__;
  {
    int shape[1];
    shape[0] = 1+(n-1)*incx;
    rblapack_y_out__ = na_make_object(NA_SCOMPLEX, 1, shape, cNArray);
  }
  y_out__ = NA_PTR_TYPE(rblapack_y_out__, complex*);
  MEMCPY(y_out__, y, complex, NA_TOTAL(rblapack_y));
  rblapack_y = rblapack_y_out__;
  y = y_out__;
  {
    int shape[1];
    shape[0] = 1+(n-1)*incx;
    rblapack_z_out__ = na_make_object(NA_SCOMPLEX, 1, shape, cNArray);
  }
  z_out__ = NA_PTR_TYPE(rblapack_z_out__, complex*);
  MEMCPY(z_out__, z, complex, NA_TOTAL(rblapack_z));
  rblapack_z = rblapack_z_out__;
  z = z_out__;

  clar2v_(&n, x, y, z, &incx, c, s, &incc);

  return rb_ary_new3(3, rblapack_x, rblapack_y, rblapack_z);
}

void
init_lapack_clar2v(VALUE mLapack){
  rb_define_module_function(mLapack, "clar2v", rblapack_clar2v, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}

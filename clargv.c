#include "rb_lapack.h"

extern VOID clargv_(integer *n, complex *x, integer *incx, complex *y, integer *incy, real *c, integer *incc);

static VALUE sHelp, sUsage;

static VALUE
rblapack_clargv(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_n;
  integer n; 
  VALUE rblapack_x;
  complex *x; 
  VALUE rblapack_incx;
  integer incx; 
  VALUE rblapack_y;
  complex *y; 
  VALUE rblapack_incy;
  integer incy; 
  VALUE rblapack_incc;
  integer incc; 
  VALUE rblapack_c;
  real *c; 
  VALUE rblapack_x_out__;
  complex *x_out__;
  VALUE rblapack_y_out__;
  complex *y_out__;


  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  c, x, y = NumRu::Lapack.clargv( n, x, incx, y, incy, incc, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE CLARGV( N, X, INCX, Y, INCY, C, INCC )\n\n*  Purpose\n*  =======\n*\n*  CLARGV generates a vector of complex plane rotations with real\n*  cosines, determined by elements of the complex vectors x and y.\n*  For i = 1,2,...,n\n*\n*     (        c(i)   s(i) ) ( x(i) ) = ( r(i) )\n*     ( -conjg(s(i))  c(i) ) ( y(i) ) = (   0  )\n*\n*     where c(i)**2 + ABS(s(i))**2 = 1\n*\n*  The following conventions are used (these are the same as in CLARTG,\n*  but differ from the BLAS1 routine CROTG):\n*     If y(i)=0, then c(i)=1 and s(i)=0.\n*     If x(i)=0, then c(i)=0 and s(i) is chosen so that r(i) is real.\n*\n\n*  Arguments\n*  =========\n*\n*  N       (input) INTEGER\n*          The number of plane rotations to be generated.\n*\n*  X       (input/output) COMPLEX array, dimension (1+(N-1)*INCX)\n*          On entry, the vector x.\n*          On exit, x(i) is overwritten by r(i), for i = 1,...,n.\n*\n*  INCX    (input) INTEGER\n*          The increment between elements of X. INCX > 0.\n*\n*  Y       (input/output) COMPLEX array, dimension (1+(N-1)*INCY)\n*          On entry, the vector y.\n*          On exit, the sines of the plane rotations.\n*\n*  INCY    (input) INTEGER\n*          The increment between elements of Y. INCY > 0.\n*\n*  C       (output) REAL array, dimension (1+(N-1)*INCC)\n*          The cosines of the plane rotations.\n*\n*  INCC    (input) INTEGER\n*          The increment between elements of C. INCC > 0.\n*\n\n*  Further Details\n*  ======= =======\n*\n*  6-6-96 - Modified with a new algorithm by W. Kahan and J. Demmel\n*\n*  This version has a few statements commented out for thread safety\n*  (machine parameters are computed on each entry). 10 feb 03, SJH.\n*\n*  =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  c, x, y = NumRu::Lapack.clargv( n, x, incx, y, incy, incc, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 6)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 6)", argc);
  rblapack_n = argv[0];
  rblapack_x = argv[1];
  rblapack_incx = argv[2];
  rblapack_y = argv[3];
  rblapack_incy = argv[4];
  rblapack_incc = argv[5];
  if (rb_options != Qnil) {
  }

  incy = NUM2INT(rblapack_incy);
  incc = NUM2INT(rblapack_incc);
  n = NUM2INT(rblapack_n);
  incx = NUM2INT(rblapack_incx);
  if (!NA_IsNArray(rblapack_y))
    rb_raise(rb_eArgError, "y (4th argument) must be NArray");
  if (NA_RANK(rblapack_y) != 1)
    rb_raise(rb_eArgError, "rank of y (4th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_y) != (1+(n-1)*incy))
    rb_raise(rb_eRuntimeError, "shape 0 of y must be %d", 1+(n-1)*incy);
  if (NA_TYPE(rblapack_y) != NA_SCOMPLEX)
    rblapack_y = na_change_type(rblapack_y, NA_SCOMPLEX);
  y = NA_PTR_TYPE(rblapack_y, complex*);
  if (!NA_IsNArray(rblapack_x))
    rb_raise(rb_eArgError, "x (2th argument) must be NArray");
  if (NA_RANK(rblapack_x) != 1)
    rb_raise(rb_eArgError, "rank of x (2th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_x) != (1+(n-1)*incx))
    rb_raise(rb_eRuntimeError, "shape 0 of x must be %d", 1+(n-1)*incx);
  if (NA_TYPE(rblapack_x) != NA_SCOMPLEX)
    rblapack_x = na_change_type(rblapack_x, NA_SCOMPLEX);
  x = NA_PTR_TYPE(rblapack_x, complex*);
  {
    int shape[1];
    shape[0] = 1+(n-1)*incc;
    rblapack_c = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  c = NA_PTR_TYPE(rblapack_c, real*);
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
    shape[0] = 1+(n-1)*incy;
    rblapack_y_out__ = na_make_object(NA_SCOMPLEX, 1, shape, cNArray);
  }
  y_out__ = NA_PTR_TYPE(rblapack_y_out__, complex*);
  MEMCPY(y_out__, y, complex, NA_TOTAL(rblapack_y));
  rblapack_y = rblapack_y_out__;
  y = y_out__;

  clargv_(&n, x, &incx, y, &incy, c, &incc);

  return rb_ary_new3(3, rblapack_c, rblapack_x, rblapack_y);
}

void
init_lapack_clargv(VALUE mLapack){
  rb_define_module_function(mLapack, "clargv", rblapack_clargv, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}

#include "rb_lapack.h"

extern VOID dlarnv_(integer *idist, integer *iseed, integer *n, doublereal *x);

static VALUE
rb_dlarnv(int argc, VALUE *argv, VALUE self){
  VALUE rb_idist;
  integer idist; 
  VALUE rb_iseed;
  integer *iseed; 
  VALUE rb_n;
  integer n; 
  VALUE rb_x;
  doublereal *x; 
  VALUE rb_iseed_out__;
  integer *iseed_out__;


  if (argc == 0) {
    printf("%s\n", "USAGE:\n  x, iseed = NumRu::Lapack.dlarnv( idist, iseed, n)\n    or\n  NumRu::Lapack.dlarnv  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE DLARNV( IDIST, ISEED, N, X )\n\n*  Purpose\n*  =======\n*\n*  DLARNV returns a vector of n random real numbers from a uniform or\n*  normal distribution.\n*\n\n*  Arguments\n*  =========\n*\n*  IDIST   (input) INTEGER\n*          Specifies the distribution of the random numbers:\n*          = 1:  uniform (0,1)\n*          = 2:  uniform (-1,1)\n*          = 3:  normal (0,1)\n*\n*  ISEED   (input/output) INTEGER array, dimension (4)\n*          On entry, the seed of the random number generator; the array\n*          elements must be between 0 and 4095, and ISEED(4) must be\n*          odd.\n*          On exit, the seed is updated.\n*\n*  N       (input) INTEGER\n*          The number of random numbers to be generated.\n*\n*  X       (output) DOUBLE PRECISION array, dimension (N)\n*          The generated random numbers.\n*\n\n*  Further Details\n*  ===============\n*\n*  This routine calls the auxiliary routine DLARUV to generate random\n*  real numbers from a uniform (0,1) distribution, in batches of up to\n*  128 using vectorisable code. The Box-Muller method is used to\n*  transform numbers from a uniform to a normal distribution.\n*\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 3)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 3)", argc);
  rb_idist = argv[0];
  rb_iseed = argv[1];
  rb_n = argv[2];

  n = NUM2INT(rb_n);
  idist = NUM2INT(rb_idist);
  if (!NA_IsNArray(rb_iseed))
    rb_raise(rb_eArgError, "iseed (2th argument) must be NArray");
  if (NA_RANK(rb_iseed) != 1)
    rb_raise(rb_eArgError, "rank of iseed (2th argument) must be %d", 1);
  if (NA_SHAPE0(rb_iseed) != (4))
    rb_raise(rb_eRuntimeError, "shape 0 of iseed must be %d", 4);
  if (NA_TYPE(rb_iseed) != NA_LINT)
    rb_iseed = na_change_type(rb_iseed, NA_LINT);
  iseed = NA_PTR_TYPE(rb_iseed, integer*);
  {
    int shape[1];
    shape[0] = MAX(1,n);
    rb_x = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  x = NA_PTR_TYPE(rb_x, doublereal*);
  {
    int shape[1];
    shape[0] = 4;
    rb_iseed_out__ = na_make_object(NA_LINT, 1, shape, cNArray);
  }
  iseed_out__ = NA_PTR_TYPE(rb_iseed_out__, integer*);
  MEMCPY(iseed_out__, iseed, integer, NA_TOTAL(rb_iseed));
  rb_iseed = rb_iseed_out__;
  iseed = iseed_out__;

  dlarnv_(&idist, iseed, &n, x);

  return rb_ary_new3(2, rb_x, rb_iseed);
}

void
init_lapack_dlarnv(VALUE mLapack){
  rb_define_module_function(mLapack, "dlarnv", rb_dlarnv, -1);
}

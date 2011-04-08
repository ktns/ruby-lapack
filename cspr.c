#include "rb_lapack.h"

extern VOID cspr_(char *uplo, integer *n, complex *alpha, complex *x, integer *incx, complex *ap);

static VALUE sHelp, sUsage;

static VALUE
rblapack_cspr(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_uplo;
  char uplo; 
  VALUE rblapack_n;
  integer n; 
  VALUE rblapack_alpha;
  complex alpha; 
  VALUE rblapack_x;
  complex *x; 
  VALUE rblapack_incx;
  integer incx; 
  VALUE rblapack_ap;
  complex *ap; 
  VALUE rblapack_ap_out__;
  complex *ap_out__;


  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  ap = NumRu::Lapack.cspr( uplo, n, alpha, x, incx, ap, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE CSPR( UPLO, N, ALPHA, X, INCX, AP )\n\n*  Purpose\n*  =======\n*\n*  CSPR    performs the symmetric rank 1 operation\n*\n*     A := alpha*x*conjg( x' ) + A,\n*\n*  where alpha is a complex scalar, x is an n element vector and A is an\n*  n by n symmetric matrix, supplied in packed form.\n*\n\n*  Arguments\n*  ==========\n*\n*  UPLO     (input) CHARACTER*1\n*           On entry, UPLO specifies whether the upper or lower\n*           triangular part of the matrix A is supplied in the packed\n*           array AP as follows:\n*\n*              UPLO = 'U' or 'u'   The upper triangular part of A is\n*                                  supplied in AP.\n*\n*              UPLO = 'L' or 'l'   The lower triangular part of A is\n*                                  supplied in AP.\n*\n*           Unchanged on exit.\n*\n*  N        (input) INTEGER\n*           On entry, N specifies the order of the matrix A.\n*           N must be at least zero.\n*           Unchanged on exit.\n*\n*  ALPHA    (input) COMPLEX\n*           On entry, ALPHA specifies the scalar alpha.\n*           Unchanged on exit.\n*\n*  X        (input) COMPLEX array, dimension at least\n*           ( 1 + ( N - 1 )*abs( INCX ) ).\n*           Before entry, the incremented array X must contain the N-\n*           element vector x.\n*           Unchanged on exit.\n*\n*  INCX     (input) INTEGER\n*           On entry, INCX specifies the increment for the elements of\n*           X. INCX must not be zero.\n*           Unchanged on exit.\n*\n*  AP       (input/output) COMPLEX array, dimension at least\n*           ( ( N*( N + 1 ) )/2 ).\n*           Before entry, with  UPLO = 'U' or 'u', the array AP must\n*           contain the upper triangular part of the symmetric matrix\n*           packed sequentially, column by column, so that AP( 1 )\n*           contains a( 1, 1 ), AP( 2 ) and AP( 3 ) contain a( 1, 2 )\n*           and a( 2, 2 ) respectively, and so on. On exit, the array\n*           AP is overwritten by the upper triangular part of the\n*           updated matrix.\n*           Before entry, with UPLO = 'L' or 'l', the array AP must\n*           contain the lower triangular part of the symmetric matrix\n*           packed sequentially, column by column, so that AP( 1 )\n*           contains a( 1, 1 ), AP( 2 ) and AP( 3 ) contain a( 2, 1 )\n*           and a( 3, 1 ) respectively, and so on. On exit, the array\n*           AP is overwritten by the lower triangular part of the\n*           updated matrix.\n*           Note that the imaginary parts of the diagonal elements need\n*           not be set, they are assumed to be zero, and on exit they\n*           are set to zero.\n*\n\n* =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  ap = NumRu::Lapack.cspr( uplo, n, alpha, x, incx, ap, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 6)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 6)", argc);
  rblapack_uplo = argv[0];
  rblapack_n = argv[1];
  rblapack_alpha = argv[2];
  rblapack_x = argv[3];
  rblapack_incx = argv[4];
  rblapack_ap = argv[5];
  if (rb_options != Qnil) {
  }

  uplo = StringValueCStr(rblapack_uplo)[0];
  n = NUM2INT(rblapack_n);
  alpha.r = (real)NUM2DBL(rb_funcall(rblapack_alpha, rb_intern("real"), 0));
  alpha.i = (real)NUM2DBL(rb_funcall(rblapack_alpha, rb_intern("imag"), 0));
  incx = NUM2INT(rblapack_incx);
  if (!NA_IsNArray(rblapack_ap))
    rb_raise(rb_eArgError, "ap (6th argument) must be NArray");
  if (NA_RANK(rblapack_ap) != 1)
    rb_raise(rb_eArgError, "rank of ap (6th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_ap) != (( n*( n + 1 ) )/2))
    rb_raise(rb_eRuntimeError, "shape 0 of ap must be %d", ( n*( n + 1 ) )/2);
  if (NA_TYPE(rblapack_ap) != NA_SCOMPLEX)
    rblapack_ap = na_change_type(rblapack_ap, NA_SCOMPLEX);
  ap = NA_PTR_TYPE(rblapack_ap, complex*);
  if (!NA_IsNArray(rblapack_x))
    rb_raise(rb_eArgError, "x (4th argument) must be NArray");
  if (NA_RANK(rblapack_x) != 1)
    rb_raise(rb_eArgError, "rank of x (4th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_x) != (1 + ( n - 1 )*abs( incx )))
    rb_raise(rb_eRuntimeError, "shape 0 of x must be %d", 1 + ( n - 1 )*abs( incx ));
  if (NA_TYPE(rblapack_x) != NA_SCOMPLEX)
    rblapack_x = na_change_type(rblapack_x, NA_SCOMPLEX);
  x = NA_PTR_TYPE(rblapack_x, complex*);
  {
    int shape[1];
    shape[0] = ( n*( n + 1 ) )/2;
    rblapack_ap_out__ = na_make_object(NA_SCOMPLEX, 1, shape, cNArray);
  }
  ap_out__ = NA_PTR_TYPE(rblapack_ap_out__, complex*);
  MEMCPY(ap_out__, ap, complex, NA_TOTAL(rblapack_ap));
  rblapack_ap = rblapack_ap_out__;
  ap = ap_out__;

  cspr_(&uplo, &n, &alpha, x, &incx, ap);

  return rblapack_ap;
}

void
init_lapack_cspr(VALUE mLapack){
  rb_define_module_function(mLapack, "cspr", rblapack_cspr, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}

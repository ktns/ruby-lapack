#include "rb_lapack.h"

extern VOID zlaqhp_(char *uplo, integer *n, doublecomplex *ap, doublereal *s, doublereal *scond, doublereal *amax, char *equed);

static VALUE sHelp, sUsage;

static VALUE
rblapack_zlaqhp(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_uplo;
  char uplo; 
  VALUE rblapack_ap;
  doublecomplex *ap; 
  VALUE rblapack_s;
  doublereal *s; 
  VALUE rblapack_scond;
  doublereal scond; 
  VALUE rblapack_amax;
  doublereal amax; 
  VALUE rblapack_equed;
  char equed; 
  VALUE rblapack_ap_out__;
  doublecomplex *ap_out__;

  integer n;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  equed, ap = NumRu::Lapack.zlaqhp( uplo, ap, s, scond, amax, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE ZLAQHP( UPLO, N, AP, S, SCOND, AMAX, EQUED )\n\n*  Purpose\n*  =======\n*\n*  ZLAQHP equilibrates a Hermitian matrix A using the scaling factors\n*  in the vector S.\n*\n\n*  Arguments\n*  =========\n*\n*  UPLO    (input) CHARACTER*1\n*          Specifies whether the upper or lower triangular part of the\n*          Hermitian matrix A is stored.\n*          = 'U':  Upper triangular\n*          = 'L':  Lower triangular\n*\n*  N       (input) INTEGER\n*          The order of the matrix A.  N >= 0.\n*\n*  AP      (input/output) COMPLEX*16 array, dimension (N*(N+1)/2)\n*          On entry, the upper or lower triangle of the Hermitian matrix\n*          A, packed columnwise in a linear array.  The j-th column of A\n*          is stored in the array AP as follows:\n*          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;\n*          if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = A(i,j) for j<=i<=n.\n*\n*          On exit, the equilibrated matrix:  diag(S) * A * diag(S), in\n*          the same storage format as A.\n*\n*  S       (input) DOUBLE PRECISION array, dimension (N)\n*          The scale factors for A.\n*\n*  SCOND   (input) DOUBLE PRECISION\n*          Ratio of the smallest S(i) to the largest S(i).\n*\n*  AMAX    (input) DOUBLE PRECISION\n*          Absolute value of largest matrix entry.\n*\n*  EQUED   (output) CHARACTER*1\n*          Specifies whether or not equilibration was done.\n*          = 'N':  No equilibration.\n*          = 'Y':  Equilibration was done, i.e., A has been replaced by\n*                  diag(S) * A * diag(S).\n*\n*  Internal Parameters\n*  ===================\n*\n*  THRESH is a threshold value used to decide if scaling should be done\n*  based on the ratio of the scaling factors.  If SCOND < THRESH,\n*  scaling is done.\n*\n*  LARGE and SMALL are threshold values used to decide if scaling should\n*  be done based on the absolute size of the largest matrix element.\n*  If AMAX > LARGE or AMAX < SMALL, scaling is done.\n*\n\n*  =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  equed, ap = NumRu::Lapack.zlaqhp( uplo, ap, s, scond, amax, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 5)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 5)", argc);
  rblapack_uplo = argv[0];
  rblapack_ap = argv[1];
  rblapack_s = argv[2];
  rblapack_scond = argv[3];
  rblapack_amax = argv[4];
  if (rb_options != Qnil) {
  }

  scond = NUM2DBL(rblapack_scond);
  uplo = StringValueCStr(rblapack_uplo)[0];
  if (!NA_IsNArray(rblapack_s))
    rb_raise(rb_eArgError, "s (3th argument) must be NArray");
  if (NA_RANK(rblapack_s) != 1)
    rb_raise(rb_eArgError, "rank of s (3th argument) must be %d", 1);
  n = NA_SHAPE0(rblapack_s);
  if (NA_TYPE(rblapack_s) != NA_DFLOAT)
    rblapack_s = na_change_type(rblapack_s, NA_DFLOAT);
  s = NA_PTR_TYPE(rblapack_s, doublereal*);
  amax = NUM2DBL(rblapack_amax);
  if (!NA_IsNArray(rblapack_ap))
    rb_raise(rb_eArgError, "ap (2th argument) must be NArray");
  if (NA_RANK(rblapack_ap) != 1)
    rb_raise(rb_eArgError, "rank of ap (2th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_ap) != (n*(n+1)/2))
    rb_raise(rb_eRuntimeError, "shape 0 of ap must be %d", n*(n+1)/2);
  if (NA_TYPE(rblapack_ap) != NA_DCOMPLEX)
    rblapack_ap = na_change_type(rblapack_ap, NA_DCOMPLEX);
  ap = NA_PTR_TYPE(rblapack_ap, doublecomplex*);
  {
    int shape[1];
    shape[0] = n*(n+1)/2;
    rblapack_ap_out__ = na_make_object(NA_DCOMPLEX, 1, shape, cNArray);
  }
  ap_out__ = NA_PTR_TYPE(rblapack_ap_out__, doublecomplex*);
  MEMCPY(ap_out__, ap, doublecomplex, NA_TOTAL(rblapack_ap));
  rblapack_ap = rblapack_ap_out__;
  ap = ap_out__;

  zlaqhp_(&uplo, &n, ap, s, &scond, &amax, &equed);

  rblapack_equed = rb_str_new(&equed,1);
  return rb_ary_new3(2, rblapack_equed, rblapack_ap);
}

void
init_lapack_zlaqhp(VALUE mLapack){
  rb_define_module_function(mLapack, "zlaqhp", rblapack_zlaqhp, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}

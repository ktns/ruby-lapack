#include "rb_lapack.h"

extern VOID sppequ_(char *uplo, integer *n, real *ap, real *s, real *scond, real *amax, integer *info);

static VALUE sHelp, sUsage;

static VALUE
rblapack_sppequ(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_uplo;
  char uplo; 
  VALUE rblapack_ap;
  real *ap; 
  VALUE rblapack_s;
  real *s; 
  VALUE rblapack_scond;
  real scond; 
  VALUE rblapack_amax;
  real amax; 
  VALUE rblapack_info;
  integer info; 

  integer ldap;
  integer n;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  s, scond, amax, info = NumRu::Lapack.sppequ( uplo, ap, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE SPPEQU( UPLO, N, AP, S, SCOND, AMAX, INFO )\n\n*  Purpose\n*  =======\n*\n*  SPPEQU computes row and column scalings intended to equilibrate a\n*  symmetric positive definite matrix A in packed storage and reduce\n*  its condition number (with respect to the two-norm).  S contains the\n*  scale factors, S(i)=1/sqrt(A(i,i)), chosen so that the scaled matrix\n*  B with elements B(i,j)=S(i)*A(i,j)*S(j) has ones on the diagonal.\n*  This choice of S puts the condition number of B within a factor N of\n*  the smallest possible condition number over all possible diagonal\n*  scalings.\n*\n\n*  Arguments\n*  =========\n*\n*  UPLO    (input) CHARACTER*1\n*          = 'U':  Upper triangle of A is stored;\n*          = 'L':  Lower triangle of A is stored.\n*\n*  N       (input) INTEGER\n*          The order of the matrix A.  N >= 0.\n*\n*  AP      (input) REAL array, dimension (N*(N+1)/2)\n*          The upper or lower triangle of the symmetric matrix A, packed\n*          columnwise in a linear array.  The j-th column of A is stored\n*          in the array AP as follows:\n*          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;\n*          if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = A(i,j) for j<=i<=n.\n*\n*  S       (output) REAL array, dimension (N)\n*          If INFO = 0, S contains the scale factors for A.\n*\n*  SCOND   (output) REAL\n*          If INFO = 0, S contains the ratio of the smallest S(i) to\n*          the largest S(i).  If SCOND >= 0.1 and AMAX is neither too\n*          large nor too small, it is not worth scaling by S.\n*\n*  AMAX    (output) REAL\n*          Absolute value of largest matrix element.  If AMAX is very\n*          close to overflow or very close to underflow, the matrix\n*          should be scaled.\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit\n*          < 0:  if INFO = -i, the i-th argument had an illegal value\n*          > 0:  if INFO = i, the i-th diagonal element is nonpositive.\n*\n\n*  =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  s, scond, amax, info = NumRu::Lapack.sppequ( uplo, ap, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 2)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 2)", argc);
  rblapack_uplo = argv[0];
  rblapack_ap = argv[1];
  if (rb_options != Qnil) {
  }

  uplo = StringValueCStr(rblapack_uplo)[0];
  if (!NA_IsNArray(rblapack_ap))
    rb_raise(rb_eArgError, "ap (2th argument) must be NArray");
  if (NA_RANK(rblapack_ap) != 1)
    rb_raise(rb_eArgError, "rank of ap (2th argument) must be %d", 1);
  ldap = NA_SHAPE0(rblapack_ap);
  if (NA_TYPE(rblapack_ap) != NA_SFLOAT)
    rblapack_ap = na_change_type(rblapack_ap, NA_SFLOAT);
  ap = NA_PTR_TYPE(rblapack_ap, real*);
  n = ((int)sqrtf(ldap*8+1.0f)-1)/2;
  {
    int shape[1];
    shape[0] = n;
    rblapack_s = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  s = NA_PTR_TYPE(rblapack_s, real*);

  sppequ_(&uplo, &n, ap, s, &scond, &amax, &info);

  rblapack_scond = rb_float_new((double)scond);
  rblapack_amax = rb_float_new((double)amax);
  rblapack_info = INT2NUM(info);
  return rb_ary_new3(4, rblapack_s, rblapack_scond, rblapack_amax, rblapack_info);
}

void
init_lapack_sppequ(VALUE mLapack){
  rb_define_module_function(mLapack, "sppequ", rblapack_sppequ, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}

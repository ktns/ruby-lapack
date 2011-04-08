#include "rb_lapack.h"

extern VOID sopgtr_(char *uplo, integer *n, real *ap, real *tau, real *q, integer *ldq, real *work, integer *info);

static VALUE sHelp, sUsage;

static VALUE
rblapack_sopgtr(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_uplo;
  char uplo; 
  VALUE rblapack_ap;
  real *ap; 
  VALUE rblapack_tau;
  real *tau; 
  VALUE rblapack_q;
  real *q; 
  VALUE rblapack_info;
  integer info; 
  real *work;

  integer ldap;
  integer ldtau;
  integer ldq;
  integer n;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  q, info = NumRu::Lapack.sopgtr( uplo, ap, tau, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE SOPGTR( UPLO, N, AP, TAU, Q, LDQ, WORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  SOPGTR generates a real orthogonal matrix Q which is defined as the\n*  product of n-1 elementary reflectors H(i) of order n, as returned by\n*  SSPTRD using packed storage:\n*\n*  if UPLO = 'U', Q = H(n-1) . . . H(2) H(1),\n*\n*  if UPLO = 'L', Q = H(1) H(2) . . . H(n-1).\n*\n\n*  Arguments\n*  =========\n*\n*  UPLO    (input) CHARACTER*1\n*          = 'U': Upper triangular packed storage used in previous\n*                 call to SSPTRD;\n*          = 'L': Lower triangular packed storage used in previous\n*                 call to SSPTRD.\n*\n*  N       (input) INTEGER\n*          The order of the matrix Q. N >= 0.\n*\n*  AP      (input) REAL array, dimension (N*(N+1)/2)\n*          The vectors which define the elementary reflectors, as\n*          returned by SSPTRD.\n*\n*  TAU     (input) REAL array, dimension (N-1)\n*          TAU(i) must contain the scalar factor of the elementary\n*          reflector H(i), as returned by SSPTRD.\n*\n*  Q       (output) REAL array, dimension (LDQ,N)\n*          The N-by-N orthogonal matrix Q.\n*\n*  LDQ     (input) INTEGER\n*          The leading dimension of the array Q. LDQ >= max(1,N).\n*\n*  WORK    (workspace) REAL array, dimension (N-1)\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit\n*          < 0:  if INFO = -i, the i-th argument had an illegal value\n*\n\n*  =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  q, info = NumRu::Lapack.sopgtr( uplo, ap, tau, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 3)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 3)", argc);
  rblapack_uplo = argv[0];
  rblapack_ap = argv[1];
  rblapack_tau = argv[2];
  if (rb_options != Qnil) {
  }

  uplo = StringValueCStr(rblapack_uplo)[0];
  if (!NA_IsNArray(rblapack_tau))
    rb_raise(rb_eArgError, "tau (3th argument) must be NArray");
  if (NA_RANK(rblapack_tau) != 1)
    rb_raise(rb_eArgError, "rank of tau (3th argument) must be %d", 1);
  ldtau = NA_SHAPE0(rblapack_tau);
  if (NA_TYPE(rblapack_tau) != NA_SFLOAT)
    rblapack_tau = na_change_type(rblapack_tau, NA_SFLOAT);
  tau = NA_PTR_TYPE(rblapack_tau, real*);
  if (!NA_IsNArray(rblapack_ap))
    rb_raise(rb_eArgError, "ap (2th argument) must be NArray");
  if (NA_RANK(rblapack_ap) != 1)
    rb_raise(rb_eArgError, "rank of ap (2th argument) must be %d", 1);
  ldap = NA_SHAPE0(rblapack_ap);
  if (NA_TYPE(rblapack_ap) != NA_SFLOAT)
    rblapack_ap = na_change_type(rblapack_ap, NA_SFLOAT);
  ap = NA_PTR_TYPE(rblapack_ap, real*);
  n = ldtau+1;
  ldq = MAX(1,n);
  {
    int shape[2];
    shape[0] = ldq;
    shape[1] = n;
    rblapack_q = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  q = NA_PTR_TYPE(rblapack_q, real*);
  work = ALLOC_N(real, (n-1));

  sopgtr_(&uplo, &n, ap, tau, q, &ldq, work, &info);

  free(work);
  rblapack_info = INT2NUM(info);
  return rb_ary_new3(2, rblapack_q, rblapack_info);
}

void
init_lapack_sopgtr(VALUE mLapack){
  rb_define_module_function(mLapack, "sopgtr", rblapack_sopgtr, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}

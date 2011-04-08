#include "rb_lapack.h"

extern VOID sorg2r_(integer *m, integer *n, integer *k, real *a, integer *lda, real *tau, real *work, integer *info);

static VALUE sHelp, sUsage;

static VALUE
rblapack_sorg2r(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_m;
  integer m; 
  VALUE rblapack_a;
  real *a; 
  VALUE rblapack_tau;
  real *tau; 
  VALUE rblapack_info;
  integer info; 
  VALUE rblapack_a_out__;
  real *a_out__;
  real *work;

  integer lda;
  integer n;
  integer k;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  info, a = NumRu::Lapack.sorg2r( m, a, tau, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE SORG2R( M, N, K, A, LDA, TAU, WORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  SORG2R generates an m by n real matrix Q with orthonormal columns,\n*  which is defined as the first n columns of a product of k elementary\n*  reflectors of order m\n*\n*        Q  =  H(1) H(2) . . . H(k)\n*\n*  as returned by SGEQRF.\n*\n\n*  Arguments\n*  =========\n*\n*  M       (input) INTEGER\n*          The number of rows of the matrix Q. M >= 0.\n*\n*  N       (input) INTEGER\n*          The number of columns of the matrix Q. M >= N >= 0.\n*\n*  K       (input) INTEGER\n*          The number of elementary reflectors whose product defines the\n*          matrix Q. N >= K >= 0.\n*\n*  A       (input/output) REAL array, dimension (LDA,N)\n*          On entry, the i-th column must contain the vector which\n*          defines the elementary reflector H(i), for i = 1,2,...,k, as\n*          returned by SGEQRF in the first k columns of its array\n*          argument A.\n*          On exit, the m-by-n matrix Q.\n*\n*  LDA     (input) INTEGER\n*          The first dimension of the array A. LDA >= max(1,M).\n*\n*  TAU     (input) REAL array, dimension (K)\n*          TAU(i) must contain the scalar factor of the elementary\n*          reflector H(i), as returned by SGEQRF.\n*\n*  WORK    (workspace) REAL array, dimension (N)\n*\n*  INFO    (output) INTEGER\n*          = 0: successful exit\n*          < 0: if INFO = -i, the i-th argument has an illegal value\n*\n\n*  =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  info, a = NumRu::Lapack.sorg2r( m, a, tau, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 3)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 3)", argc);
  rblapack_m = argv[0];
  rblapack_a = argv[1];
  rblapack_tau = argv[2];
  if (rb_options != Qnil) {
  }

  if (!NA_IsNArray(rblapack_a))
    rb_raise(rb_eArgError, "a (2th argument) must be NArray");
  if (NA_RANK(rblapack_a) != 2)
    rb_raise(rb_eArgError, "rank of a (2th argument) must be %d", 2);
  n = NA_SHAPE1(rblapack_a);
  lda = NA_SHAPE0(rblapack_a);
  if (NA_TYPE(rblapack_a) != NA_SFLOAT)
    rblapack_a = na_change_type(rblapack_a, NA_SFLOAT);
  a = NA_PTR_TYPE(rblapack_a, real*);
  m = NUM2INT(rblapack_m);
  if (!NA_IsNArray(rblapack_tau))
    rb_raise(rb_eArgError, "tau (3th argument) must be NArray");
  if (NA_RANK(rblapack_tau) != 1)
    rb_raise(rb_eArgError, "rank of tau (3th argument) must be %d", 1);
  k = NA_SHAPE0(rblapack_tau);
  if (NA_TYPE(rblapack_tau) != NA_SFLOAT)
    rblapack_tau = na_change_type(rblapack_tau, NA_SFLOAT);
  tau = NA_PTR_TYPE(rblapack_tau, real*);
  {
    int shape[2];
    shape[0] = lda;
    shape[1] = n;
    rblapack_a_out__ = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  a_out__ = NA_PTR_TYPE(rblapack_a_out__, real*);
  MEMCPY(a_out__, a, real, NA_TOTAL(rblapack_a));
  rblapack_a = rblapack_a_out__;
  a = a_out__;
  work = ALLOC_N(real, (n));

  sorg2r_(&m, &n, &k, a, &lda, tau, work, &info);

  free(work);
  rblapack_info = INT2NUM(info);
  return rb_ary_new3(2, rblapack_info, rblapack_a);
}

void
init_lapack_sorg2r(VALUE mLapack){
  rb_define_module_function(mLapack, "sorg2r", rblapack_sorg2r, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}

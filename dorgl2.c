#include "rb_lapack.h"

extern VOID dorgl2_(integer *m, integer *n, integer *k, doublereal *a, integer *lda, doublereal *tau, doublereal *work, integer *info);

static VALUE sHelp, sUsage;

static VALUE
rblapack_dorgl2(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_a;
  doublereal *a; 
  VALUE rblapack_tau;
  doublereal *tau; 
  VALUE rblapack_info;
  integer info; 
  VALUE rblapack_a_out__;
  doublereal *a_out__;
  doublereal *work;

  integer lda;
  integer n;
  integer k;
  integer m;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  info, a = NumRu::Lapack.dorgl2( a, tau, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE DORGL2( M, N, K, A, LDA, TAU, WORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  DORGL2 generates an m by n real matrix Q with orthonormal rows,\n*  which is defined as the first m rows of a product of k elementary\n*  reflectors of order n\n*\n*        Q  =  H(k) . . . H(2) H(1)\n*\n*  as returned by DGELQF.\n*\n\n*  Arguments\n*  =========\n*\n*  M       (input) INTEGER\n*          The number of rows of the matrix Q. M >= 0.\n*\n*  N       (input) INTEGER\n*          The number of columns of the matrix Q. N >= M.\n*\n*  K       (input) INTEGER\n*          The number of elementary reflectors whose product defines the\n*          matrix Q. M >= K >= 0.\n*\n*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)\n*          On entry, the i-th row must contain the vector which defines\n*          the elementary reflector H(i), for i = 1,2,...,k, as returned\n*          by DGELQF in the first k rows of its array argument A.\n*          On exit, the m-by-n matrix Q.\n*\n*  LDA     (input) INTEGER\n*          The first dimension of the array A. LDA >= max(1,M).\n*\n*  TAU     (input) DOUBLE PRECISION array, dimension (K)\n*          TAU(i) must contain the scalar factor of the elementary\n*          reflector H(i), as returned by DGELQF.\n*\n*  WORK    (workspace) DOUBLE PRECISION array, dimension (M)\n*\n*  INFO    (output) INTEGER\n*          = 0: successful exit\n*          < 0: if INFO = -i, the i-th argument has an illegal value\n*\n\n*  =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  info, a = NumRu::Lapack.dorgl2( a, tau, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 2)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 2)", argc);
  rblapack_a = argv[0];
  rblapack_tau = argv[1];
  if (rb_options != Qnil) {
  }

  if (!NA_IsNArray(rblapack_a))
    rb_raise(rb_eArgError, "a (1th argument) must be NArray");
  if (NA_RANK(rblapack_a) != 2)
    rb_raise(rb_eArgError, "rank of a (1th argument) must be %d", 2);
  n = NA_SHAPE1(rblapack_a);
  lda = NA_SHAPE0(rblapack_a);
  if (NA_TYPE(rblapack_a) != NA_DFLOAT)
    rblapack_a = na_change_type(rblapack_a, NA_DFLOAT);
  a = NA_PTR_TYPE(rblapack_a, doublereal*);
  if (!NA_IsNArray(rblapack_tau))
    rb_raise(rb_eArgError, "tau (2th argument) must be NArray");
  if (NA_RANK(rblapack_tau) != 1)
    rb_raise(rb_eArgError, "rank of tau (2th argument) must be %d", 1);
  k = NA_SHAPE0(rblapack_tau);
  if (NA_TYPE(rblapack_tau) != NA_DFLOAT)
    rblapack_tau = na_change_type(rblapack_tau, NA_DFLOAT);
  tau = NA_PTR_TYPE(rblapack_tau, doublereal*);
  m = lda;
  {
    int shape[2];
    shape[0] = lda;
    shape[1] = n;
    rblapack_a_out__ = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  a_out__ = NA_PTR_TYPE(rblapack_a_out__, doublereal*);
  MEMCPY(a_out__, a, doublereal, NA_TOTAL(rblapack_a));
  rblapack_a = rblapack_a_out__;
  a = a_out__;
  work = ALLOC_N(doublereal, (m));

  dorgl2_(&m, &n, &k, a, &lda, tau, work, &info);

  free(work);
  rblapack_info = INT2NUM(info);
  return rb_ary_new3(2, rblapack_info, rblapack_a);
}

void
init_lapack_dorgl2(VALUE mLapack){
  rb_define_module_function(mLapack, "dorgl2", rblapack_dorgl2, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}

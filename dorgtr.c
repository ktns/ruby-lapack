#include "rb_lapack.h"

extern VOID dorgtr_(char *uplo, integer *n, doublereal *a, integer *lda, doublereal *tau, doublereal *work, integer *lwork, integer *info);

static VALUE sHelp, sUsage;

static VALUE
rblapack_dorgtr(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_uplo;
  char uplo; 
  VALUE rblapack_a;
  doublereal *a; 
  VALUE rblapack_tau;
  doublereal *tau; 
  VALUE rblapack_lwork;
  integer lwork; 
  VALUE rblapack_work;
  doublereal *work; 
  VALUE rblapack_info;
  integer info; 
  VALUE rblapack_a_out__;
  doublereal *a_out__;

  integer lda;
  integer n;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  work, info, a = NumRu::Lapack.dorgtr( uplo, a, tau, lwork, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE DORGTR( UPLO, N, A, LDA, TAU, WORK, LWORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  DORGTR generates a real orthogonal matrix Q which is defined as the\n*  product of n-1 elementary reflectors of order N, as returned by\n*  DSYTRD:\n*\n*  if UPLO = 'U', Q = H(n-1) . . . H(2) H(1),\n*\n*  if UPLO = 'L', Q = H(1) H(2) . . . H(n-1).\n*\n\n*  Arguments\n*  =========\n*\n*  UPLO    (input) CHARACTER*1\n*          = 'U': Upper triangle of A contains elementary reflectors\n*                 from DSYTRD;\n*          = 'L': Lower triangle of A contains elementary reflectors\n*                 from DSYTRD.\n*\n*  N       (input) INTEGER\n*          The order of the matrix Q. N >= 0.\n*\n*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)\n*          On entry, the vectors which define the elementary reflectors,\n*          as returned by DSYTRD.\n*          On exit, the N-by-N orthogonal matrix Q.\n*\n*  LDA     (input) INTEGER\n*          The leading dimension of the array A. LDA >= max(1,N).\n*\n*  TAU     (input) DOUBLE PRECISION array, dimension (N-1)\n*          TAU(i) must contain the scalar factor of the elementary\n*          reflector H(i), as returned by DSYTRD.\n*\n*  WORK    (workspace/output) DOUBLE PRECISION array, dimension (MAX(1,LWORK))\n*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.\n*\n*  LWORK   (input) INTEGER\n*          The dimension of the array WORK. LWORK >= max(1,N-1).\n*          For optimum performance LWORK >= (N-1)*NB, where NB is\n*          the optimal blocksize.\n*\n*          If LWORK = -1, then a workspace query is assumed; the routine\n*          only calculates the optimal size of the WORK array, returns\n*          this value as the first entry of the WORK array, and no error\n*          message related to LWORK is issued by XERBLA.\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit\n*          < 0:  if INFO = -i, the i-th argument had an illegal value\n*\n\n*  =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  work, info, a = NumRu::Lapack.dorgtr( uplo, a, tau, lwork, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 4)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 4)", argc);
  rblapack_uplo = argv[0];
  rblapack_a = argv[1];
  rblapack_tau = argv[2];
  rblapack_lwork = argv[3];
  if (rb_options != Qnil) {
  }

  if (!NA_IsNArray(rblapack_a))
    rb_raise(rb_eArgError, "a (2th argument) must be NArray");
  if (NA_RANK(rblapack_a) != 2)
    rb_raise(rb_eArgError, "rank of a (2th argument) must be %d", 2);
  n = NA_SHAPE1(rblapack_a);
  lda = NA_SHAPE0(rblapack_a);
  if (NA_TYPE(rblapack_a) != NA_DFLOAT)
    rblapack_a = na_change_type(rblapack_a, NA_DFLOAT);
  a = NA_PTR_TYPE(rblapack_a, doublereal*);
  lwork = NUM2INT(rblapack_lwork);
  uplo = StringValueCStr(rblapack_uplo)[0];
  if (!NA_IsNArray(rblapack_tau))
    rb_raise(rb_eArgError, "tau (3th argument) must be NArray");
  if (NA_RANK(rblapack_tau) != 1)
    rb_raise(rb_eArgError, "rank of tau (3th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_tau) != (n-1))
    rb_raise(rb_eRuntimeError, "shape 0 of tau must be %d", n-1);
  if (NA_TYPE(rblapack_tau) != NA_DFLOAT)
    rblapack_tau = na_change_type(rblapack_tau, NA_DFLOAT);
  tau = NA_PTR_TYPE(rblapack_tau, doublereal*);
  {
    int shape[1];
    shape[0] = MAX(1,lwork);
    rblapack_work = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  work = NA_PTR_TYPE(rblapack_work, doublereal*);
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

  dorgtr_(&uplo, &n, a, &lda, tau, work, &lwork, &info);

  rblapack_info = INT2NUM(info);
  return rb_ary_new3(3, rblapack_work, rblapack_info, rblapack_a);
}

void
init_lapack_dorgtr(VALUE mLapack){
  rb_define_module_function(mLapack, "dorgtr", rblapack_dorgtr, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}

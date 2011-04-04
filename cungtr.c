#include "rb_lapack.h"

extern VOID cungtr_(char *uplo, integer *n, complex *a, integer *lda, complex *tau, complex *work, integer *lwork, integer *info);

static VALUE
rb_cungtr(int argc, VALUE *argv, VALUE self){
  VALUE rb_uplo;
  char uplo; 
  VALUE rb_a;
  complex *a; 
  VALUE rb_tau;
  complex *tau; 
  VALUE rb_lwork;
  integer lwork; 
  VALUE rb_work;
  complex *work; 
  VALUE rb_info;
  integer info; 
  VALUE rb_a_out__;
  complex *a_out__;

  integer lda;
  integer n;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  work, info, a = NumRu::Lapack.cungtr( uplo, a, tau, lwork)\n    or\n  NumRu::Lapack.cungtr  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE CUNGTR( UPLO, N, A, LDA, TAU, WORK, LWORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  CUNGTR generates a complex unitary matrix Q which is defined as the\n*  product of n-1 elementary reflectors of order N, as returned by\n*  CHETRD:\n*\n*  if UPLO = 'U', Q = H(n-1) . . . H(2) H(1),\n*\n*  if UPLO = 'L', Q = H(1) H(2) . . . H(n-1).\n*\n\n*  Arguments\n*  =========\n*\n*  UPLO    (input) CHARACTER*1\n*          = 'U': Upper triangle of A contains elementary reflectors\n*                 from CHETRD;\n*          = 'L': Lower triangle of A contains elementary reflectors\n*                 from CHETRD.\n*\n*  N       (input) INTEGER\n*          The order of the matrix Q. N >= 0.\n*\n*  A       (input/output) COMPLEX array, dimension (LDA,N)\n*          On entry, the vectors which define the elementary reflectors,\n*          as returned by CHETRD.\n*          On exit, the N-by-N unitary matrix Q.\n*\n*  LDA     (input) INTEGER\n*          The leading dimension of the array A. LDA >= N.\n*\n*  TAU     (input) COMPLEX array, dimension (N-1)\n*          TAU(i) must contain the scalar factor of the elementary\n*          reflector H(i), as returned by CHETRD.\n*\n*  WORK    (workspace/output) COMPLEX array, dimension (MAX(1,LWORK))\n*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.\n*\n*  LWORK   (input) INTEGER\n*          The dimension of the array WORK. LWORK >= N-1.\n*          For optimum performance LWORK >= (N-1)*NB, where NB is\n*          the optimal blocksize.\n*\n*          If LWORK = -1, then a workspace query is assumed; the routine\n*          only calculates the optimal size of the WORK array, returns\n*          this value as the first entry of the WORK array, and no error\n*          message related to LWORK is issued by XERBLA.\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit\n*          < 0:  if INFO = -i, the i-th argument had an illegal value\n*\n\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 4)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 4)", argc);
  rb_uplo = argv[0];
  rb_a = argv[1];
  rb_tau = argv[2];
  rb_lwork = argv[3];

  if (!NA_IsNArray(rb_a))
    rb_raise(rb_eArgError, "a (2th argument) must be NArray");
  if (NA_RANK(rb_a) != 2)
    rb_raise(rb_eArgError, "rank of a (2th argument) must be %d", 2);
  n = NA_SHAPE1(rb_a);
  lda = NA_SHAPE0(rb_a);
  if (NA_TYPE(rb_a) != NA_SCOMPLEX)
    rb_a = na_change_type(rb_a, NA_SCOMPLEX);
  a = NA_PTR_TYPE(rb_a, complex*);
  lwork = NUM2INT(rb_lwork);
  uplo = StringValueCStr(rb_uplo)[0];
  if (!NA_IsNArray(rb_tau))
    rb_raise(rb_eArgError, "tau (3th argument) must be NArray");
  if (NA_RANK(rb_tau) != 1)
    rb_raise(rb_eArgError, "rank of tau (3th argument) must be %d", 1);
  if (NA_SHAPE0(rb_tau) != (n-1))
    rb_raise(rb_eRuntimeError, "shape 0 of tau must be %d", n-1);
  if (NA_TYPE(rb_tau) != NA_SCOMPLEX)
    rb_tau = na_change_type(rb_tau, NA_SCOMPLEX);
  tau = NA_PTR_TYPE(rb_tau, complex*);
  {
    int shape[1];
    shape[0] = MAX(1,lwork);
    rb_work = na_make_object(NA_SCOMPLEX, 1, shape, cNArray);
  }
  work = NA_PTR_TYPE(rb_work, complex*);
  {
    int shape[2];
    shape[0] = lda;
    shape[1] = n;
    rb_a_out__ = na_make_object(NA_SCOMPLEX, 2, shape, cNArray);
  }
  a_out__ = NA_PTR_TYPE(rb_a_out__, complex*);
  MEMCPY(a_out__, a, complex, NA_TOTAL(rb_a));
  rb_a = rb_a_out__;
  a = a_out__;

  cungtr_(&uplo, &n, a, &lda, tau, work, &lwork, &info);

  rb_info = INT2NUM(info);
  return rb_ary_new3(3, rb_work, rb_info, rb_a);
}

void
init_lapack_cungtr(VALUE mLapack){
  rb_define_module_function(mLapack, "cungtr", rb_cungtr, -1);
}

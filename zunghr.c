#include "rb_lapack.h"

extern VOID zunghr_(integer *n, integer *ilo, integer *ihi, doublecomplex *a, integer *lda, doublecomplex *tau, doublecomplex *work, integer *lwork, integer *info);

static VALUE
rb_zunghr(int argc, VALUE *argv, VALUE self){
  VALUE rb_ilo;
  integer ilo; 
  VALUE rb_ihi;
  integer ihi; 
  VALUE rb_a;
  doublecomplex *a; 
  VALUE rb_tau;
  doublecomplex *tau; 
  VALUE rb_lwork;
  integer lwork; 
  VALUE rb_work;
  doublecomplex *work; 
  VALUE rb_info;
  integer info; 
  VALUE rb_a_out__;
  doublecomplex *a_out__;

  integer lda;
  integer n;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  work, info, a = NumRu::Lapack.zunghr( ilo, ihi, a, tau, lwork)\n    or\n  NumRu::Lapack.zunghr  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE ZUNGHR( N, ILO, IHI, A, LDA, TAU, WORK, LWORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  ZUNGHR generates a complex unitary matrix Q which is defined as the\n*  product of IHI-ILO elementary reflectors of order N, as returned by\n*  ZGEHRD:\n*\n*  Q = H(ilo) H(ilo+1) . . . H(ihi-1).\n*\n\n*  Arguments\n*  =========\n*\n*  N       (input) INTEGER\n*          The order of the matrix Q. N >= 0.\n*\n*  ILO     (input) INTEGER\n*  IHI     (input) INTEGER\n*          ILO and IHI must have the same values as in the previous call\n*          of ZGEHRD. Q is equal to the unit matrix except in the\n*          submatrix Q(ilo+1:ihi,ilo+1:ihi).\n*          1 <= ILO <= IHI <= N, if N > 0; ILO=1 and IHI=0, if N=0.\n*\n*  A       (input/output) COMPLEX*16 array, dimension (LDA,N)\n*          On entry, the vectors which define the elementary reflectors,\n*          as returned by ZGEHRD.\n*          On exit, the N-by-N unitary matrix Q.\n*\n*  LDA     (input) INTEGER\n*          The leading dimension of the array A. LDA >= max(1,N).\n*\n*  TAU     (input) COMPLEX*16 array, dimension (N-1)\n*          TAU(i) must contain the scalar factor of the elementary\n*          reflector H(i), as returned by ZGEHRD.\n*\n*  WORK    (workspace/output) COMPLEX*16 array, dimension (MAX(1,LWORK))\n*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.\n*\n*  LWORK   (input) INTEGER\n*          The dimension of the array WORK. LWORK >= IHI-ILO.\n*          For optimum performance LWORK >= (IHI-ILO)*NB, where NB is\n*          the optimal blocksize.\n*\n*          If LWORK = -1, then a workspace query is assumed; the routine\n*          only calculates the optimal size of the WORK array, returns\n*          this value as the first entry of the WORK array, and no error\n*          message related to LWORK is issued by XERBLA.\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit\n*          < 0:  if INFO = -i, the i-th argument had an illegal value\n*\n\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 5)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 5)", argc);
  rb_ilo = argv[0];
  rb_ihi = argv[1];
  rb_a = argv[2];
  rb_tau = argv[3];
  rb_lwork = argv[4];

  if (!NA_IsNArray(rb_a))
    rb_raise(rb_eArgError, "a (3th argument) must be NArray");
  if (NA_RANK(rb_a) != 2)
    rb_raise(rb_eArgError, "rank of a (3th argument) must be %d", 2);
  n = NA_SHAPE1(rb_a);
  lda = NA_SHAPE0(rb_a);
  if (NA_TYPE(rb_a) != NA_DCOMPLEX)
    rb_a = na_change_type(rb_a, NA_DCOMPLEX);
  a = NA_PTR_TYPE(rb_a, doublecomplex*);
  ilo = NUM2INT(rb_ilo);
  ihi = NUM2INT(rb_ihi);
  lwork = NUM2INT(rb_lwork);
  if (!NA_IsNArray(rb_tau))
    rb_raise(rb_eArgError, "tau (4th argument) must be NArray");
  if (NA_RANK(rb_tau) != 1)
    rb_raise(rb_eArgError, "rank of tau (4th argument) must be %d", 1);
  if (NA_SHAPE0(rb_tau) != (n-1))
    rb_raise(rb_eRuntimeError, "shape 0 of tau must be %d", n-1);
  if (NA_TYPE(rb_tau) != NA_DCOMPLEX)
    rb_tau = na_change_type(rb_tau, NA_DCOMPLEX);
  tau = NA_PTR_TYPE(rb_tau, doublecomplex*);
  {
    int shape[1];
    shape[0] = MAX(1,lwork);
    rb_work = na_make_object(NA_DCOMPLEX, 1, shape, cNArray);
  }
  work = NA_PTR_TYPE(rb_work, doublecomplex*);
  {
    int shape[2];
    shape[0] = lda;
    shape[1] = n;
    rb_a_out__ = na_make_object(NA_DCOMPLEX, 2, shape, cNArray);
  }
  a_out__ = NA_PTR_TYPE(rb_a_out__, doublecomplex*);
  MEMCPY(a_out__, a, doublecomplex, NA_TOTAL(rb_a));
  rb_a = rb_a_out__;
  a = a_out__;

  zunghr_(&n, &ilo, &ihi, a, &lda, tau, work, &lwork, &info);

  rb_info = INT2NUM(info);
  return rb_ary_new3(3, rb_work, rb_info, rb_a);
}

void
init_lapack_zunghr(VALUE mLapack){
  rb_define_module_function(mLapack, "zunghr", rb_zunghr, -1);
}

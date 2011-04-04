#include "rb_lapack.h"

extern VOID cungr2_(integer *m, integer *n, integer *k, complex *a, integer *lda, complex *tau, complex *work, integer *info);

static VALUE
rb_cungr2(int argc, VALUE *argv, VALUE self){
  VALUE rb_a;
  complex *a; 
  VALUE rb_tau;
  complex *tau; 
  VALUE rb_info;
  integer info; 
  VALUE rb_a_out__;
  complex *a_out__;
  complex *work;

  integer lda;
  integer n;
  integer k;
  integer m;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  info, a = NumRu::Lapack.cungr2( a, tau)\n    or\n  NumRu::Lapack.cungr2  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE CUNGR2( M, N, K, A, LDA, TAU, WORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  CUNGR2 generates an m by n complex matrix Q with orthonormal rows,\n*  which is defined as the last m rows of a product of k elementary\n*  reflectors of order n\n*\n*        Q  =  H(1)' H(2)' . . . H(k)'\n*\n*  as returned by CGERQF.\n*\n\n*  Arguments\n*  =========\n*\n*  M       (input) INTEGER\n*          The number of rows of the matrix Q. M >= 0.\n*\n*  N       (input) INTEGER\n*          The number of columns of the matrix Q. N >= M.\n*\n*  K       (input) INTEGER\n*          The number of elementary reflectors whose product defines the\n*          matrix Q. M >= K >= 0.\n*\n*  A       (input/output) COMPLEX array, dimension (LDA,N)\n*          On entry, the (m-k+i)-th row must contain the vector which\n*          defines the elementary reflector H(i), for i = 1,2,...,k, as\n*          returned by CGERQF in the last k rows of its array argument\n*          A.\n*          On exit, the m-by-n matrix Q.\n*\n*  LDA     (input) INTEGER\n*          The first dimension of the array A. LDA >= max(1,M).\n*\n*  TAU     (input) COMPLEX array, dimension (K)\n*          TAU(i) must contain the scalar factor of the elementary\n*          reflector H(i), as returned by CGERQF.\n*\n*  WORK    (workspace) COMPLEX array, dimension (M)\n*\n*  INFO    (output) INTEGER\n*          = 0: successful exit\n*          < 0: if INFO = -i, the i-th argument has an illegal value\n*\n\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 2)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 2)", argc);
  rb_a = argv[0];
  rb_tau = argv[1];

  if (!NA_IsNArray(rb_a))
    rb_raise(rb_eArgError, "a (1th argument) must be NArray");
  if (NA_RANK(rb_a) != 2)
    rb_raise(rb_eArgError, "rank of a (1th argument) must be %d", 2);
  n = NA_SHAPE1(rb_a);
  lda = NA_SHAPE0(rb_a);
  if (NA_TYPE(rb_a) != NA_SCOMPLEX)
    rb_a = na_change_type(rb_a, NA_SCOMPLEX);
  a = NA_PTR_TYPE(rb_a, complex*);
  if (!NA_IsNArray(rb_tau))
    rb_raise(rb_eArgError, "tau (2th argument) must be NArray");
  if (NA_RANK(rb_tau) != 1)
    rb_raise(rb_eArgError, "rank of tau (2th argument) must be %d", 1);
  k = NA_SHAPE0(rb_tau);
  if (NA_TYPE(rb_tau) != NA_SCOMPLEX)
    rb_tau = na_change_type(rb_tau, NA_SCOMPLEX);
  tau = NA_PTR_TYPE(rb_tau, complex*);
  m = lda;
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
  work = ALLOC_N(complex, (m));

  cungr2_(&m, &n, &k, a, &lda, tau, work, &info);

  free(work);
  rb_info = INT2NUM(info);
  return rb_ary_new3(2, rb_info, rb_a);
}

void
init_lapack_cungr2(VALUE mLapack){
  rb_define_module_function(mLapack, "cungr2", rb_cungr2, -1);
}

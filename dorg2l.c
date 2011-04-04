#include "rb_lapack.h"

extern VOID dorg2l_(integer *m, integer *n, integer *k, doublereal *a, integer *lda, doublereal *tau, doublereal *work, integer *info);

static VALUE
rb_dorg2l(int argc, VALUE *argv, VALUE self){
  VALUE rb_m;
  integer m; 
  VALUE rb_a;
  doublereal *a; 
  VALUE rb_tau;
  doublereal *tau; 
  VALUE rb_info;
  integer info; 
  VALUE rb_a_out__;
  doublereal *a_out__;
  doublereal *work;

  integer lda;
  integer n;
  integer k;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  info, a = NumRu::Lapack.dorg2l( m, a, tau)\n    or\n  NumRu::Lapack.dorg2l  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE DORG2L( M, N, K, A, LDA, TAU, WORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  DORG2L generates an m by n real matrix Q with orthonormal columns,\n*  which is defined as the last n columns of a product of k elementary\n*  reflectors of order m\n*\n*        Q  =  H(k) . . . H(2) H(1)\n*\n*  as returned by DGEQLF.\n*\n\n*  Arguments\n*  =========\n*\n*  M       (input) INTEGER\n*          The number of rows of the matrix Q. M >= 0.\n*\n*  N       (input) INTEGER\n*          The number of columns of the matrix Q. M >= N >= 0.\n*\n*  K       (input) INTEGER\n*          The number of elementary reflectors whose product defines the\n*          matrix Q. N >= K >= 0.\n*\n*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)\n*          On entry, the (n-k+i)-th column must contain the vector which\n*          defines the elementary reflector H(i), for i = 1,2,...,k, as\n*          returned by DGEQLF in the last k columns of its array\n*          argument A.\n*          On exit, the m by n matrix Q.\n*\n*  LDA     (input) INTEGER\n*          The first dimension of the array A. LDA >= max(1,M).\n*\n*  TAU     (input) DOUBLE PRECISION array, dimension (K)\n*          TAU(i) must contain the scalar factor of the elementary\n*          reflector H(i), as returned by DGEQLF.\n*\n*  WORK    (workspace) DOUBLE PRECISION array, dimension (N)\n*\n*  INFO    (output) INTEGER\n*          = 0: successful exit\n*          < 0: if INFO = -i, the i-th argument has an illegal value\n*\n\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 3)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 3)", argc);
  rb_m = argv[0];
  rb_a = argv[1];
  rb_tau = argv[2];

  if (!NA_IsNArray(rb_a))
    rb_raise(rb_eArgError, "a (2th argument) must be NArray");
  if (NA_RANK(rb_a) != 2)
    rb_raise(rb_eArgError, "rank of a (2th argument) must be %d", 2);
  n = NA_SHAPE1(rb_a);
  lda = NA_SHAPE0(rb_a);
  if (NA_TYPE(rb_a) != NA_DFLOAT)
    rb_a = na_change_type(rb_a, NA_DFLOAT);
  a = NA_PTR_TYPE(rb_a, doublereal*);
  m = NUM2INT(rb_m);
  if (!NA_IsNArray(rb_tau))
    rb_raise(rb_eArgError, "tau (3th argument) must be NArray");
  if (NA_RANK(rb_tau) != 1)
    rb_raise(rb_eArgError, "rank of tau (3th argument) must be %d", 1);
  k = NA_SHAPE0(rb_tau);
  if (NA_TYPE(rb_tau) != NA_DFLOAT)
    rb_tau = na_change_type(rb_tau, NA_DFLOAT);
  tau = NA_PTR_TYPE(rb_tau, doublereal*);
  {
    int shape[2];
    shape[0] = lda;
    shape[1] = n;
    rb_a_out__ = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  a_out__ = NA_PTR_TYPE(rb_a_out__, doublereal*);
  MEMCPY(a_out__, a, doublereal, NA_TOTAL(rb_a));
  rb_a = rb_a_out__;
  a = a_out__;
  work = ALLOC_N(doublereal, (n));

  dorg2l_(&m, &n, &k, a, &lda, tau, work, &info);

  free(work);
  rb_info = INT2NUM(info);
  return rb_ary_new3(2, rb_info, rb_a);
}

void
init_lapack_dorg2l(VALUE mLapack){
  rb_define_module_function(mLapack, "dorg2l", rb_dorg2l, -1);
}

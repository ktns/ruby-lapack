#include "rb_lapack.h"

extern VOID dtzrqf_(integer *m, integer *n, doublereal *a, integer *lda, doublereal *tau, integer *info);

static VALUE
rb_dtzrqf(int argc, VALUE *argv, VALUE self){
  VALUE rb_a;
  doublereal *a; 
  VALUE rb_tau;
  doublereal *tau; 
  VALUE rb_info;
  integer info; 
  VALUE rb_a_out__;
  doublereal *a_out__;

  integer lda;
  integer n;
  integer m;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  tau, info, a = NumRu::Lapack.dtzrqf( a)\n    or\n  NumRu::Lapack.dtzrqf  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE DTZRQF( M, N, A, LDA, TAU, INFO )\n\n*  Purpose\n*  =======\n*\n*  This routine is deprecated and has been replaced by routine DTZRZF.\n*\n*  DTZRQF reduces the M-by-N ( M<=N ) real upper trapezoidal matrix A\n*  to upper triangular form by means of orthogonal transformations.\n*\n*  The upper trapezoidal matrix A is factored as\n*\n*     A = ( R  0 ) * Z,\n*\n*  where Z is an N-by-N orthogonal matrix and R is an M-by-M upper\n*  triangular matrix.\n*\n\n*  Arguments\n*  =========\n*\n*  M       (input) INTEGER\n*          The number of rows of the matrix A.  M >= 0.\n*\n*  N       (input) INTEGER\n*          The number of columns of the matrix A.  N >= M.\n*\n*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)\n*          On entry, the leading M-by-N upper trapezoidal part of the\n*          array A must contain the matrix to be factorized.\n*          On exit, the leading M-by-M upper triangular part of A\n*          contains the upper triangular matrix R, and elements M+1 to\n*          N of the first M rows of A, with the array TAU, represent the\n*          orthogonal matrix Z as a product of M elementary reflectors.\n*\n*  LDA     (input) INTEGER\n*          The leading dimension of the array A.  LDA >= max(1,M).\n*\n*  TAU     (output) DOUBLE PRECISION array, dimension (M)\n*          The scalar factors of the elementary reflectors.\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit\n*          < 0:  if INFO = -i, the i-th argument had an illegal value\n*\n\n*  Further Details\n*  ===============\n*\n*  The factorization is obtained by Householder's method.  The kth\n*  transformation matrix, Z( k ), which is used to introduce zeros into\n*  the ( m - k + 1 )th row of A, is given in the form\n*\n*     Z( k ) = ( I     0   ),\n*              ( 0  T( k ) )\n*\n*  where\n*\n*     T( k ) = I - tau*u( k )*u( k )',   u( k ) = (   1    ),\n*                                                 (   0    )\n*                                                 ( z( k ) )\n*\n*  tau is a scalar and z( k ) is an ( n - m ) element vector.\n*  tau and z( k ) are chosen to annihilate the elements of the kth row\n*  of X.\n*\n*  The scalar tau is returned in the kth element of TAU and the vector\n*  u( k ) in the kth row of A, such that the elements of z( k ) are\n*  in  a( k, m + 1 ), ..., a( k, n ). The elements of R are returned in\n*  the upper triangular part of A.\n*\n*  Z is given by\n*\n*     Z =  Z( 1 ) * Z( 2 ) * ... * Z( m ).\n*\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 1)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 1)", argc);
  rb_a = argv[0];

  if (!NA_IsNArray(rb_a))
    rb_raise(rb_eArgError, "a (1th argument) must be NArray");
  if (NA_RANK(rb_a) != 2)
    rb_raise(rb_eArgError, "rank of a (1th argument) must be %d", 2);
  n = NA_SHAPE1(rb_a);
  lda = NA_SHAPE0(rb_a);
  if (NA_TYPE(rb_a) != NA_DFLOAT)
    rb_a = na_change_type(rb_a, NA_DFLOAT);
  a = NA_PTR_TYPE(rb_a, doublereal*);
  m = lda;
  {
    int shape[1];
    shape[0] = m;
    rb_tau = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
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

  dtzrqf_(&m, &n, a, &lda, tau, &info);

  rb_info = INT2NUM(info);
  return rb_ary_new3(3, rb_tau, rb_info, rb_a);
}

void
init_lapack_dtzrqf(VALUE mLapack){
  rb_define_module_function(mLapack, "dtzrqf", rb_dtzrqf, -1);
}

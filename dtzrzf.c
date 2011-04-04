#include "rb_lapack.h"

extern VOID dtzrzf_(integer *m, integer *n, doublereal *a, integer *lda, doublereal *tau, doublereal *work, integer *lwork, integer *info);

static VALUE
rb_dtzrzf(int argc, VALUE *argv, VALUE self){
  VALUE rb_a;
  doublereal *a; 
  VALUE rb_lwork;
  integer lwork; 
  VALUE rb_tau;
  doublereal *tau; 
  VALUE rb_work;
  doublereal *work; 
  VALUE rb_info;
  integer info; 
  VALUE rb_a_out__;
  doublereal *a_out__;

  integer lda;
  integer n;
  integer m;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  tau, work, info, a = NumRu::Lapack.dtzrzf( a, lwork)\n    or\n  NumRu::Lapack.dtzrzf  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE DTZRZF( M, N, A, LDA, TAU, WORK, LWORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  DTZRZF reduces the M-by-N ( M<=N ) real upper trapezoidal matrix A\n*  to upper triangular form by means of orthogonal transformations.\n*\n*  The upper trapezoidal matrix A is factored as\n*\n*     A = ( R  0 ) * Z,\n*\n*  where Z is an N-by-N orthogonal matrix and R is an M-by-M upper\n*  triangular matrix.\n*\n\n*  Arguments\n*  =========\n*\n*  M       (input) INTEGER\n*          The number of rows of the matrix A.  M >= 0.\n*\n*  N       (input) INTEGER\n*          The number of columns of the matrix A.  N >= M.\n*\n*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)\n*          On entry, the leading M-by-N upper trapezoidal part of the\n*          array A must contain the matrix to be factorized.\n*          On exit, the leading M-by-M upper triangular part of A\n*          contains the upper triangular matrix R, and elements M+1 to\n*          N of the first M rows of A, with the array TAU, represent the\n*          orthogonal matrix Z as a product of M elementary reflectors.\n*\n*  LDA     (input) INTEGER\n*          The leading dimension of the array A.  LDA >= max(1,M).\n*\n*  TAU     (output) DOUBLE PRECISION array, dimension (M)\n*          The scalar factors of the elementary reflectors.\n*\n*  WORK    (workspace/output) DOUBLE PRECISION array, dimension (MAX(1,LWORK))\n*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.\n*\n*  LWORK   (input) INTEGER\n*          The dimension of the array WORK.  LWORK >= max(1,M).\n*          For optimum performance LWORK >= M*NB, where NB is\n*          the optimal blocksize.\n*\n*          If LWORK = -1, then a workspace query is assumed; the routine\n*          only calculates the optimal size of the WORK array, returns\n*          this value as the first entry of the WORK array, and no error\n*          message related to LWORK is issued by XERBLA.\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit\n*          < 0:  if INFO = -i, the i-th argument had an illegal value\n*\n\n*  Further Details\n*  ===============\n*\n*  Based on contributions by\n*    A. Petitet, Computer Science Dept., Univ. of Tenn., Knoxville, USA\n*\n*  The factorization is obtained by Householder's method.  The kth\n*  transformation matrix, Z( k ), which is used to introduce zeros into\n*  the ( m - k + 1 )th row of A, is given in the form\n*\n*     Z( k ) = ( I     0   ),\n*              ( 0  T( k ) )\n*\n*  where\n*\n*     T( k ) = I - tau*u( k )*u( k )',   u( k ) = (   1    ),\n*                                                 (   0    )\n*                                                 ( z( k ) )\n*\n*  tau is a scalar and z( k ) is an ( n - m ) element vector.\n*  tau and z( k ) are chosen to annihilate the elements of the kth row\n*  of X.\n*\n*  The scalar tau is returned in the kth element of TAU and the vector\n*  u( k ) in the kth row of A, such that the elements of z( k ) are\n*  in  a( k, m + 1 ), ..., a( k, n ). The elements of R are returned in\n*  the upper triangular part of A.\n*\n*  Z is given by\n*\n*     Z =  Z( 1 ) * Z( 2 ) * ... * Z( m ).\n*\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 2)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 2)", argc);
  rb_a = argv[0];
  rb_lwork = argv[1];

  if (!NA_IsNArray(rb_a))
    rb_raise(rb_eArgError, "a (1th argument) must be NArray");
  if (NA_RANK(rb_a) != 2)
    rb_raise(rb_eArgError, "rank of a (1th argument) must be %d", 2);
  n = NA_SHAPE1(rb_a);
  lda = NA_SHAPE0(rb_a);
  if (NA_TYPE(rb_a) != NA_DFLOAT)
    rb_a = na_change_type(rb_a, NA_DFLOAT);
  a = NA_PTR_TYPE(rb_a, doublereal*);
  lwork = NUM2INT(rb_lwork);
  m = lda;
  {
    int shape[1];
    shape[0] = m;
    rb_tau = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  tau = NA_PTR_TYPE(rb_tau, doublereal*);
  {
    int shape[1];
    shape[0] = MAX(1,lwork);
    rb_work = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  work = NA_PTR_TYPE(rb_work, doublereal*);
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

  dtzrzf_(&m, &n, a, &lda, tau, work, &lwork, &info);

  rb_info = INT2NUM(info);
  return rb_ary_new3(4, rb_tau, rb_work, rb_info, rb_a);
}

void
init_lapack_dtzrzf(VALUE mLapack){
  rb_define_module_function(mLapack, "dtzrzf", rb_dtzrzf, -1);
}

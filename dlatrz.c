#include "rb_lapack.h"

static VALUE
rb_dlatrz(int argc, VALUE *argv, VALUE self){
  VALUE rb_l;
  integer l; 
  VALUE rb_a;
  doublereal *a; 
  VALUE rb_tau;
  doublereal *tau; 
  VALUE rb_a_out__;
  doublereal *a_out__;
  doublereal *work;

  integer lda;
  integer n;
  integer m;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  tau, a = NumRu::Lapack.dlatrz( l, a)\n    or\n  NumRu::Lapack.dlatrz  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE DLATRZ( M, N, L, A, LDA, TAU, WORK )\n\n*  Purpose\n*  =======\n*\n*  DLATRZ factors the M-by-(M+L) real upper trapezoidal matrix\n*  [ A1 A2 ] = [ A(1:M,1:M) A(1:M,N-L+1:N) ] as ( R  0 ) * Z, by means\n*  of orthogonal transformations.  Z is an (M+L)-by-(M+L) orthogonal\n*  matrix and, R and A1 are M-by-M upper triangular matrices.\n*\n\n*  Arguments\n*  =========\n*\n*  M       (input) INTEGER\n*          The number of rows of the matrix A.  M >= 0.\n*\n*  N       (input) INTEGER\n*          The number of columns of the matrix A.  N >= 0.\n*\n*  L       (input) INTEGER\n*          The number of columns of the matrix A containing the\n*          meaningful part of the Householder vectors. N-M >= L >= 0.\n*\n*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)\n*          On entry, the leading M-by-N upper trapezoidal part of the\n*          array A must contain the matrix to be factorized.\n*          On exit, the leading M-by-M upper triangular part of A\n*          contains the upper triangular matrix R, and elements N-L+1 to\n*          N of the first M rows of A, with the array TAU, represent the\n*          orthogonal matrix Z as a product of M elementary reflectors.\n*\n*  LDA     (input) INTEGER\n*          The leading dimension of the array A.  LDA >= max(1,M).\n*\n*  TAU     (output) DOUBLE PRECISION array, dimension (M)\n*          The scalar factors of the elementary reflectors.\n*\n*  WORK    (workspace) DOUBLE PRECISION array, dimension (M)\n*\n\n*  Further Details\n*  ===============\n*\n*  Based on contributions by\n*    A. Petitet, Computer Science Dept., Univ. of Tenn., Knoxville, USA\n*\n*  The factorization is obtained by Householder's method.  The kth\n*  transformation matrix, Z( k ), which is used to introduce zeros into\n*  the ( m - k + 1 )th row of A, is given in the form\n*\n*     Z( k ) = ( I     0   ),\n*              ( 0  T( k ) )\n*\n*  where\n*\n*     T( k ) = I - tau*u( k )*u( k )',   u( k ) = (   1    ),\n*                                                 (   0    )\n*                                                 ( z( k ) )\n*\n*  tau is a scalar and z( k ) is an l element vector. tau and z( k )\n*  are chosen to annihilate the elements of the kth row of A2.\n*\n*  The scalar tau is returned in the kth element of TAU and the vector\n*  u( k ) in the kth row of A2, such that the elements of z( k ) are\n*  in  a( k, l + 1 ), ..., a( k, n ). The elements of R are returned in\n*  the upper triangular part of A1.\n*\n*  Z is given by\n*\n*     Z =  Z( 1 ) * Z( 2 ) * ... * Z( m ).\n*\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 2)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 2)", argc);
  rb_l = argv[0];
  rb_a = argv[1];

  l = NUM2INT(rb_l);
  if (!NA_IsNArray(rb_a))
    rb_raise(rb_eArgError, "a (2th argument) must be NArray");
  if (NA_RANK(rb_a) != 2)
    rb_raise(rb_eArgError, "rank of a (2th argument) must be %d", 2);
  lda = NA_SHAPE0(rb_a);
  n = NA_SHAPE1(rb_a);
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
  work = ALLOC_N(doublereal, (m));

  dlatrz_(&m, &n, &l, a, &lda, tau, work);

  free(work);
  return rb_ary_new3(2, rb_tau, rb_a);
}

void
init_lapack_dlatrz(VALUE mLapack){
  rb_define_module_function(mLapack, "dlatrz", rb_dlatrz, -1);
}

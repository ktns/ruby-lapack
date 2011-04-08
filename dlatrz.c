#include "rb_lapack.h"

extern VOID dlatrz_(integer *m, integer *n, integer *l, doublereal *a, integer *lda, doublereal *tau, doublereal *work);

static VALUE sHelp, sUsage;

static VALUE
rblapack_dlatrz(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_l;
  integer l; 
  VALUE rblapack_a;
  doublereal *a; 
  VALUE rblapack_tau;
  doublereal *tau; 
  VALUE rblapack_a_out__;
  doublereal *a_out__;
  doublereal *work;

  integer lda;
  integer n;
  integer m;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  tau, a = NumRu::Lapack.dlatrz( l, a, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE DLATRZ( M, N, L, A, LDA, TAU, WORK )\n\n*  Purpose\n*  =======\n*\n*  DLATRZ factors the M-by-(M+L) real upper trapezoidal matrix\n*  [ A1 A2 ] = [ A(1:M,1:M) A(1:M,N-L+1:N) ] as ( R  0 ) * Z, by means\n*  of orthogonal transformations.  Z is an (M+L)-by-(M+L) orthogonal\n*  matrix and, R and A1 are M-by-M upper triangular matrices.\n*\n\n*  Arguments\n*  =========\n*\n*  M       (input) INTEGER\n*          The number of rows of the matrix A.  M >= 0.\n*\n*  N       (input) INTEGER\n*          The number of columns of the matrix A.  N >= 0.\n*\n*  L       (input) INTEGER\n*          The number of columns of the matrix A containing the\n*          meaningful part of the Householder vectors. N-M >= L >= 0.\n*\n*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)\n*          On entry, the leading M-by-N upper trapezoidal part of the\n*          array A must contain the matrix to be factorized.\n*          On exit, the leading M-by-M upper triangular part of A\n*          contains the upper triangular matrix R, and elements N-L+1 to\n*          N of the first M rows of A, with the array TAU, represent the\n*          orthogonal matrix Z as a product of M elementary reflectors.\n*\n*  LDA     (input) INTEGER\n*          The leading dimension of the array A.  LDA >= max(1,M).\n*\n*  TAU     (output) DOUBLE PRECISION array, dimension (M)\n*          The scalar factors of the elementary reflectors.\n*\n*  WORK    (workspace) DOUBLE PRECISION array, dimension (M)\n*\n\n*  Further Details\n*  ===============\n*\n*  Based on contributions by\n*    A. Petitet, Computer Science Dept., Univ. of Tenn., Knoxville, USA\n*\n*  The factorization is obtained by Householder's method.  The kth\n*  transformation matrix, Z( k ), which is used to introduce zeros into\n*  the ( m - k + 1 )th row of A, is given in the form\n*\n*     Z( k ) = ( I     0   ),\n*              ( 0  T( k ) )\n*\n*  where\n*\n*     T( k ) = I - tau*u( k )*u( k )',   u( k ) = (   1    ),\n*                                                 (   0    )\n*                                                 ( z( k ) )\n*\n*  tau is a scalar and z( k ) is an l element vector. tau and z( k )\n*  are chosen to annihilate the elements of the kth row of A2.\n*\n*  The scalar tau is returned in the kth element of TAU and the vector\n*  u( k ) in the kth row of A2, such that the elements of z( k ) are\n*  in  a( k, l + 1 ), ..., a( k, n ). The elements of R are returned in\n*  the upper triangular part of A1.\n*\n*  Z is given by\n*\n*     Z =  Z( 1 ) * Z( 2 ) * ... * Z( m ).\n*\n*  =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  tau, a = NumRu::Lapack.dlatrz( l, a, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 2)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 2)", argc);
  rblapack_l = argv[0];
  rblapack_a = argv[1];
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
  l = NUM2INT(rblapack_l);
  m = lda;
  {
    int shape[1];
    shape[0] = m;
    rblapack_tau = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  tau = NA_PTR_TYPE(rblapack_tau, doublereal*);
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

  dlatrz_(&m, &n, &l, a, &lda, tau, work);

  free(work);
  return rb_ary_new3(2, rblapack_tau, rblapack_a);
}

void
init_lapack_dlatrz(VALUE mLapack){
  rb_define_module_function(mLapack, "dlatrz", rblapack_dlatrz, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}

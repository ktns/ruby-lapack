#include "rb_lapack.h"

extern VOID dlaswp_(integer *n, doublereal *a, integer *lda, integer *k1, integer *k2, integer *ipiv, integer *incx);

static VALUE
rb_dlaswp(int argc, VALUE *argv, VALUE self){
  VALUE rb_a;
  doublereal *a; 
  VALUE rb_k1;
  integer k1; 
  VALUE rb_k2;
  integer k2; 
  VALUE rb_ipiv;
  integer *ipiv; 
  VALUE rb_incx;
  integer incx; 
  VALUE rb_a_out__;
  doublereal *a_out__;

  integer lda;
  integer n;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  a = NumRu::Lapack.dlaswp( a, k1, k2, ipiv, incx)\n    or\n  NumRu::Lapack.dlaswp  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE DLASWP( N, A, LDA, K1, K2, IPIV, INCX )\n\n*  Purpose\n*  =======\n*\n*  DLASWP performs a series of row interchanges on the matrix A.\n*  One row interchange is initiated for each of rows K1 through K2 of A.\n*\n\n*  Arguments\n*  =========\n*\n*  N       (input) INTEGER\n*          The number of columns of the matrix A.\n*\n*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)\n*          On entry, the matrix of column dimension N to which the row\n*          interchanges will be applied.\n*          On exit, the permuted matrix.\n*\n*  LDA     (input) INTEGER\n*          The leading dimension of the array A.\n*\n*  K1      (input) INTEGER\n*          The first element of IPIV for which a row interchange will\n*          be done.\n*\n*  K2      (input) INTEGER\n*          The last element of IPIV for which a row interchange will\n*          be done.\n*\n*  IPIV    (input) INTEGER array, dimension (K2*abs(INCX))\n*          The vector of pivot indices.  Only the elements in positions\n*          K1 through K2 of IPIV are accessed.\n*          IPIV(K) = L implies rows K and L are to be interchanged.\n*\n*  INCX    (input) INTEGER\n*          The increment between successive values of IPIV.  If IPIV\n*          is negative, the pivots are applied in reverse order.\n*\n\n*  Further Details\n*  ===============\n*\n*  Modified by\n*   R. C. Whaley, Computer Science Dept., Univ. of Tenn., Knoxville, USA\n*\n* =====================================================================\n*\n*     .. Local Scalars ..\n      INTEGER            I, I1, I2, INC, IP, IX, IX0, J, K, N32\n      DOUBLE PRECISION   TEMP\n*     ..\n\n");
    return Qnil;
  }
  if (argc != 5)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 5)", argc);
  rb_a = argv[0];
  rb_k1 = argv[1];
  rb_k2 = argv[2];
  rb_ipiv = argv[3];
  rb_incx = argv[4];

  k2 = NUM2INT(rb_k2);
  if (!NA_IsNArray(rb_a))
    rb_raise(rb_eArgError, "a (1th argument) must be NArray");
  if (NA_RANK(rb_a) != 2)
    rb_raise(rb_eArgError, "rank of a (1th argument) must be %d", 2);
  n = NA_SHAPE1(rb_a);
  lda = NA_SHAPE0(rb_a);
  if (NA_TYPE(rb_a) != NA_DFLOAT)
    rb_a = na_change_type(rb_a, NA_DFLOAT);
  a = NA_PTR_TYPE(rb_a, doublereal*);
  k1 = NUM2INT(rb_k1);
  incx = NUM2INT(rb_incx);
  if (!NA_IsNArray(rb_ipiv))
    rb_raise(rb_eArgError, "ipiv (4th argument) must be NArray");
  if (NA_RANK(rb_ipiv) != 1)
    rb_raise(rb_eArgError, "rank of ipiv (4th argument) must be %d", 1);
  if (NA_SHAPE0(rb_ipiv) != (k2*abs(incx)))
    rb_raise(rb_eRuntimeError, "shape 0 of ipiv must be %d", k2*abs(incx));
  if (NA_TYPE(rb_ipiv) != NA_LINT)
    rb_ipiv = na_change_type(rb_ipiv, NA_LINT);
  ipiv = NA_PTR_TYPE(rb_ipiv, integer*);
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

  dlaswp_(&n, a, &lda, &k1, &k2, ipiv, &incx);

  return rb_a;
}

void
init_lapack_dlaswp(VALUE mLapack){
  rb_define_module_function(mLapack, "dlaswp", rb_dlaswp, -1);
}

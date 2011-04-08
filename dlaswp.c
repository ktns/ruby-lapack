#include "rb_lapack.h"

extern VOID dlaswp_(integer *n, doublereal *a, integer *lda, integer *k1, integer *k2, integer *ipiv, integer *incx);

static VALUE sHelp, sUsage;

static VALUE
rblapack_dlaswp(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_a;
  doublereal *a; 
  VALUE rblapack_k1;
  integer k1; 
  VALUE rblapack_k2;
  integer k2; 
  VALUE rblapack_ipiv;
  integer *ipiv; 
  VALUE rblapack_incx;
  integer incx; 
  VALUE rblapack_a_out__;
  doublereal *a_out__;

  integer lda;
  integer n;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  a = NumRu::Lapack.dlaswp( a, k1, k2, ipiv, incx, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE DLASWP( N, A, LDA, K1, K2, IPIV, INCX )\n\n*  Purpose\n*  =======\n*\n*  DLASWP performs a series of row interchanges on the matrix A.\n*  One row interchange is initiated for each of rows K1 through K2 of A.\n*\n\n*  Arguments\n*  =========\n*\n*  N       (input) INTEGER\n*          The number of columns of the matrix A.\n*\n*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)\n*          On entry, the matrix of column dimension N to which the row\n*          interchanges will be applied.\n*          On exit, the permuted matrix.\n*\n*  LDA     (input) INTEGER\n*          The leading dimension of the array A.\n*\n*  K1      (input) INTEGER\n*          The first element of IPIV for which a row interchange will\n*          be done.\n*\n*  K2      (input) INTEGER\n*          The last element of IPIV for which a row interchange will\n*          be done.\n*\n*  IPIV    (input) INTEGER array, dimension (K2*abs(INCX))\n*          The vector of pivot indices.  Only the elements in positions\n*          K1 through K2 of IPIV are accessed.\n*          IPIV(K) = L implies rows K and L are to be interchanged.\n*\n*  INCX    (input) INTEGER\n*          The increment between successive values of IPIV.  If IPIV\n*          is negative, the pivots are applied in reverse order.\n*\n\n*  Further Details\n*  ===============\n*\n*  Modified by\n*   R. C. Whaley, Computer Science Dept., Univ. of Tenn., Knoxville, USA\n*\n* =====================================================================\n*\n*     .. Local Scalars ..\n      INTEGER            I, I1, I2, INC, IP, IX, IX0, J, K, N32\n      DOUBLE PRECISION   TEMP\n*     ..\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  a = NumRu::Lapack.dlaswp( a, k1, k2, ipiv, incx, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 5)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 5)", argc);
  rblapack_a = argv[0];
  rblapack_k1 = argv[1];
  rblapack_k2 = argv[2];
  rblapack_ipiv = argv[3];
  rblapack_incx = argv[4];
  if (rb_options != Qnil) {
  }

  k2 = NUM2INT(rblapack_k2);
  if (!NA_IsNArray(rblapack_a))
    rb_raise(rb_eArgError, "a (1th argument) must be NArray");
  if (NA_RANK(rblapack_a) != 2)
    rb_raise(rb_eArgError, "rank of a (1th argument) must be %d", 2);
  n = NA_SHAPE1(rblapack_a);
  lda = NA_SHAPE0(rblapack_a);
  if (NA_TYPE(rblapack_a) != NA_DFLOAT)
    rblapack_a = na_change_type(rblapack_a, NA_DFLOAT);
  a = NA_PTR_TYPE(rblapack_a, doublereal*);
  k1 = NUM2INT(rblapack_k1);
  incx = NUM2INT(rblapack_incx);
  if (!NA_IsNArray(rblapack_ipiv))
    rb_raise(rb_eArgError, "ipiv (4th argument) must be NArray");
  if (NA_RANK(rblapack_ipiv) != 1)
    rb_raise(rb_eArgError, "rank of ipiv (4th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_ipiv) != (k2*abs(incx)))
    rb_raise(rb_eRuntimeError, "shape 0 of ipiv must be %d", k2*abs(incx));
  if (NA_TYPE(rblapack_ipiv) != NA_LINT)
    rblapack_ipiv = na_change_type(rblapack_ipiv, NA_LINT);
  ipiv = NA_PTR_TYPE(rblapack_ipiv, integer*);
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

  dlaswp_(&n, a, &lda, &k1, &k2, ipiv, &incx);

  return rblapack_a;
}

void
init_lapack_dlaswp(VALUE mLapack){
  rb_define_module_function(mLapack, "dlaswp", rblapack_dlaswp, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}

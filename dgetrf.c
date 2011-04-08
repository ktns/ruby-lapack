#include "rb_lapack.h"

extern VOID dgetrf_(integer *m, integer *n, doublereal *a, integer *lda, integer *ipiv, integer *info);

static VALUE sHelp, sUsage;

static VALUE
rblapack_dgetrf(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_m;
  integer m; 
  VALUE rblapack_a;
  doublereal *a; 
  VALUE rblapack_ipiv;
  integer *ipiv; 
  VALUE rblapack_info;
  integer info; 
  VALUE rblapack_a_out__;
  doublereal *a_out__;

  integer lda;
  integer n;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  ipiv, info, a = NumRu::Lapack.dgetrf( m, a, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE DGETRF( M, N, A, LDA, IPIV, INFO )\n\n*  Purpose\n*  =======\n*\n*  DGETRF computes an LU factorization of a general M-by-N matrix A\n*  using partial pivoting with row interchanges.\n*\n*  The factorization has the form\n*     A = P * L * U\n*  where P is a permutation matrix, L is lower triangular with unit\n*  diagonal elements (lower trapezoidal if m > n), and U is upper\n*  triangular (upper trapezoidal if m < n).\n*\n*  This is the right-looking Level 3 BLAS version of the algorithm.\n*\n\n*  Arguments\n*  =========\n*\n*  M       (input) INTEGER\n*          The number of rows of the matrix A.  M >= 0.\n*\n*  N       (input) INTEGER\n*          The number of columns of the matrix A.  N >= 0.\n*\n*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)\n*          On entry, the M-by-N matrix to be factored.\n*          On exit, the factors L and U from the factorization\n*          A = P*L*U; the unit diagonal elements of L are not stored.\n*\n*  LDA     (input) INTEGER\n*          The leading dimension of the array A.  LDA >= max(1,M).\n*\n*  IPIV    (output) INTEGER array, dimension (min(M,N))\n*          The pivot indices; for 1 <= i <= min(M,N), row i of the\n*          matrix was interchanged with row IPIV(i).\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit\n*          < 0:  if INFO = -i, the i-th argument had an illegal value\n*          > 0:  if INFO = i, U(i,i) is exactly zero. The factorization\n*                has been completed, but the factor U is exactly\n*                singular, and division by zero will occur if it is used\n*                to solve a system of equations.\n*\n\n*  =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  ipiv, info, a = NumRu::Lapack.dgetrf( m, a, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 2)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 2)", argc);
  rblapack_m = argv[0];
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
  m = NUM2INT(rblapack_m);
  {
    int shape[1];
    shape[0] = MIN(m,n);
    rblapack_ipiv = na_make_object(NA_LINT, 1, shape, cNArray);
  }
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

  dgetrf_(&m, &n, a, &lda, ipiv, &info);

  rblapack_info = INT2NUM(info);
  return rb_ary_new3(3, rblapack_ipiv, rblapack_info, rblapack_a);
}

void
init_lapack_dgetrf(VALUE mLapack){
  rb_define_module_function(mLapack, "dgetrf", rblapack_dgetrf, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}

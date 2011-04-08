#include "rb_lapack.h"

extern VOID slagtf_(integer *n, real *a, real *lambda, real *b, real *c, real *tol, real *d, integer *in, integer *info);

static VALUE sHelp, sUsage;

static VALUE
rblapack_slagtf(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_a;
  real *a; 
  VALUE rblapack_lambda;
  real lambda; 
  VALUE rblapack_b;
  real *b; 
  VALUE rblapack_c;
  real *c; 
  VALUE rblapack_tol;
  real tol; 
  VALUE rblapack_d;
  real *d; 
  VALUE rblapack_in;
  integer *in; 
  VALUE rblapack_info;
  integer info; 
  VALUE rblapack_a_out__;
  real *a_out__;
  VALUE rblapack_b_out__;
  real *b_out__;
  VALUE rblapack_c_out__;
  real *c_out__;

  integer n;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  d, in, info, a, b, c = NumRu::Lapack.slagtf( a, lambda, b, c, tol, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE SLAGTF( N, A, LAMBDA, B, C, TOL, D, IN, INFO )\n\n*  Purpose\n*  =======\n*\n*  SLAGTF factorizes the matrix (T - lambda*I), where T is an n by n\n*  tridiagonal matrix and lambda is a scalar, as\n*\n*     T - lambda*I = PLU,\n*\n*  where P is a permutation matrix, L is a unit lower tridiagonal matrix\n*  with at most one non-zero sub-diagonal elements per column and U is\n*  an upper triangular matrix with at most two non-zero super-diagonal\n*  elements per column.\n*\n*  The factorization is obtained by Gaussian elimination with partial\n*  pivoting and implicit row scaling.\n*\n*  The parameter LAMBDA is included in the routine so that SLAGTF may\n*  be used, in conjunction with SLAGTS, to obtain eigenvectors of T by\n*  inverse iteration.\n*\n\n*  Arguments\n*  =========\n*\n*  N       (input) INTEGER\n*          The order of the matrix T.\n*\n*  A       (input/output) REAL array, dimension (N)\n*          On entry, A must contain the diagonal elements of T.\n*\n*          On exit, A is overwritten by the n diagonal elements of the\n*          upper triangular matrix U of the factorization of T.\n*\n*  LAMBDA  (input) REAL\n*          On entry, the scalar lambda.\n*\n*  B       (input/output) REAL array, dimension (N-1)\n*          On entry, B must contain the (n-1) super-diagonal elements of\n*          T.\n*\n*          On exit, B is overwritten by the (n-1) super-diagonal\n*          elements of the matrix U of the factorization of T.\n*\n*  C       (input/output) REAL array, dimension (N-1)\n*          On entry, C must contain the (n-1) sub-diagonal elements of\n*          T.\n*\n*          On exit, C is overwritten by the (n-1) sub-diagonal elements\n*          of the matrix L of the factorization of T.\n*\n*  TOL     (input) REAL\n*          On entry, a relative tolerance used to indicate whether or\n*          not the matrix (T - lambda*I) is nearly singular. TOL should\n*          normally be chose as approximately the largest relative error\n*          in the elements of T. For example, if the elements of T are\n*          correct to about 4 significant figures, then TOL should be\n*          set to about 5*10**(-4). If TOL is supplied as less than eps,\n*          where eps is the relative machine precision, then the value\n*          eps is used in place of TOL.\n*\n*  D       (output) REAL array, dimension (N-2)\n*          On exit, D is overwritten by the (n-2) second super-diagonal\n*          elements of the matrix U of the factorization of T.\n*\n*  IN      (output) INTEGER array, dimension (N)\n*          On exit, IN contains details of the permutation matrix P. If\n*          an interchange occurred at the kth step of the elimination,\n*          then IN(k) = 1, otherwise IN(k) = 0. The element IN(n)\n*          returns the smallest positive integer j such that\n*\n*             abs( u(j,j) ).le. norm( (T - lambda*I)(j) )*TOL,\n*\n*          where norm( A(j) ) denotes the sum of the absolute values of\n*          the jth row of the matrix A. If no such j exists then IN(n)\n*          is returned as zero. If IN(n) is returned as positive, then a\n*          diagonal element of U is small, indicating that\n*          (T - lambda*I) is singular or nearly singular,\n*\n*  INFO    (output) INTEGER\n*          = 0   : successful exit\n*          .lt. 0: if INFO = -k, the kth argument had an illegal value\n*\n\n* =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  d, in, info, a, b, c = NumRu::Lapack.slagtf( a, lambda, b, c, tol, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 5)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 5)", argc);
  rblapack_a = argv[0];
  rblapack_lambda = argv[1];
  rblapack_b = argv[2];
  rblapack_c = argv[3];
  rblapack_tol = argv[4];
  if (rb_options != Qnil) {
  }

  if (!NA_IsNArray(rblapack_a))
    rb_raise(rb_eArgError, "a (1th argument) must be NArray");
  if (NA_RANK(rblapack_a) != 1)
    rb_raise(rb_eArgError, "rank of a (1th argument) must be %d", 1);
  n = NA_SHAPE0(rblapack_a);
  if (NA_TYPE(rblapack_a) != NA_SFLOAT)
    rblapack_a = na_change_type(rblapack_a, NA_SFLOAT);
  a = NA_PTR_TYPE(rblapack_a, real*);
  tol = (real)NUM2DBL(rblapack_tol);
  lambda = (real)NUM2DBL(rblapack_lambda);
  if (!NA_IsNArray(rblapack_b))
    rb_raise(rb_eArgError, "b (3th argument) must be NArray");
  if (NA_RANK(rblapack_b) != 1)
    rb_raise(rb_eArgError, "rank of b (3th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_b) != (n-1))
    rb_raise(rb_eRuntimeError, "shape 0 of b must be %d", n-1);
  if (NA_TYPE(rblapack_b) != NA_SFLOAT)
    rblapack_b = na_change_type(rblapack_b, NA_SFLOAT);
  b = NA_PTR_TYPE(rblapack_b, real*);
  if (!NA_IsNArray(rblapack_c))
    rb_raise(rb_eArgError, "c (4th argument) must be NArray");
  if (NA_RANK(rblapack_c) != 1)
    rb_raise(rb_eArgError, "rank of c (4th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_c) != (n-1))
    rb_raise(rb_eRuntimeError, "shape 0 of c must be %d", n-1);
  if (NA_TYPE(rblapack_c) != NA_SFLOAT)
    rblapack_c = na_change_type(rblapack_c, NA_SFLOAT);
  c = NA_PTR_TYPE(rblapack_c, real*);
  {
    int shape[1];
    shape[0] = n-2;
    rblapack_d = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  d = NA_PTR_TYPE(rblapack_d, real*);
  {
    int shape[1];
    shape[0] = n;
    rblapack_in = na_make_object(NA_LINT, 1, shape, cNArray);
  }
  in = NA_PTR_TYPE(rblapack_in, integer*);
  {
    int shape[1];
    shape[0] = n;
    rblapack_a_out__ = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  a_out__ = NA_PTR_TYPE(rblapack_a_out__, real*);
  MEMCPY(a_out__, a, real, NA_TOTAL(rblapack_a));
  rblapack_a = rblapack_a_out__;
  a = a_out__;
  {
    int shape[1];
    shape[0] = n-1;
    rblapack_b_out__ = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  b_out__ = NA_PTR_TYPE(rblapack_b_out__, real*);
  MEMCPY(b_out__, b, real, NA_TOTAL(rblapack_b));
  rblapack_b = rblapack_b_out__;
  b = b_out__;
  {
    int shape[1];
    shape[0] = n-1;
    rblapack_c_out__ = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  c_out__ = NA_PTR_TYPE(rblapack_c_out__, real*);
  MEMCPY(c_out__, c, real, NA_TOTAL(rblapack_c));
  rblapack_c = rblapack_c_out__;
  c = c_out__;

  slagtf_(&n, a, &lambda, b, c, &tol, d, in, &info);

  rblapack_info = INT2NUM(info);
  return rb_ary_new3(6, rblapack_d, rblapack_in, rblapack_info, rblapack_a, rblapack_b, rblapack_c);
}

void
init_lapack_slagtf(VALUE mLapack){
  rb_define_module_function(mLapack, "slagtf", rblapack_slagtf, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}

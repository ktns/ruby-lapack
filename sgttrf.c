#include "rb_lapack.h"

extern VOID sgttrf_(integer *n, real *dl, real *d, real *du, real *du2, integer *ipiv, integer *info);

static VALUE sHelp, sUsage;

static VALUE
rblapack_sgttrf(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_dl;
  real *dl; 
  VALUE rblapack_d;
  real *d; 
  VALUE rblapack_du;
  real *du; 
  VALUE rblapack_du2;
  real *du2; 
  VALUE rblapack_ipiv;
  integer *ipiv; 
  VALUE rblapack_info;
  integer info; 
  VALUE rblapack_dl_out__;
  real *dl_out__;
  VALUE rblapack_d_out__;
  real *d_out__;
  VALUE rblapack_du_out__;
  real *du_out__;

  integer n;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  du2, ipiv, info, dl, d, du = NumRu::Lapack.sgttrf( dl, d, du, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE SGTTRF( N, DL, D, DU, DU2, IPIV, INFO )\n\n*  Purpose\n*  =======\n*\n*  SGTTRF computes an LU factorization of a real tridiagonal matrix A\n*  using elimination with partial pivoting and row interchanges.\n*\n*  The factorization has the form\n*     A = L * U\n*  where L is a product of permutation and unit lower bidiagonal\n*  matrices and U is upper triangular with nonzeros in only the main\n*  diagonal and first two superdiagonals.\n*\n\n*  Arguments\n*  =========\n*\n*  N       (input) INTEGER\n*          The order of the matrix A.\n*\n*  DL      (input/output) REAL array, dimension (N-1)\n*          On entry, DL must contain the (n-1) sub-diagonal elements of\n*          A.\n*\n*          On exit, DL is overwritten by the (n-1) multipliers that\n*          define the matrix L from the LU factorization of A.\n*\n*  D       (input/output) REAL array, dimension (N)\n*          On entry, D must contain the diagonal elements of A.\n*\n*          On exit, D is overwritten by the n diagonal elements of the\n*          upper triangular matrix U from the LU factorization of A.\n*\n*  DU      (input/output) REAL array, dimension (N-1)\n*          On entry, DU must contain the (n-1) super-diagonal elements\n*          of A.\n*\n*          On exit, DU is overwritten by the (n-1) elements of the first\n*          super-diagonal of U.\n*\n*  DU2     (output) REAL array, dimension (N-2)\n*          On exit, DU2 is overwritten by the (n-2) elements of the\n*          second super-diagonal of U.\n*\n*  IPIV    (output) INTEGER array, dimension (N)\n*          The pivot indices; for 1 <= i <= n, row i of the matrix was\n*          interchanged with row IPIV(i).  IPIV(i) will always be either\n*          i or i+1; IPIV(i) = i indicates a row interchange was not\n*          required.\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit\n*          < 0:  if INFO = -k, the k-th argument had an illegal value\n*          > 0:  if INFO = k, U(k,k) is exactly zero. The factorization\n*                has been completed, but the factor U is exactly\n*                singular, and division by zero will occur if it is used\n*                to solve a system of equations.\n*\n\n*  =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  du2, ipiv, info, dl, d, du = NumRu::Lapack.sgttrf( dl, d, du, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 3)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 3)", argc);
  rblapack_dl = argv[0];
  rblapack_d = argv[1];
  rblapack_du = argv[2];
  if (rb_options != Qnil) {
  }

  if (!NA_IsNArray(rblapack_d))
    rb_raise(rb_eArgError, "d (2th argument) must be NArray");
  if (NA_RANK(rblapack_d) != 1)
    rb_raise(rb_eArgError, "rank of d (2th argument) must be %d", 1);
  n = NA_SHAPE0(rblapack_d);
  if (NA_TYPE(rblapack_d) != NA_SFLOAT)
    rblapack_d = na_change_type(rblapack_d, NA_SFLOAT);
  d = NA_PTR_TYPE(rblapack_d, real*);
  if (!NA_IsNArray(rblapack_du))
    rb_raise(rb_eArgError, "du (3th argument) must be NArray");
  if (NA_RANK(rblapack_du) != 1)
    rb_raise(rb_eArgError, "rank of du (3th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_du) != (n-1))
    rb_raise(rb_eRuntimeError, "shape 0 of du must be %d", n-1);
  if (NA_TYPE(rblapack_du) != NA_SFLOAT)
    rblapack_du = na_change_type(rblapack_du, NA_SFLOAT);
  du = NA_PTR_TYPE(rblapack_du, real*);
  if (!NA_IsNArray(rblapack_dl))
    rb_raise(rb_eArgError, "dl (1th argument) must be NArray");
  if (NA_RANK(rblapack_dl) != 1)
    rb_raise(rb_eArgError, "rank of dl (1th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_dl) != (n-1))
    rb_raise(rb_eRuntimeError, "shape 0 of dl must be %d", n-1);
  if (NA_TYPE(rblapack_dl) != NA_SFLOAT)
    rblapack_dl = na_change_type(rblapack_dl, NA_SFLOAT);
  dl = NA_PTR_TYPE(rblapack_dl, real*);
  {
    int shape[1];
    shape[0] = n-2;
    rblapack_du2 = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  du2 = NA_PTR_TYPE(rblapack_du2, real*);
  {
    int shape[1];
    shape[0] = n;
    rblapack_ipiv = na_make_object(NA_LINT, 1, shape, cNArray);
  }
  ipiv = NA_PTR_TYPE(rblapack_ipiv, integer*);
  {
    int shape[1];
    shape[0] = n-1;
    rblapack_dl_out__ = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  dl_out__ = NA_PTR_TYPE(rblapack_dl_out__, real*);
  MEMCPY(dl_out__, dl, real, NA_TOTAL(rblapack_dl));
  rblapack_dl = rblapack_dl_out__;
  dl = dl_out__;
  {
    int shape[1];
    shape[0] = n;
    rblapack_d_out__ = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  d_out__ = NA_PTR_TYPE(rblapack_d_out__, real*);
  MEMCPY(d_out__, d, real, NA_TOTAL(rblapack_d));
  rblapack_d = rblapack_d_out__;
  d = d_out__;
  {
    int shape[1];
    shape[0] = n-1;
    rblapack_du_out__ = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  du_out__ = NA_PTR_TYPE(rblapack_du_out__, real*);
  MEMCPY(du_out__, du, real, NA_TOTAL(rblapack_du));
  rblapack_du = rblapack_du_out__;
  du = du_out__;

  sgttrf_(&n, dl, d, du, du2, ipiv, &info);

  rblapack_info = INT2NUM(info);
  return rb_ary_new3(6, rblapack_du2, rblapack_ipiv, rblapack_info, rblapack_dl, rblapack_d, rblapack_du);
}

void
init_lapack_sgttrf(VALUE mLapack){
  rb_define_module_function(mLapack, "sgttrf", rblapack_sgttrf, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}

#include "rb_lapack.h"

static VALUE
rb_cgttrf(int argc, VALUE *argv, VALUE self){
  VALUE rb_dl;
  complex *dl; 
  VALUE rb_d;
  complex *d; 
  VALUE rb_du;
  complex *du; 
  VALUE rb_du2;
  complex *du2; 
  VALUE rb_ipiv;
  integer *ipiv; 
  VALUE rb_info;
  integer info; 
  VALUE rb_dl_out__;
  complex *dl_out__;
  VALUE rb_d_out__;
  complex *d_out__;
  VALUE rb_du_out__;
  complex *du_out__;

  integer n;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  du2, ipiv, info, dl, d, du = NumRu::Lapack.cgttrf( dl, d, du)\n    or\n  NumRu::Lapack.cgttrf  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE CGTTRF( N, DL, D, DU, DU2, IPIV, INFO )\n\n*  Purpose\n*  =======\n*\n*  CGTTRF computes an LU factorization of a complex tridiagonal matrix A\n*  using elimination with partial pivoting and row interchanges.\n*\n*  The factorization has the form\n*     A = L * U\n*  where L is a product of permutation and unit lower bidiagonal\n*  matrices and U is upper triangular with nonzeros in only the main\n*  diagonal and first two superdiagonals.\n*\n\n*  Arguments\n*  =========\n*\n*  N       (input) INTEGER\n*          The order of the matrix A.\n*\n*  DL      (input/output) COMPLEX array, dimension (N-1)\n*          On entry, DL must contain the (n-1) sub-diagonal elements of\n*          A.\n*\n*          On exit, DL is overwritten by the (n-1) multipliers that\n*          define the matrix L from the LU factorization of A.\n*\n*  D       (input/output) COMPLEX array, dimension (N)\n*          On entry, D must contain the diagonal elements of A.\n*\n*          On exit, D is overwritten by the n diagonal elements of the\n*          upper triangular matrix U from the LU factorization of A.\n*\n*  DU      (input/output) COMPLEX array, dimension (N-1)\n*          On entry, DU must contain the (n-1) super-diagonal elements\n*          of A.\n*\n*          On exit, DU is overwritten by the (n-1) elements of the first\n*          super-diagonal of U.\n*\n*  DU2     (output) COMPLEX array, dimension (N-2)\n*          On exit, DU2 is overwritten by the (n-2) elements of the\n*          second super-diagonal of U.\n*\n*  IPIV    (output) INTEGER array, dimension (N)\n*          The pivot indices; for 1 <= i <= n, row i of the matrix was\n*          interchanged with row IPIV(i).  IPIV(i) will always be either\n*          i or i+1; IPIV(i) = i indicates a row interchange was not\n*          required.\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit\n*          < 0:  if INFO = -k, the k-th argument had an illegal value\n*          > 0:  if INFO = k, U(k,k) is exactly zero. The factorization\n*                has been completed, but the factor U is exactly\n*                singular, and division by zero will occur if it is used\n*                to solve a system of equations.\n*\n\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 3)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 3)", argc);
  rb_dl = argv[0];
  rb_d = argv[1];
  rb_du = argv[2];

  if (!NA_IsNArray(rb_d))
    rb_raise(rb_eArgError, "d (1th argument) must be NArray");
  if (NA_RANK(rb_d) != 1)
    rb_raise(rb_eArgError, "rank of d (1th argument) must be %d", 1);
  n = NA_SHAPE0(rb_d);
  if (NA_TYPE(rb_d) != NA_SCOMPLEX)
    rb_d = na_change_type(rb_d, NA_SCOMPLEX);
  d = NA_PTR_TYPE(rb_d, complex*);
  if (!NA_IsNArray(rb_dl))
    rb_raise(rb_eArgError, "dl (2th argument) must be NArray");
  if (NA_RANK(rb_dl) != 1)
    rb_raise(rb_eArgError, "rank of dl (2th argument) must be %d", 1);
  if (NA_SHAPE0(rb_dl) != (n-1))
    rb_raise(rb_eRuntimeError, "shape 0 of dl must be %d", n-1);
  if (NA_TYPE(rb_dl) != NA_SCOMPLEX)
    rb_dl = na_change_type(rb_dl, NA_SCOMPLEX);
  dl = NA_PTR_TYPE(rb_dl, complex*);
  if (!NA_IsNArray(rb_du))
    rb_raise(rb_eArgError, "du (3th argument) must be NArray");
  if (NA_RANK(rb_du) != 1)
    rb_raise(rb_eArgError, "rank of du (3th argument) must be %d", 1);
  if (NA_SHAPE0(rb_du) != (n-1))
    rb_raise(rb_eRuntimeError, "shape 0 of du must be %d", n-1);
  if (NA_TYPE(rb_du) != NA_SCOMPLEX)
    rb_du = na_change_type(rb_du, NA_SCOMPLEX);
  du = NA_PTR_TYPE(rb_du, complex*);
  {
    int shape[1];
    shape[0] = n-2;
    rb_du2 = na_make_object(NA_SCOMPLEX, 1, shape, cNArray);
  }
  du2 = NA_PTR_TYPE(rb_du2, complex*);
  {
    int shape[1];
    shape[0] = DIM_LEN(n);
    rb_ipiv = na_make_object(NA_LINT, 1, shape, cNArray);
  }
  ipiv = NA_PTR_TYPE(rb_ipiv, integer*);
  {
    int shape[1];
    shape[0] = n-1;
    rb_dl_out__ = na_make_object(NA_SCOMPLEX, 1, shape, cNArray);
  }
  dl_out__ = NA_PTR_TYPE(rb_dl_out__, complex*);
  MEMCPY(dl_out__, dl, complex, NA_TOTAL(rb_dl));
  rb_dl = rb_dl_out__;
  dl = dl_out__;
  {
    int shape[1];
    shape[0] = n;
    rb_d_out__ = na_make_object(NA_SCOMPLEX, 1, shape, cNArray);
  }
  d_out__ = NA_PTR_TYPE(rb_d_out__, complex*);
  MEMCPY(d_out__, d, complex, NA_TOTAL(rb_d));
  rb_d = rb_d_out__;
  d = d_out__;
  {
    int shape[1];
    shape[0] = n-1;
    rb_du_out__ = na_make_object(NA_SCOMPLEX, 1, shape, cNArray);
  }
  du_out__ = NA_PTR_TYPE(rb_du_out__, complex*);
  MEMCPY(du_out__, du, complex, NA_TOTAL(rb_du));
  rb_du = rb_du_out__;
  du = du_out__;

  cgttrf_(&n, dl, d, du, du2, ipiv, &info);

  rb_info = INT2NUM(info);
  return rb_ary_new3(6, rb_du2, rb_ipiv, rb_info, rb_dl, rb_d, rb_du);
}

void
init_lapack_cgttrf(VALUE mLapack){
  rb_define_module_function(mLapack, "cgttrf", rb_cgttrf, -1);
}

#include "rb_lapack.h"

extern VOID sgtsv_(integer *n, integer *nrhs, real *dl, real *d, real *du, real *b, integer *ldb, integer *info);

static VALUE
rb_sgtsv(int argc, VALUE *argv, VALUE self){
  VALUE rb_dl;
  real *dl; 
  VALUE rb_d;
  real *d; 
  VALUE rb_du;
  real *du; 
  VALUE rb_b;
  real *b; 
  VALUE rb_info;
  integer info; 
  VALUE rb_dl_out__;
  real *dl_out__;
  VALUE rb_d_out__;
  real *d_out__;
  VALUE rb_du_out__;
  real *du_out__;
  VALUE rb_b_out__;
  real *b_out__;

  integer n;
  integer ldb;
  integer nrhs;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  info, dl, d, du, b = NumRu::Lapack.sgtsv( dl, d, du, b)\n    or\n  NumRu::Lapack.sgtsv  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE SGTSV( N, NRHS, DL, D, DU, B, LDB, INFO )\n\n*  Purpose\n*  =======\n*\n*  SGTSV  solves the equation\n*\n*     A*X = B,\n*\n*  where A is an n by n tridiagonal matrix, by Gaussian elimination with\n*  partial pivoting.\n*\n*  Note that the equation  A'*X = B  may be solved by interchanging the\n*  order of the arguments DU and DL.\n*\n\n*  Arguments\n*  =========\n*\n*  N       (input) INTEGER\n*          The order of the matrix A.  N >= 0.\n*\n*  NRHS    (input) INTEGER\n*          The number of right hand sides, i.e., the number of columns\n*          of the matrix B.  NRHS >= 0.\n*\n*  DL      (input/output) REAL array, dimension (N-1)\n*          On entry, DL must contain the (n-1) sub-diagonal elements of\n*          A.\n*\n*          On exit, DL is overwritten by the (n-2) elements of the\n*          second super-diagonal of the upper triangular matrix U from\n*          the LU factorization of A, in DL(1), ..., DL(n-2).\n*\n*  D       (input/output) REAL array, dimension (N)\n*          On entry, D must contain the diagonal elements of A.\n*\n*          On exit, D is overwritten by the n diagonal elements of U.\n*\n*  DU      (input/output) REAL array, dimension (N-1)\n*          On entry, DU must contain the (n-1) super-diagonal elements\n*          of A.\n*\n*          On exit, DU is overwritten by the (n-1) elements of the first\n*          super-diagonal of U.\n*\n*  B       (input/output) REAL array, dimension (LDB,NRHS)\n*          On entry, the N by NRHS matrix of right hand side matrix B.\n*          On exit, if INFO = 0, the N by NRHS solution matrix X.\n*\n*  LDB     (input) INTEGER\n*          The leading dimension of the array B.  LDB >= max(1,N).\n*\n*  INFO    (output) INTEGER\n*          = 0: successful exit\n*          < 0: if INFO = -i, the i-th argument had an illegal value\n*          > 0: if INFO = i, U(i,i) is exactly zero, and the solution\n*               has not been computed.  The factorization has not been\n*               completed unless i = N.\n*\n\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 4)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 4)", argc);
  rb_dl = argv[0];
  rb_d = argv[1];
  rb_du = argv[2];
  rb_b = argv[3];

  if (!NA_IsNArray(rb_b))
    rb_raise(rb_eArgError, "b (4th argument) must be NArray");
  if (NA_RANK(rb_b) != 2)
    rb_raise(rb_eArgError, "rank of b (4th argument) must be %d", 2);
  nrhs = NA_SHAPE1(rb_b);
  ldb = NA_SHAPE0(rb_b);
  if (NA_TYPE(rb_b) != NA_SFLOAT)
    rb_b = na_change_type(rb_b, NA_SFLOAT);
  b = NA_PTR_TYPE(rb_b, real*);
  if (!NA_IsNArray(rb_d))
    rb_raise(rb_eArgError, "d (2th argument) must be NArray");
  if (NA_RANK(rb_d) != 1)
    rb_raise(rb_eArgError, "rank of d (2th argument) must be %d", 1);
  n = NA_SHAPE0(rb_d);
  if (NA_TYPE(rb_d) != NA_SFLOAT)
    rb_d = na_change_type(rb_d, NA_SFLOAT);
  d = NA_PTR_TYPE(rb_d, real*);
  if (!NA_IsNArray(rb_du))
    rb_raise(rb_eArgError, "du (3th argument) must be NArray");
  if (NA_RANK(rb_du) != 1)
    rb_raise(rb_eArgError, "rank of du (3th argument) must be %d", 1);
  if (NA_SHAPE0(rb_du) != (n-1))
    rb_raise(rb_eRuntimeError, "shape 0 of du must be %d", n-1);
  if (NA_TYPE(rb_du) != NA_SFLOAT)
    rb_du = na_change_type(rb_du, NA_SFLOAT);
  du = NA_PTR_TYPE(rb_du, real*);
  if (!NA_IsNArray(rb_dl))
    rb_raise(rb_eArgError, "dl (1th argument) must be NArray");
  if (NA_RANK(rb_dl) != 1)
    rb_raise(rb_eArgError, "rank of dl (1th argument) must be %d", 1);
  if (NA_SHAPE0(rb_dl) != (n-1))
    rb_raise(rb_eRuntimeError, "shape 0 of dl must be %d", n-1);
  if (NA_TYPE(rb_dl) != NA_SFLOAT)
    rb_dl = na_change_type(rb_dl, NA_SFLOAT);
  dl = NA_PTR_TYPE(rb_dl, real*);
  {
    int shape[1];
    shape[0] = n-1;
    rb_dl_out__ = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  dl_out__ = NA_PTR_TYPE(rb_dl_out__, real*);
  MEMCPY(dl_out__, dl, real, NA_TOTAL(rb_dl));
  rb_dl = rb_dl_out__;
  dl = dl_out__;
  {
    int shape[1];
    shape[0] = n;
    rb_d_out__ = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  d_out__ = NA_PTR_TYPE(rb_d_out__, real*);
  MEMCPY(d_out__, d, real, NA_TOTAL(rb_d));
  rb_d = rb_d_out__;
  d = d_out__;
  {
    int shape[1];
    shape[0] = n-1;
    rb_du_out__ = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  du_out__ = NA_PTR_TYPE(rb_du_out__, real*);
  MEMCPY(du_out__, du, real, NA_TOTAL(rb_du));
  rb_du = rb_du_out__;
  du = du_out__;
  {
    int shape[2];
    shape[0] = ldb;
    shape[1] = nrhs;
    rb_b_out__ = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  b_out__ = NA_PTR_TYPE(rb_b_out__, real*);
  MEMCPY(b_out__, b, real, NA_TOTAL(rb_b));
  rb_b = rb_b_out__;
  b = b_out__;

  sgtsv_(&n, &nrhs, dl, d, du, b, &ldb, &info);

  rb_info = INT2NUM(info);
  return rb_ary_new3(5, rb_info, rb_dl, rb_d, rb_du, rb_b);
}

void
init_lapack_sgtsv(VALUE mLapack){
  rb_define_module_function(mLapack, "sgtsv", rb_sgtsv, -1);
}

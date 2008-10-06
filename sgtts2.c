#include "rb_lapack.h"

static VALUE
rb_sgtts2(int argc, VALUE *argv, VALUE self){
  VALUE rb_itrans;
  integer itrans; 
  VALUE rb_dl;
  real *dl; 
  VALUE rb_d;
  real *d; 
  VALUE rb_du;
  real *du; 
  VALUE rb_du2;
  real *du2; 
  VALUE rb_ipiv;
  integer *ipiv; 
  VALUE rb_b;
  real *b; 
  VALUE rb_b_out__;
  real *b_out__;

  integer n;
  integer ldb;
  integer nrhs;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  b = NumRu::Lapack.sgtts2( itrans, dl, d, du, du2, ipiv, b)\n    or\n  NumRu::Lapack.sgtts2  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE SGTTS2( ITRANS, N, NRHS, DL, D, DU, DU2, IPIV, B, LDB )\n\n*  Purpose\n*  =======\n*\n*  SGTTS2 solves one of the systems of equations\n*     A*X = B  or  A'*X = B,\n*  with a tridiagonal matrix A using the LU factorization computed\n*  by SGTTRF.\n*\n\n*  Arguments\n*  =========\n*\n*  ITRANS  (input) INTEGER\n*          Specifies the form of the system of equations.\n*          = 0:  A * X = B  (No transpose)\n*          = 1:  A'* X = B  (Transpose)\n*          = 2:  A'* X = B  (Conjugate transpose = Transpose)\n*\n*  N       (input) INTEGER\n*          The order of the matrix A.\n*\n*  NRHS    (input) INTEGER\n*          The number of right hand sides, i.e., the number of columns\n*          of the matrix B.  NRHS >= 0.\n*\n*  DL      (input) REAL array, dimension (N-1)\n*          The (n-1) multipliers that define the matrix L from the\n*          LU factorization of A.\n*\n*  D       (input) REAL array, dimension (N)\n*          The n diagonal elements of the upper triangular matrix U from\n*          the LU factorization of A.\n*\n*  DU      (input) REAL array, dimension (N-1)\n*          The (n-1) elements of the first super-diagonal of U.\n*\n*  DU2     (input) REAL array, dimension (N-2)\n*          The (n-2) elements of the second super-diagonal of U.\n*\n*  IPIV    (input) INTEGER array, dimension (N)\n*          The pivot indices; for 1 <= i <= n, row i of the matrix was\n*          interchanged with row IPIV(i).  IPIV(i) will always be either\n*          i or i+1; IPIV(i) = i indicates a row interchange was not\n*          required.\n*\n*  B       (input/output) REAL array, dimension (LDB,NRHS)\n*          On entry, the matrix of right hand side vectors B.\n*          On exit, B is overwritten by the solution vectors X.\n*\n*  LDB     (input) INTEGER\n*          The leading dimension of the array B.  LDB >= max(1,N).\n*\n\n*  =====================================================================\n*\n*     .. Local Scalars ..\n      INTEGER            I, IP, J\n      REAL               TEMP\n*     ..\n\n");
    return Qnil;
  }
  if (argc != 7)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 7)", argc);
  rb_itrans = argv[0];
  rb_dl = argv[1];
  rb_d = argv[2];
  rb_du = argv[3];
  rb_du2 = argv[4];
  rb_ipiv = argv[5];
  rb_b = argv[6];

  itrans = NUM2INT(rb_itrans);
  if (!NA_IsNArray(rb_ipiv))
    rb_raise(rb_eArgError, "ipiv (2th argument) must be NArray");
  if (NA_RANK(rb_ipiv) != 1)
    rb_raise(rb_eArgError, "rank of ipiv (2th argument) must be %d", 1);
  n = NA_SHAPE0(rb_ipiv);
  if (NA_TYPE(rb_ipiv) != NA_LINT)
    rb_ipiv = na_change_type(rb_ipiv, NA_LINT);
  ipiv = NA_PTR_TYPE(rb_ipiv, integer*);
  if (!NA_IsNArray(rb_d))
    rb_raise(rb_eArgError, "d (3th argument) must be NArray");
  if (NA_RANK(rb_d) != 1)
    rb_raise(rb_eArgError, "rank of d (3th argument) must be %d", 1);
  if (NA_SHAPE0(rb_d) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of d must be the same as shape 0 of ipiv");
  if (NA_TYPE(rb_d) != NA_SFLOAT)
    rb_d = na_change_type(rb_d, NA_SFLOAT);
  d = NA_PTR_TYPE(rb_d, real*);
  if (!NA_IsNArray(rb_dl))
    rb_raise(rb_eArgError, "dl (4th argument) must be NArray");
  if (NA_RANK(rb_dl) != 1)
    rb_raise(rb_eArgError, "rank of dl (4th argument) must be %d", 1);
  if (NA_SHAPE0(rb_dl) != (n-1))
    rb_raise(rb_eRuntimeError, "shape 0 of dl must be %d", n-1);
  if (NA_TYPE(rb_dl) != NA_SFLOAT)
    rb_dl = na_change_type(rb_dl, NA_SFLOAT);
  dl = NA_PTR_TYPE(rb_dl, real*);
  if (!NA_IsNArray(rb_du2))
    rb_raise(rb_eArgError, "du2 (5th argument) must be NArray");
  if (NA_RANK(rb_du2) != 1)
    rb_raise(rb_eArgError, "rank of du2 (5th argument) must be %d", 1);
  if (NA_SHAPE0(rb_du2) != (n-2))
    rb_raise(rb_eRuntimeError, "shape 0 of du2 must be %d", n-2);
  if (NA_TYPE(rb_du2) != NA_SFLOAT)
    rb_du2 = na_change_type(rb_du2, NA_SFLOAT);
  du2 = NA_PTR_TYPE(rb_du2, real*);
  if (!NA_IsNArray(rb_du))
    rb_raise(rb_eArgError, "du (6th argument) must be NArray");
  if (NA_RANK(rb_du) != 1)
    rb_raise(rb_eArgError, "rank of du (6th argument) must be %d", 1);
  if (NA_SHAPE0(rb_du) != (n-1))
    rb_raise(rb_eRuntimeError, "shape 0 of du must be %d", n-1);
  if (NA_TYPE(rb_du) != NA_SFLOAT)
    rb_du = na_change_type(rb_du, NA_SFLOAT);
  du = NA_PTR_TYPE(rb_du, real*);
  if (!NA_IsNArray(rb_b))
    rb_raise(rb_eArgError, "b (7th argument) must be NArray");
  if (NA_RANK(rb_b) != 2)
    rb_raise(rb_eArgError, "rank of b (7th argument) must be %d", 2);
  ldb = NA_SHAPE0(rb_b);
  nrhs = NA_SHAPE1(rb_b);
  if (NA_TYPE(rb_b) != NA_SFLOAT)
    rb_b = na_change_type(rb_b, NA_SFLOAT);
  b = NA_PTR_TYPE(rb_b, real*);
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

  sgtts2_(&itrans, &n, &nrhs, dl, d, du, du2, ipiv, b, &ldb);

  return rb_b;
}

void
init_lapack_sgtts2(VALUE mLapack){
  rb_define_module_function(mLapack, "sgtts2", rb_sgtts2, -1);
}

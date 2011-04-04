#include "rb_lapack.h"

extern VOID dgttrs_(char *trans, integer *n, integer *nrhs, doublereal *dl, doublereal *d, doublereal *du, doublereal *du2, integer *ipiv, doublereal *b, integer *ldb, integer *info);

static VALUE
rb_dgttrs(int argc, VALUE *argv, VALUE self){
  VALUE rb_trans;
  char trans; 
  VALUE rb_dl;
  doublereal *dl; 
  VALUE rb_d;
  doublereal *d; 
  VALUE rb_du;
  doublereal *du; 
  VALUE rb_du2;
  doublereal *du2; 
  VALUE rb_ipiv;
  integer *ipiv; 
  VALUE rb_b;
  doublereal *b; 
  VALUE rb_info;
  integer info; 
  VALUE rb_b_out__;
  doublereal *b_out__;

  integer n;
  integer ldb;
  integer nrhs;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  info, b = NumRu::Lapack.dgttrs( trans, dl, d, du, du2, ipiv, b)\n    or\n  NumRu::Lapack.dgttrs  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE DGTTRS( TRANS, N, NRHS, DL, D, DU, DU2, IPIV, B, LDB, INFO )\n\n*  Purpose\n*  =======\n*\n*  DGTTRS solves one of the systems of equations\n*     A*X = B  or  A'*X = B,\n*  with a tridiagonal matrix A using the LU factorization computed\n*  by DGTTRF.\n*\n\n*  Arguments\n*  =========\n*\n*  TRANS   (input) CHARACTER*1\n*          Specifies the form of the system of equations.\n*          = 'N':  A * X = B  (No transpose)\n*          = 'T':  A'* X = B  (Transpose)\n*          = 'C':  A'* X = B  (Conjugate transpose = Transpose)\n*\n*  N       (input) INTEGER\n*          The order of the matrix A.\n*\n*  NRHS    (input) INTEGER\n*          The number of right hand sides, i.e., the number of columns\n*          of the matrix B.  NRHS >= 0.\n*\n*  DL      (input) DOUBLE PRECISION array, dimension (N-1)\n*          The (n-1) multipliers that define the matrix L from the\n*          LU factorization of A.\n*\n*  D       (input) DOUBLE PRECISION array, dimension (N)\n*          The n diagonal elements of the upper triangular matrix U from\n*          the LU factorization of A.\n*\n*  DU      (input) DOUBLE PRECISION array, dimension (N-1)\n*          The (n-1) elements of the first super-diagonal of U.\n*\n*  DU2     (input) DOUBLE PRECISION array, dimension (N-2)\n*          The (n-2) elements of the second super-diagonal of U.\n*\n*  IPIV    (input) INTEGER array, dimension (N)\n*          The pivot indices; for 1 <= i <= n, row i of the matrix was\n*          interchanged with row IPIV(i).  IPIV(i) will always be either\n*          i or i+1; IPIV(i) = i indicates a row interchange was not\n*          required.\n*\n*  B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)\n*          On entry, the matrix of right hand side vectors B.\n*          On exit, B is overwritten by the solution vectors X.\n*\n*  LDB     (input) INTEGER\n*          The leading dimension of the array B.  LDB >= max(1,N).\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit\n*          < 0:  if INFO = -i, the i-th argument had an illegal value\n*\n\n*  =====================================================================\n*\n*     .. Local Scalars ..\n      LOGICAL            NOTRAN\n      INTEGER            ITRANS, J, JB, NB\n*     ..\n*     .. External Functions ..\n      INTEGER            ILAENV\n      EXTERNAL           ILAENV\n*     ..\n*     .. External Subroutines ..\n      EXTERNAL           DGTTS2, XERBLA\n*     ..\n*     .. Intrinsic Functions ..\n      INTRINSIC          MAX, MIN\n*     ..\n\n");
    return Qnil;
  }
  if (argc != 7)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 7)", argc);
  rb_trans = argv[0];
  rb_dl = argv[1];
  rb_d = argv[2];
  rb_du = argv[3];
  rb_du2 = argv[4];
  rb_ipiv = argv[5];
  rb_b = argv[6];

  if (!NA_IsNArray(rb_ipiv))
    rb_raise(rb_eArgError, "ipiv (6th argument) must be NArray");
  if (NA_RANK(rb_ipiv) != 1)
    rb_raise(rb_eArgError, "rank of ipiv (6th argument) must be %d", 1);
  n = NA_SHAPE0(rb_ipiv);
  if (NA_TYPE(rb_ipiv) != NA_LINT)
    rb_ipiv = na_change_type(rb_ipiv, NA_LINT);
  ipiv = NA_PTR_TYPE(rb_ipiv, integer*);
  trans = StringValueCStr(rb_trans)[0];
  if (!NA_IsNArray(rb_b))
    rb_raise(rb_eArgError, "b (7th argument) must be NArray");
  if (NA_RANK(rb_b) != 2)
    rb_raise(rb_eArgError, "rank of b (7th argument) must be %d", 2);
  nrhs = NA_SHAPE1(rb_b);
  ldb = NA_SHAPE0(rb_b);
  if (NA_TYPE(rb_b) != NA_DFLOAT)
    rb_b = na_change_type(rb_b, NA_DFLOAT);
  b = NA_PTR_TYPE(rb_b, doublereal*);
  if (!NA_IsNArray(rb_d))
    rb_raise(rb_eArgError, "d (3th argument) must be NArray");
  if (NA_RANK(rb_d) != 1)
    rb_raise(rb_eArgError, "rank of d (3th argument) must be %d", 1);
  if (NA_SHAPE0(rb_d) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of d must be the same as shape 0 of ipiv");
  if (NA_TYPE(rb_d) != NA_DFLOAT)
    rb_d = na_change_type(rb_d, NA_DFLOAT);
  d = NA_PTR_TYPE(rb_d, doublereal*);
  if (!NA_IsNArray(rb_du))
    rb_raise(rb_eArgError, "du (4th argument) must be NArray");
  if (NA_RANK(rb_du) != 1)
    rb_raise(rb_eArgError, "rank of du (4th argument) must be %d", 1);
  if (NA_SHAPE0(rb_du) != (n-1))
    rb_raise(rb_eRuntimeError, "shape 0 of du must be %d", n-1);
  if (NA_TYPE(rb_du) != NA_DFLOAT)
    rb_du = na_change_type(rb_du, NA_DFLOAT);
  du = NA_PTR_TYPE(rb_du, doublereal*);
  if (!NA_IsNArray(rb_dl))
    rb_raise(rb_eArgError, "dl (2th argument) must be NArray");
  if (NA_RANK(rb_dl) != 1)
    rb_raise(rb_eArgError, "rank of dl (2th argument) must be %d", 1);
  if (NA_SHAPE0(rb_dl) != (n-1))
    rb_raise(rb_eRuntimeError, "shape 0 of dl must be %d", n-1);
  if (NA_TYPE(rb_dl) != NA_DFLOAT)
    rb_dl = na_change_type(rb_dl, NA_DFLOAT);
  dl = NA_PTR_TYPE(rb_dl, doublereal*);
  if (!NA_IsNArray(rb_du2))
    rb_raise(rb_eArgError, "du2 (5th argument) must be NArray");
  if (NA_RANK(rb_du2) != 1)
    rb_raise(rb_eArgError, "rank of du2 (5th argument) must be %d", 1);
  if (NA_SHAPE0(rb_du2) != (n-2))
    rb_raise(rb_eRuntimeError, "shape 0 of du2 must be %d", n-2);
  if (NA_TYPE(rb_du2) != NA_DFLOAT)
    rb_du2 = na_change_type(rb_du2, NA_DFLOAT);
  du2 = NA_PTR_TYPE(rb_du2, doublereal*);
  {
    int shape[2];
    shape[0] = ldb;
    shape[1] = nrhs;
    rb_b_out__ = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  b_out__ = NA_PTR_TYPE(rb_b_out__, doublereal*);
  MEMCPY(b_out__, b, doublereal, NA_TOTAL(rb_b));
  rb_b = rb_b_out__;
  b = b_out__;

  dgttrs_(&trans, &n, &nrhs, dl, d, du, du2, ipiv, b, &ldb, &info);

  rb_info = INT2NUM(info);
  return rb_ary_new3(2, rb_info, rb_b);
}

void
init_lapack_dgttrs(VALUE mLapack){
  rb_define_module_function(mLapack, "dgttrs", rb_dgttrs, -1);
}

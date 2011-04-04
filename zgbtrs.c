#include "rb_lapack.h"

extern VOID zgbtrs_(char *trans, integer *n, integer *kl, integer *ku, integer *nrhs, doublecomplex *ab, integer *ldab, integer *ipiv, doublecomplex *b, integer *ldb, integer *info);

static VALUE
rb_zgbtrs(int argc, VALUE *argv, VALUE self){
  VALUE rb_trans;
  char trans; 
  VALUE rb_kl;
  integer kl; 
  VALUE rb_ku;
  integer ku; 
  VALUE rb_ab;
  doublecomplex *ab; 
  VALUE rb_ipiv;
  integer *ipiv; 
  VALUE rb_b;
  doublecomplex *b; 
  VALUE rb_info;
  integer info; 
  VALUE rb_b_out__;
  doublecomplex *b_out__;

  integer ldab;
  integer n;
  integer ldb;
  integer nrhs;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  info, b = NumRu::Lapack.zgbtrs( trans, kl, ku, ab, ipiv, b)\n    or\n  NumRu::Lapack.zgbtrs  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE ZGBTRS( TRANS, N, KL, KU, NRHS, AB, LDAB, IPIV, B, LDB, INFO )\n\n*  Purpose\n*  =======\n*\n*  ZGBTRS solves a system of linear equations\n*     A * X = B,  A**T * X = B,  or  A**H * X = B\n*  with a general band matrix A using the LU factorization computed\n*  by ZGBTRF.\n*\n\n*  Arguments\n*  =========\n*\n*  TRANS   (input) CHARACTER*1\n*          Specifies the form of the system of equations.\n*          = 'N':  A * X = B     (No transpose)\n*          = 'T':  A**T * X = B  (Transpose)\n*          = 'C':  A**H * X = B  (Conjugate transpose)\n*\n*  N       (input) INTEGER\n*          The order of the matrix A.  N >= 0.\n*\n*  KL      (input) INTEGER\n*          The number of subdiagonals within the band of A.  KL >= 0.\n*\n*  KU      (input) INTEGER\n*          The number of superdiagonals within the band of A.  KU >= 0.\n*\n*  NRHS    (input) INTEGER\n*          The number of right hand sides, i.e., the number of columns\n*          of the matrix B.  NRHS >= 0.\n*\n*  AB      (input) COMPLEX*16 array, dimension (LDAB,N)\n*          Details of the LU factorization of the band matrix A, as\n*          computed by ZGBTRF.  U is stored as an upper triangular band\n*          matrix with KL+KU superdiagonals in rows 1 to KL+KU+1, and\n*          the multipliers used during the factorization are stored in\n*          rows KL+KU+2 to 2*KL+KU+1.\n*\n*  LDAB    (input) INTEGER\n*          The leading dimension of the array AB.  LDAB >= 2*KL+KU+1.\n*\n*  IPIV    (input) INTEGER array, dimension (N)\n*          The pivot indices; for 1 <= i <= N, row i of the matrix was\n*          interchanged with row IPIV(i).\n*\n*  B       (input/output) COMPLEX*16 array, dimension (LDB,NRHS)\n*          On entry, the right hand side matrix B.\n*          On exit, the solution matrix X.\n*\n*  LDB     (input) INTEGER\n*          The leading dimension of the array B.  LDB >= max(1,N).\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit\n*          < 0:  if INFO = -i, the i-th argument had an illegal value\n*\n\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 6)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 6)", argc);
  rb_trans = argv[0];
  rb_kl = argv[1];
  rb_ku = argv[2];
  rb_ab = argv[3];
  rb_ipiv = argv[4];
  rb_b = argv[5];

  if (!NA_IsNArray(rb_ipiv))
    rb_raise(rb_eArgError, "ipiv (5th argument) must be NArray");
  if (NA_RANK(rb_ipiv) != 1)
    rb_raise(rb_eArgError, "rank of ipiv (5th argument) must be %d", 1);
  n = NA_SHAPE0(rb_ipiv);
  if (NA_TYPE(rb_ipiv) != NA_LINT)
    rb_ipiv = na_change_type(rb_ipiv, NA_LINT);
  ipiv = NA_PTR_TYPE(rb_ipiv, integer*);
  if (!NA_IsNArray(rb_ab))
    rb_raise(rb_eArgError, "ab (4th argument) must be NArray");
  if (NA_RANK(rb_ab) != 2)
    rb_raise(rb_eArgError, "rank of ab (4th argument) must be %d", 2);
  if (NA_SHAPE1(rb_ab) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of ab must be the same as shape 0 of ipiv");
  ldab = NA_SHAPE0(rb_ab);
  if (NA_TYPE(rb_ab) != NA_DCOMPLEX)
    rb_ab = na_change_type(rb_ab, NA_DCOMPLEX);
  ab = NA_PTR_TYPE(rb_ab, doublecomplex*);
  kl = NUM2INT(rb_kl);
  if (!NA_IsNArray(rb_b))
    rb_raise(rb_eArgError, "b (6th argument) must be NArray");
  if (NA_RANK(rb_b) != 2)
    rb_raise(rb_eArgError, "rank of b (6th argument) must be %d", 2);
  nrhs = NA_SHAPE1(rb_b);
  ldb = NA_SHAPE0(rb_b);
  if (NA_TYPE(rb_b) != NA_DCOMPLEX)
    rb_b = na_change_type(rb_b, NA_DCOMPLEX);
  b = NA_PTR_TYPE(rb_b, doublecomplex*);
  ku = NUM2INT(rb_ku);
  trans = StringValueCStr(rb_trans)[0];
  {
    int shape[2];
    shape[0] = ldb;
    shape[1] = nrhs;
    rb_b_out__ = na_make_object(NA_DCOMPLEX, 2, shape, cNArray);
  }
  b_out__ = NA_PTR_TYPE(rb_b_out__, doublecomplex*);
  MEMCPY(b_out__, b, doublecomplex, NA_TOTAL(rb_b));
  rb_b = rb_b_out__;
  b = b_out__;

  zgbtrs_(&trans, &n, &kl, &ku, &nrhs, ab, &ldab, ipiv, b, &ldb, &info);

  rb_info = INT2NUM(info);
  return rb_ary_new3(2, rb_info, rb_b);
}

void
init_lapack_zgbtrs(VALUE mLapack){
  rb_define_module_function(mLapack, "zgbtrs", rb_zgbtrs, -1);
}

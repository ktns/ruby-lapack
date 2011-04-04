#include "rb_lapack.h"

extern VOID cgbtrs_(char *trans, integer *n, integer *kl, integer *ku, integer *nrhs, complex *ab, integer *ldab, integer *ipiv, complex *b, integer *ldb, integer *info);

static VALUE
rb_cgbtrs(int argc, VALUE *argv, VALUE self){
  VALUE rb_trans;
  char trans; 
  VALUE rb_kl;
  integer kl; 
  VALUE rb_ku;
  integer ku; 
  VALUE rb_ab;
  complex *ab; 
  VALUE rb_ipiv;
  integer *ipiv; 
  VALUE rb_b;
  complex *b; 
  VALUE rb_info;
  integer info; 
  VALUE rb_b_out__;
  complex *b_out__;

  integer ldab;
  integer n;
  integer ldb;
  integer nrhs;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  info, b = NumRu::Lapack.cgbtrs( trans, kl, ku, ab, ipiv, b)\n    or\n  NumRu::Lapack.cgbtrs  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE CGBTRS( TRANS, N, KL, KU, NRHS, AB, LDAB, IPIV, B, LDB, INFO )\n\n*  Purpose\n*  =======\n*\n*  CGBTRS solves a system of linear equations\n*     A * X = B,  A**T * X = B,  or  A**H * X = B\n*  with a general band matrix A using the LU factorization computed\n*  by CGBTRF.\n*\n\n*  Arguments\n*  =========\n*\n*  TRANS   (input) CHARACTER*1\n*          Specifies the form of the system of equations.\n*          = 'N':  A * X = B     (No transpose)\n*          = 'T':  A**T * X = B  (Transpose)\n*          = 'C':  A**H * X = B  (Conjugate transpose)\n*\n*  N       (input) INTEGER\n*          The order of the matrix A.  N >= 0.\n*\n*  KL      (input) INTEGER\n*          The number of subdiagonals within the band of A.  KL >= 0.\n*\n*  KU      (input) INTEGER\n*          The number of superdiagonals within the band of A.  KU >= 0.\n*\n*  NRHS    (input) INTEGER\n*          The number of right hand sides, i.e., the number of columns\n*          of the matrix B.  NRHS >= 0.\n*\n*  AB      (input) COMPLEX array, dimension (LDAB,N)\n*          Details of the LU factorization of the band matrix A, as\n*          computed by CGBTRF.  U is stored as an upper triangular band\n*          matrix with KL+KU superdiagonals in rows 1 to KL+KU+1, and\n*          the multipliers used during the factorization are stored in\n*          rows KL+KU+2 to 2*KL+KU+1.\n*\n*  LDAB    (input) INTEGER\n*          The leading dimension of the array AB.  LDAB >= 2*KL+KU+1.\n*\n*  IPIV    (input) INTEGER array, dimension (N)\n*          The pivot indices; for 1 <= i <= N, row i of the matrix was\n*          interchanged with row IPIV(i).\n*\n*  B       (input/output) COMPLEX array, dimension (LDB,NRHS)\n*          On entry, the right hand side matrix B.\n*          On exit, the solution matrix X.\n*\n*  LDB     (input) INTEGER\n*          The leading dimension of the array B.  LDB >= max(1,N).\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit\n*          < 0:  if INFO = -i, the i-th argument had an illegal value\n*\n\n*  =====================================================================\n*\n\n");
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
  if (NA_TYPE(rb_ab) != NA_SCOMPLEX)
    rb_ab = na_change_type(rb_ab, NA_SCOMPLEX);
  ab = NA_PTR_TYPE(rb_ab, complex*);
  kl = NUM2INT(rb_kl);
  if (!NA_IsNArray(rb_b))
    rb_raise(rb_eArgError, "b (6th argument) must be NArray");
  if (NA_RANK(rb_b) != 2)
    rb_raise(rb_eArgError, "rank of b (6th argument) must be %d", 2);
  nrhs = NA_SHAPE1(rb_b);
  ldb = NA_SHAPE0(rb_b);
  if (NA_TYPE(rb_b) != NA_SCOMPLEX)
    rb_b = na_change_type(rb_b, NA_SCOMPLEX);
  b = NA_PTR_TYPE(rb_b, complex*);
  ku = NUM2INT(rb_ku);
  trans = StringValueCStr(rb_trans)[0];
  {
    int shape[2];
    shape[0] = ldb;
    shape[1] = nrhs;
    rb_b_out__ = na_make_object(NA_SCOMPLEX, 2, shape, cNArray);
  }
  b_out__ = NA_PTR_TYPE(rb_b_out__, complex*);
  MEMCPY(b_out__, b, complex, NA_TOTAL(rb_b));
  rb_b = rb_b_out__;
  b = b_out__;

  cgbtrs_(&trans, &n, &kl, &ku, &nrhs, ab, &ldab, ipiv, b, &ldb, &info);

  rb_info = INT2NUM(info);
  return rb_ary_new3(2, rb_info, rb_b);
}

void
init_lapack_cgbtrs(VALUE mLapack){
  rb_define_module_function(mLapack, "cgbtrs", rb_cgbtrs, -1);
}

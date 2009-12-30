#include "rb_lapack.h"

static VALUE
rb_cgbsv(int argc, VALUE *argv, VALUE self){
  VALUE rb_kl;
  integer kl; 
  VALUE rb_ku;
  integer ku; 
  VALUE rb_ab;
  complex *ab; 
  VALUE rb_b;
  complex *b; 
  VALUE rb_ipiv;
  integer *ipiv; 
  VALUE rb_info;
  integer info; 
  VALUE rb_ab_out__;
  complex *ab_out__;
  VALUE rb_b_out__;
  complex *b_out__;

  integer ldab;
  integer n;
  integer ldb;
  integer nrhs;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  ipiv, info, ab, b = NumRu::Lapack.cgbsv( kl, ku, ab, b)\n    or\n  NumRu::Lapack.cgbsv  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE CGBSV( N, KL, KU, NRHS, AB, LDAB, IPIV, B, LDB, INFO )\n\n*  Purpose\n*  =======\n*\n*  CGBSV computes the solution to a complex system of linear equations\n*  A * X = B, where A is a band matrix of order N with KL subdiagonals\n*  and KU superdiagonals, and X and B are N-by-NRHS matrices.\n*\n*  The LU decomposition with partial pivoting and row interchanges is\n*  used to factor A as A = L * U, where L is a product of permutation\n*  and unit lower triangular matrices with KL subdiagonals, and U is\n*  upper triangular with KL+KU superdiagonals.  The factored form of A\n*  is then used to solve the system of equations A * X = B.\n*\n\n*  Arguments\n*  =========\n*\n*  N       (input) INTEGER\n*          The number of linear equations, i.e., the order of the\n*          matrix A.  N >= 0.\n*\n*  KL      (input) INTEGER\n*          The number of subdiagonals within the band of A.  KL >= 0.\n*\n*  KU      (input) INTEGER\n*          The number of superdiagonals within the band of A.  KU >= 0.\n*\n*  NRHS    (input) INTEGER\n*          The number of right hand sides, i.e., the number of columns\n*          of the matrix B.  NRHS >= 0.\n*\n*  AB      (input/output) COMPLEX array, dimension (LDAB,N)\n*          On entry, the matrix A in band storage, in rows KL+1 to\n*          2*KL+KU+1; rows 1 to KL of the array need not be set.\n*          The j-th column of A is stored in the j-th column of the\n*          array AB as follows:\n*          AB(KL+KU+1+i-j,j) = A(i,j) for max(1,j-KU)<=i<=min(N,j+KL)\n*          On exit, details of the factorization: U is stored as an\n*          upper triangular band matrix with KL+KU superdiagonals in\n*          rows 1 to KL+KU+1, and the multipliers used during the\n*          factorization are stored in rows KL+KU+2 to 2*KL+KU+1.\n*          See below for further details.\n*\n*  LDAB    (input) INTEGER\n*          The leading dimension of the array AB.  LDAB >= 2*KL+KU+1.\n*\n*  IPIV    (output) INTEGER array, dimension (N)\n*          The pivot indices that define the permutation matrix P;\n*          row i of the matrix was interchanged with row IPIV(i).\n*\n*  B       (input/output) COMPLEX array, dimension (LDB,NRHS)\n*          On entry, the N-by-NRHS right hand side matrix B.\n*          On exit, if INFO = 0, the N-by-NRHS solution matrix X.\n*\n*  LDB     (input) INTEGER\n*          The leading dimension of the array B.  LDB >= max(1,N).\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit\n*          < 0:  if INFO = -i, the i-th argument had an illegal value\n*          > 0:  if INFO = i, U(i,i) is exactly zero.  The factorization\n*                has been completed, but the factor U is exactly\n*                singular, and the solution has not been computed.\n*\n\n*  Further Details\n*  ===============\n*\n*  The band storage scheme is illustrated by the following example, when\n*  M = N = 6, KL = 2, KU = 1:\n*\n*  On entry:                       On exit:\n*\n*      *    *    *    +    +    +       *    *    *   u14  u25  u36\n*      *    *    +    +    +    +       *    *   u13  u24  u35  u46\n*      *   a12  a23  a34  a45  a56      *   u12  u23  u34  u45  u56\n*     a11  a22  a33  a44  a55  a66     u11  u22  u33  u44  u55  u66\n*     a21  a32  a43  a54  a65   *      m21  m32  m43  m54  m65   *\n*     a31  a42  a53  a64   *    *      m31  m42  m53  m64   *    *\n*\n*  Array elements marked * are not used by the routine; elements marked\n*  + need not be set on entry, but are required by the routine to store\n*  elements of U because of fill-in resulting from the row interchanges.\n*\n*  =====================================================================\n*\n*     .. External Subroutines ..\n      EXTERNAL           CGBTRF, CGBTRS, XERBLA\n*     ..\n*     .. Intrinsic Functions ..\n      INTRINSIC          MAX\n*     ..\n\n");
    return Qnil;
  }
  if (argc != 4)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 4)", argc);
  rb_kl = argv[0];
  rb_ku = argv[1];
  rb_ab = argv[2];
  rb_b = argv[3];

  kl = NUM2INT(rb_kl);
  ku = NUM2INT(rb_ku);
  if (!NA_IsNArray(rb_ab))
    rb_raise(rb_eArgError, "ab (3th argument) must be NArray");
  if (NA_RANK(rb_ab) != 2)
    rb_raise(rb_eArgError, "rank of ab (3th argument) must be %d", 2);
  ldab = NA_SHAPE0(rb_ab);
  n = NA_SHAPE1(rb_ab);
  if (NA_TYPE(rb_ab) != NA_SCOMPLEX)
    rb_ab = na_change_type(rb_ab, NA_SCOMPLEX);
  ab = NA_PTR_TYPE(rb_ab, complex*);
  if (!NA_IsNArray(rb_b))
    rb_raise(rb_eArgError, "b (4th argument) must be NArray");
  if (NA_RANK(rb_b) != 2)
    rb_raise(rb_eArgError, "rank of b (4th argument) must be %d", 2);
  ldb = NA_SHAPE0(rb_b);
  nrhs = NA_SHAPE1(rb_b);
  if (NA_TYPE(rb_b) != NA_SCOMPLEX)
    rb_b = na_change_type(rb_b, NA_SCOMPLEX);
  b = NA_PTR_TYPE(rb_b, complex*);
  {
    int shape[1];
    shape[0] = DIM_LEN(n);
    rb_ipiv = na_make_object(NA_LINT, 1, shape, cNArray);
  }
  ipiv = NA_PTR_TYPE(rb_ipiv, integer*);
  {
    int shape[2];
    shape[0] = ldab;
    shape[1] = n;
    rb_ab_out__ = na_make_object(NA_SCOMPLEX, 2, shape, cNArray);
  }
  ab_out__ = NA_PTR_TYPE(rb_ab_out__, complex*);
  MEMCPY(ab_out__, ab, complex, NA_TOTAL(rb_ab));
  rb_ab = rb_ab_out__;
  ab = ab_out__;
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

  cgbsv_(&n, &kl, &ku, &nrhs, ab, &ldab, ipiv, b, &ldb, &info);

  rb_info = INT2NUM(info);
  return rb_ary_new3(4, rb_ipiv, rb_info, rb_ab, rb_b);
}

void
init_lapack_cgbsv(VALUE mLapack){
  rb_define_module_function(mLapack, "cgbsv", rb_cgbsv, -1);
}

#include "rb_lapack.h"

static VALUE
rb_cpttrs(int argc, VALUE *argv, VALUE self){
  VALUE rb_uplo;
  char uplo; 
  VALUE rb_d;
  real *d; 
  VALUE rb_e;
  complex *e; 
  VALUE rb_b;
  complex *b; 
  VALUE rb_info;
  integer info; 
  VALUE rb_b_out__;
  complex *b_out__;

  integer n;
  integer ldb;
  integer nrhs;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  info, b = NumRu::Lapack.cpttrs( uplo, d, e, b)\n    or\n  NumRu::Lapack.cpttrs  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE CPTTRS( UPLO, N, NRHS, D, E, B, LDB, INFO )\n\n*  Purpose\n*  =======\n*\n*  CPTTRS solves a tridiagonal system of the form\n*     A * X = B\n*  using the factorization A = U'*D*U or A = L*D*L' computed by CPTTRF.\n*  D is a diagonal matrix specified in the vector D, U (or L) is a unit\n*  bidiagonal matrix whose superdiagonal (subdiagonal) is specified in\n*  the vector E, and X and B are N by NRHS matrices.\n*\n\n*  Arguments\n*  =========\n*\n*  UPLO    (input) CHARACTER*1\n*          Specifies the form of the factorization and whether the\n*          vector E is the superdiagonal of the upper bidiagonal factor\n*          U or the subdiagonal of the lower bidiagonal factor L.\n*          = 'U':  A = U'*D*U, E is the superdiagonal of U\n*          = 'L':  A = L*D*L', E is the subdiagonal of L\n*\n*  N       (input) INTEGER\n*          The order of the tridiagonal matrix A.  N >= 0.\n*\n*  NRHS    (input) INTEGER\n*          The number of right hand sides, i.e., the number of columns\n*          of the matrix B.  NRHS >= 0.\n*\n*  D       (input) REAL array, dimension (N)\n*          The n diagonal elements of the diagonal matrix D from the\n*          factorization A = U'*D*U or A = L*D*L'.\n*\n*  E       (input) COMPLEX array, dimension (N-1)\n*          If UPLO = 'U', the (n-1) superdiagonal elements of the unit\n*          bidiagonal factor U from the factorization A = U'*D*U.\n*          If UPLO = 'L', the (n-1) subdiagonal elements of the unit\n*          bidiagonal factor L from the factorization A = L*D*L'.\n*\n*  B       (input/output) REAL array, dimension (LDB,NRHS)\n*          On entry, the right hand side vectors B for the system of\n*          linear equations.\n*          On exit, the solution vectors, X.\n*\n*  LDB     (input) INTEGER\n*          The leading dimension of the array B.  LDB >= max(1,N).\n*\n*  INFO    (output) INTEGER\n*          = 0: successful exit\n*          < 0: if INFO = -k, the k-th argument had an illegal value\n*\n\n*  =====================================================================\n*\n*     .. Local Scalars ..\n      LOGICAL            UPPER\n      INTEGER            IUPLO, J, JB, NB\n*     ..\n*     .. External Functions ..\n      INTEGER            ILAENV\n      EXTERNAL           ILAENV\n*     ..\n*     .. External Subroutines ..\n      EXTERNAL           CPTTS2, XERBLA\n*     ..\n*     .. Intrinsic Functions ..\n      INTRINSIC          MAX, MIN\n*     ..\n\n");
    return Qnil;
  }
  if (argc != 4)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 4)", argc);
  rb_uplo = argv[0];
  rb_d = argv[1];
  rb_e = argv[2];
  rb_b = argv[3];

  uplo = StringValueCStr(rb_uplo)[0];
  if (!NA_IsNArray(rb_d))
    rb_raise(rb_eArgError, "d (2th argument) must be NArray");
  if (NA_RANK(rb_d) != 1)
    rb_raise(rb_eArgError, "rank of d (2th argument) must be %d", 1);
  n = NA_SHAPE0(rb_d);
  if (NA_TYPE(rb_d) != NA_SFLOAT)
    rb_d = na_change_type(rb_d, NA_SFLOAT);
  d = NA_PTR_TYPE(rb_d, real*);
  if (!NA_IsNArray(rb_e))
    rb_raise(rb_eArgError, "e (3th argument) must be NArray");
  if (NA_RANK(rb_e) != 1)
    rb_raise(rb_eArgError, "rank of e (3th argument) must be %d", 1);
  if (NA_SHAPE0(rb_e) != (n-1))
    rb_raise(rb_eRuntimeError, "shape 0 of e must be %d", n-1);
  if (NA_TYPE(rb_e) != NA_SCOMPLEX)
    rb_e = na_change_type(rb_e, NA_SCOMPLEX);
  e = NA_PTR_TYPE(rb_e, complex*);
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
    int shape[2];
    shape[0] = ldb;
    shape[1] = nrhs;
    rb_b_out__ = na_make_object(NA_SCOMPLEX, 2, shape, cNArray);
  }
  b_out__ = NA_PTR_TYPE(rb_b_out__, complex*);
  MEMCPY(b_out__, b, complex, NA_TOTAL(rb_b));
  rb_b = rb_b_out__;
  b = b_out__;

  cpttrs_(&uplo, &n, &nrhs, d, e, b, &ldb, &info);

  rb_info = INT2NUM(info);
  return rb_ary_new3(2, rb_info, rb_b);
}

void
init_lapack_cpttrs(VALUE mLapack){
  rb_define_module_function(mLapack, "cpttrs", rb_cpttrs, -1);
}

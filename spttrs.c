#include "rb_lapack.h"

extern VOID spttrs_(integer *n, integer *nrhs, real *d, real *e, real *b, integer *ldb, integer *info);

static VALUE
rb_spttrs(int argc, VALUE *argv, VALUE self){
  VALUE rb_d;
  real *d; 
  VALUE rb_e;
  real *e; 
  VALUE rb_b;
  real *b; 
  VALUE rb_info;
  integer info; 
  VALUE rb_b_out__;
  real *b_out__;

  integer n;
  integer ldb;
  integer nrhs;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  info, b = NumRu::Lapack.spttrs( d, e, b)\n    or\n  NumRu::Lapack.spttrs  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE SPTTRS( N, NRHS, D, E, B, LDB, INFO )\n\n*  Purpose\n*  =======\n*\n*  SPTTRS solves a tridiagonal system of the form\n*     A * X = B\n*  using the L*D*L' factorization of A computed by SPTTRF.  D is a\n*  diagonal matrix specified in the vector D, L is a unit bidiagonal\n*  matrix whose subdiagonal is specified in the vector E, and X and B\n*  are N by NRHS matrices.\n*\n\n*  Arguments\n*  =========\n*\n*  N       (input) INTEGER\n*          The order of the tridiagonal matrix A.  N >= 0.\n*\n*  NRHS    (input) INTEGER\n*          The number of right hand sides, i.e., the number of columns\n*          of the matrix B.  NRHS >= 0.\n*\n*  D       (input) REAL array, dimension (N)\n*          The n diagonal elements of the diagonal matrix D from the\n*          L*D*L' factorization of A.\n*\n*  E       (input) REAL array, dimension (N-1)\n*          The (n-1) subdiagonal elements of the unit bidiagonal factor\n*          L from the L*D*L' factorization of A.  E can also be regarded\n*          as the superdiagonal of the unit bidiagonal factor U from the\n*          factorization A = U'*D*U.\n*\n*  B       (input/output) REAL array, dimension (LDB,NRHS)\n*          On entry, the right hand side vectors B for the system of\n*          linear equations.\n*          On exit, the solution vectors, X.\n*\n*  LDB     (input) INTEGER\n*          The leading dimension of the array B.  LDB >= max(1,N).\n*\n*  INFO    (output) INTEGER\n*          = 0: successful exit\n*          < 0: if INFO = -k, the k-th argument had an illegal value\n*\n\n*  =====================================================================\n*\n*     .. Local Scalars ..\n      INTEGER            J, JB, NB\n*     ..\n*     .. External Functions ..\n      INTEGER            ILAENV\n      EXTERNAL           ILAENV\n*     ..\n*     .. External Subroutines ..\n      EXTERNAL           SPTTS2, XERBLA\n*     ..\n*     .. Intrinsic Functions ..\n      INTRINSIC          MAX, MIN\n*     ..\n\n");
    return Qnil;
  }
  if (argc != 3)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 3)", argc);
  rb_d = argv[0];
  rb_e = argv[1];
  rb_b = argv[2];

  if (!NA_IsNArray(rb_b))
    rb_raise(rb_eArgError, "b (3th argument) must be NArray");
  if (NA_RANK(rb_b) != 2)
    rb_raise(rb_eArgError, "rank of b (3th argument) must be %d", 2);
  nrhs = NA_SHAPE1(rb_b);
  ldb = NA_SHAPE0(rb_b);
  if (NA_TYPE(rb_b) != NA_SFLOAT)
    rb_b = na_change_type(rb_b, NA_SFLOAT);
  b = NA_PTR_TYPE(rb_b, real*);
  if (!NA_IsNArray(rb_d))
    rb_raise(rb_eArgError, "d (1th argument) must be NArray");
  if (NA_RANK(rb_d) != 1)
    rb_raise(rb_eArgError, "rank of d (1th argument) must be %d", 1);
  n = NA_SHAPE0(rb_d);
  if (NA_TYPE(rb_d) != NA_SFLOAT)
    rb_d = na_change_type(rb_d, NA_SFLOAT);
  d = NA_PTR_TYPE(rb_d, real*);
  if (!NA_IsNArray(rb_e))
    rb_raise(rb_eArgError, "e (2th argument) must be NArray");
  if (NA_RANK(rb_e) != 1)
    rb_raise(rb_eArgError, "rank of e (2th argument) must be %d", 1);
  if (NA_SHAPE0(rb_e) != (n-1))
    rb_raise(rb_eRuntimeError, "shape 0 of e must be %d", n-1);
  if (NA_TYPE(rb_e) != NA_SFLOAT)
    rb_e = na_change_type(rb_e, NA_SFLOAT);
  e = NA_PTR_TYPE(rb_e, real*);
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

  spttrs_(&n, &nrhs, d, e, b, &ldb, &info);

  rb_info = INT2NUM(info);
  return rb_ary_new3(2, rb_info, rb_b);
}

void
init_lapack_spttrs(VALUE mLapack){
  rb_define_module_function(mLapack, "spttrs", rb_spttrs, -1);
}

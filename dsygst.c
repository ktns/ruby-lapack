#include "rb_lapack.h"

extern VOID dsygst_(integer *itype, char *uplo, integer *n, doublereal *a, integer *lda, doublereal *b, integer *ldb, integer *info);

static VALUE
rb_dsygst(int argc, VALUE *argv, VALUE self){
  VALUE rb_itype;
  integer itype; 
  VALUE rb_uplo;
  char uplo; 
  VALUE rb_a;
  doublereal *a; 
  VALUE rb_b;
  doublereal *b; 
  VALUE rb_info;
  integer info; 
  VALUE rb_a_out__;
  doublereal *a_out__;

  integer lda;
  integer n;
  integer ldb;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  info, a = NumRu::Lapack.dsygst( itype, uplo, a, b)\n    or\n  NumRu::Lapack.dsygst  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE DSYGST( ITYPE, UPLO, N, A, LDA, B, LDB, INFO )\n\n*  Purpose\n*  =======\n*\n*  DSYGST reduces a real symmetric-definite generalized eigenproblem\n*  to standard form.\n*\n*  If ITYPE = 1, the problem is A*x = lambda*B*x,\n*  and A is overwritten by inv(U**T)*A*inv(U) or inv(L)*A*inv(L**T)\n*\n*  If ITYPE = 2 or 3, the problem is A*B*x = lambda*x or\n*  B*A*x = lambda*x, and A is overwritten by U*A*U**T or L**T*A*L.\n*\n*  B must have been previously factorized as U**T*U or L*L**T by DPOTRF.\n*\n\n*  Arguments\n*  =========\n*\n*  ITYPE   (input) INTEGER\n*          = 1: compute inv(U**T)*A*inv(U) or inv(L)*A*inv(L**T);\n*          = 2 or 3: compute U*A*U**T or L**T*A*L.\n*\n*  UPLO    (input) CHARACTER*1\n*          = 'U':  Upper triangle of A is stored and B is factored as\n*                  U**T*U;\n*          = 'L':  Lower triangle of A is stored and B is factored as\n*                  L*L**T.\n*\n*  N       (input) INTEGER\n*          The order of the matrices A and B.  N >= 0.\n*\n*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)\n*          On entry, the symmetric matrix A.  If UPLO = 'U', the leading\n*          N-by-N upper triangular part of A contains the upper\n*          triangular part of the matrix A, and the strictly lower\n*          triangular part of A is not referenced.  If UPLO = 'L', the\n*          leading N-by-N lower triangular part of A contains the lower\n*          triangular part of the matrix A, and the strictly upper\n*          triangular part of A is not referenced.\n*\n*          On exit, if INFO = 0, the transformed matrix, stored in the\n*          same format as A.\n*\n*  LDA     (input) INTEGER\n*          The leading dimension of the array A.  LDA >= max(1,N).\n*\n*  B       (input) DOUBLE PRECISION array, dimension (LDB,N)\n*          The triangular factor from the Cholesky factorization of B,\n*          as returned by DPOTRF.\n*\n*  LDB     (input) INTEGER\n*          The leading dimension of the array B.  LDB >= max(1,N).\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit\n*          < 0:  if INFO = -i, the i-th argument had an illegal value\n*\n\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 4)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 4)", argc);
  rb_itype = argv[0];
  rb_uplo = argv[1];
  rb_a = argv[2];
  rb_b = argv[3];

  if (!NA_IsNArray(rb_a))
    rb_raise(rb_eArgError, "a (3th argument) must be NArray");
  if (NA_RANK(rb_a) != 2)
    rb_raise(rb_eArgError, "rank of a (3th argument) must be %d", 2);
  n = NA_SHAPE1(rb_a);
  lda = NA_SHAPE0(rb_a);
  if (NA_TYPE(rb_a) != NA_DFLOAT)
    rb_a = na_change_type(rb_a, NA_DFLOAT);
  a = NA_PTR_TYPE(rb_a, doublereal*);
  if (!NA_IsNArray(rb_b))
    rb_raise(rb_eArgError, "b (4th argument) must be NArray");
  if (NA_RANK(rb_b) != 2)
    rb_raise(rb_eArgError, "rank of b (4th argument) must be %d", 2);
  if (NA_SHAPE1(rb_b) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of b must be the same as shape 1 of a");
  ldb = NA_SHAPE0(rb_b);
  if (NA_TYPE(rb_b) != NA_DFLOAT)
    rb_b = na_change_type(rb_b, NA_DFLOAT);
  b = NA_PTR_TYPE(rb_b, doublereal*);
  itype = NUM2INT(rb_itype);
  uplo = StringValueCStr(rb_uplo)[0];
  {
    int shape[2];
    shape[0] = lda;
    shape[1] = n;
    rb_a_out__ = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  a_out__ = NA_PTR_TYPE(rb_a_out__, doublereal*);
  MEMCPY(a_out__, a, doublereal, NA_TOTAL(rb_a));
  rb_a = rb_a_out__;
  a = a_out__;

  dsygst_(&itype, &uplo, &n, a, &lda, b, &ldb, &info);

  rb_info = INT2NUM(info);
  return rb_ary_new3(2, rb_info, rb_a);
}

void
init_lapack_dsygst(VALUE mLapack){
  rb_define_module_function(mLapack, "dsygst", rb_dsygst, -1);
}

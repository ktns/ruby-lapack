#include "rb_lapack.h"

static VALUE
rb_claset(int argc, VALUE *argv, VALUE self){
  VALUE rb_uplo;
  char uplo; 
  VALUE rb_m;
  integer m; 
  VALUE rb_alpha;
  complex alpha; 
  VALUE rb_beta;
  complex beta; 
  VALUE rb_a;
  complex *a; 
  VALUE rb_a_out__;
  complex *a_out__;

  integer lda;
  integer n;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  a = NumRu::Lapack.claset( uplo, m, alpha, beta, a)\n    or\n  NumRu::Lapack.claset  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE CLASET( UPLO, M, N, ALPHA, BETA, A, LDA )\n\n*  Purpose\n*  =======\n*\n*  CLASET initializes a 2-D array A to BETA on the diagonal and\n*  ALPHA on the offdiagonals.\n*\n\n*  Arguments\n*  =========\n*\n*  UPLO    (input) CHARACTER*1\n*          Specifies the part of the matrix A to be set.\n*          = 'U':      Upper triangular part is set. The lower triangle\n*                      is unchanged.\n*          = 'L':      Lower triangular part is set. The upper triangle\n*                      is unchanged.\n*          Otherwise:  All of the matrix A is set.\n*\n*  M       (input) INTEGER\n*          On entry, M specifies the number of rows of A.\n*\n*  N       (input) INTEGER\n*          On entry, N specifies the number of columns of A.\n*\n*  ALPHA   (input) COMPLEX\n*          All the offdiagonal array elements are set to ALPHA.\n*\n*  BETA    (input) COMPLEX\n*          All the diagonal array elements are set to BETA.\n*\n*  A       (input/output) COMPLEX array, dimension (LDA,N)\n*          On entry, the m by n matrix A.\n*          On exit, A(i,j) = ALPHA, 1 <= i <= m, 1 <= j <= n, i.ne.j;\n*                   A(i,i) = BETA , 1 <= i <= min(m,n)\n*\n*  LDA     (input) INTEGER\n*          The leading dimension of the array A.  LDA >= max(1,M).\n*\n\n*  =====================================================================\n*\n*     .. Local Scalars ..\n      INTEGER            I, J\n*     ..\n*     .. External Functions ..\n      LOGICAL            LSAME\n      EXTERNAL           LSAME\n*     ..\n*     .. Intrinsic Functions ..\n      INTRINSIC          MIN\n*     ..\n\n");
    return Qnil;
  }
  if (argc != 5)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 5)", argc);
  rb_uplo = argv[0];
  rb_m = argv[1];
  rb_alpha = argv[2];
  rb_beta = argv[3];
  rb_a = argv[4];

  uplo = StringValueCStr(rb_uplo)[0];
  m = NUM2INT(rb_m);
  alpha.r = (real)NUM2DBL(rb_funcall(rb_alpha, rb_intern("real"), 0));
  alpha.i = (real)NUM2DBL(rb_funcall(rb_alpha, rb_intern("imag"), 0));
  beta.r = (real)NUM2DBL(rb_funcall(rb_beta, rb_intern("real"), 0));
  beta.i = (real)NUM2DBL(rb_funcall(rb_beta, rb_intern("imag"), 0));
  if (!NA_IsNArray(rb_a))
    rb_raise(rb_eArgError, "a (5th argument) must be NArray");
  if (NA_RANK(rb_a) != 2)
    rb_raise(rb_eArgError, "rank of a (5th argument) must be %d", 2);
  lda = NA_SHAPE0(rb_a);
  n = NA_SHAPE1(rb_a);
  if (NA_TYPE(rb_a) != NA_SCOMPLEX)
    rb_a = na_change_type(rb_a, NA_SCOMPLEX);
  a = NA_PTR_TYPE(rb_a, complex*);
  {
    int shape[2];
    shape[0] = lda;
    shape[1] = n;
    rb_a_out__ = na_make_object(NA_SCOMPLEX, 2, shape, cNArray);
  }
  a_out__ = NA_PTR_TYPE(rb_a_out__, complex*);
  MEMCPY(a_out__, a, complex, NA_TOTAL(rb_a));
  rb_a = rb_a_out__;
  a = a_out__;

  claset_(&uplo, &m, &n, &alpha, &beta, a, &lda);

  return rb_a;
}

void
init_lapack_claset(VALUE mLapack){
  rb_define_module_function(mLapack, "claset", rb_claset, -1);
}

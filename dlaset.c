#include "rb_lapack.h"

static VALUE
rb_dlaset(int argc, VALUE *argv, VALUE self){
  VALUE rb_uplo;
  char uplo; 
  VALUE rb_m;
  integer m; 
  VALUE rb_alpha;
  doublereal alpha; 
  VALUE rb_beta;
  doublereal beta; 
  VALUE rb_a;
  doublereal *a; 
  VALUE rb_a_out__;
  doublereal *a_out__;

  integer lda;
  integer n;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  a = NumRu::Lapack.dlaset( uplo, m, alpha, beta, a)\n    or\n  NumRu::Lapack.dlaset  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE DLASET( UPLO, M, N, ALPHA, BETA, A, LDA )\n\n*  Purpose\n*  =======\n*\n*  DLASET initializes an m-by-n matrix A to BETA on the diagonal and\n*  ALPHA on the offdiagonals.\n*\n\n*  Arguments\n*  =========\n*\n*  UPLO    (input) CHARACTER*1\n*          Specifies the part of the matrix A to be set.\n*          = 'U':      Upper triangular part is set; the strictly lower\n*                      triangular part of A is not changed.\n*          = 'L':      Lower triangular part is set; the strictly upper\n*                      triangular part of A is not changed.\n*          Otherwise:  All of the matrix A is set.\n*\n*  M       (input) INTEGER\n*          The number of rows of the matrix A.  M >= 0.\n*\n*  N       (input) INTEGER\n*          The number of columns of the matrix A.  N >= 0.\n*\n*  ALPHA   (input) DOUBLE PRECISION\n*          The constant to which the offdiagonal elements are to be set.\n*\n*  BETA    (input) DOUBLE PRECISION\n*          The constant to which the diagonal elements are to be set.\n*\n*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)\n*          On exit, the leading m-by-n submatrix of A is set as follows:\n*\n*          if UPLO = 'U', A(i,j) = ALPHA, 1<=i<=j-1, 1<=j<=n,\n*          if UPLO = 'L', A(i,j) = ALPHA, j+1<=i<=m, 1<=j<=n,\n*          otherwise,     A(i,j) = ALPHA, 1<=i<=m, 1<=j<=n, i.ne.j,\n*\n*          and, for all UPLO, A(i,i) = BETA, 1<=i<=min(m,n).\n*\n*  LDA     (input) INTEGER\n*          The leading dimension of the array A.  LDA >= max(1,M).\n*\n\n* =====================================================================\n*\n*     .. Local Scalars ..\n      INTEGER            I, J\n*     ..\n*     .. External Functions ..\n      LOGICAL            LSAME\n      EXTERNAL           LSAME\n*     ..\n*     .. Intrinsic Functions ..\n      INTRINSIC          MIN\n*     ..\n\n");
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
  alpha = NUM2DBL(rb_alpha);
  beta = NUM2DBL(rb_beta);
  if (!NA_IsNArray(rb_a))
    rb_raise(rb_eArgError, "a (5th argument) must be NArray");
  if (NA_RANK(rb_a) != 2)
    rb_raise(rb_eArgError, "rank of a (5th argument) must be %d", 2);
  lda = NA_SHAPE0(rb_a);
  n = NA_SHAPE1(rb_a);
  if (NA_TYPE(rb_a) != NA_DFLOAT)
    rb_a = na_change_type(rb_a, NA_DFLOAT);
  a = NA_PTR_TYPE(rb_a, doublereal*);
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

  dlaset_(&uplo, &m, &n, &alpha, &beta, a, &lda);

  return rb_a;
}

void
init_lapack_dlaset(VALUE mLapack){
  rb_define_module_function(mLapack, "dlaset", rb_dlaset, -1);
}

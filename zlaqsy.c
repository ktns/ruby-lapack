#include "rb_lapack.h"

extern VOID zlaqsy_(char *uplo, integer *n, doublecomplex *a, integer *lda, doublereal *s, doublereal *scond, doublereal *amax, char *equed);

static VALUE
rb_zlaqsy(int argc, VALUE *argv, VALUE self){
  VALUE rb_uplo;
  char uplo; 
  VALUE rb_a;
  doublecomplex *a; 
  VALUE rb_s;
  doublereal *s; 
  VALUE rb_scond;
  doublereal scond; 
  VALUE rb_amax;
  doublereal amax; 
  VALUE rb_equed;
  char equed; 
  VALUE rb_a_out__;
  doublecomplex *a_out__;

  integer lda;
  integer n;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  equed, a = NumRu::Lapack.zlaqsy( uplo, a, s, scond, amax)\n    or\n  NumRu::Lapack.zlaqsy  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE ZLAQSY( UPLO, N, A, LDA, S, SCOND, AMAX, EQUED )\n\n*  Purpose\n*  =======\n*\n*  ZLAQSY equilibrates a symmetric matrix A using the scaling factors\n*  in the vector S.\n*\n\n*  Arguments\n*  =========\n*\n*  UPLO    (input) CHARACTER*1\n*          Specifies whether the upper or lower triangular part of the\n*          symmetric matrix A is stored.\n*          = 'U':  Upper triangular\n*          = 'L':  Lower triangular\n*\n*  N       (input) INTEGER\n*          The order of the matrix A.  N >= 0.\n*\n*  A       (input/output) COMPLEX*16 array, dimension (LDA,N)\n*          On entry, the symmetric matrix A.  If UPLO = 'U', the leading\n*          n by n upper triangular part of A contains the upper\n*          triangular part of the matrix A, and the strictly lower\n*          triangular part of A is not referenced.  If UPLO = 'L', the\n*          leading n by n lower triangular part of A contains the lower\n*          triangular part of the matrix A, and the strictly upper\n*          triangular part of A is not referenced.\n*\n*          On exit, if EQUED = 'Y', the equilibrated matrix:\n*          diag(S) * A * diag(S).\n*\n*  LDA     (input) INTEGER\n*          The leading dimension of the array A.  LDA >= max(N,1).\n*\n*  S       (input) DOUBLE PRECISION array, dimension (N)\n*          The scale factors for A.\n*\n*  SCOND   (input) DOUBLE PRECISION\n*          Ratio of the smallest S(i) to the largest S(i).\n*\n*  AMAX    (input) DOUBLE PRECISION\n*          Absolute value of largest matrix entry.\n*\n*  EQUED   (output) CHARACTER*1\n*          Specifies whether or not equilibration was done.\n*          = 'N':  No equilibration.\n*          = 'Y':  Equilibration was done, i.e., A has been replaced by\n*                  diag(S) * A * diag(S).\n*\n*  Internal Parameters\n*  ===================\n*\n*  THRESH is a threshold value used to decide if scaling should be done\n*  based on the ratio of the scaling factors.  If SCOND < THRESH,\n*  scaling is done.\n*\n*  LARGE and SMALL are threshold values used to decide if scaling should\n*  be done based on the absolute size of the largest matrix element.\n*  If AMAX > LARGE or AMAX < SMALL, scaling is done.\n*\n\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 5)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 5)", argc);
  rb_uplo = argv[0];
  rb_a = argv[1];
  rb_s = argv[2];
  rb_scond = argv[3];
  rb_amax = argv[4];

  scond = NUM2DBL(rb_scond);
  if (!NA_IsNArray(rb_a))
    rb_raise(rb_eArgError, "a (2th argument) must be NArray");
  if (NA_RANK(rb_a) != 2)
    rb_raise(rb_eArgError, "rank of a (2th argument) must be %d", 2);
  n = NA_SHAPE1(rb_a);
  lda = NA_SHAPE0(rb_a);
  if (NA_TYPE(rb_a) != NA_DCOMPLEX)
    rb_a = na_change_type(rb_a, NA_DCOMPLEX);
  a = NA_PTR_TYPE(rb_a, doublecomplex*);
  amax = NUM2DBL(rb_amax);
  if (!NA_IsNArray(rb_s))
    rb_raise(rb_eArgError, "s (3th argument) must be NArray");
  if (NA_RANK(rb_s) != 1)
    rb_raise(rb_eArgError, "rank of s (3th argument) must be %d", 1);
  if (NA_SHAPE0(rb_s) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of s must be the same as shape 1 of a");
  if (NA_TYPE(rb_s) != NA_DFLOAT)
    rb_s = na_change_type(rb_s, NA_DFLOAT);
  s = NA_PTR_TYPE(rb_s, doublereal*);
  uplo = StringValueCStr(rb_uplo)[0];
  {
    int shape[2];
    shape[0] = lda;
    shape[1] = n;
    rb_a_out__ = na_make_object(NA_DCOMPLEX, 2, shape, cNArray);
  }
  a_out__ = NA_PTR_TYPE(rb_a_out__, doublecomplex*);
  MEMCPY(a_out__, a, doublecomplex, NA_TOTAL(rb_a));
  rb_a = rb_a_out__;
  a = a_out__;

  zlaqsy_(&uplo, &n, a, &lda, s, &scond, &amax, &equed);

  rb_equed = rb_str_new(&equed,1);
  return rb_ary_new3(2, rb_equed, rb_a);
}

void
init_lapack_zlaqsy(VALUE mLapack){
  rb_define_module_function(mLapack, "zlaqsy", rb_zlaqsy, -1);
}

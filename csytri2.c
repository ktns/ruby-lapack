#include "rb_lapack.h"

extern VOID csytri2_(char *uplo, integer *n, complex *a, integer *lda, integer *ipiv, complex *work, integer *lwork, integer *info);

static VALUE
rb_csytri2(int argc, VALUE *argv, VALUE self){
  VALUE rb_uplo;
  char uplo; 
  VALUE rb_a;
  complex *a; 
  VALUE rb_ipiv;
  integer *ipiv; 
  VALUE rb_lwork;
  integer lwork; 
  VALUE rb_info;
  integer info; 
  VALUE rb_a_out__;
  complex *a_out__;
  integer c__1;
  integer nb;
  integer c__m1;
  complex *work;

  integer lda;
  integer n;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  info, a = NumRu::Lapack.csytri2( uplo, a, ipiv, lwork)\n    or\n  NumRu::Lapack.csytri2  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE CSYTRI2( UPLO, N, A, LDA, IPIV, WORK, LWORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  CSYTRI2 computes the inverse of a complex symmetric indefinite matrix\n*  A using the factorization A = U*D*U**T or A = L*D*L**T computed by\n*  CSYTRF. CSYTRI2 sets the LEADING DIMENSION of the workspace\n*  before calling CSYTRI2X that actually computes the inverse.\n*\n\n*  Arguments\n*  =========\n*\n*  UPLO    (input) CHARACTER*1\n*          Specifies whether the details of the factorization are stored\n*          as an upper or lower triangular matrix.\n*          = 'U':  Upper triangular, form is A = U*D*U**T;\n*          = 'L':  Lower triangular, form is A = L*D*L**T.\n*\n*  N       (input) INTEGER\n*          The order of the matrix A.  N >= 0.\n*\n*  A       (input/output) COMPLEX array, dimension (LDA,N)\n*          On entry, the NB diagonal matrix D and the multipliers\n*          used to obtain the factor U or L as computed by CSYTRF.\n*\n*          On exit, if INFO = 0, the (symmetric) inverse of the original\n*          matrix.  If UPLO = 'U', the upper triangular part of the\n*          inverse is formed and the part of A below the diagonal is not\n*          referenced; if UPLO = 'L' the lower triangular part of the\n*          inverse is formed and the part of A above the diagonal is\n*          not referenced.\n*\n*  LDA     (input) INTEGER\n*          The leading dimension of the array A.  LDA >= max(1,N).\n*\n*  IPIV    (input) INTEGER array, dimension (N)\n*          Details of the interchanges and the NB structure of D\n*          as determined by CSYTRF.\n*\n*  WORK    (workspace) COMPLEX array, dimension (N+NB+1)*(NB+3)\n*\n*  LWORK   (input) INTEGER\n*          The dimension of the array WORK.\n*          WORK is size >= (N+NB+1)*(NB+3)\n*          If LDWORK = -1, then a workspace query is assumed; the routine\n*           calculates:\n*              - the optimal size of the WORK array, returns\n*          this value as the first entry of the WORK array,\n*              - and no error message related to LDWORK is issued by XERBLA.\n*\n*  INFO    (output) INTEGER\n*          = 0: successful exit\n*          < 0: if INFO = -i, the i-th argument had an illegal value\n*          > 0: if INFO = i, D(i,i) = 0; the matrix is singular and its\n*               inverse could not be computed.\n*\n\n*  =====================================================================\n*\n*     .. Local Scalars ..\n      LOGICAL            UPPER, LQUERY\n      INTEGER            MINSIZE, NBMAX\n*     ..\n*     .. External Functions ..\n      LOGICAL            LSAME\n      INTEGER            ILAENV\n      EXTERNAL           LSAME, ILAENV\n*     ..\n*     .. External Subroutines ..\n      EXTERNAL           CSYTRI2X\n*     ..\n\n");
    return Qnil;
  }
  if (argc != 4)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 4)", argc);
  rb_uplo = argv[0];
  rb_a = argv[1];
  rb_ipiv = argv[2];
  rb_lwork = argv[3];

  if (!NA_IsNArray(rb_ipiv))
    rb_raise(rb_eArgError, "ipiv (3th argument) must be NArray");
  if (NA_RANK(rb_ipiv) != 1)
    rb_raise(rb_eArgError, "rank of ipiv (3th argument) must be %d", 1);
  n = NA_SHAPE0(rb_ipiv);
  if (NA_TYPE(rb_ipiv) != NA_LINT)
    rb_ipiv = na_change_type(rb_ipiv, NA_LINT);
  ipiv = NA_PTR_TYPE(rb_ipiv, integer*);
  if (!NA_IsNArray(rb_a))
    rb_raise(rb_eArgError, "a (2th argument) must be NArray");
  if (NA_RANK(rb_a) != 2)
    rb_raise(rb_eArgError, "rank of a (2th argument) must be %d", 2);
  if (NA_SHAPE1(rb_a) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of a must be the same as shape 0 of ipiv");
  lda = NA_SHAPE0(rb_a);
  if (NA_TYPE(rb_a) != NA_SCOMPLEX)
    rb_a = na_change_type(rb_a, NA_SCOMPLEX);
  a = NA_PTR_TYPE(rb_a, complex*);
  c__1 = 1;
  uplo = StringValueCStr(rb_uplo)[0];
  c__m1 = -1;
  lwork = NUM2INT(rb_lwork);
  nb = ilaenv_(&c__1, "CSYTRF", &uplo, &n, &c__m1, &c__m1, &c__m1);
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
  work = ALLOC_N(complex, ((n+nb+1)*(nb+3)));

  csytri2_(&uplo, &n, a, &lda, ipiv, work, &lwork, &info);

  free(work);
  rb_info = INT2NUM(info);
  return rb_ary_new3(2, rb_info, rb_a);
}

void
init_lapack_csytri2(VALUE mLapack){
  rb_define_module_function(mLapack, "csytri2", rb_csytri2, -1);
}

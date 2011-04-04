#include "rb_lapack.h"

extern VOID dlat2s_(char *uplo, integer *n, doublereal *a, integer *lda, real *sa, integer *ldsa, integer *info);

static VALUE
rb_dlat2s(int argc, VALUE *argv, VALUE self){
  VALUE rb_uplo;
  char uplo; 
  VALUE rb_a;
  doublereal *a; 
  VALUE rb_sa;
  real *sa; 
  VALUE rb_info;
  integer info; 

  integer lda;
  integer n;
  integer ldsa;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  sa, info = NumRu::Lapack.dlat2s( uplo, a)\n    or\n  NumRu::Lapack.dlat2s  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE DLAT2S( UPLO, N, A, LDA, SA, LDSA, INFO )\n\n*  Purpose\n*  =======\n*\n*  DLAT2S converts a DOUBLE PRECISION triangular matrix, SA, to a SINGLE\n*  PRECISION triangular matrix, A.\n*\n*  RMAX is the overflow for the SINGLE PRECISION arithmetic\n*  DLAS2S checks that all the entries of A are between -RMAX and\n*  RMAX. If not the convertion is aborted and a flag is raised.\n*\n*  This is an auxiliary routine so there is no argument checking.\n*\n\n*  Arguments\n*  =========\n*\n*  UPLO    (input) CHARACTER*1\n*          = 'U':  A is upper triangular;\n*          = 'L':  A is lower triangular.\n*\n*  N       (input) INTEGER\n*          The number of rows and columns of the matrix A.  N >= 0.\n*\n*  A       (input) DOUBLE PRECISION array, dimension (LDA,N)\n*          On entry, the N-by-N triangular coefficient matrix A.\n*\n*  LDA     (input) INTEGER\n*          The leading dimension of the array A.  LDA >= max(1,N).\n*\n*  SA      (output) REAL array, dimension (LDSA,N)\n*          Only the UPLO part of SA is referenced.  On exit, if INFO=0,\n*          the N-by-N coefficient matrix SA; if INFO>0, the content of\n*          the UPLO part of SA is unspecified.\n*\n*  LDSA    (input) INTEGER\n*          The leading dimension of the array SA.  LDSA >= max(1,M).\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit.\n*          = 1:  an entry of the matrix A is greater than the SINGLE\n*                PRECISION overflow threshold, in this case, the content\n*                of the UPLO part of SA in exit is unspecified.\n*\n*  =========\n*\n*     .. Local Scalars ..\n      INTEGER            I, J\n      DOUBLE PRECISION   RMAX\n      LOGICAL            UPPER\n*     ..\n*     .. External Functions ..\n      REAL               SLAMCH\n      LOGICAL            LSAME\n      EXTERNAL           SLAMCH, LSAME\n*     ..\n\n");
    return Qnil;
  }
  if (argc != 2)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 2)", argc);
  rb_uplo = argv[0];
  rb_a = argv[1];

  if (!NA_IsNArray(rb_a))
    rb_raise(rb_eArgError, "a (2th argument) must be NArray");
  if (NA_RANK(rb_a) != 2)
    rb_raise(rb_eArgError, "rank of a (2th argument) must be %d", 2);
  n = NA_SHAPE1(rb_a);
  lda = NA_SHAPE0(rb_a);
  if (NA_TYPE(rb_a) != NA_DFLOAT)
    rb_a = na_change_type(rb_a, NA_DFLOAT);
  a = NA_PTR_TYPE(rb_a, doublereal*);
  uplo = StringValueCStr(rb_uplo)[0];
  ldsa = MAX(1,n);
  {
    int shape[2];
    shape[0] = ldsa;
    shape[1] = n;
    rb_sa = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  sa = NA_PTR_TYPE(rb_sa, real*);

  dlat2s_(&uplo, &n, a, &lda, sa, &ldsa, &info);

  rb_info = INT2NUM(info);
  return rb_ary_new3(2, rb_sa, rb_info);
}

void
init_lapack_dlat2s(VALUE mLapack){
  rb_define_module_function(mLapack, "dlat2s", rb_dlat2s, -1);
}

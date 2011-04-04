#include "rb_lapack.h"

extern VOID zlag2c_(integer *m, integer *n, doublecomplex *a, integer *lda, complex *sa, integer *ldsa, integer *info);

static VALUE
rb_zlag2c(int argc, VALUE *argv, VALUE self){
  VALUE rb_m;
  integer m; 
  VALUE rb_a;
  doublecomplex *a; 
  VALUE rb_sa;
  complex *sa; 
  VALUE rb_info;
  integer info; 

  integer lda;
  integer n;
  integer ldsa;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  sa, info = NumRu::Lapack.zlag2c( m, a)\n    or\n  NumRu::Lapack.zlag2c  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE ZLAG2C( M, N, A, LDA, SA, LDSA, INFO )\n\n*  Purpose\n*  =======\n*\n*  ZLAG2C converts a COMPLEX*16 matrix, SA, to a COMPLEX matrix, A.\n*\n*  RMAX is the overflow for the SINGLE PRECISION arithmetic\n*  ZLAG2C checks that all the entries of A are between -RMAX and\n*  RMAX. If not the convertion is aborted and a flag is raised.\n*\n*  This is an auxiliary routine so there is no argument checking.\n*\n\n*  Arguments\n*  =========\n*\n*  M       (input) INTEGER\n*          The number of lines of the matrix A.  M >= 0.\n*\n*  N       (input) INTEGER\n*          The number of columns of the matrix A.  N >= 0.\n*\n*  A       (input) COMPLEX*16 array, dimension (LDA,N)\n*          On entry, the M-by-N coefficient matrix A.\n*\n*  LDA     (input) INTEGER\n*          The leading dimension of the array A.  LDA >= max(1,M).\n*\n*  SA      (output) COMPLEX array, dimension (LDSA,N)\n*          On exit, if INFO=0, the M-by-N coefficient matrix SA; if\n*          INFO>0, the content of SA is unspecified.\n*\n*  LDSA    (input) INTEGER\n*          The leading dimension of the array SA.  LDSA >= max(1,M).\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit.\n*          = 1:  an entry of the matrix A is greater than the SINGLE\n*                PRECISION overflow threshold, in this case, the content\n*                of SA in exit is unspecified.\n*\n*  =========\n*\n*     .. Local Scalars ..\n      INTEGER            I, J\n      DOUBLE PRECISION   RMAX\n*     ..\n*     .. Intrinsic Functions ..\n      INTRINSIC          DBLE, DIMAG\n*     ..\n*     .. External Functions ..\n      REAL               SLAMCH\n      EXTERNAL           SLAMCH\n*     ..\n\n");
    return Qnil;
  }
  if (argc != 2)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 2)", argc);
  rb_m = argv[0];
  rb_a = argv[1];

  if (!NA_IsNArray(rb_a))
    rb_raise(rb_eArgError, "a (2th argument) must be NArray");
  if (NA_RANK(rb_a) != 2)
    rb_raise(rb_eArgError, "rank of a (2th argument) must be %d", 2);
  n = NA_SHAPE1(rb_a);
  lda = NA_SHAPE0(rb_a);
  if (NA_TYPE(rb_a) != NA_DCOMPLEX)
    rb_a = na_change_type(rb_a, NA_DCOMPLEX);
  a = NA_PTR_TYPE(rb_a, doublecomplex*);
  m = NUM2INT(rb_m);
  ldsa = MAX(1,m);
  {
    int shape[2];
    shape[0] = ldsa;
    shape[1] = n;
    rb_sa = na_make_object(NA_SCOMPLEX, 2, shape, cNArray);
  }
  sa = NA_PTR_TYPE(rb_sa, complex*);

  zlag2c_(&m, &n, a, &lda, sa, &ldsa, &info);

  rb_info = INT2NUM(info);
  return rb_ary_new3(2, rb_sa, rb_info);
}

void
init_lapack_zlag2c(VALUE mLapack){
  rb_define_module_function(mLapack, "zlag2c", rb_zlag2c, -1);
}

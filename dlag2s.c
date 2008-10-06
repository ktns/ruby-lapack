#include "rb_lapack.h"

static VALUE
rb_dlag2s(int argc, VALUE *argv, VALUE self){
  VALUE rb_m;
  integer m; 
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
    printf("%s\n", "USAGE:\n  sa, info = NumRu::Lapack.dlag2s( m, a)\n    or\n  NumRu::Lapack.dlag2s  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE DLAG2S( M, N, A, LDA, SA, LDSA, INFO)\n\n*  Purpose\n*  =======\n*\n*  DLAG2S converts a DOUBLE PRECISION matrix, SA, to a SINGLE\n*  PRECISION matrix, A.\n*\n*  RMAX is the overflow for the SINGLE PRECISION arithmetic\n*  DLAG2S checks that all the entries of A are between -RMAX and\n*  RMAX. If not the convertion is aborted and a flag is raised.\n*\n*  This is a helper routine so there is no argument checking.\n*\n\n*  Arguments\n*  =========\n*\n*  M       (input) INTEGER\n*          The number of lines of the matrix A.  M >= 0.\n*\n*  N       (input) INTEGER\n*          The number of columns of the matrix A.  N >= 0.\n*\n*  A       (input) DOUBLE PRECISION array, dimension (LDA,N)\n*          On entry, the M-by-N coefficient matrix A.\n*\n*  LDA     (input) INTEGER\n*          The leading dimension of the array A.  LDA >= max(1,M).\n*\n*  SA      (output) REAL array, dimension (LDSA,N)\n*          On exit, if INFO=0, the M-by-N coefficient matrix SA.\n*\n*  LDSA    (input) INTEGER\n*          The leading dimension of the array SA.  LDSA >= max(1,M).\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit\n*          > 0:  if INFO = k, the (i,j) entry of the matrix A has\n*                overflowed when moving from DOUBLE PRECISION to SINGLE\n*                k is given by k = (i-1)*LDA+j\n*\n*  =========\n*\n*     .. Local Scalars ..\n      INTEGER I,J\n      DOUBLE PRECISION RMAX\n*     ..\n*     .. External Functions ..\n      REAL SLAMCH\n      EXTERNAL SLAMCH\n*     ..\n\n");
    return Qnil;
  }
  if (argc != 2)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 2)", argc);
  rb_m = argv[0];
  rb_a = argv[1];

  m = NUM2INT(rb_m);
  if (!NA_IsNArray(rb_a))
    rb_raise(rb_eArgError, "a (2th argument) must be NArray");
  if (NA_RANK(rb_a) != 2)
    rb_raise(rb_eArgError, "rank of a (2th argument) must be %d", 2);
  lda = NA_SHAPE0(rb_a);
  n = NA_SHAPE1(rb_a);
  if (NA_TYPE(rb_a) != NA_DFLOAT)
    rb_a = na_change_type(rb_a, NA_DFLOAT);
  a = NA_PTR_TYPE(rb_a, doublereal*);
  ldsa = MAX(1,m);
  {
    int shape[2];
    shape[0] = ldsa;
    shape[1] = n;
    rb_sa = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  sa = NA_PTR_TYPE(rb_sa, real*);

  dlag2s_(&m, &n, a, &lda, sa, &ldsa, &info);

  rb_info = INT2NUM(info);
  return rb_ary_new3(2, rb_sa, rb_info);
}

void
init_lapack_dlag2s(VALUE mLapack){
  rb_define_module_function(mLapack, "dlag2s", rb_dlag2s, -1);
}

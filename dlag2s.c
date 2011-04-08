#include "rb_lapack.h"

extern VOID dlag2s_(integer *m, integer *n, doublereal *a, integer *lda, real *sa, integer *ldsa, integer *info);

static VALUE sHelp, sUsage;

static VALUE
rblapack_dlag2s(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_m;
  integer m; 
  VALUE rblapack_a;
  doublereal *a; 
  VALUE rblapack_sa;
  real *sa; 
  VALUE rblapack_info;
  integer info; 

  integer lda;
  integer n;
  integer ldsa;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  sa, info = NumRu::Lapack.dlag2s( m, a, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE DLAG2S( M, N, A, LDA, SA, LDSA, INFO )\n\n*  Purpose\n*  =======\n*\n*  DLAG2S converts a DOUBLE PRECISION matrix, SA, to a SINGLE\n*  PRECISION matrix, A.\n*\n*  RMAX is the overflow for the SINGLE PRECISION arithmetic\n*  DLAG2S checks that all the entries of A are between -RMAX and\n*  RMAX. If not the convertion is aborted and a flag is raised.\n*\n*  This is an auxiliary routine so there is no argument checking.\n*\n\n*  Arguments\n*  =========\n*\n*  M       (input) INTEGER\n*          The number of lines of the matrix A.  M >= 0.\n*\n*  N       (input) INTEGER\n*          The number of columns of the matrix A.  N >= 0.\n*\n*  A       (input) DOUBLE PRECISION array, dimension (LDA,N)\n*          On entry, the M-by-N coefficient matrix A.\n*\n*  LDA     (input) INTEGER\n*          The leading dimension of the array A.  LDA >= max(1,M).\n*\n*  SA      (output) REAL array, dimension (LDSA,N)\n*          On exit, if INFO=0, the M-by-N coefficient matrix SA; if\n*          INFO>0, the content of SA is unspecified.\n*\n*  LDSA    (input) INTEGER\n*          The leading dimension of the array SA.  LDSA >= max(1,M).\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit.\n*          = 1:  an entry of the matrix A is greater than the SINGLE\n*                PRECISION overflow threshold, in this case, the content\n*                of SA in exit is unspecified.\n*\n*  =========\n*\n*     .. Local Scalars ..\n      INTEGER            I, J\n      DOUBLE PRECISION   RMAX\n*     ..\n*     .. External Functions ..\n      REAL               SLAMCH\n      EXTERNAL           SLAMCH\n*     ..\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  sa, info = NumRu::Lapack.dlag2s( m, a, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 2)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 2)", argc);
  rblapack_m = argv[0];
  rblapack_a = argv[1];
  if (rb_options != Qnil) {
  }

  if (!NA_IsNArray(rblapack_a))
    rb_raise(rb_eArgError, "a (2th argument) must be NArray");
  if (NA_RANK(rblapack_a) != 2)
    rb_raise(rb_eArgError, "rank of a (2th argument) must be %d", 2);
  n = NA_SHAPE1(rblapack_a);
  lda = NA_SHAPE0(rblapack_a);
  if (NA_TYPE(rblapack_a) != NA_DFLOAT)
    rblapack_a = na_change_type(rblapack_a, NA_DFLOAT);
  a = NA_PTR_TYPE(rblapack_a, doublereal*);
  m = NUM2INT(rblapack_m);
  ldsa = MAX(1,m);
  {
    int shape[2];
    shape[0] = ldsa;
    shape[1] = n;
    rblapack_sa = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  sa = NA_PTR_TYPE(rblapack_sa, real*);

  dlag2s_(&m, &n, a, &lda, sa, &ldsa, &info);

  rblapack_info = INT2NUM(info);
  return rb_ary_new3(2, rblapack_sa, rblapack_info);
}

void
init_lapack_dlag2s(VALUE mLapack){
  rb_define_module_function(mLapack, "dlag2s", rblapack_dlag2s, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}

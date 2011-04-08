#include "rb_lapack.h"

extern VOID slag2d_(integer *m, integer *n, real *sa, integer *ldsa, doublereal *a, integer *lda, integer *info);

static VALUE sHelp, sUsage;

static VALUE
rblapack_slag2d(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_m;
  integer m; 
  VALUE rblapack_sa;
  real *sa; 
  VALUE rblapack_a;
  doublereal *a; 
  VALUE rblapack_info;
  integer info; 

  integer ldsa;
  integer n;
  integer lda;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  a, info = NumRu::Lapack.slag2d( m, sa, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE SLAG2D( M, N, SA, LDSA, A, LDA, INFO )\n\n*  Purpose\n*  =======\n*\n*  SLAG2D converts a SINGLE PRECISION matrix, SA, to a DOUBLE\n*  PRECISION matrix, A.\n*\n*  Note that while it is possible to overflow while converting\n*  from double to single, it is not possible to overflow when\n*  converting from single to double.\n*\n*  This is an auxiliary routine so there is no argument checking.\n*\n\n*  Arguments\n*  =========\n*\n*  M       (input) INTEGER\n*          The number of lines of the matrix A.  M >= 0.\n*\n*  N       (input) INTEGER\n*          The number of columns of the matrix A.  N >= 0.\n*\n*  SA      (input) REAL array, dimension (LDSA,N)\n*          On entry, the M-by-N coefficient matrix SA.\n*\n*  LDSA    (input) INTEGER\n*          The leading dimension of the array SA.  LDSA >= max(1,M).\n*\n*  A       (output) DOUBLE PRECISION array, dimension (LDA,N)\n*          On exit, the M-by-N coefficient matrix A.\n*\n*  LDA     (input) INTEGER\n*          The leading dimension of the array A.  LDA >= max(1,M).\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit\n*  =========\n*\n*     .. Local Scalars ..\n      INTEGER            I, J\n*     ..\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  a, info = NumRu::Lapack.slag2d( m, sa, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 2)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 2)", argc);
  rblapack_m = argv[0];
  rblapack_sa = argv[1];
  if (rb_options != Qnil) {
  }

  m = NUM2INT(rblapack_m);
  if (!NA_IsNArray(rblapack_sa))
    rb_raise(rb_eArgError, "sa (2th argument) must be NArray");
  if (NA_RANK(rblapack_sa) != 2)
    rb_raise(rb_eArgError, "rank of sa (2th argument) must be %d", 2);
  n = NA_SHAPE1(rblapack_sa);
  ldsa = NA_SHAPE0(rblapack_sa);
  if (NA_TYPE(rblapack_sa) != NA_SFLOAT)
    rblapack_sa = na_change_type(rblapack_sa, NA_SFLOAT);
  sa = NA_PTR_TYPE(rblapack_sa, real*);
  lda = MAX(1,m);
  {
    int shape[2];
    shape[0] = lda;
    shape[1] = n;
    rblapack_a = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  a = NA_PTR_TYPE(rblapack_a, doublereal*);

  slag2d_(&m, &n, sa, &ldsa, a, &lda, &info);

  rblapack_info = INT2NUM(info);
  return rb_ary_new3(2, rblapack_a, rblapack_info);
}

void
init_lapack_slag2d(VALUE mLapack){
  rb_define_module_function(mLapack, "slag2d", rblapack_slag2d, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}

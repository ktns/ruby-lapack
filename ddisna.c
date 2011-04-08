#include "rb_lapack.h"

extern VOID ddisna_(char *job, integer *m, integer *n, doublereal *d, doublereal *sep, integer *info);

static VALUE sHelp, sUsage;

static VALUE
rblapack_ddisna(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_job;
  char job; 
  VALUE rblapack_n;
  integer n; 
  VALUE rblapack_d;
  doublereal *d; 
  VALUE rblapack_sep;
  doublereal *sep; 
  VALUE rblapack_info;
  integer info; 

  integer m;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  sep, info = NumRu::Lapack.ddisna( job, n, d, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE DDISNA( JOB, M, N, D, SEP, INFO )\n\n*  Purpose\n*  =======\n*\n*  DDISNA computes the reciprocal condition numbers for the eigenvectors\n*  of a real symmetric or complex Hermitian matrix or for the left or\n*  right singular vectors of a general m-by-n matrix. The reciprocal\n*  condition number is the 'gap' between the corresponding eigenvalue or\n*  singular value and the nearest other one.\n*\n*  The bound on the error, measured by angle in radians, in the I-th\n*  computed vector is given by\n*\n*         DLAMCH( 'E' ) * ( ANORM / SEP( I ) )\n*\n*  where ANORM = 2-norm(A) = max( abs( D(j) ) ).  SEP(I) is not allowed\n*  to be smaller than DLAMCH( 'E' )*ANORM in order to limit the size of\n*  the error bound.\n*\n*  DDISNA may also be used to compute error bounds for eigenvectors of\n*  the generalized symmetric definite eigenproblem.\n*\n\n*  Arguments\n*  =========\n*\n*  JOB     (input) CHARACTER*1\n*          Specifies for which problem the reciprocal condition numbers\n*          should be computed:\n*          = 'E':  the eigenvectors of a symmetric/Hermitian matrix;\n*          = 'L':  the left singular vectors of a general matrix;\n*          = 'R':  the right singular vectors of a general matrix.\n*\n*  M       (input) INTEGER\n*          The number of rows of the matrix. M >= 0.\n*\n*  N       (input) INTEGER\n*          If JOB = 'L' or 'R', the number of columns of the matrix,\n*          in which case N >= 0. Ignored if JOB = 'E'.\n*\n*  D       (input) DOUBLE PRECISION array, dimension (M) if JOB = 'E'\n*                              dimension (min(M,N)) if JOB = 'L' or 'R'\n*          The eigenvalues (if JOB = 'E') or singular values (if JOB =\n*          'L' or 'R') of the matrix, in either increasing or decreasing\n*          order. If singular values, they must be non-negative.\n*\n*  SEP     (output) DOUBLE PRECISION array, dimension (M) if JOB = 'E'\n*                               dimension (min(M,N)) if JOB = 'L' or 'R'\n*          The reciprocal condition numbers of the vectors.\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit.\n*          < 0:  if INFO = -i, the i-th argument had an illegal value.\n*\n\n*  =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  sep, info = NumRu::Lapack.ddisna( job, n, d, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 3)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 3)", argc);
  rblapack_job = argv[0];
  rblapack_n = argv[1];
  rblapack_d = argv[2];
  if (rb_options != Qnil) {
  }

  n = NUM2INT(rblapack_n);
  if (!NA_IsNArray(rblapack_d))
    rb_raise(rb_eArgError, "d (3th argument) must be NArray");
  if (NA_RANK(rblapack_d) != 1)
    rb_raise(rb_eArgError, "rank of d (3th argument) must be %d", 1);
  m = NA_SHAPE0(rblapack_d);
  if (NA_TYPE(rblapack_d) != NA_DFLOAT)
    rblapack_d = na_change_type(rblapack_d, NA_DFLOAT);
  d = NA_PTR_TYPE(rblapack_d, doublereal*);
  job = StringValueCStr(rblapack_job)[0];
  {
    int shape[1];
    shape[0] = lsame_(&job,"E") ? m : ((lsame_(&job,"L")) || (lsame_(&job,"R"))) ? MIN(m,n) : 0;
    rblapack_sep = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  sep = NA_PTR_TYPE(rblapack_sep, doublereal*);

  ddisna_(&job, &m, &n, d, sep, &info);

  rblapack_info = INT2NUM(info);
  return rb_ary_new3(2, rblapack_sep, rblapack_info);
}

void
init_lapack_ddisna(VALUE mLapack){
  rb_define_module_function(mLapack, "ddisna", rblapack_ddisna, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}

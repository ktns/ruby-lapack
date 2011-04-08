#include "rb_lapack.h"

extern VOID claqge_(integer *m, integer *n, complex *a, integer *lda, real *r, real *c, real *rowcnd, real *colcnd, real *amax, char *equed);

static VALUE sHelp, sUsage;

static VALUE
rblapack_claqge(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_a;
  complex *a; 
  VALUE rblapack_r;
  real *r; 
  VALUE rblapack_c;
  real *c; 
  VALUE rblapack_rowcnd;
  real rowcnd; 
  VALUE rblapack_colcnd;
  real colcnd; 
  VALUE rblapack_amax;
  real amax; 
  VALUE rblapack_equed;
  char equed; 
  VALUE rblapack_a_out__;
  complex *a_out__;

  integer lda;
  integer n;
  integer m;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  equed, a = NumRu::Lapack.claqge( a, r, c, rowcnd, colcnd, amax, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE CLAQGE( M, N, A, LDA, R, C, ROWCND, COLCND, AMAX, EQUED )\n\n*  Purpose\n*  =======\n*\n*  CLAQGE equilibrates a general M by N matrix A using the row and\n*  column scaling factors in the vectors R and C.\n*\n\n*  Arguments\n*  =========\n*\n*  M       (input) INTEGER\n*          The number of rows of the matrix A.  M >= 0.\n*\n*  N       (input) INTEGER\n*          The number of columns of the matrix A.  N >= 0.\n*\n*  A       (input/output) COMPLEX array, dimension (LDA,N)\n*          On entry, the M by N matrix A.\n*          On exit, the equilibrated matrix.  See EQUED for the form of\n*          the equilibrated matrix.\n*\n*  LDA     (input) INTEGER\n*          The leading dimension of the array A.  LDA >= max(M,1).\n*\n*  R       (input) REAL array, dimension (M)\n*          The row scale factors for A.\n*\n*  C       (input) REAL array, dimension (N)\n*          The column scale factors for A.\n*\n*  ROWCND  (input) REAL\n*          Ratio of the smallest R(i) to the largest R(i).\n*\n*  COLCND  (input) REAL\n*          Ratio of the smallest C(i) to the largest C(i).\n*\n*  AMAX    (input) REAL\n*          Absolute value of largest matrix entry.\n*\n*  EQUED   (output) CHARACTER*1\n*          Specifies the form of equilibration that was done.\n*          = 'N':  No equilibration\n*          = 'R':  Row equilibration, i.e., A has been premultiplied by\n*                  diag(R).\n*          = 'C':  Column equilibration, i.e., A has been postmultiplied\n*                  by diag(C).\n*          = 'B':  Both row and column equilibration, i.e., A has been\n*                  replaced by diag(R) * A * diag(C).\n*\n*  Internal Parameters\n*  ===================\n*\n*  THRESH is a threshold value used to decide if row or column scaling\n*  should be done based on the ratio of the row or column scaling\n*  factors.  If ROWCND < THRESH, row scaling is done, and if\n*  COLCND < THRESH, column scaling is done.\n*\n*  LARGE and SMALL are threshold values used to decide if row scaling\n*  should be done based on the absolute size of the largest matrix\n*  element.  If AMAX > LARGE or AMAX < SMALL, row scaling is done.\n*\n\n*  =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  equed, a = NumRu::Lapack.claqge( a, r, c, rowcnd, colcnd, amax, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 6)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 6)", argc);
  rblapack_a = argv[0];
  rblapack_r = argv[1];
  rblapack_c = argv[2];
  rblapack_rowcnd = argv[3];
  rblapack_colcnd = argv[4];
  rblapack_amax = argv[5];
  if (rb_options != Qnil) {
  }

  if (!NA_IsNArray(rblapack_a))
    rb_raise(rb_eArgError, "a (1th argument) must be NArray");
  if (NA_RANK(rblapack_a) != 2)
    rb_raise(rb_eArgError, "rank of a (1th argument) must be %d", 2);
  n = NA_SHAPE1(rblapack_a);
  lda = NA_SHAPE0(rblapack_a);
  if (NA_TYPE(rblapack_a) != NA_SCOMPLEX)
    rblapack_a = na_change_type(rblapack_a, NA_SCOMPLEX);
  a = NA_PTR_TYPE(rblapack_a, complex*);
  if (!NA_IsNArray(rblapack_c))
    rb_raise(rb_eArgError, "c (3th argument) must be NArray");
  if (NA_RANK(rblapack_c) != 1)
    rb_raise(rb_eArgError, "rank of c (3th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_c) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of c must be the same as shape 1 of a");
  if (NA_TYPE(rblapack_c) != NA_SFLOAT)
    rblapack_c = na_change_type(rblapack_c, NA_SFLOAT);
  c = NA_PTR_TYPE(rblapack_c, real*);
  amax = (real)NUM2DBL(rblapack_amax);
  colcnd = (real)NUM2DBL(rblapack_colcnd);
  if (!NA_IsNArray(rblapack_r))
    rb_raise(rb_eArgError, "r (2th argument) must be NArray");
  if (NA_RANK(rblapack_r) != 1)
    rb_raise(rb_eArgError, "rank of r (2th argument) must be %d", 1);
  m = NA_SHAPE0(rblapack_r);
  if (NA_TYPE(rblapack_r) != NA_SFLOAT)
    rblapack_r = na_change_type(rblapack_r, NA_SFLOAT);
  r = NA_PTR_TYPE(rblapack_r, real*);
  rowcnd = (real)NUM2DBL(rblapack_rowcnd);
  {
    int shape[2];
    shape[0] = lda;
    shape[1] = n;
    rblapack_a_out__ = na_make_object(NA_SCOMPLEX, 2, shape, cNArray);
  }
  a_out__ = NA_PTR_TYPE(rblapack_a_out__, complex*);
  MEMCPY(a_out__, a, complex, NA_TOTAL(rblapack_a));
  rblapack_a = rblapack_a_out__;
  a = a_out__;

  claqge_(&m, &n, a, &lda, r, c, &rowcnd, &colcnd, &amax, &equed);

  rblapack_equed = rb_str_new(&equed,1);
  return rb_ary_new3(2, rblapack_equed, rblapack_a);
}

void
init_lapack_claqge(VALUE mLapack){
  rb_define_module_function(mLapack, "claqge", rblapack_claqge, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}

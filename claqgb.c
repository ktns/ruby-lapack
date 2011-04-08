#include "rb_lapack.h"

extern VOID claqgb_(integer *m, integer *n, integer *kl, integer *ku, complex *ab, integer *ldab, real *r, real *c, real *rowcnd, real *colcnd, real *amax, char *equed);

static VALUE sHelp, sUsage;

static VALUE
rblapack_claqgb(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_kl;
  integer kl; 
  VALUE rblapack_ku;
  integer ku; 
  VALUE rblapack_ab;
  complex *ab; 
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
  VALUE rblapack_ab_out__;
  complex *ab_out__;

  integer ldab;
  integer n;
  integer m;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  equed, ab = NumRu::Lapack.claqgb( kl, ku, ab, r, c, rowcnd, colcnd, amax, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE CLAQGB( M, N, KL, KU, AB, LDAB, R, C, ROWCND, COLCND, AMAX, EQUED )\n\n*  Purpose\n*  =======\n*\n*  CLAQGB equilibrates a general M by N band matrix A with KL\n*  subdiagonals and KU superdiagonals using the row and scaling factors\n*  in the vectors R and C.\n*\n\n*  Arguments\n*  =========\n*\n*  M       (input) INTEGER\n*          The number of rows of the matrix A.  M >= 0.\n*\n*  N       (input) INTEGER\n*          The number of columns of the matrix A.  N >= 0.\n*\n*  KL      (input) INTEGER\n*          The number of subdiagonals within the band of A.  KL >= 0.\n*\n*  KU      (input) INTEGER\n*          The number of superdiagonals within the band of A.  KU >= 0.\n*\n*  AB      (input/output) COMPLEX array, dimension (LDAB,N)\n*          On entry, the matrix A in band storage, in rows 1 to KL+KU+1.\n*          The j-th column of A is stored in the j-th column of the\n*          array AB as follows:\n*          AB(ku+1+i-j,j) = A(i,j) for max(1,j-ku)<=i<=min(m,j+kl)\n*\n*          On exit, the equilibrated matrix, in the same storage format\n*          as A.  See EQUED for the form of the equilibrated matrix.\n*\n*  LDAB    (input) INTEGER\n*          The leading dimension of the array AB.  LDA >= KL+KU+1.\n*\n*  R       (input) REAL array, dimension (M)\n*          The row scale factors for A.\n*\n*  C       (input) REAL array, dimension (N)\n*          The column scale factors for A.\n*\n*  ROWCND  (input) REAL\n*          Ratio of the smallest R(i) to the largest R(i).\n*\n*  COLCND  (input) REAL\n*          Ratio of the smallest C(i) to the largest C(i).\n*\n*  AMAX    (input) REAL\n*          Absolute value of largest matrix entry.\n*\n*  EQUED   (output) CHARACTER*1\n*          Specifies the form of equilibration that was done.\n*          = 'N':  No equilibration\n*          = 'R':  Row equilibration, i.e., A has been premultiplied by\n*                  diag(R).\n*          = 'C':  Column equilibration, i.e., A has been postmultiplied\n*                  by diag(C).\n*          = 'B':  Both row and column equilibration, i.e., A has been\n*                  replaced by diag(R) * A * diag(C).\n*\n*  Internal Parameters\n*  ===================\n*\n*  THRESH is a threshold value used to decide if row or column scaling\n*  should be done based on the ratio of the row or column scaling\n*  factors.  If ROWCND < THRESH, row scaling is done, and if\n*  COLCND < THRESH, column scaling is done.\n*\n*  LARGE and SMALL are threshold values used to decide if row scaling\n*  should be done based on the absolute size of the largest matrix\n*  element.  If AMAX > LARGE or AMAX < SMALL, row scaling is done.\n*\n\n*  =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  equed, ab = NumRu::Lapack.claqgb( kl, ku, ab, r, c, rowcnd, colcnd, amax, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 8)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 8)", argc);
  rblapack_kl = argv[0];
  rblapack_ku = argv[1];
  rblapack_ab = argv[2];
  rblapack_r = argv[3];
  rblapack_c = argv[4];
  rblapack_rowcnd = argv[5];
  rblapack_colcnd = argv[6];
  rblapack_amax = argv[7];
  if (rb_options != Qnil) {
  }

  if (!NA_IsNArray(rblapack_ab))
    rb_raise(rb_eArgError, "ab (3th argument) must be NArray");
  if (NA_RANK(rblapack_ab) != 2)
    rb_raise(rb_eArgError, "rank of ab (3th argument) must be %d", 2);
  n = NA_SHAPE1(rblapack_ab);
  ldab = NA_SHAPE0(rblapack_ab);
  if (NA_TYPE(rblapack_ab) != NA_SCOMPLEX)
    rblapack_ab = na_change_type(rblapack_ab, NA_SCOMPLEX);
  ab = NA_PTR_TYPE(rblapack_ab, complex*);
  kl = NUM2INT(rblapack_kl);
  if (!NA_IsNArray(rblapack_c))
    rb_raise(rb_eArgError, "c (5th argument) must be NArray");
  if (NA_RANK(rblapack_c) != 1)
    rb_raise(rb_eArgError, "rank of c (5th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_c) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of c must be the same as shape 1 of ab");
  if (NA_TYPE(rblapack_c) != NA_SFLOAT)
    rblapack_c = na_change_type(rblapack_c, NA_SFLOAT);
  c = NA_PTR_TYPE(rblapack_c, real*);
  amax = (real)NUM2DBL(rblapack_amax);
  colcnd = (real)NUM2DBL(rblapack_colcnd);
  if (!NA_IsNArray(rblapack_r))
    rb_raise(rb_eArgError, "r (4th argument) must be NArray");
  if (NA_RANK(rblapack_r) != 1)
    rb_raise(rb_eArgError, "rank of r (4th argument) must be %d", 1);
  m = NA_SHAPE0(rblapack_r);
  if (NA_TYPE(rblapack_r) != NA_SFLOAT)
    rblapack_r = na_change_type(rblapack_r, NA_SFLOAT);
  r = NA_PTR_TYPE(rblapack_r, real*);
  rowcnd = (real)NUM2DBL(rblapack_rowcnd);
  ku = NUM2INT(rblapack_ku);
  {
    int shape[2];
    shape[0] = ldab;
    shape[1] = n;
    rblapack_ab_out__ = na_make_object(NA_SCOMPLEX, 2, shape, cNArray);
  }
  ab_out__ = NA_PTR_TYPE(rblapack_ab_out__, complex*);
  MEMCPY(ab_out__, ab, complex, NA_TOTAL(rblapack_ab));
  rblapack_ab = rblapack_ab_out__;
  ab = ab_out__;

  claqgb_(&m, &n, &kl, &ku, ab, &ldab, r, c, &rowcnd, &colcnd, &amax, &equed);

  rblapack_equed = rb_str_new(&equed,1);
  return rb_ary_new3(2, rblapack_equed, rblapack_ab);
}

void
init_lapack_claqgb(VALUE mLapack){
  rb_define_module_function(mLapack, "claqgb", rblapack_claqgb, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}

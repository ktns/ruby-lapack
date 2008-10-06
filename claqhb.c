#include "rb_lapack.h"

static VALUE
rb_claqhb(int argc, VALUE *argv, VALUE self){
  VALUE rb_uplo;
  char uplo; 
  VALUE rb_kd;
  integer kd; 
  VALUE rb_ab;
  complex *ab; 
  VALUE rb_scond;
  real scond; 
  VALUE rb_amax;
  real amax; 
  VALUE rb_s;
  real *s; 
  VALUE rb_equed;
  char equed; 
  VALUE rb_ab_out__;
  complex *ab_out__;

  integer ldab;
  integer n;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  s, equed, ab = NumRu::Lapack.claqhb( uplo, kd, ab, scond, amax)\n    or\n  NumRu::Lapack.claqhb  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE CLAQHB( UPLO, N, KD, AB, LDAB, S, SCOND, AMAX, EQUED )\n\n*  Purpose\n*  =======\n*\n*  CLAQHB equilibrates an Hermitian band matrix A using the scaling\n*  factors in the vector S.\n*\n\n*  Arguments\n*  =========\n*\n*  UPLO    (input) CHARACTER*1\n*          Specifies whether the upper or lower triangular part of the\n*          symmetric matrix A is stored.\n*          = 'U':  Upper triangular\n*          = 'L':  Lower triangular\n*\n*  N       (input) INTEGER\n*          The order of the matrix A.  N >= 0.\n*\n*  KD      (input) INTEGER\n*          The number of super-diagonals of the matrix A if UPLO = 'U',\n*          or the number of sub-diagonals if UPLO = 'L'.  KD >= 0.\n*\n*  AB      (input/output) COMPLEX array, dimension (LDAB,N)\n*          On entry, the upper or lower triangle of the symmetric band\n*          matrix A, stored in the first KD+1 rows of the array.  The\n*          j-th column of A is stored in the j-th column of the array AB\n*          as follows:\n*          if UPLO = 'U', AB(kd+1+i-j,j) = A(i,j) for max(1,j-kd)<=i<=j;\n*          if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=min(n,j+kd).\n*\n*          On exit, if INFO = 0, the triangular factor U or L from the\n*          Cholesky factorization A = U'*U or A = L*L' of the band\n*          matrix A, in the same storage format as A.\n*\n*  LDAB    (input) INTEGER\n*          The leading dimension of the array AB.  LDAB >= KD+1.\n*\n*  S       (output) REAL array, dimension (N)\n*          The scale factors for A.\n*\n*  SCOND   (input) REAL\n*          Ratio of the smallest S(i) to the largest S(i).\n*\n*  AMAX    (input) REAL\n*          Absolute value of largest matrix entry.\n*\n*  EQUED   (output) CHARACTER*1\n*          Specifies whether or not equilibration was done.\n*          = 'N':  No equilibration.\n*          = 'Y':  Equilibration was done, i.e., A has been replaced by\n*                  diag(S) * A * diag(S).\n*\n*  Internal Parameters\n*  ===================\n*\n*  THRESH is a threshold value used to decide if scaling should be done\n*  based on the ratio of the scaling factors.  If SCOND < THRESH,\n*  scaling is done.\n*\n*  LARGE and SMALL are threshold values used to decide if scaling should\n*  be done based on the absolute size of the largest matrix element.\n*  If AMAX > LARGE or AMAX < SMALL, scaling is done.\n*\n\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 5)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 5)", argc);
  rb_uplo = argv[0];
  rb_kd = argv[1];
  rb_ab = argv[2];
  rb_scond = argv[3];
  rb_amax = argv[4];

  uplo = StringValueCStr(rb_uplo)[0];
  kd = NUM2INT(rb_kd);
  scond = (real)NUM2DBL(rb_scond);
  amax = (real)NUM2DBL(rb_amax);
  if (!NA_IsNArray(rb_ab))
    rb_raise(rb_eArgError, "ab (3th argument) must be NArray");
  if (NA_RANK(rb_ab) != 2)
    rb_raise(rb_eArgError, "rank of ab (3th argument) must be %d", 2);
  ldab = NA_SHAPE0(rb_ab);
  n = NA_SHAPE1(rb_ab);
  if (NA_TYPE(rb_ab) != NA_SCOMPLEX)
    rb_ab = na_change_type(rb_ab, NA_SCOMPLEX);
  ab = NA_PTR_TYPE(rb_ab, complex*);
  {
    int shape[1];
    shape[0] = n;
    rb_s = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  s = NA_PTR_TYPE(rb_s, real*);
  {
    int shape[2];
    shape[0] = ldab;
    shape[1] = n;
    rb_ab_out__ = na_make_object(NA_SCOMPLEX, 2, shape, cNArray);
  }
  ab_out__ = NA_PTR_TYPE(rb_ab_out__, complex*);
  MEMCPY(ab_out__, ab, complex, NA_TOTAL(rb_ab));
  rb_ab = rb_ab_out__;
  ab = ab_out__;

  claqhb_(&uplo, &n, &kd, ab, &ldab, s, &scond, &amax, &equed);

  rb_equed = rb_str_new(&equed,1);
  return rb_ary_new3(3, rb_s, rb_equed, rb_ab);
}

void
init_lapack_claqhb(VALUE mLapack){
  rb_define_module_function(mLapack, "claqhb", rb_claqhb, -1);
}

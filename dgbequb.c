#include "rb_lapack.h"

extern VOID dgbequb_(integer *m, integer *n, integer *kl, integer *ku, doublereal *ab, integer *ldab, doublereal *r, doublereal *c, doublereal *rowcnd, doublereal *colcnd, doublereal *amax, integer *info);

static VALUE sHelp, sUsage;

static VALUE
rblapack_dgbequb(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_kl;
  integer kl; 
  VALUE rblapack_ku;
  integer ku; 
  VALUE rblapack_ab;
  doublereal *ab; 
  VALUE rblapack_r;
  doublereal *r; 
  VALUE rblapack_c;
  doublereal *c; 
  VALUE rblapack_rowcnd;
  doublereal rowcnd; 
  VALUE rblapack_colcnd;
  doublereal colcnd; 
  VALUE rblapack_amax;
  doublereal amax; 
  VALUE rblapack_info;
  integer info; 

  integer ldab;
  integer n;
  integer m;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  r, c, rowcnd, colcnd, amax, info = NumRu::Lapack.dgbequb( kl, ku, ab, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE DGBEQUB( M, N, KL, KU, AB, LDAB, R, C, ROWCND, COLCND, AMAX, INFO )\n\n*  Purpose\n*  =======\n*\n*  DGBEQUB computes row and column scalings intended to equilibrate an\n*  M-by-N matrix A and reduce its condition number.  R returns the row\n*  scale factors and C the column scale factors, chosen to try to make\n*  the largest element in each row and column of the matrix B with\n*  elements B(i,j)=R(i)*A(i,j)*C(j) have an absolute value of at most\n*  the radix.\n*\n*  R(i) and C(j) are restricted to be a power of the radix between\n*  SMLNUM = smallest safe number and BIGNUM = largest safe number.  Use\n*  of these scaling factors is not guaranteed to reduce the condition\n*  number of A but works well in practice.\n*\n*  This routine differs from DGEEQU by restricting the scaling factors\n*  to a power of the radix.  Baring over- and underflow, scaling by\n*  these factors introduces no additional rounding errors.  However, the\n*  scaled entries' magnitured are no longer approximately 1 but lie\n*  between sqrt(radix) and 1/sqrt(radix).\n*\n\n*  Arguments\n*  =========\n*\n*  M       (input) INTEGER\n*          The number of rows of the matrix A.  M >= 0.\n*\n*  N       (input) INTEGER\n*          The number of columns of the matrix A.  N >= 0.\n*\n*  KL      (input) INTEGER\n*          The number of subdiagonals within the band of A.  KL >= 0.\n*\n*  KU      (input) INTEGER\n*          The number of superdiagonals within the band of A.  KU >= 0.\n*\n*  AB      (input) DOUBLE PRECISION array, dimension (LDAB,N)\n*          On entry, the matrix A in band storage, in rows 1 to KL+KU+1.\n*          The j-th column of A is stored in the j-th column of the\n*          array AB as follows:\n*          AB(KU+1+i-j,j) = A(i,j) for max(1,j-KU)<=i<=min(N,j+kl)\n*\n*  LDAB    (input) INTEGER\n*          The leading dimension of the array A.  LDAB >= max(1,M).\n*\n*  R       (output) DOUBLE PRECISION array, dimension (M)\n*          If INFO = 0 or INFO > M, R contains the row scale factors\n*          for A.\n*\n*  C       (output) DOUBLE PRECISION array, dimension (N)\n*          If INFO = 0,  C contains the column scale factors for A.\n*\n*  ROWCND  (output) DOUBLE PRECISION\n*          If INFO = 0 or INFO > M, ROWCND contains the ratio of the\n*          smallest R(i) to the largest R(i).  If ROWCND >= 0.1 and\n*          AMAX is neither too large nor too small, it is not worth\n*          scaling by R.\n*\n*  COLCND  (output) DOUBLE PRECISION\n*          If INFO = 0, COLCND contains the ratio of the smallest\n*          C(i) to the largest C(i).  If COLCND >= 0.1, it is not\n*          worth scaling by C.\n*\n*  AMAX    (output) DOUBLE PRECISION\n*          Absolute value of largest matrix element.  If AMAX is very\n*          close to overflow or very close to underflow, the matrix\n*          should be scaled.\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit\n*          < 0:  if INFO = -i, the i-th argument had an illegal value\n*          > 0:  if INFO = i,  and i is\n*                <= M:  the i-th row of A is exactly zero\n*                >  M:  the (i-M)-th column of A is exactly zero\n*\n\n*  =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  r, c, rowcnd, colcnd, amax, info = NumRu::Lapack.dgbequb( kl, ku, ab, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 3)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 3)", argc);
  rblapack_kl = argv[0];
  rblapack_ku = argv[1];
  rblapack_ab = argv[2];
  if (rb_options != Qnil) {
  }

  if (!NA_IsNArray(rblapack_ab))
    rb_raise(rb_eArgError, "ab (3th argument) must be NArray");
  if (NA_RANK(rblapack_ab) != 2)
    rb_raise(rb_eArgError, "rank of ab (3th argument) must be %d", 2);
  n = NA_SHAPE1(rblapack_ab);
  ldab = NA_SHAPE0(rblapack_ab);
  if (ldab != (m))
    rb_raise(rb_eRuntimeError, "shape 0 of ab must be %d", m);
  m = ldab;
  if (NA_TYPE(rblapack_ab) != NA_DFLOAT)
    rblapack_ab = na_change_type(rblapack_ab, NA_DFLOAT);
  ab = NA_PTR_TYPE(rblapack_ab, doublereal*);
  kl = NUM2INT(rblapack_kl);
  ku = NUM2INT(rblapack_ku);
  ldab = m;
  {
    int shape[1];
    shape[0] = m;
    rblapack_r = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  r = NA_PTR_TYPE(rblapack_r, doublereal*);
  {
    int shape[1];
    shape[0] = n;
    rblapack_c = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  c = NA_PTR_TYPE(rblapack_c, doublereal*);

  dgbequb_(&m, &n, &kl, &ku, ab, &ldab, r, c, &rowcnd, &colcnd, &amax, &info);

  rblapack_rowcnd = rb_float_new((double)rowcnd);
  rblapack_colcnd = rb_float_new((double)colcnd);
  rblapack_amax = rb_float_new((double)amax);
  rblapack_info = INT2NUM(info);
  return rb_ary_new3(6, rblapack_r, rblapack_c, rblapack_rowcnd, rblapack_colcnd, rblapack_amax, rblapack_info);
}

void
init_lapack_dgbequb(VALUE mLapack){
  rb_define_module_function(mLapack, "dgbequb", rblapack_dgbequb, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}

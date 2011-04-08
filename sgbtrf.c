#include "rb_lapack.h"

extern VOID sgbtrf_(integer *m, integer *n, integer *kl, integer *ku, real *ab, integer *ldab, integer *ipiv, integer *info);

static VALUE sHelp, sUsage;

static VALUE
rblapack_sgbtrf(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_m;
  integer m; 
  VALUE rblapack_kl;
  integer kl; 
  VALUE rblapack_ku;
  integer ku; 
  VALUE rblapack_ab;
  real *ab; 
  VALUE rblapack_ipiv;
  integer *ipiv; 
  VALUE rblapack_info;
  integer info; 
  VALUE rblapack_ab_out__;
  real *ab_out__;

  integer ldab;
  integer n;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  ipiv, info, ab = NumRu::Lapack.sgbtrf( m, kl, ku, ab, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE SGBTRF( M, N, KL, KU, AB, LDAB, IPIV, INFO )\n\n*  Purpose\n*  =======\n*\n*  SGBTRF computes an LU factorization of a real m-by-n band matrix A\n*  using partial pivoting with row interchanges.\n*\n*  This is the blocked version of the algorithm, calling Level 3 BLAS.\n*\n\n*  Arguments\n*  =========\n*\n*  M       (input) INTEGER\n*          The number of rows of the matrix A.  M >= 0.\n*\n*  N       (input) INTEGER\n*          The number of columns of the matrix A.  N >= 0.\n*\n*  KL      (input) INTEGER\n*          The number of subdiagonals within the band of A.  KL >= 0.\n*\n*  KU      (input) INTEGER\n*          The number of superdiagonals within the band of A.  KU >= 0.\n*\n*  AB      (input/output) REAL array, dimension (LDAB,N)\n*          On entry, the matrix A in band storage, in rows KL+1 to\n*          2*KL+KU+1; rows 1 to KL of the array need not be set.\n*          The j-th column of A is stored in the j-th column of the\n*          array AB as follows:\n*          AB(kl+ku+1+i-j,j) = A(i,j) for max(1,j-ku)<=i<=min(m,j+kl)\n*\n*          On exit, details of the factorization: U is stored as an\n*          upper triangular band matrix with KL+KU superdiagonals in\n*          rows 1 to KL+KU+1, and the multipliers used during the\n*          factorization are stored in rows KL+KU+2 to 2*KL+KU+1.\n*          See below for further details.\n*\n*  LDAB    (input) INTEGER\n*          The leading dimension of the array AB.  LDAB >= 2*KL+KU+1.\n*\n*  IPIV    (output) INTEGER array, dimension (min(M,N))\n*          The pivot indices; for 1 <= i <= min(M,N), row i of the\n*          matrix was interchanged with row IPIV(i).\n*\n*  INFO    (output) INTEGER\n*          = 0: successful exit\n*          < 0: if INFO = -i, the i-th argument had an illegal value\n*          > 0: if INFO = +i, U(i,i) is exactly zero. The factorization\n*               has been completed, but the factor U is exactly\n*               singular, and division by zero will occur if it is used\n*               to solve a system of equations.\n*\n\n*  Further Details\n*  ===============\n*\n*  The band storage scheme is illustrated by the following example, when\n*  M = N = 6, KL = 2, KU = 1:\n*\n*  On entry:                       On exit:\n*\n*      *    *    *    +    +    +       *    *    *   u14  u25  u36\n*      *    *    +    +    +    +       *    *   u13  u24  u35  u46\n*      *   a12  a23  a34  a45  a56      *   u12  u23  u34  u45  u56\n*     a11  a22  a33  a44  a55  a66     u11  u22  u33  u44  u55  u66\n*     a21  a32  a43  a54  a65   *      m21  m32  m43  m54  m65   *\n*     a31  a42  a53  a64   *    *      m31  m42  m53  m64   *    *\n*\n*  Array elements marked * are not used by the routine; elements marked\n*  + need not be set on entry, but are required by the routine to store\n*  elements of U because of fill-in resulting from the row interchanges.\n*\n*  =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  ipiv, info, ab = NumRu::Lapack.sgbtrf( m, kl, ku, ab, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 4)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 4)", argc);
  rblapack_m = argv[0];
  rblapack_kl = argv[1];
  rblapack_ku = argv[2];
  rblapack_ab = argv[3];
  if (rb_options != Qnil) {
  }

  if (!NA_IsNArray(rblapack_ab))
    rb_raise(rb_eArgError, "ab (4th argument) must be NArray");
  if (NA_RANK(rblapack_ab) != 2)
    rb_raise(rb_eArgError, "rank of ab (4th argument) must be %d", 2);
  n = NA_SHAPE1(rblapack_ab);
  ldab = NA_SHAPE0(rblapack_ab);
  if (NA_TYPE(rblapack_ab) != NA_SFLOAT)
    rblapack_ab = na_change_type(rblapack_ab, NA_SFLOAT);
  ab = NA_PTR_TYPE(rblapack_ab, real*);
  kl = NUM2INT(rblapack_kl);
  m = NUM2INT(rblapack_m);
  ku = NUM2INT(rblapack_ku);
  {
    int shape[1];
    shape[0] = MIN(m,n);
    rblapack_ipiv = na_make_object(NA_LINT, 1, shape, cNArray);
  }
  ipiv = NA_PTR_TYPE(rblapack_ipiv, integer*);
  {
    int shape[2];
    shape[0] = ldab;
    shape[1] = n;
    rblapack_ab_out__ = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  ab_out__ = NA_PTR_TYPE(rblapack_ab_out__, real*);
  MEMCPY(ab_out__, ab, real, NA_TOTAL(rblapack_ab));
  rblapack_ab = rblapack_ab_out__;
  ab = ab_out__;

  sgbtrf_(&m, &n, &kl, &ku, ab, &ldab, ipiv, &info);

  rblapack_info = INT2NUM(info);
  return rb_ary_new3(3, rblapack_ipiv, rblapack_info, rblapack_ab);
}

void
init_lapack_sgbtrf(VALUE mLapack){
  rb_define_module_function(mLapack, "sgbtrf", rblapack_sgbtrf, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}

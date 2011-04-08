#include "rb_lapack.h"

extern VOID ssbgst_(char *vect, char *uplo, integer *n, integer *ka, integer *kb, real *ab, integer *ldab, real *bb, integer *ldbb, real *x, integer *ldx, real *work, integer *info);

static VALUE sHelp, sUsage;

static VALUE
rblapack_ssbgst(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_vect;
  char vect; 
  VALUE rblapack_uplo;
  char uplo; 
  VALUE rblapack_ka;
  integer ka; 
  VALUE rblapack_kb;
  integer kb; 
  VALUE rblapack_ab;
  real *ab; 
  VALUE rblapack_bb;
  real *bb; 
  VALUE rblapack_x;
  real *x; 
  VALUE rblapack_info;
  integer info; 
  VALUE rblapack_ab_out__;
  real *ab_out__;
  real *work;

  integer ldab;
  integer n;
  integer ldbb;
  integer ldx;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  x, info, ab = NumRu::Lapack.ssbgst( vect, uplo, ka, kb, ab, bb, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE SSBGST( VECT, UPLO, N, KA, KB, AB, LDAB, BB, LDBB, X, LDX, WORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  SSBGST reduces a real symmetric-definite banded generalized\n*  eigenproblem  A*x = lambda*B*x  to standard form  C*y = lambda*y,\n*  such that C has the same bandwidth as A.\n*\n*  B must have been previously factorized as S**T*S by SPBSTF, using a\n*  split Cholesky factorization. A is overwritten by C = X**T*A*X, where\n*  X = S**(-1)*Q and Q is an orthogonal matrix chosen to preserve the\n*  bandwidth of A.\n*\n\n*  Arguments\n*  =========\n*\n*  VECT    (input) CHARACTER*1\n*          = 'N':  do not form the transformation matrix X;\n*          = 'V':  form X.\n*\n*  UPLO    (input) CHARACTER*1\n*          = 'U':  Upper triangle of A is stored;\n*          = 'L':  Lower triangle of A is stored.\n*\n*  N       (input) INTEGER\n*          The order of the matrices A and B.  N >= 0.\n*\n*  KA      (input) INTEGER\n*          The number of superdiagonals of the matrix A if UPLO = 'U',\n*          or the number of subdiagonals if UPLO = 'L'.  KA >= 0.\n*\n*  KB      (input) INTEGER\n*          The number of superdiagonals of the matrix B if UPLO = 'U',\n*          or the number of subdiagonals if UPLO = 'L'.  KA >= KB >= 0.\n*\n*  AB      (input/output) REAL array, dimension (LDAB,N)\n*          On entry, the upper or lower triangle of the symmetric band\n*          matrix A, stored in the first ka+1 rows of the array.  The\n*          j-th column of A is stored in the j-th column of the array AB\n*          as follows:\n*          if UPLO = 'U', AB(ka+1+i-j,j) = A(i,j) for max(1,j-ka)<=i<=j;\n*          if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=min(n,j+ka).\n*\n*          On exit, the transformed matrix X**T*A*X, stored in the same\n*          format as A.\n*\n*  LDAB    (input) INTEGER\n*          The leading dimension of the array AB.  LDAB >= KA+1.\n*\n*  BB      (input) REAL array, dimension (LDBB,N)\n*          The banded factor S from the split Cholesky factorization of\n*          B, as returned by SPBSTF, stored in the first KB+1 rows of\n*          the array.\n*\n*  LDBB    (input) INTEGER\n*          The leading dimension of the array BB.  LDBB >= KB+1.\n*\n*  X       (output) REAL array, dimension (LDX,N)\n*          If VECT = 'V', the n-by-n matrix X.\n*          If VECT = 'N', the array X is not referenced.\n*\n*  LDX     (input) INTEGER\n*          The leading dimension of the array X.\n*          LDX >= max(1,N) if VECT = 'V'; LDX >= 1 otherwise.\n*\n*  WORK    (workspace) REAL array, dimension (2*N)\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit\n*          < 0:  if INFO = -i, the i-th argument had an illegal value.\n*\n\n*  =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  x, info, ab = NumRu::Lapack.ssbgst( vect, uplo, ka, kb, ab, bb, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 6)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 6)", argc);
  rblapack_vect = argv[0];
  rblapack_uplo = argv[1];
  rblapack_ka = argv[2];
  rblapack_kb = argv[3];
  rblapack_ab = argv[4];
  rblapack_bb = argv[5];
  if (rb_options != Qnil) {
  }

  uplo = StringValueCStr(rblapack_uplo)[0];
  if (!NA_IsNArray(rblapack_bb))
    rb_raise(rb_eArgError, "bb (6th argument) must be NArray");
  if (NA_RANK(rblapack_bb) != 2)
    rb_raise(rb_eArgError, "rank of bb (6th argument) must be %d", 2);
  n = NA_SHAPE1(rblapack_bb);
  ldbb = NA_SHAPE0(rblapack_bb);
  if (NA_TYPE(rblapack_bb) != NA_SFLOAT)
    rblapack_bb = na_change_type(rblapack_bb, NA_SFLOAT);
  bb = NA_PTR_TYPE(rblapack_bb, real*);
  ka = NUM2INT(rblapack_ka);
  vect = StringValueCStr(rblapack_vect)[0];
  kb = NUM2INT(rblapack_kb);
  if (!NA_IsNArray(rblapack_ab))
    rb_raise(rb_eArgError, "ab (5th argument) must be NArray");
  if (NA_RANK(rblapack_ab) != 2)
    rb_raise(rb_eArgError, "rank of ab (5th argument) must be %d", 2);
  if (NA_SHAPE1(rblapack_ab) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of ab must be the same as shape 1 of bb");
  ldab = NA_SHAPE0(rblapack_ab);
  if (NA_TYPE(rblapack_ab) != NA_SFLOAT)
    rblapack_ab = na_change_type(rblapack_ab, NA_SFLOAT);
  ab = NA_PTR_TYPE(rblapack_ab, real*);
  ldx = lsame_(&vect,"V") ? MAX(1,n) : 1;
  {
    int shape[2];
    shape[0] = ldx;
    shape[1] = n;
    rblapack_x = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  x = NA_PTR_TYPE(rblapack_x, real*);
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
  work = ALLOC_N(real, (2*n));

  ssbgst_(&vect, &uplo, &n, &ka, &kb, ab, &ldab, bb, &ldbb, x, &ldx, work, &info);

  free(work);
  rblapack_info = INT2NUM(info);
  return rb_ary_new3(3, rblapack_x, rblapack_info, rblapack_ab);
}

void
init_lapack_ssbgst(VALUE mLapack){
  rb_define_module_function(mLapack, "ssbgst", rblapack_ssbgst, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}

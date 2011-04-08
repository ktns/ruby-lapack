#include "rb_lapack.h"

extern VOID dsyswapr_(char *uplo, integer *n, doublereal *a, integer *i1, integer *i2);

static VALUE sHelp, sUsage;

static VALUE
rblapack_dsyswapr(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_uplo;
  char uplo; 
  VALUE rblapack_a;
  doublereal *a; 
  VALUE rblapack_i1;
  integer i1; 
  VALUE rblapack_i2;
  integer i2; 
  VALUE rblapack_a_out__;
  doublereal *a_out__;

  integer lda;
  integer n;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  a = NumRu::Lapack.dsyswapr( uplo, a, i1, i2, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE DSYSWAPR( UPLO, N, A, I1, I2)\n\n*  Purpose\n*  =======\n*\n*  DSYSWAPR applies an elementary permutation on the rows and the columns of\n*  a symmetric matrix.\n*\n\n*  Arguments\n*  =========\n*\n*  UPLO    (input) CHARACTER*1\n*          Specifies whether the details of the factorization are stored\n*          as an upper or lower triangular matrix.\n*          = 'U':  Upper triangular, form is A = U*D*U**T;\n*          = 'L':  Lower triangular, form is A = L*D*L**T.\n*\n*  N       (input) INTEGER\n*          The order of the matrix A.  N >= 0.\n*\n*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)\n*          On entry, the NB diagonal matrix D and the multipliers\n*          used to obtain the factor U or L as computed by DSYTRF.\n*\n*          On exit, if INFO = 0, the (symmetric) inverse of the original\n*          matrix.  If UPLO = 'U', the upper triangular part of the\n*          inverse is formed and the part of A below the diagonal is not\n*          referenced; if UPLO = 'L' the lower triangular part of the\n*          inverse is formed and the part of A above the diagonal is\n*          not referenced.\n*\n*  I1      (input) INTEGER\n*          Index of the first row to swap\n*\n*  I2      (input) INTEGER\n*          Index of the second row to swap\n*\n\n*  =====================================================================\n*\n*     ..\n*     .. Local Scalars ..\n      LOGICAL            UPPER\n      INTEGER            I\n      DOUBLE PRECISION   TMP\n*\n*     .. External Functions ..\n      LOGICAL            LSAME\n      EXTERNAL           LSAME\n*     ..\n*     .. External Subroutines ..\n      EXTERNAL         DSWAP\n*     ..\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  a = NumRu::Lapack.dsyswapr( uplo, a, i1, i2, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 4)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 4)", argc);
  rblapack_uplo = argv[0];
  rblapack_a = argv[1];
  rblapack_i1 = argv[2];
  rblapack_i2 = argv[3];
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
  i1 = NUM2INT(rblapack_i1);
  i2 = NUM2INT(rblapack_i2);
  uplo = StringValueCStr(rblapack_uplo)[0];
  {
    int shape[2];
    shape[0] = lda;
    shape[1] = n;
    rblapack_a_out__ = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  a_out__ = NA_PTR_TYPE(rblapack_a_out__, doublereal*);
  MEMCPY(a_out__, a, doublereal, NA_TOTAL(rblapack_a));
  rblapack_a = rblapack_a_out__;
  a = a_out__;

  dsyswapr_(&uplo, &n, a, &i1, &i2);

  return rblapack_a;
}

void
init_lapack_dsyswapr(VALUE mLapack){
  rb_define_module_function(mLapack, "dsyswapr", rblapack_dsyswapr, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}

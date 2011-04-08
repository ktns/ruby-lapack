#include "rb_lapack.h"

extern VOID csyswapr_(char *uplo, integer *n, complex *a, integer *i1, integer *i2);

static VALUE sHelp, sUsage;

static VALUE
rblapack_csyswapr(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_uplo;
  char uplo; 
  VALUE rblapack_a;
  complex *a; 
  VALUE rblapack_i1;
  integer i1; 
  VALUE rblapack_i2;
  integer i2; 
  VALUE rblapack_a_out__;
  complex *a_out__;

  integer lda;
  integer n;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  a = NumRu::Lapack.csyswapr( uplo, a, i1, i2, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE CSYSWAPR( UPLO, N, A, I1, I2)\n\n*  Purpose\n*  =======\n*\n*  CSYSWAPR applies an elementary permutation on the rows and the columns of\n*  a symmetric matrix.\n*\n\n*  Arguments\n*  =========\n*\n*  UPLO    (input) CHARACTER*1\n*          Specifies whether the details of the factorization are stored\n*          as an upper or lower triangular matrix.\n*          = 'U':  Upper triangular, form is A = U*D*U**T;\n*          = 'L':  Lower triangular, form is A = L*D*L**T.\n*\n*  N       (input) INTEGER\n*          The order of the matrix A.  N >= 0.\n*\n*  A       (input/output) COMPLEX array, dimension (LDA,N)\n*          On entry, the NB diagonal matrix D and the multipliers\n*          used to obtain the factor U or L as computed by CSYTRF.\n*\n*          On exit, if INFO = 0, the (symmetric) inverse of the original\n*          matrix.  If UPLO = 'U', the upper triangular part of the\n*          inverse is formed and the part of A below the diagonal is not\n*          referenced; if UPLO = 'L' the lower triangular part of the\n*          inverse is formed and the part of A above the diagonal is\n*          not referenced.\n*\n*  I1      (input) INTEGER\n*          Index of the first row to swap\n*\n*  I2      (input) INTEGER\n*          Index of the second row to swap\n*\n\n*  =====================================================================\n*\n*     ..\n*     .. Local Scalars ..\n      LOGICAL            UPPER\n      INTEGER            I\n      COMPLEX            TMP\n*\n*     .. External Functions ..\n      LOGICAL            LSAME\n      EXTERNAL           LSAME\n*     ..\n*     .. External Subroutines ..\n      EXTERNAL           CSWAP\n*     ..\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  a = NumRu::Lapack.csyswapr( uplo, a, i1, i2, [:usage => usage, :help => help])\n");
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
  if (NA_TYPE(rblapack_a) != NA_SCOMPLEX)
    rblapack_a = na_change_type(rblapack_a, NA_SCOMPLEX);
  a = NA_PTR_TYPE(rblapack_a, complex*);
  i1 = NUM2INT(rblapack_i1);
  i2 = NUM2INT(rblapack_i2);
  uplo = StringValueCStr(rblapack_uplo)[0];
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

  csyswapr_(&uplo, &n, a, &i1, &i2);

  return rblapack_a;
}

void
init_lapack_csyswapr(VALUE mLapack){
  rb_define_module_function(mLapack, "csyswapr", rblapack_csyswapr, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}

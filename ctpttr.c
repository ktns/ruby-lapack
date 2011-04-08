#include "rb_lapack.h"

extern VOID ctpttr_(char *uplo, integer *n, complex *ap, complex *a, integer *lda, integer *info);

static VALUE sHelp, sUsage;

static VALUE
rblapack_ctpttr(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_uplo;
  char uplo; 
  VALUE rblapack_ap;
  complex *ap; 
  VALUE rblapack_a;
  complex *a; 
  VALUE rblapack_info;
  integer info; 

  integer ldap;
  integer lda;
  integer n;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  a, info = NumRu::Lapack.ctpttr( uplo, ap, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE CTPTTR( UPLO, N, AP, A, LDA, INFO )\n\n*  Purpose\n*  =======\n*\n*  CTPTTR copies a triangular matrix A from standard packed format (TP)\n*  to standard full format (TR).\n*\n\n*  Arguments\n*  =========\n*\n*  UPLO    (input) CHARACTER*1\n*          = 'U':  A is upper triangular.\n*          = 'L':  A is lower triangular.\n*\n*  N       (input) INTEGER\n*          The order of the matrix A. N >= 0.\n*\n*  AP      (input) COMPLEX array, dimension ( N*(N+1)/2 ),\n*          On entry, the upper or lower triangular matrix A, packed\n*          columnwise in a linear array. The j-th column of A is stored\n*          in the array AP as follows:\n*          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;\n*          if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = A(i,j) for j<=i<=n.\n*\n*  A       (output) COMPLEX array, dimension ( LDA, N )\n*          On exit, the triangular matrix A.  If UPLO = 'U', the leading\n*          N-by-N upper triangular part of A contains the upper\n*          triangular part of the matrix A, and the strictly lower\n*          triangular part of A is not referenced.  If UPLO = 'L', the\n*          leading N-by-N lower triangular part of A contains the lower\n*          triangular part of the matrix A, and the strictly upper\n*          triangular part of A is not referenced.\n*\n*  LDA     (input) INTEGER\n*          The leading dimension of the array A.  LDA >= max(1,N).\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit\n*          < 0:  if INFO = -i, the i-th argument had an illegal value\n*\n\n*  =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  a, info = NumRu::Lapack.ctpttr( uplo, ap, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 2)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 2)", argc);
  rblapack_uplo = argv[0];
  rblapack_ap = argv[1];
  if (rb_options != Qnil) {
  }

  uplo = StringValueCStr(rblapack_uplo)[0];
  if (!NA_IsNArray(rblapack_ap))
    rb_raise(rb_eArgError, "ap (2th argument) must be NArray");
  if (NA_RANK(rblapack_ap) != 1)
    rb_raise(rb_eArgError, "rank of ap (2th argument) must be %d", 1);
  ldap = NA_SHAPE0(rblapack_ap);
  if (NA_TYPE(rblapack_ap) != NA_SCOMPLEX)
    rblapack_ap = na_change_type(rblapack_ap, NA_SCOMPLEX);
  ap = NA_PTR_TYPE(rblapack_ap, complex*);
  n = ((int)sqrtf(ldap*8+1.0f)-1)/2;
  lda = MAX(1,n);
  {
    int shape[2];
    shape[0] = lda;
    shape[1] = n;
    rblapack_a = na_make_object(NA_SCOMPLEX, 2, shape, cNArray);
  }
  a = NA_PTR_TYPE(rblapack_a, complex*);

  ctpttr_(&uplo, &n, ap, a, &lda, &info);

  rblapack_info = INT2NUM(info);
  return rb_ary_new3(2, rblapack_a, rblapack_info);
}

void
init_lapack_ctpttr(VALUE mLapack){
  rb_define_module_function(mLapack, "ctpttr", rblapack_ctpttr, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}

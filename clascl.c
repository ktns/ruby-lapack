#include "rb_lapack.h"

extern VOID clascl_(char *type, integer *kl, integer *ku, real *cfrom, real *cto, integer *m, integer *n, complex *a, integer *lda, integer *info);

static VALUE sHelp, sUsage;

static VALUE
rblapack_clascl(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_type;
  char type; 
  VALUE rblapack_kl;
  integer kl; 
  VALUE rblapack_ku;
  integer ku; 
  VALUE rblapack_cfrom;
  real cfrom; 
  VALUE rblapack_cto;
  real cto; 
  VALUE rblapack_m;
  integer m; 
  VALUE rblapack_a;
  complex *a; 
  VALUE rblapack_info;
  integer info; 
  VALUE rblapack_a_out__;
  complex *a_out__;

  integer lda;
  integer n;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  info, a = NumRu::Lapack.clascl( type, kl, ku, cfrom, cto, m, a, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE CLASCL( TYPE, KL, KU, CFROM, CTO, M, N, A, LDA, INFO )\n\n*  Purpose\n*  =======\n*\n*  CLASCL multiplies the M by N complex matrix A by the real scalar\n*  CTO/CFROM.  This is done without over/underflow as long as the final\n*  result CTO*A(I,J)/CFROM does not over/underflow. TYPE specifies that\n*  A may be full, upper triangular, lower triangular, upper Hessenberg,\n*  or banded.\n*\n\n*  Arguments\n*  =========\n*\n*  TYPE    (input) CHARACTER*1\n*          TYPE indices the storage type of the input matrix.\n*          = 'G':  A is a full matrix.\n*          = 'L':  A is a lower triangular matrix.\n*          = 'U':  A is an upper triangular matrix.\n*          = 'H':  A is an upper Hessenberg matrix.\n*          = 'B':  A is a symmetric band matrix with lower bandwidth KL\n*                  and upper bandwidth KU and with the only the lower\n*                  half stored.\n*          = 'Q':  A is a symmetric band matrix with lower bandwidth KL\n*                  and upper bandwidth KU and with the only the upper\n*                  half stored.\n*          = 'Z':  A is a band matrix with lower bandwidth KL and upper\n*                  bandwidth KU. See CGBTRF for storage details.\n*\n*  KL      (input) INTEGER\n*          The lower bandwidth of A.  Referenced only if TYPE = 'B',\n*          'Q' or 'Z'.\n*\n*  KU      (input) INTEGER\n*          The upper bandwidth of A.  Referenced only if TYPE = 'B',\n*          'Q' or 'Z'.\n*\n*  CFROM   (input) REAL\n*  CTO     (input) REAL\n*          The matrix A is multiplied by CTO/CFROM. A(I,J) is computed\n*          without over/underflow if the final result CTO*A(I,J)/CFROM\n*          can be represented without over/underflow.  CFROM must be\n*          nonzero.\n*\n*  M       (input) INTEGER\n*          The number of rows of the matrix A.  M >= 0.\n*\n*  N       (input) INTEGER\n*          The number of columns of the matrix A.  N >= 0.\n*\n*  A       (input/output) COMPLEX array, dimension (LDA,N)\n*          The matrix to be multiplied by CTO/CFROM.  See TYPE for the\n*          storage type.\n*\n*  LDA     (input) INTEGER\n*          The leading dimension of the array A.  LDA >= max(1,M).\n*\n*  INFO    (output) INTEGER\n*          0  - successful exit\n*          <0 - if INFO = -i, the i-th argument had an illegal value.\n*\n\n*  =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  info, a = NumRu::Lapack.clascl( type, kl, ku, cfrom, cto, m, a, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 7)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 7)", argc);
  rblapack_type = argv[0];
  rblapack_kl = argv[1];
  rblapack_ku = argv[2];
  rblapack_cfrom = argv[3];
  rblapack_cto = argv[4];
  rblapack_m = argv[5];
  rblapack_a = argv[6];
  if (rb_options != Qnil) {
  }

  if (!NA_IsNArray(rblapack_a))
    rb_raise(rb_eArgError, "a (7th argument) must be NArray");
  if (NA_RANK(rblapack_a) != 2)
    rb_raise(rb_eArgError, "rank of a (7th argument) must be %d", 2);
  n = NA_SHAPE1(rblapack_a);
  lda = NA_SHAPE0(rblapack_a);
  if (NA_TYPE(rblapack_a) != NA_SCOMPLEX)
    rblapack_a = na_change_type(rblapack_a, NA_SCOMPLEX);
  a = NA_PTR_TYPE(rblapack_a, complex*);
  kl = NUM2INT(rblapack_kl);
  m = NUM2INT(rblapack_m);
  cfrom = (real)NUM2DBL(rblapack_cfrom);
  type = StringValueCStr(rblapack_type)[0];
  cto = (real)NUM2DBL(rblapack_cto);
  ku = NUM2INT(rblapack_ku);
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

  clascl_(&type, &kl, &ku, &cfrom, &cto, &m, &n, a, &lda, &info);

  rblapack_info = INT2NUM(info);
  return rb_ary_new3(2, rblapack_info, rblapack_a);
}

void
init_lapack_clascl(VALUE mLapack){
  rb_define_module_function(mLapack, "clascl", rblapack_clascl, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}

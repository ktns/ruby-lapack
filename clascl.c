#include "rb_lapack.h"

static VALUE
rb_clascl(int argc, VALUE *argv, VALUE self){
  VALUE rb_type;
  char type; 
  VALUE rb_kl;
  integer kl; 
  VALUE rb_ku;
  integer ku; 
  VALUE rb_cfrom;
  real cfrom; 
  VALUE rb_cto;
  real cto; 
  VALUE rb_m;
  integer m; 
  VALUE rb_a;
  complex *a; 
  VALUE rb_info;
  integer info; 
  VALUE rb_a_out__;
  complex *a_out__;

  integer lda;
  integer n;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  info, a = NumRu::Lapack.clascl( type, kl, ku, cfrom, cto, m, a)\n    or\n  NumRu::Lapack.clascl  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE CLASCL( TYPE, KL, KU, CFROM, CTO, M, N, A, LDA, INFO )\n\n*  Purpose\n*  =======\n*\n*  CLASCL multiplies the M by N complex matrix A by the real scalar\n*  CTO/CFROM.  This is done without over/underflow as long as the final\n*  result CTO*A(I,J)/CFROM does not over/underflow. TYPE specifies that\n*  A may be full, upper triangular, lower triangular, upper Hessenberg,\n*  or banded.\n*\n\n*  Arguments\n*  =========\n*\n*  TYPE    (input) CHARACTER*1\n*          TYPE indices the storage type of the input matrix.\n*          = 'G':  A is a full matrix.\n*          = 'L':  A is a lower triangular matrix.\n*          = 'U':  A is an upper triangular matrix.\n*          = 'H':  A is an upper Hessenberg matrix.\n*          = 'B':  A is a symmetric band matrix with lower bandwidth KL\n*                  and upper bandwidth KU and with the only the lower\n*                  half stored.\n*          = 'Q':  A is a symmetric band matrix with lower bandwidth KL\n*                  and upper bandwidth KU and with the only the upper\n*                  half stored.\n*          = 'Z':  A is a band matrix with lower bandwidth KL and upper\n*                  bandwidth KU.\n*\n*  KL      (input) INTEGER\n*          The lower bandwidth of A.  Referenced only if TYPE = 'B',\n*          'Q' or 'Z'.\n*\n*  KU      (input) INTEGER\n*          The upper bandwidth of A.  Referenced only if TYPE = 'B',\n*          'Q' or 'Z'.\n*\n*  CFROM   (input) REAL\n*  CTO     (input) REAL\n*          The matrix A is multiplied by CTO/CFROM. A(I,J) is computed\n*          without over/underflow if the final result CTO*A(I,J)/CFROM\n*          can be represented without over/underflow.  CFROM must be\n*          nonzero.\n*\n*  M       (input) INTEGER\n*          The number of rows of the matrix A.  M >= 0.\n*\n*  N       (input) INTEGER\n*          The number of columns of the matrix A.  N >= 0.\n*\n*  A       (input/output) COMPLEX array, dimension (LDA,N)\n*          The matrix to be multiplied by CTO/CFROM.  See TYPE for the\n*          storage type.\n*\n*  LDA     (input) INTEGER\n*          The leading dimension of the array A.  LDA >= max(1,M).\n*\n*  INFO    (output) INTEGER\n*          0  - successful exit\n*          <0 - if INFO = -i, the i-th argument had an illegal value.\n*\n\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 7)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 7)", argc);
  rb_type = argv[0];
  rb_kl = argv[1];
  rb_ku = argv[2];
  rb_cfrom = argv[3];
  rb_cto = argv[4];
  rb_m = argv[5];
  rb_a = argv[6];

  type = StringValueCStr(rb_type)[0];
  kl = NUM2INT(rb_kl);
  ku = NUM2INT(rb_ku);
  cfrom = (real)NUM2DBL(rb_cfrom);
  cto = (real)NUM2DBL(rb_cto);
  m = NUM2INT(rb_m);
  if (!NA_IsNArray(rb_a))
    rb_raise(rb_eArgError, "a (7th argument) must be NArray");
  if (NA_RANK(rb_a) != 2)
    rb_raise(rb_eArgError, "rank of a (7th argument) must be %d", 2);
  lda = NA_SHAPE0(rb_a);
  n = NA_SHAPE1(rb_a);
  if (NA_TYPE(rb_a) != NA_SCOMPLEX)
    rb_a = na_change_type(rb_a, NA_SCOMPLEX);
  a = NA_PTR_TYPE(rb_a, complex*);
  {
    int shape[2];
    shape[0] = lda;
    shape[1] = n;
    rb_a_out__ = na_make_object(NA_SCOMPLEX, 2, shape, cNArray);
  }
  a_out__ = NA_PTR_TYPE(rb_a_out__, complex*);
  MEMCPY(a_out__, a, complex, NA_TOTAL(rb_a));
  rb_a = rb_a_out__;
  a = a_out__;

  clascl_(&type, &kl, &ku, &cfrom, &cto, &m, &n, a, &lda, &info);

  rb_info = INT2NUM(info);
  return rb_ary_new3(2, rb_info, rb_a);
}

void
init_lapack_clascl(VALUE mLapack){
  rb_define_module_function(mLapack, "clascl", rb_clascl, -1);
}

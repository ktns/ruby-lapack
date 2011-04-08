#include "rb_lapack.h"

extern VOID claqr1_(integer *n, complex *h, integer *ldh, complex *s1, complex *s2, complex *v);

static VALUE sHelp, sUsage;

static VALUE
rblapack_claqr1(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_h;
  complex *h; 
  VALUE rblapack_s1;
  complex s1; 
  VALUE rblapack_s2;
  complex s2; 
  VALUE rblapack_v;
  complex *v; 

  integer ldh;
  integer n;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  v = NumRu::Lapack.claqr1( h, s1, s2, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE CLAQR1( N, H, LDH, S1, S2, V )\n\n*       Given a 2-by-2 or 3-by-3 matrix H, CLAQR1 sets v to a\n*       scalar multiple of the first column of the product\n*\n*       (*)  K = (H - s1*I)*(H - s2*I)\n*\n*       scaling to avoid overflows and most underflows.\n*\n*       This is useful for starting double implicit shift bulges\n*       in the QR algorithm.\n*\n*\n\n*       N      (input) integer\n*              Order of the matrix H. N must be either 2 or 3.\n*\n*       H      (input) COMPLEX array of dimension (LDH,N)\n*              The 2-by-2 or 3-by-3 matrix H in (*).\n*\n*       LDH    (input) integer\n*              The leading dimension of H as declared in\n*              the calling procedure.  LDH.GE.N\n*\n*       S1     (input) COMPLEX\n*       S2     S1 and S2 are the shifts defining K in (*) above.\n*\n*       V      (output) COMPLEX array of dimension N\n*              A scalar multiple of the first column of the\n*              matrix K in (*).\n*\n\n*     ================================================================\n*     Based on contributions by\n*        Karen Braman and Ralph Byers, Department of Mathematics,\n*        University of Kansas, USA\n*\n*     ================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  v = NumRu::Lapack.claqr1( h, s1, s2, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 3)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 3)", argc);
  rblapack_h = argv[0];
  rblapack_s1 = argv[1];
  rblapack_s2 = argv[2];
  if (rb_options != Qnil) {
  }

  s1.r = (real)NUM2DBL(rb_funcall(rblapack_s1, rb_intern("real"), 0));
  s1.i = (real)NUM2DBL(rb_funcall(rblapack_s1, rb_intern("imag"), 0));
  if (!NA_IsNArray(rblapack_h))
    rb_raise(rb_eArgError, "h (1th argument) must be NArray");
  if (NA_RANK(rblapack_h) != 2)
    rb_raise(rb_eArgError, "rank of h (1th argument) must be %d", 2);
  n = NA_SHAPE1(rblapack_h);
  ldh = NA_SHAPE0(rblapack_h);
  if (NA_TYPE(rblapack_h) != NA_SCOMPLEX)
    rblapack_h = na_change_type(rblapack_h, NA_SCOMPLEX);
  h = NA_PTR_TYPE(rblapack_h, complex*);
  s2.r = (real)NUM2DBL(rb_funcall(rblapack_s2, rb_intern("real"), 0));
  s2.i = (real)NUM2DBL(rb_funcall(rblapack_s2, rb_intern("imag"), 0));
  {
    int shape[1];
    shape[0] = n;
    rblapack_v = na_make_object(NA_SCOMPLEX, 1, shape, cNArray);
  }
  v = NA_PTR_TYPE(rblapack_v, complex*);

  claqr1_(&n, h, &ldh, &s1, &s2, v);

  return rblapack_v;
}

void
init_lapack_claqr1(VALUE mLapack){
  rb_define_module_function(mLapack, "claqr1", rblapack_claqr1, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}

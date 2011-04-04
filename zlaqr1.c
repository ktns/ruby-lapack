#include "rb_lapack.h"

extern VOID zlaqr1_(integer *n, doublecomplex *h, integer *ldh, doublecomplex *s1, doublecomplex *s2, doublecomplex *v);

static VALUE
rb_zlaqr1(int argc, VALUE *argv, VALUE self){
  VALUE rb_h;
  doublecomplex *h; 
  VALUE rb_s1;
  doublecomplex s1; 
  VALUE rb_s2;
  doublecomplex s2; 
  VALUE rb_v;
  doublecomplex *v; 

  integer ldh;
  integer n;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  v = NumRu::Lapack.zlaqr1( h, s1, s2)\n    or\n  NumRu::Lapack.zlaqr1  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE ZLAQR1( N, H, LDH, S1, S2, V )\n\n*       Given a 2-by-2 or 3-by-3 matrix H, ZLAQR1 sets v to a\n*       scalar multiple of the first column of the product\n*\n*       (*)  K = (H - s1*I)*(H - s2*I)\n*\n*       scaling to avoid overflows and most underflows.\n*\n*       This is useful for starting double implicit shift bulges\n*       in the QR algorithm.\n*\n*\n\n*       N      (input) integer\n*              Order of the matrix H. N must be either 2 or 3.\n*\n*       H      (input) COMPLEX*16 array of dimension (LDH,N)\n*              The 2-by-2 or 3-by-3 matrix H in (*).\n*\n*       LDH    (input) integer\n*              The leading dimension of H as declared in\n*              the calling procedure.  LDH.GE.N\n*\n*       S1     (input) COMPLEX*16\n*       S2     S1 and S2 are the shifts defining K in (*) above.\n*\n*       V      (output) COMPLEX*16 array of dimension N\n*              A scalar multiple of the first column of the\n*              matrix K in (*).\n*\n\n*     ================================================================\n*     Based on contributions by\n*        Karen Braman and Ralph Byers, Department of Mathematics,\n*        University of Kansas, USA\n*\n*     ================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 3)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 3)", argc);
  rb_h = argv[0];
  rb_s1 = argv[1];
  rb_s2 = argv[2];

  s1.r = NUM2DBL(rb_funcall(rb_s1, rb_intern("real"), 0));
  s1.i = NUM2DBL(rb_funcall(rb_s1, rb_intern("imag"), 0));
  if (!NA_IsNArray(rb_h))
    rb_raise(rb_eArgError, "h (1th argument) must be NArray");
  if (NA_RANK(rb_h) != 2)
    rb_raise(rb_eArgError, "rank of h (1th argument) must be %d", 2);
  n = NA_SHAPE1(rb_h);
  ldh = NA_SHAPE0(rb_h);
  if (NA_TYPE(rb_h) != NA_DCOMPLEX)
    rb_h = na_change_type(rb_h, NA_DCOMPLEX);
  h = NA_PTR_TYPE(rb_h, doublecomplex*);
  s2.r = NUM2DBL(rb_funcall(rb_s2, rb_intern("real"), 0));
  s2.i = NUM2DBL(rb_funcall(rb_s2, rb_intern("imag"), 0));
  {
    int shape[1];
    shape[0] = n;
    rb_v = na_make_object(NA_DCOMPLEX, 1, shape, cNArray);
  }
  v = NA_PTR_TYPE(rb_v, doublecomplex*);

  zlaqr1_(&n, h, &ldh, &s1, &s2, v);

  return rb_v;
}

void
init_lapack_zlaqr1(VALUE mLapack){
  rb_define_module_function(mLapack, "zlaqr1", rb_zlaqr1, -1);
}

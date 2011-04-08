#include "rb_lapack.h"

extern VOID dlaqr1_(integer *n, doublereal *h, integer *ldh, doublereal *sr1, doublereal *si1, doublereal *sr2, doublereal *si2, doublereal *v);

static VALUE sHelp, sUsage;

static VALUE
rblapack_dlaqr1(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_h;
  doublereal *h; 
  VALUE rblapack_sr1;
  doublereal sr1; 
  VALUE rblapack_si1;
  doublereal si1; 
  VALUE rblapack_sr2;
  doublereal sr2; 
  VALUE rblapack_si2;
  doublereal si2; 
  VALUE rblapack_v;
  doublereal *v; 

  integer ldh;
  integer n;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  v = NumRu::Lapack.dlaqr1( h, sr1, si1, sr2, si2, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE DLAQR1( N, H, LDH, SR1, SI1, SR2, SI2, V )\n\n*       Given a 2-by-2 or 3-by-3 matrix H, DLAQR1 sets v to a\n*       scalar multiple of the first column of the product\n*\n*       (*)  K = (H - (sr1 + i*si1)*I)*(H - (sr2 + i*si2)*I)\n*\n*       scaling to avoid overflows and most underflows. It\n*       is assumed that either\n*\n*               1) sr1 = sr2 and si1 = -si2\n*           or\n*               2) si1 = si2 = 0.\n*\n*       This is useful for starting double implicit shift bulges\n*       in the QR algorithm.\n*\n*\n\n*       N      (input) integer\n*              Order of the matrix H. N must be either 2 or 3.\n*\n*       H      (input) DOUBLE PRECISION array of dimension (LDH,N)\n*              The 2-by-2 or 3-by-3 matrix H in (*).\n*\n*       LDH    (input) integer\n*              The leading dimension of H as declared in\n*              the calling procedure.  LDH.GE.N\n*\n*       SR1    (input) DOUBLE PRECISION\n*       SI1    The shifts in (*).\n*       SR2\n*       SI2\n*\n*       V      (output) DOUBLE PRECISION array of dimension N\n*              A scalar multiple of the first column of the\n*              matrix K in (*).\n*\n\n*     ================================================================\n*     Based on contributions by\n*        Karen Braman and Ralph Byers, Department of Mathematics,\n*        University of Kansas, USA\n*\n*     ================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  v = NumRu::Lapack.dlaqr1( h, sr1, si1, sr2, si2, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 5)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 5)", argc);
  rblapack_h = argv[0];
  rblapack_sr1 = argv[1];
  rblapack_si1 = argv[2];
  rblapack_sr2 = argv[3];
  rblapack_si2 = argv[4];
  if (rb_options != Qnil) {
  }

  si1 = NUM2DBL(rblapack_si1);
  si2 = NUM2DBL(rblapack_si2);
  sr1 = NUM2DBL(rblapack_sr1);
  sr2 = NUM2DBL(rblapack_sr2);
  if (!NA_IsNArray(rblapack_h))
    rb_raise(rb_eArgError, "h (1th argument) must be NArray");
  if (NA_RANK(rblapack_h) != 2)
    rb_raise(rb_eArgError, "rank of h (1th argument) must be %d", 2);
  n = NA_SHAPE1(rblapack_h);
  ldh = NA_SHAPE0(rblapack_h);
  if (NA_TYPE(rblapack_h) != NA_DFLOAT)
    rblapack_h = na_change_type(rblapack_h, NA_DFLOAT);
  h = NA_PTR_TYPE(rblapack_h, doublereal*);
  {
    int shape[1];
    shape[0] = n;
    rblapack_v = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  v = NA_PTR_TYPE(rblapack_v, doublereal*);

  dlaqr1_(&n, h, &ldh, &sr1, &si1, &sr2, &si2, v);

  return rblapack_v;
}

void
init_lapack_dlaqr1(VALUE mLapack){
  rb_define_module_function(mLapack, "dlaqr1", rblapack_dlaqr1, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}

#include "rb_lapack.h"

extern VOID zlaesy_(doublecomplex *a, doublecomplex *b, doublecomplex *c, doublecomplex *rt1, doublecomplex *rt2, doublecomplex *evscal, doublecomplex *cs1, doublecomplex *sn1);

static VALUE sHelp, sUsage;

static VALUE
rblapack_zlaesy(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_a;
  doublecomplex a; 
  VALUE rblapack_b;
  doublecomplex b; 
  VALUE rblapack_c;
  doublecomplex c; 
  VALUE rblapack_rt1;
  doublecomplex rt1; 
  VALUE rblapack_rt2;
  doublecomplex rt2; 
  VALUE rblapack_evscal;
  doublecomplex evscal; 
  VALUE rblapack_cs1;
  doublecomplex cs1; 
  VALUE rblapack_sn1;
  doublecomplex sn1; 


  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  rt1, rt2, evscal, cs1, sn1 = NumRu::Lapack.zlaesy( a, b, c, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE ZLAESY( A, B, C, RT1, RT2, EVSCAL, CS1, SN1 )\n\n*  Purpose\n*  =======\n*\n*  ZLAESY computes the eigendecomposition of a 2-by-2 symmetric matrix\n*     ( ( A, B );( B, C ) )\n*  provided the norm of the matrix of eigenvectors is larger than\n*  some threshold value.\n*\n*  RT1 is the eigenvalue of larger absolute value, and RT2 of\n*  smaller absolute value.  If the eigenvectors are computed, then\n*  on return ( CS1, SN1 ) is the unit eigenvector for RT1, hence\n*\n*  [  CS1     SN1   ] . [ A  B ] . [ CS1    -SN1   ] = [ RT1  0  ]\n*  [ -SN1     CS1   ]   [ B  C ]   [ SN1     CS1   ]   [  0  RT2 ]\n*\n\n*  Arguments\n*  =========\n*\n*  A       (input) COMPLEX*16\n*          The ( 1, 1 ) element of input matrix.\n*\n*  B       (input) COMPLEX*16\n*          The ( 1, 2 ) element of input matrix.  The ( 2, 1 ) element\n*          is also given by B, since the 2-by-2 matrix is symmetric.\n*\n*  C       (input) COMPLEX*16\n*          The ( 2, 2 ) element of input matrix.\n*\n*  RT1     (output) COMPLEX*16\n*          The eigenvalue of larger modulus.\n*\n*  RT2     (output) COMPLEX*16\n*          The eigenvalue of smaller modulus.\n*\n*  EVSCAL  (output) COMPLEX*16\n*          The complex value by which the eigenvector matrix was scaled\n*          to make it orthonormal.  If EVSCAL is zero, the eigenvectors\n*          were not computed.  This means one of two things:  the 2-by-2\n*          matrix could not be diagonalized, or the norm of the matrix\n*          of eigenvectors before scaling was larger than the threshold\n*          value THRESH (set below).\n*\n*  CS1     (output) COMPLEX*16\n*  SN1     (output) COMPLEX*16\n*          If EVSCAL .NE. 0,  ( CS1, SN1 ) is the unit right eigenvector\n*          for RT1.\n*\n\n* =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  rt1, rt2, evscal, cs1, sn1 = NumRu::Lapack.zlaesy( a, b, c, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 3)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 3)", argc);
  rblapack_a = argv[0];
  rblapack_b = argv[1];
  rblapack_c = argv[2];
  if (rb_options != Qnil) {
  }

  a.r = NUM2DBL(rb_funcall(rblapack_a, rb_intern("real"), 0));
  a.i = NUM2DBL(rb_funcall(rblapack_a, rb_intern("imag"), 0));
  b.r = NUM2DBL(rb_funcall(rblapack_b, rb_intern("real"), 0));
  b.i = NUM2DBL(rb_funcall(rblapack_b, rb_intern("imag"), 0));
  c.r = NUM2DBL(rb_funcall(rblapack_c, rb_intern("real"), 0));
  c.i = NUM2DBL(rb_funcall(rblapack_c, rb_intern("imag"), 0));

  zlaesy_(&a, &b, &c, &rt1, &rt2, &evscal, &cs1, &sn1);

  rblapack_rt1 = rb_funcall(rb_gv_get("Complex"), rb_intern("new"), 2, rb_float_new((double)(rt1.r)), rb_float_new((double)(rt1.i)));
  rblapack_rt2 = rb_funcall(rb_gv_get("Complex"), rb_intern("new"), 2, rb_float_new((double)(rt2.r)), rb_float_new((double)(rt2.i)));
  rblapack_evscal = rb_funcall(rb_gv_get("Complex"), rb_intern("new"), 2, rb_float_new((double)(evscal.r)), rb_float_new((double)(evscal.i)));
  rblapack_cs1 = rb_funcall(rb_gv_get("Complex"), rb_intern("new"), 2, rb_float_new((double)(cs1.r)), rb_float_new((double)(cs1.i)));
  rblapack_sn1 = rb_funcall(rb_gv_get("Complex"), rb_intern("new"), 2, rb_float_new((double)(sn1.r)), rb_float_new((double)(sn1.i)));
  return rb_ary_new3(5, rblapack_rt1, rblapack_rt2, rblapack_evscal, rblapack_cs1, rblapack_sn1);
}

void
init_lapack_zlaesy(VALUE mLapack){
  rb_define_module_function(mLapack, "zlaesy", rblapack_zlaesy, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}

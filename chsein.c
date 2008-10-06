#include "rb_lapack.h"

static VALUE
rb_chsein(int argc, VALUE *argv, VALUE self){
  VALUE rb_side;
  char side; 
  VALUE rb_eigsrc;
  char eigsrc; 
  VALUE rb_initv;
  char initv; 
  VALUE rb_select;
  logical *select; 
  VALUE rb_h;
  complex *h; 
  VALUE rb_w;
  complex *w; 
  VALUE rb_vl;
  complex *vl; 
  VALUE rb_vr;
  complex *vr; 
  VALUE rb_m;
  integer m; 
  VALUE rb_ifaill;
  integer *ifaill; 
  VALUE rb_ifailr;
  integer *ifailr; 
  VALUE rb_info;
  integer info; 
  VALUE rb_w_out__;
  complex *w_out__;
  VALUE rb_vl_out__;
  complex *vl_out__;
  VALUE rb_vr_out__;
  complex *vr_out__;
  complex *work;
  real *rwork;

  integer n;
  integer ldh;
  integer ldvl;
  integer mm;
  integer ldvr;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  m, ifaill, ifailr, info, w, vl, vr = NumRu::Lapack.chsein( side, eigsrc, initv, select, h, w, vl, vr)\n    or\n  NumRu::Lapack.chsein  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE CHSEIN( SIDE, EIGSRC, INITV, SELECT, N, H, LDH, W, VL, LDVL, VR, LDVR, MM, M, WORK, RWORK, IFAILL, IFAILR, INFO )\n\n*  Purpose\n*  =======\n*\n*  CHSEIN uses inverse iteration to find specified right and/or left\n*  eigenvectors of a complex upper Hessenberg matrix H.\n*\n*  The right eigenvector x and the left eigenvector y of the matrix H\n*  corresponding to an eigenvalue w are defined by:\n*\n*               H * x = w * x,     y**h * H = w * y**h\n*\n*  where y**h denotes the conjugate transpose of the vector y.\n*\n\n*  Arguments\n*  =========\n*\n*  SIDE    (input) CHARACTER*1\n*          = 'R': compute right eigenvectors only;\n*          = 'L': compute left eigenvectors only;\n*          = 'B': compute both right and left eigenvectors.\n*\n*  EIGSRC  (input) CHARACTER*1\n*          Specifies the source of eigenvalues supplied in W:\n*          = 'Q': the eigenvalues were found using CHSEQR; thus, if\n*                 H has zero subdiagonal elements, and so is\n*                 block-triangular, then the j-th eigenvalue can be\n*                 assumed to be an eigenvalue of the block containing\n*                 the j-th row/column.  This property allows CHSEIN to\n*                 perform inverse iteration on just one diagonal block.\n*          = 'N': no assumptions are made on the correspondence\n*                 between eigenvalues and diagonal blocks.  In this\n*                 case, CHSEIN must always perform inverse iteration\n*                 using the whole matrix H.\n*\n*  INITV   (input) CHARACTER*1\n*          = 'N': no initial vectors are supplied;\n*          = 'U': user-supplied initial vectors are stored in the arrays\n*                 VL and/or VR.\n*\n*  SELECT  (input) LOGICAL array, dimension (N)\n*          Specifies the eigenvectors to be computed. To select the\n*          eigenvector corresponding to the eigenvalue W(j),\n*          SELECT(j) must be set to .TRUE..\n*\n*  N       (input) INTEGER\n*          The order of the matrix H.  N >= 0.\n*\n*  H       (input) COMPLEX array, dimension (LDH,N)\n*          The upper Hessenberg matrix H.\n*\n*  LDH     (input) INTEGER\n*          The leading dimension of the array H.  LDH >= max(1,N).\n*\n*  W       (input/output) COMPLEX array, dimension (N)\n*          On entry, the eigenvalues of H.\n*          On exit, the real parts of W may have been altered since\n*          close eigenvalues are perturbed slightly in searching for\n*          independent eigenvectors.\n*\n*  VL      (input/output) COMPLEX array, dimension (LDVL,MM)\n*          On entry, if INITV = 'U' and SIDE = 'L' or 'B', VL must\n*          contain starting vectors for the inverse iteration for the\n*          left eigenvectors; the starting vector for each eigenvector\n*          must be in the same column in which the eigenvector will be\n*          stored.\n*          On exit, if SIDE = 'L' or 'B', the left eigenvectors\n*          specified by SELECT will be stored consecutively in the\n*          columns of VL, in the same order as their eigenvalues.\n*          If SIDE = 'R', VL is not referenced.\n*\n*  LDVL    (input) INTEGER\n*          The leading dimension of the array VL.\n*          LDVL >= max(1,N) if SIDE = 'L' or 'B'; LDVL >= 1 otherwise.\n*\n*  VR      (input/output) COMPLEX array, dimension (LDVR,MM)\n*          On entry, if INITV = 'U' and SIDE = 'R' or 'B', VR must\n*          contain starting vectors for the inverse iteration for the\n*          right eigenvectors; the starting vector for each eigenvector\n*          must be in the same column in which the eigenvector will be\n*          stored.\n*          On exit, if SIDE = 'R' or 'B', the right eigenvectors\n*          specified by SELECT will be stored consecutively in the\n*          columns of VR, in the same order as their eigenvalues.\n*          If SIDE = 'L', VR is not referenced.\n*\n*  LDVR    (input) INTEGER\n*          The leading dimension of the array VR.\n*          LDVR >= max(1,N) if SIDE = 'R' or 'B'; LDVR >= 1 otherwise.\n*\n*  MM      (input) INTEGER\n*          The number of columns in the arrays VL and/or VR. MM >= M.\n*\n*  M       (output) INTEGER\n*          The number of columns in the arrays VL and/or VR required to\n*          store the eigenvectors (= the number of .TRUE. elements in\n*          SELECT).\n*\n*  WORK    (workspace) COMPLEX array, dimension (N*N)\n*\n*  RWORK   (workspace) REAL array, dimension (N)\n*\n*  IFAILL  (output) INTEGER array, dimension (MM)\n*          If SIDE = 'L' or 'B', IFAILL(i) = j > 0 if the left\n*          eigenvector in the i-th column of VL (corresponding to the\n*          eigenvalue w(j)) failed to converge; IFAILL(i) = 0 if the\n*          eigenvector converged satisfactorily.\n*          If SIDE = 'R', IFAILL is not referenced.\n*\n*  IFAILR  (output) INTEGER array, dimension (MM)\n*          If SIDE = 'R' or 'B', IFAILR(i) = j > 0 if the right\n*          eigenvector in the i-th column of VR (corresponding to the\n*          eigenvalue w(j)) failed to converge; IFAILR(i) = 0 if the\n*          eigenvector converged satisfactorily.\n*          If SIDE = 'L', IFAILR is not referenced.\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit\n*          < 0:  if INFO = -i, the i-th argument had an illegal value\n*          > 0:  if INFO = i, i is the number of eigenvectors which\n*                failed to converge; see IFAILL and IFAILR for further\n*                details.\n*\n\n*  Further Details\n*  ===============\n*\n*  Each eigenvector is normalized so that the element of largest\n*  magnitude has magnitude 1; here the magnitude of a complex number\n*  (x,y) is taken to be |x|+|y|.\n*\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 8)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 8)", argc);
  rb_side = argv[0];
  rb_eigsrc = argv[1];
  rb_initv = argv[2];
  rb_select = argv[3];
  rb_h = argv[4];
  rb_w = argv[5];
  rb_vl = argv[6];
  rb_vr = argv[7];

  side = StringValueCStr(rb_side)[0];
  eigsrc = StringValueCStr(rb_eigsrc)[0];
  initv = StringValueCStr(rb_initv)[0];
  if (!NA_IsNArray(rb_select))
    rb_raise(rb_eArgError, "select (4th argument) must be NArray");
  if (NA_RANK(rb_select) != 1)
    rb_raise(rb_eArgError, "rank of select (4th argument) must be %d", 1);
  n = NA_SHAPE0(rb_select);
  if (NA_TYPE(rb_select) != NA_LINT)
    rb_select = na_change_type(rb_select, NA_LINT);
  select = NA_PTR_TYPE(rb_select, logical*);
  if (!NA_IsNArray(rb_h))
    rb_raise(rb_eArgError, "h (5th argument) must be NArray");
  if (NA_RANK(rb_h) != 2)
    rb_raise(rb_eArgError, "rank of h (5th argument) must be %d", 2);
  ldh = NA_SHAPE0(rb_h);
  if (NA_SHAPE1(rb_h) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of h must be the same as shape 0 of select");
  if (NA_TYPE(rb_h) != NA_SCOMPLEX)
    rb_h = na_change_type(rb_h, NA_SCOMPLEX);
  h = NA_PTR_TYPE(rb_h, complex*);
  if (!NA_IsNArray(rb_w))
    rb_raise(rb_eArgError, "w (6th argument) must be NArray");
  if (NA_RANK(rb_w) != 1)
    rb_raise(rb_eArgError, "rank of w (6th argument) must be %d", 1);
  if (NA_SHAPE0(rb_w) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of w must be the same as shape 0 of select");
  if (NA_TYPE(rb_w) != NA_SCOMPLEX)
    rb_w = na_change_type(rb_w, NA_SCOMPLEX);
  w = NA_PTR_TYPE(rb_w, complex*);
  if (!NA_IsNArray(rb_vl))
    rb_raise(rb_eArgError, "vl (7th argument) must be NArray");
  if (NA_RANK(rb_vl) != 2)
    rb_raise(rb_eArgError, "rank of vl (7th argument) must be %d", 2);
  ldvl = NA_SHAPE0(rb_vl);
  mm = NA_SHAPE1(rb_vl);
  if (NA_TYPE(rb_vl) != NA_SCOMPLEX)
    rb_vl = na_change_type(rb_vl, NA_SCOMPLEX);
  vl = NA_PTR_TYPE(rb_vl, complex*);
  if (!NA_IsNArray(rb_vr))
    rb_raise(rb_eArgError, "vr (8th argument) must be NArray");
  if (NA_RANK(rb_vr) != 2)
    rb_raise(rb_eArgError, "rank of vr (8th argument) must be %d", 2);
  ldvr = NA_SHAPE0(rb_vr);
  if (NA_SHAPE1(rb_vr) != mm)
    rb_raise(rb_eRuntimeError, "shape 1 of vr must be the same as shape 1 of vl");
  if (NA_TYPE(rb_vr) != NA_SCOMPLEX)
    rb_vr = na_change_type(rb_vr, NA_SCOMPLEX);
  vr = NA_PTR_TYPE(rb_vr, complex*);
  {
    int shape[1];
    shape[0] = mm;
    rb_ifaill = na_make_object(NA_LINT, 1, shape, cNArray);
  }
  ifaill = NA_PTR_TYPE(rb_ifaill, integer*);
  {
    int shape[1];
    shape[0] = mm;
    rb_ifailr = na_make_object(NA_LINT, 1, shape, cNArray);
  }
  ifailr = NA_PTR_TYPE(rb_ifailr, integer*);
  {
    int shape[1];
    shape[0] = n;
    rb_w_out__ = na_make_object(NA_SCOMPLEX, 1, shape, cNArray);
  }
  w_out__ = NA_PTR_TYPE(rb_w_out__, complex*);
  MEMCPY(w_out__, w, complex, NA_TOTAL(rb_w));
  rb_w = rb_w_out__;
  w = w_out__;
  {
    int shape[2];
    shape[0] = ldvl;
    shape[1] = mm;
    rb_vl_out__ = na_make_object(NA_SCOMPLEX, 2, shape, cNArray);
  }
  vl_out__ = NA_PTR_TYPE(rb_vl_out__, complex*);
  MEMCPY(vl_out__, vl, complex, NA_TOTAL(rb_vl));
  rb_vl = rb_vl_out__;
  vl = vl_out__;
  {
    int shape[2];
    shape[0] = ldvr;
    shape[1] = mm;
    rb_vr_out__ = na_make_object(NA_SCOMPLEX, 2, shape, cNArray);
  }
  vr_out__ = NA_PTR_TYPE(rb_vr_out__, complex*);
  MEMCPY(vr_out__, vr, complex, NA_TOTAL(rb_vr));
  rb_vr = rb_vr_out__;
  vr = vr_out__;
  work = ALLOC_N(complex, (n*n));
  rwork = ALLOC_N(real, (n));

  chsein_(&side, &eigsrc, &initv, select, &n, h, &ldh, w, vl, &ldvl, vr, &ldvr, &mm, &m, work, rwork, ifaill, ifailr, &info);

  free(work);
  free(rwork);
  rb_m = INT2NUM(m);
  rb_info = INT2NUM(info);
  return rb_ary_new3(7, rb_m, rb_ifaill, rb_ifailr, rb_info, rb_w, rb_vl, rb_vr);
}

void
init_lapack_chsein(VALUE mLapack){
  rb_define_module_function(mLapack, "chsein", rb_chsein, -1);
}

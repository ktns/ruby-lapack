#include "rb_lapack.h"

extern VOID zhsein_(char *side, char *eigsrc, char *initv, logical *select, integer *n, doublecomplex *h, integer *ldh, doublecomplex *w, doublecomplex *vl, integer *ldvl, doublecomplex *vr, integer *ldvr, integer *mm, integer *m, doublecomplex *work, doublereal *rwork, integer *ifaill, integer *ifailr, integer *info);

static VALUE sHelp, sUsage;

static VALUE
rblapack_zhsein(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_side;
  char side; 
  VALUE rblapack_eigsrc;
  char eigsrc; 
  VALUE rblapack_initv;
  char initv; 
  VALUE rblapack_select;
  logical *select; 
  VALUE rblapack_h;
  doublecomplex *h; 
  VALUE rblapack_w;
  doublecomplex *w; 
  VALUE rblapack_vl;
  doublecomplex *vl; 
  VALUE rblapack_vr;
  doublecomplex *vr; 
  VALUE rblapack_m;
  integer m; 
  VALUE rblapack_ifaill;
  integer *ifaill; 
  VALUE rblapack_ifailr;
  integer *ifailr; 
  VALUE rblapack_info;
  integer info; 
  VALUE rblapack_w_out__;
  doublecomplex *w_out__;
  VALUE rblapack_vl_out__;
  doublecomplex *vl_out__;
  VALUE rblapack_vr_out__;
  doublecomplex *vr_out__;
  doublecomplex *work;
  doublereal *rwork;

  integer n;
  integer ldh;
  integer ldvl;
  integer mm;
  integer ldvr;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  m, ifaill, ifailr, info, w, vl, vr = NumRu::Lapack.zhsein( side, eigsrc, initv, select, h, w, vl, vr, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE ZHSEIN( SIDE, EIGSRC, INITV, SELECT, N, H, LDH, W, VL, LDVL, VR, LDVR, MM, M, WORK, RWORK, IFAILL, IFAILR, INFO )\n\n*  Purpose\n*  =======\n*\n*  ZHSEIN uses inverse iteration to find specified right and/or left\n*  eigenvectors of a complex upper Hessenberg matrix H.\n*\n*  The right eigenvector x and the left eigenvector y of the matrix H\n*  corresponding to an eigenvalue w are defined by:\n*\n*               H * x = w * x,     y**h * H = w * y**h\n*\n*  where y**h denotes the conjugate transpose of the vector y.\n*\n\n*  Arguments\n*  =========\n*\n*  SIDE    (input) CHARACTER*1\n*          = 'R': compute right eigenvectors only;\n*          = 'L': compute left eigenvectors only;\n*          = 'B': compute both right and left eigenvectors.\n*\n*  EIGSRC  (input) CHARACTER*1\n*          Specifies the source of eigenvalues supplied in W:\n*          = 'Q': the eigenvalues were found using ZHSEQR; thus, if\n*                 H has zero subdiagonal elements, and so is\n*                 block-triangular, then the j-th eigenvalue can be\n*                 assumed to be an eigenvalue of the block containing\n*                 the j-th row/column.  This property allows ZHSEIN to\n*                 perform inverse iteration on just one diagonal block.\n*          = 'N': no assumptions are made on the correspondence\n*                 between eigenvalues and diagonal blocks.  In this\n*                 case, ZHSEIN must always perform inverse iteration\n*                 using the whole matrix H.\n*\n*  INITV   (input) CHARACTER*1\n*          = 'N': no initial vectors are supplied;\n*          = 'U': user-supplied initial vectors are stored in the arrays\n*                 VL and/or VR.\n*\n*  SELECT  (input) LOGICAL array, dimension (N)\n*          Specifies the eigenvectors to be computed. To select the\n*          eigenvector corresponding to the eigenvalue W(j),\n*          SELECT(j) must be set to .TRUE..\n*\n*  N       (input) INTEGER\n*          The order of the matrix H.  N >= 0.\n*\n*  H       (input) COMPLEX*16 array, dimension (LDH,N)\n*          The upper Hessenberg matrix H.\n*\n*  LDH     (input) INTEGER\n*          The leading dimension of the array H.  LDH >= max(1,N).\n*\n*  W       (input/output) COMPLEX*16 array, dimension (N)\n*          On entry, the eigenvalues of H.\n*          On exit, the real parts of W may have been altered since\n*          close eigenvalues are perturbed slightly in searching for\n*          independent eigenvectors.\n*\n*  VL      (input/output) COMPLEX*16 array, dimension (LDVL,MM)\n*          On entry, if INITV = 'U' and SIDE = 'L' or 'B', VL must\n*          contain starting vectors for the inverse iteration for the\n*          left eigenvectors; the starting vector for each eigenvector\n*          must be in the same column in which the eigenvector will be\n*          stored.\n*          On exit, if SIDE = 'L' or 'B', the left eigenvectors\n*          specified by SELECT will be stored consecutively in the\n*          columns of VL, in the same order as their eigenvalues.\n*          If SIDE = 'R', VL is not referenced.\n*\n*  LDVL    (input) INTEGER\n*          The leading dimension of the array VL.\n*          LDVL >= max(1,N) if SIDE = 'L' or 'B'; LDVL >= 1 otherwise.\n*\n*  VR      (input/output) COMPLEX*16 array, dimension (LDVR,MM)\n*          On entry, if INITV = 'U' and SIDE = 'R' or 'B', VR must\n*          contain starting vectors for the inverse iteration for the\n*          right eigenvectors; the starting vector for each eigenvector\n*          must be in the same column in which the eigenvector will be\n*          stored.\n*          On exit, if SIDE = 'R' or 'B', the right eigenvectors\n*          specified by SELECT will be stored consecutively in the\n*          columns of VR, in the same order as their eigenvalues.\n*          If SIDE = 'L', VR is not referenced.\n*\n*  LDVR    (input) INTEGER\n*          The leading dimension of the array VR.\n*          LDVR >= max(1,N) if SIDE = 'R' or 'B'; LDVR >= 1 otherwise.\n*\n*  MM      (input) INTEGER\n*          The number of columns in the arrays VL and/or VR. MM >= M.\n*\n*  M       (output) INTEGER\n*          The number of columns in the arrays VL and/or VR required to\n*          store the eigenvectors (= the number of .TRUE. elements in\n*          SELECT).\n*\n*  WORK    (workspace) COMPLEX*16 array, dimension (N*N)\n*\n*  RWORK   (workspace) DOUBLE PRECISION array, dimension (N)\n*\n*  IFAILL  (output) INTEGER array, dimension (MM)\n*          If SIDE = 'L' or 'B', IFAILL(i) = j > 0 if the left\n*          eigenvector in the i-th column of VL (corresponding to the\n*          eigenvalue w(j)) failed to converge; IFAILL(i) = 0 if the\n*          eigenvector converged satisfactorily.\n*          If SIDE = 'R', IFAILL is not referenced.\n*\n*  IFAILR  (output) INTEGER array, dimension (MM)\n*          If SIDE = 'R' or 'B', IFAILR(i) = j > 0 if the right\n*          eigenvector in the i-th column of VR (corresponding to the\n*          eigenvalue w(j)) failed to converge; IFAILR(i) = 0 if the\n*          eigenvector converged satisfactorily.\n*          If SIDE = 'L', IFAILR is not referenced.\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit\n*          < 0:  if INFO = -i, the i-th argument had an illegal value\n*          > 0:  if INFO = i, i is the number of eigenvectors which\n*                failed to converge; see IFAILL and IFAILR for further\n*                details.\n*\n\n*  Further Details\n*  ===============\n*\n*  Each eigenvector is normalized so that the element of largest\n*  magnitude has magnitude 1; here the magnitude of a complex number\n*  (x,y) is taken to be |x|+|y|.\n*\n*  =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  m, ifaill, ifailr, info, w, vl, vr = NumRu::Lapack.zhsein( side, eigsrc, initv, select, h, w, vl, vr, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 8)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 8)", argc);
  rblapack_side = argv[0];
  rblapack_eigsrc = argv[1];
  rblapack_initv = argv[2];
  rblapack_select = argv[3];
  rblapack_h = argv[4];
  rblapack_w = argv[5];
  rblapack_vl = argv[6];
  rblapack_vr = argv[7];
  if (rb_options != Qnil) {
  }

  if (!NA_IsNArray(rblapack_vl))
    rb_raise(rb_eArgError, "vl (7th argument) must be NArray");
  if (NA_RANK(rblapack_vl) != 2)
    rb_raise(rb_eArgError, "rank of vl (7th argument) must be %d", 2);
  mm = NA_SHAPE1(rblapack_vl);
  ldvl = NA_SHAPE0(rblapack_vl);
  if (NA_TYPE(rblapack_vl) != NA_DCOMPLEX)
    rblapack_vl = na_change_type(rblapack_vl, NA_DCOMPLEX);
  vl = NA_PTR_TYPE(rblapack_vl, doublecomplex*);
  if (!NA_IsNArray(rblapack_w))
    rb_raise(rb_eArgError, "w (6th argument) must be NArray");
  if (NA_RANK(rblapack_w) != 1)
    rb_raise(rb_eArgError, "rank of w (6th argument) must be %d", 1);
  n = NA_SHAPE0(rblapack_w);
  if (NA_TYPE(rblapack_w) != NA_DCOMPLEX)
    rblapack_w = na_change_type(rblapack_w, NA_DCOMPLEX);
  w = NA_PTR_TYPE(rblapack_w, doublecomplex*);
  side = StringValueCStr(rblapack_side)[0];
  eigsrc = StringValueCStr(rblapack_eigsrc)[0];
  if (!NA_IsNArray(rblapack_vr))
    rb_raise(rb_eArgError, "vr (8th argument) must be NArray");
  if (NA_RANK(rblapack_vr) != 2)
    rb_raise(rb_eArgError, "rank of vr (8th argument) must be %d", 2);
  if (NA_SHAPE1(rblapack_vr) != mm)
    rb_raise(rb_eRuntimeError, "shape 1 of vr must be the same as shape 1 of vl");
  ldvr = NA_SHAPE0(rblapack_vr);
  if (NA_TYPE(rblapack_vr) != NA_DCOMPLEX)
    rblapack_vr = na_change_type(rblapack_vr, NA_DCOMPLEX);
  vr = NA_PTR_TYPE(rblapack_vr, doublecomplex*);
  initv = StringValueCStr(rblapack_initv)[0];
  if (!NA_IsNArray(rblapack_h))
    rb_raise(rb_eArgError, "h (5th argument) must be NArray");
  if (NA_RANK(rblapack_h) != 2)
    rb_raise(rb_eArgError, "rank of h (5th argument) must be %d", 2);
  if (NA_SHAPE1(rblapack_h) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of h must be the same as shape 0 of w");
  ldh = NA_SHAPE0(rblapack_h);
  if (NA_TYPE(rblapack_h) != NA_DCOMPLEX)
    rblapack_h = na_change_type(rblapack_h, NA_DCOMPLEX);
  h = NA_PTR_TYPE(rblapack_h, doublecomplex*);
  if (!NA_IsNArray(rblapack_select))
    rb_raise(rb_eArgError, "select (4th argument) must be NArray");
  if (NA_RANK(rblapack_select) != 1)
    rb_raise(rb_eArgError, "rank of select (4th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_select) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of select must be the same as shape 0 of w");
  if (NA_TYPE(rblapack_select) != NA_LINT)
    rblapack_select = na_change_type(rblapack_select, NA_LINT);
  select = NA_PTR_TYPE(rblapack_select, logical*);
  {
    int shape[1];
    shape[0] = mm;
    rblapack_ifaill = na_make_object(NA_LINT, 1, shape, cNArray);
  }
  ifaill = NA_PTR_TYPE(rblapack_ifaill, integer*);
  {
    int shape[1];
    shape[0] = mm;
    rblapack_ifailr = na_make_object(NA_LINT, 1, shape, cNArray);
  }
  ifailr = NA_PTR_TYPE(rblapack_ifailr, integer*);
  {
    int shape[1];
    shape[0] = n;
    rblapack_w_out__ = na_make_object(NA_DCOMPLEX, 1, shape, cNArray);
  }
  w_out__ = NA_PTR_TYPE(rblapack_w_out__, doublecomplex*);
  MEMCPY(w_out__, w, doublecomplex, NA_TOTAL(rblapack_w));
  rblapack_w = rblapack_w_out__;
  w = w_out__;
  {
    int shape[2];
    shape[0] = ldvl;
    shape[1] = mm;
    rblapack_vl_out__ = na_make_object(NA_DCOMPLEX, 2, shape, cNArray);
  }
  vl_out__ = NA_PTR_TYPE(rblapack_vl_out__, doublecomplex*);
  MEMCPY(vl_out__, vl, doublecomplex, NA_TOTAL(rblapack_vl));
  rblapack_vl = rblapack_vl_out__;
  vl = vl_out__;
  {
    int shape[2];
    shape[0] = ldvr;
    shape[1] = mm;
    rblapack_vr_out__ = na_make_object(NA_DCOMPLEX, 2, shape, cNArray);
  }
  vr_out__ = NA_PTR_TYPE(rblapack_vr_out__, doublecomplex*);
  MEMCPY(vr_out__, vr, doublecomplex, NA_TOTAL(rblapack_vr));
  rblapack_vr = rblapack_vr_out__;
  vr = vr_out__;
  work = ALLOC_N(doublecomplex, (n*n));
  rwork = ALLOC_N(doublereal, (n));

  zhsein_(&side, &eigsrc, &initv, select, &n, h, &ldh, w, vl, &ldvl, vr, &ldvr, &mm, &m, work, rwork, ifaill, ifailr, &info);

  free(work);
  free(rwork);
  rblapack_m = INT2NUM(m);
  rblapack_info = INT2NUM(info);
  return rb_ary_new3(7, rblapack_m, rblapack_ifaill, rblapack_ifailr, rblapack_info, rblapack_w, rblapack_vl, rblapack_vr);
}

void
init_lapack_zhsein(VALUE mLapack){
  rb_define_module_function(mLapack, "zhsein", rblapack_zhsein, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}

#include "rb_lapack.h"

extern VOID dhsein_(char *side, char *eigsrc, char *initv, logical *select, integer *n, doublereal *h, integer *ldh, doublereal *wr, doublereal *wi, doublereal *vl, integer *ldvl, doublereal *vr, integer *ldvr, integer *mm, integer *m, doublereal *work, integer *ifaill, integer *ifailr, integer *info);

static VALUE sHelp, sUsage;

static VALUE
rblapack_dhsein(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_side;
  char side; 
  VALUE rblapack_eigsrc;
  char eigsrc; 
  VALUE rblapack_initv;
  char initv; 
  VALUE rblapack_select;
  logical *select; 
  VALUE rblapack_h;
  doublereal *h; 
  VALUE rblapack_wr;
  doublereal *wr; 
  VALUE rblapack_wi;
  doublereal *wi; 
  VALUE rblapack_vl;
  doublereal *vl; 
  VALUE rblapack_vr;
  doublereal *vr; 
  VALUE rblapack_m;
  integer m; 
  VALUE rblapack_ifaill;
  integer *ifaill; 
  VALUE rblapack_ifailr;
  integer *ifailr; 
  VALUE rblapack_info;
  integer info; 
  VALUE rblapack_select_out__;
  logical *select_out__;
  VALUE rblapack_wr_out__;
  doublereal *wr_out__;
  VALUE rblapack_vl_out__;
  doublereal *vl_out__;
  VALUE rblapack_vr_out__;
  doublereal *vr_out__;
  doublereal *work;

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
      printf("%s\n", "USAGE:\n  m, ifaill, ifailr, info, select, wr, vl, vr = NumRu::Lapack.dhsein( side, eigsrc, initv, select, h, wr, wi, vl, vr, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE DHSEIN( SIDE, EIGSRC, INITV, SELECT, N, H, LDH, WR, WI, VL, LDVL, VR, LDVR, MM, M, WORK, IFAILL, IFAILR, INFO )\n\n*  Purpose\n*  =======\n*\n*  DHSEIN uses inverse iteration to find specified right and/or left\n*  eigenvectors of a real upper Hessenberg matrix H.\n*\n*  The right eigenvector x and the left eigenvector y of the matrix H\n*  corresponding to an eigenvalue w are defined by:\n*\n*               H * x = w * x,     y**h * H = w * y**h\n*\n*  where y**h denotes the conjugate transpose of the vector y.\n*\n\n*  Arguments\n*  =========\n*\n*  SIDE    (input) CHARACTER*1\n*          = 'R': compute right eigenvectors only;\n*          = 'L': compute left eigenvectors only;\n*          = 'B': compute both right and left eigenvectors.\n*\n*  EIGSRC  (input) CHARACTER*1\n*          Specifies the source of eigenvalues supplied in (WR,WI):\n*          = 'Q': the eigenvalues were found using DHSEQR; thus, if\n*                 H has zero subdiagonal elements, and so is\n*                 block-triangular, then the j-th eigenvalue can be\n*                 assumed to be an eigenvalue of the block containing\n*                 the j-th row/column.  This property allows DHSEIN to\n*                 perform inverse iteration on just one diagonal block.\n*          = 'N': no assumptions are made on the correspondence\n*                 between eigenvalues and diagonal blocks.  In this\n*                 case, DHSEIN must always perform inverse iteration\n*                 using the whole matrix H.\n*\n*  INITV   (input) CHARACTER*1\n*          = 'N': no initial vectors are supplied;\n*          = 'U': user-supplied initial vectors are stored in the arrays\n*                 VL and/or VR.\n*\n*  SELECT  (input/output) LOGICAL array, dimension (N)\n*          Specifies the eigenvectors to be computed. To select the\n*          real eigenvector corresponding to a real eigenvalue WR(j),\n*          SELECT(j) must be set to .TRUE.. To select the complex\n*          eigenvector corresponding to a complex eigenvalue\n*          (WR(j),WI(j)), with complex conjugate (WR(j+1),WI(j+1)),\n*          either SELECT(j) or SELECT(j+1) or both must be set to\n*          .TRUE.; then on exit SELECT(j) is .TRUE. and SELECT(j+1) is\n*          .FALSE..\n*\n*  N       (input) INTEGER\n*          The order of the matrix H.  N >= 0.\n*\n*  H       (input) DOUBLE PRECISION array, dimension (LDH,N)\n*          The upper Hessenberg matrix H.\n*\n*  LDH     (input) INTEGER\n*          The leading dimension of the array H.  LDH >= max(1,N).\n*\n*  WR      (input/output) DOUBLE PRECISION array, dimension (N)\n*  WI      (input) DOUBLE PRECISION array, dimension (N)\n*          On entry, the real and imaginary parts of the eigenvalues of\n*          H; a complex conjugate pair of eigenvalues must be stored in\n*          consecutive elements of WR and WI.\n*          On exit, WR may have been altered since close eigenvalues\n*          are perturbed slightly in searching for independent\n*          eigenvectors.\n*\n*  VL      (input/output) DOUBLE PRECISION array, dimension (LDVL,MM)\n*          On entry, if INITV = 'U' and SIDE = 'L' or 'B', VL must\n*          contain starting vectors for the inverse iteration for the\n*          left eigenvectors; the starting vector for each eigenvector\n*          must be in the same column(s) in which the eigenvector will\n*          be stored.\n*          On exit, if SIDE = 'L' or 'B', the left eigenvectors\n*          specified by SELECT will be stored consecutively in the\n*          columns of VL, in the same order as their eigenvalues. A\n*          complex eigenvector corresponding to a complex eigenvalue is\n*          stored in two consecutive columns, the first holding the real\n*          part and the second the imaginary part.\n*          If SIDE = 'R', VL is not referenced.\n*\n*  LDVL    (input) INTEGER\n*          The leading dimension of the array VL.\n*          LDVL >= max(1,N) if SIDE = 'L' or 'B'; LDVL >= 1 otherwise.\n*\n*  VR      (input/output) DOUBLE PRECISION array, dimension (LDVR,MM)\n*          On entry, if INITV = 'U' and SIDE = 'R' or 'B', VR must\n*          contain starting vectors for the inverse iteration for the\n*          right eigenvectors; the starting vector for each eigenvector\n*          must be in the same column(s) in which the eigenvector will\n*          be stored.\n*          On exit, if SIDE = 'R' or 'B', the right eigenvectors\n*          specified by SELECT will be stored consecutively in the\n*          columns of VR, in the same order as their eigenvalues. A\n*          complex eigenvector corresponding to a complex eigenvalue is\n*          stored in two consecutive columns, the first holding the real\n*          part and the second the imaginary part.\n*          If SIDE = 'L', VR is not referenced.\n*\n*  LDVR    (input) INTEGER\n*          The leading dimension of the array VR.\n*          LDVR >= max(1,N) if SIDE = 'R' or 'B'; LDVR >= 1 otherwise.\n*\n*  MM      (input) INTEGER\n*          The number of columns in the arrays VL and/or VR. MM >= M.\n*\n*  M       (output) INTEGER\n*          The number of columns in the arrays VL and/or VR required to\n*          store the eigenvectors; each selected real eigenvector\n*          occupies one column and each selected complex eigenvector\n*          occupies two columns.\n*\n*  WORK    (workspace) DOUBLE PRECISION array, dimension ((N+2)*N)\n*\n*  IFAILL  (output) INTEGER array, dimension (MM)\n*          If SIDE = 'L' or 'B', IFAILL(i) = j > 0 if the left\n*          eigenvector in the i-th column of VL (corresponding to the\n*          eigenvalue w(j)) failed to converge; IFAILL(i) = 0 if the\n*          eigenvector converged satisfactorily. If the i-th and (i+1)th\n*          columns of VL hold a complex eigenvector, then IFAILL(i) and\n*          IFAILL(i+1) are set to the same value.\n*          If SIDE = 'R', IFAILL is not referenced.\n*\n*  IFAILR  (output) INTEGER array, dimension (MM)\n*          If SIDE = 'R' or 'B', IFAILR(i) = j > 0 if the right\n*          eigenvector in the i-th column of VR (corresponding to the\n*          eigenvalue w(j)) failed to converge; IFAILR(i) = 0 if the\n*          eigenvector converged satisfactorily. If the i-th and (i+1)th\n*          columns of VR hold a complex eigenvector, then IFAILR(i) and\n*          IFAILR(i+1) are set to the same value.\n*          If SIDE = 'L', IFAILR is not referenced.\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit\n*          < 0:  if INFO = -i, the i-th argument had an illegal value\n*          > 0:  if INFO = i, i is the number of eigenvectors which\n*                failed to converge; see IFAILL and IFAILR for further\n*                details.\n*\n\n*  Further Details\n*  ===============\n*\n*  Each eigenvector is normalized so that the element of largest\n*  magnitude has magnitude 1; here the magnitude of a complex number\n*  (x,y) is taken to be |x|+|y|.\n*\n*  =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  m, ifaill, ifailr, info, select, wr, vl, vr = NumRu::Lapack.dhsein( side, eigsrc, initv, select, h, wr, wi, vl, vr, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 9)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 9)", argc);
  rblapack_side = argv[0];
  rblapack_eigsrc = argv[1];
  rblapack_initv = argv[2];
  rblapack_select = argv[3];
  rblapack_h = argv[4];
  rblapack_wr = argv[5];
  rblapack_wi = argv[6];
  rblapack_vl = argv[7];
  rblapack_vr = argv[8];
  if (rb_options != Qnil) {
  }

  if (!NA_IsNArray(rblapack_vl))
    rb_raise(rb_eArgError, "vl (8th argument) must be NArray");
  if (NA_RANK(rblapack_vl) != 2)
    rb_raise(rb_eArgError, "rank of vl (8th argument) must be %d", 2);
  mm = NA_SHAPE1(rblapack_vl);
  ldvl = NA_SHAPE0(rblapack_vl);
  if (NA_TYPE(rblapack_vl) != NA_DFLOAT)
    rblapack_vl = na_change_type(rblapack_vl, NA_DFLOAT);
  vl = NA_PTR_TYPE(rblapack_vl, doublereal*);
  side = StringValueCStr(rblapack_side)[0];
  eigsrc = StringValueCStr(rblapack_eigsrc)[0];
  if (!NA_IsNArray(rblapack_wr))
    rb_raise(rb_eArgError, "wr (6th argument) must be NArray");
  if (NA_RANK(rblapack_wr) != 1)
    rb_raise(rb_eArgError, "rank of wr (6th argument) must be %d", 1);
  n = NA_SHAPE0(rblapack_wr);
  if (NA_TYPE(rblapack_wr) != NA_DFLOAT)
    rblapack_wr = na_change_type(rblapack_wr, NA_DFLOAT);
  wr = NA_PTR_TYPE(rblapack_wr, doublereal*);
  if (!NA_IsNArray(rblapack_vr))
    rb_raise(rb_eArgError, "vr (9th argument) must be NArray");
  if (NA_RANK(rblapack_vr) != 2)
    rb_raise(rb_eArgError, "rank of vr (9th argument) must be %d", 2);
  if (NA_SHAPE1(rblapack_vr) != mm)
    rb_raise(rb_eRuntimeError, "shape 1 of vr must be the same as shape 1 of vl");
  ldvr = NA_SHAPE0(rblapack_vr);
  if (NA_TYPE(rblapack_vr) != NA_DFLOAT)
    rblapack_vr = na_change_type(rblapack_vr, NA_DFLOAT);
  vr = NA_PTR_TYPE(rblapack_vr, doublereal*);
  initv = StringValueCStr(rblapack_initv)[0];
  if (!NA_IsNArray(rblapack_wi))
    rb_raise(rb_eArgError, "wi (7th argument) must be NArray");
  if (NA_RANK(rblapack_wi) != 1)
    rb_raise(rb_eArgError, "rank of wi (7th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_wi) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of wi must be the same as shape 0 of wr");
  if (NA_TYPE(rblapack_wi) != NA_DFLOAT)
    rblapack_wi = na_change_type(rblapack_wi, NA_DFLOAT);
  wi = NA_PTR_TYPE(rblapack_wi, doublereal*);
  if (!NA_IsNArray(rblapack_h))
    rb_raise(rb_eArgError, "h (5th argument) must be NArray");
  if (NA_RANK(rblapack_h) != 2)
    rb_raise(rb_eArgError, "rank of h (5th argument) must be %d", 2);
  if (NA_SHAPE1(rblapack_h) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of h must be the same as shape 0 of wr");
  ldh = NA_SHAPE0(rblapack_h);
  if (NA_TYPE(rblapack_h) != NA_DFLOAT)
    rblapack_h = na_change_type(rblapack_h, NA_DFLOAT);
  h = NA_PTR_TYPE(rblapack_h, doublereal*);
  if (!NA_IsNArray(rblapack_select))
    rb_raise(rb_eArgError, "select (4th argument) must be NArray");
  if (NA_RANK(rblapack_select) != 1)
    rb_raise(rb_eArgError, "rank of select (4th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_select) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of select must be the same as shape 0 of wr");
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
    rblapack_select_out__ = na_make_object(NA_LINT, 1, shape, cNArray);
  }
  select_out__ = NA_PTR_TYPE(rblapack_select_out__, logical*);
  MEMCPY(select_out__, select, logical, NA_TOTAL(rblapack_select));
  rblapack_select = rblapack_select_out__;
  select = select_out__;
  {
    int shape[1];
    shape[0] = n;
    rblapack_wr_out__ = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  wr_out__ = NA_PTR_TYPE(rblapack_wr_out__, doublereal*);
  MEMCPY(wr_out__, wr, doublereal, NA_TOTAL(rblapack_wr));
  rblapack_wr = rblapack_wr_out__;
  wr = wr_out__;
  {
    int shape[2];
    shape[0] = ldvl;
    shape[1] = mm;
    rblapack_vl_out__ = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  vl_out__ = NA_PTR_TYPE(rblapack_vl_out__, doublereal*);
  MEMCPY(vl_out__, vl, doublereal, NA_TOTAL(rblapack_vl));
  rblapack_vl = rblapack_vl_out__;
  vl = vl_out__;
  {
    int shape[2];
    shape[0] = ldvr;
    shape[1] = mm;
    rblapack_vr_out__ = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  vr_out__ = NA_PTR_TYPE(rblapack_vr_out__, doublereal*);
  MEMCPY(vr_out__, vr, doublereal, NA_TOTAL(rblapack_vr));
  rblapack_vr = rblapack_vr_out__;
  vr = vr_out__;
  work = ALLOC_N(doublereal, ((n+2)*n));

  dhsein_(&side, &eigsrc, &initv, select, &n, h, &ldh, wr, wi, vl, &ldvl, vr, &ldvr, &mm, &m, work, ifaill, ifailr, &info);

  free(work);
  rblapack_m = INT2NUM(m);
  rblapack_info = INT2NUM(info);
  return rb_ary_new3(8, rblapack_m, rblapack_ifaill, rblapack_ifailr, rblapack_info, rblapack_select, rblapack_wr, rblapack_vl, rblapack_vr);
}

void
init_lapack_dhsein(VALUE mLapack){
  rb_define_module_function(mLapack, "dhsein", rblapack_dhsein, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}

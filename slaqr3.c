#include "rb_lapack.h"

extern VOID slaqr3_(logical *wantt, logical *wantz, integer *n, integer *ktop, integer *kbot, integer *nw, real *h, integer *ldh, integer *iloz, integer *ihiz, real *z, integer *ldz, integer *ns, integer *nd, real *sr, real *si, real *v, integer *ldv, integer *nh, real *t, integer *ldt, integer *nv, real *wv, integer *ldwv, real *work, integer *lwork);

static VALUE sHelp, sUsage;

static VALUE
rblapack_slaqr3(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_wantt;
  logical wantt; 
  VALUE rblapack_wantz;
  logical wantz; 
  VALUE rblapack_ktop;
  integer ktop; 
  VALUE rblapack_kbot;
  integer kbot; 
  VALUE rblapack_nw;
  integer nw; 
  VALUE rblapack_h;
  real *h; 
  VALUE rblapack_iloz;
  integer iloz; 
  VALUE rblapack_ihiz;
  integer ihiz; 
  VALUE rblapack_z;
  real *z; 
  VALUE rblapack_nh;
  integer nh; 
  VALUE rblapack_nv;
  integer nv; 
  VALUE rblapack_lwork;
  integer lwork; 
  VALUE rblapack_ns;
  integer ns; 
  VALUE rblapack_nd;
  integer nd; 
  VALUE rblapack_sr;
  real *sr; 
  VALUE rblapack_si;
  real *si; 
  VALUE rblapack_h_out__;
  real *h_out__;
  VALUE rblapack_z_out__;
  real *z_out__;
  real *v;
  real *t;
  real *wv;
  real *work;

  integer ldh;
  integer n;
  integer ldz;
  integer ldv;
  integer ldt;
  integer ldwv;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  ns, nd, sr, si, h, z = NumRu::Lapack.slaqr3( wantt, wantz, ktop, kbot, nw, h, iloz, ihiz, z, nh, nv, lwork, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE SLAQR3( WANTT, WANTZ, N, KTOP, KBOT, NW, H, LDH, ILOZ, IHIZ, Z, LDZ, NS, ND, SR, SI, V, LDV, NH, T, LDT, NV, WV, LDWV, WORK, LWORK )\n\n*     Aggressive early deflation:\n*\n*     This subroutine accepts as input an upper Hessenberg matrix\n*     H and performs an orthogonal similarity transformation\n*     designed to detect and deflate fully converged eigenvalues from\n*     a trailing principal submatrix.  On output H has been over-\n*     written by a new Hessenberg matrix that is a perturbation of\n*     an orthogonal similarity transformation of H.  It is to be\n*     hoped that the final version of H has many zero subdiagonal\n*     entries.\n*\n*     ******************************************************************\n\n*     WANTT   (input) LOGICAL\n*          If .TRUE., then the Hessenberg matrix H is fully updated\n*          so that the quasi-triangular Schur factor may be\n*          computed (in cooperation with the calling subroutine).\n*          If .FALSE., then only enough of H is updated to preserve\n*          the eigenvalues.\n*\n*     WANTZ   (input) LOGICAL\n*          If .TRUE., then the orthogonal matrix Z is updated so\n*          so that the orthogonal Schur factor may be computed\n*          (in cooperation with the calling subroutine).\n*          If .FALSE., then Z is not referenced.\n*\n*     N       (input) INTEGER\n*          The order of the matrix H and (if WANTZ is .TRUE.) the\n*          order of the orthogonal matrix Z.\n*\n*     KTOP    (input) INTEGER\n*          It is assumed that either KTOP = 1 or H(KTOP,KTOP-1)=0.\n*          KBOT and KTOP together determine an isolated block\n*          along the diagonal of the Hessenberg matrix.\n*\n*     KBOT    (input) INTEGER\n*          It is assumed without a check that either\n*          KBOT = N or H(KBOT+1,KBOT)=0.  KBOT and KTOP together\n*          determine an isolated block along the diagonal of the\n*          Hessenberg matrix.\n*\n*     NW      (input) INTEGER\n*          Deflation window size.  1 .LE. NW .LE. (KBOT-KTOP+1).\n*\n*     H       (input/output) REAL array, dimension (LDH,N)\n*          On input the initial N-by-N section of H stores the\n*          Hessenberg matrix undergoing aggressive early deflation.\n*          On output H has been transformed by an orthogonal\n*          similarity transformation, perturbed, and the returned\n*          to Hessenberg form that (it is to be hoped) has some\n*          zero subdiagonal entries.\n*\n*     LDH     (input) integer\n*          Leading dimension of H just as declared in the calling\n*          subroutine.  N .LE. LDH\n*\n*     ILOZ    (input) INTEGER\n*     IHIZ    (input) INTEGER\n*          Specify the rows of Z to which transformations must be\n*          applied if WANTZ is .TRUE.. 1 .LE. ILOZ .LE. IHIZ .LE. N.\n*\n*     Z       (input/output) REAL array, dimension (LDZ,N)\n*          IF WANTZ is .TRUE., then on output, the orthogonal\n*          similarity transformation mentioned above has been\n*          accumulated into Z(ILOZ:IHIZ,ILO:IHI) from the right.\n*          If WANTZ is .FALSE., then Z is unreferenced.\n*\n*     LDZ     (input) integer\n*          The leading dimension of Z just as declared in the\n*          calling subroutine.  1 .LE. LDZ.\n*\n*     NS      (output) integer\n*          The number of unconverged (ie approximate) eigenvalues\n*          returned in SR and SI that may be used as shifts by the\n*          calling subroutine.\n*\n*     ND      (output) integer\n*          The number of converged eigenvalues uncovered by this\n*          subroutine.\n*\n*     SR      (output) REAL array, dimension KBOT\n*     SI      (output) REAL array, dimension KBOT\n*          On output, the real and imaginary parts of approximate\n*          eigenvalues that may be used for shifts are stored in\n*          SR(KBOT-ND-NS+1) through SR(KBOT-ND) and\n*          SI(KBOT-ND-NS+1) through SI(KBOT-ND), respectively.\n*          The real and imaginary parts of converged eigenvalues\n*          are stored in SR(KBOT-ND+1) through SR(KBOT) and\n*          SI(KBOT-ND+1) through SI(KBOT), respectively.\n*\n*     V       (workspace) REAL array, dimension (LDV,NW)\n*          An NW-by-NW work array.\n*\n*     LDV     (input) integer scalar\n*          The leading dimension of V just as declared in the\n*          calling subroutine.  NW .LE. LDV\n*\n*     NH      (input) integer scalar\n*          The number of columns of T.  NH.GE.NW.\n*\n*     T       (workspace) REAL array, dimension (LDT,NW)\n*\n*     LDT     (input) integer\n*          The leading dimension of T just as declared in the\n*          calling subroutine.  NW .LE. LDT\n*\n*     NV      (input) integer\n*          The number of rows of work array WV available for\n*          workspace.  NV.GE.NW.\n*\n*     WV      (workspace) REAL array, dimension (LDWV,NW)\n*\n*     LDWV    (input) integer\n*          The leading dimension of W just as declared in the\n*          calling subroutine.  NW .LE. LDV\n*\n*     WORK    (workspace) REAL array, dimension LWORK.\n*          On exit, WORK(1) is set to an estimate of the optimal value\n*          of LWORK for the given values of N, NW, KTOP and KBOT.\n*\n*     LWORK   (input) integer\n*          The dimension of the work array WORK.  LWORK = 2*NW\n*          suffices, but greater efficiency may result from larger\n*          values of LWORK.\n*\n*          If LWORK = -1, then a workspace query is assumed; SLAQR3\n*          only estimates the optimal workspace size for the given\n*          values of N, NW, KTOP and KBOT.  The estimate is returned\n*          in WORK(1).  No error message related to LWORK is issued\n*          by XERBLA.  Neither H nor Z are accessed.\n*\n\n*     ================================================================\n*     Based on contributions by\n*        Karen Braman and Ralph Byers, Department of Mathematics,\n*        University of Kansas, USA\n*\n*     ================================================================\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  ns, nd, sr, si, h, z = NumRu::Lapack.slaqr3( wantt, wantz, ktop, kbot, nw, h, iloz, ihiz, z, nh, nv, lwork, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 12)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 12)", argc);
  rblapack_wantt = argv[0];
  rblapack_wantz = argv[1];
  rblapack_ktop = argv[2];
  rblapack_kbot = argv[3];
  rblapack_nw = argv[4];
  rblapack_h = argv[5];
  rblapack_iloz = argv[6];
  rblapack_ihiz = argv[7];
  rblapack_z = argv[8];
  rblapack_nh = argv[9];
  rblapack_nv = argv[10];
  rblapack_lwork = argv[11];
  if (rb_options != Qnil) {
  }

  ktop = NUM2INT(rblapack_ktop);
  wantz = (rblapack_wantz == Qtrue);
  nh = NUM2INT(rblapack_nh);
  if (!NA_IsNArray(rblapack_z))
    rb_raise(rb_eArgError, "z (9th argument) must be NArray");
  if (NA_RANK(rblapack_z) != 2)
    rb_raise(rb_eArgError, "rank of z (9th argument) must be %d", 2);
  n = NA_SHAPE1(rblapack_z);
  ldz = NA_SHAPE0(rblapack_z);
  if (NA_TYPE(rblapack_z) != NA_SFLOAT)
    rblapack_z = na_change_type(rblapack_z, NA_SFLOAT);
  z = NA_PTR_TYPE(rblapack_z, real*);
  kbot = NUM2INT(rblapack_kbot);
  lwork = NUM2INT(rblapack_lwork);
  iloz = NUM2INT(rblapack_iloz);
  nv = NUM2INT(rblapack_nv);
  wantt = (rblapack_wantt == Qtrue);
  ihiz = NUM2INT(rblapack_ihiz);
  if (!NA_IsNArray(rblapack_h))
    rb_raise(rb_eArgError, "h (6th argument) must be NArray");
  if (NA_RANK(rblapack_h) != 2)
    rb_raise(rb_eArgError, "rank of h (6th argument) must be %d", 2);
  if (NA_SHAPE1(rblapack_h) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of h must be the same as shape 1 of z");
  ldh = NA_SHAPE0(rblapack_h);
  if (NA_TYPE(rblapack_h) != NA_SFLOAT)
    rblapack_h = na_change_type(rblapack_h, NA_SFLOAT);
  h = NA_PTR_TYPE(rblapack_h, real*);
  nw = NUM2INT(rblapack_nw);
  ldv = nw;
  ldwv = nw;
  ldt = nw;
  {
    int shape[1];
    shape[0] = MAX(1,kbot);
    rblapack_sr = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  sr = NA_PTR_TYPE(rblapack_sr, real*);
  {
    int shape[1];
    shape[0] = MAX(1,kbot);
    rblapack_si = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  si = NA_PTR_TYPE(rblapack_si, real*);
  {
    int shape[2];
    shape[0] = ldh;
    shape[1] = n;
    rblapack_h_out__ = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  h_out__ = NA_PTR_TYPE(rblapack_h_out__, real*);
  MEMCPY(h_out__, h, real, NA_TOTAL(rblapack_h));
  rblapack_h = rblapack_h_out__;
  h = h_out__;
  {
    int shape[2];
    shape[0] = ldz;
    shape[1] = n;
    rblapack_z_out__ = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  z_out__ = NA_PTR_TYPE(rblapack_z_out__, real*);
  MEMCPY(z_out__, z, real, NA_TOTAL(rblapack_z));
  rblapack_z = rblapack_z_out__;
  z = z_out__;
  v = ALLOC_N(real, (ldv)*(MAX(1,nw)));
  t = ALLOC_N(real, (ldt)*(MAX(1,nw)));
  wv = ALLOC_N(real, (ldwv)*(MAX(1,nw)));
  work = ALLOC_N(real, (MAX(1,lwork)));

  slaqr3_(&wantt, &wantz, &n, &ktop, &kbot, &nw, h, &ldh, &iloz, &ihiz, z, &ldz, &ns, &nd, sr, si, v, &ldv, &nh, t, &ldt, &nv, wv, &ldwv, work, &lwork);

  free(v);
  free(t);
  free(wv);
  free(work);
  rblapack_ns = INT2NUM(ns);
  rblapack_nd = INT2NUM(nd);
  return rb_ary_new3(6, rblapack_ns, rblapack_nd, rblapack_sr, rblapack_si, rblapack_h, rblapack_z);
}

void
init_lapack_slaqr3(VALUE mLapack){
  rb_define_module_function(mLapack, "slaqr3", rblapack_slaqr3, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}

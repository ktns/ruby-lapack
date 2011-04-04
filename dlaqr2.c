#include "rb_lapack.h"

extern VOID dlaqr2_(logical *wantt, logical *wantz, integer *n, integer *ktop, integer *kbot, integer *nw, doublereal *h, integer *ldh, integer *iloz, integer *ihiz, doublereal *z, integer *ldz, integer *ns, integer *nd, doublereal *sr, doublereal *si, doublereal *v, integer *ldv, integer *nh, doublereal *t, integer *ldt, integer *nv, doublereal *wv, integer *ldwv, doublereal *work, integer *lwork);

static VALUE
rb_dlaqr2(int argc, VALUE *argv, VALUE self){
  VALUE rb_wantt;
  logical wantt; 
  VALUE rb_wantz;
  logical wantz; 
  VALUE rb_ktop;
  integer ktop; 
  VALUE rb_kbot;
  integer kbot; 
  VALUE rb_nw;
  integer nw; 
  VALUE rb_h;
  doublereal *h; 
  VALUE rb_iloz;
  integer iloz; 
  VALUE rb_ihiz;
  integer ihiz; 
  VALUE rb_z;
  doublereal *z; 
  VALUE rb_nh;
  integer nh; 
  VALUE rb_nv;
  integer nv; 
  VALUE rb_lwork;
  integer lwork; 
  VALUE rb_ns;
  integer ns; 
  VALUE rb_nd;
  integer nd; 
  VALUE rb_sr;
  doublereal *sr; 
  VALUE rb_si;
  doublereal *si; 
  VALUE rb_h_out__;
  doublereal *h_out__;
  VALUE rb_z_out__;
  doublereal *z_out__;
  doublereal *v;
  doublereal *t;
  doublereal *wv;
  doublereal *work;

  integer ldh;
  integer n;
  integer ldz;
  integer ldv;
  integer ldt;
  integer ldwv;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  ns, nd, sr, si, h, z = NumRu::Lapack.dlaqr2( wantt, wantz, ktop, kbot, nw, h, iloz, ihiz, z, nh, nv, lwork)\n    or\n  NumRu::Lapack.dlaqr2  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE DLAQR2( WANTT, WANTZ, N, KTOP, KBOT, NW, H, LDH, ILOZ, IHIZ, Z, LDZ, NS, ND, SR, SI, V, LDV, NH, T, LDT, NV, WV, LDWV, WORK, LWORK )\n\n*     This subroutine is identical to DLAQR3 except that it avoids\n*     recursion by calling DLAHQR instead of DLAQR4.\n*\n*\n*     ******************************************************************\n*     Aggressive early deflation:\n*\n*     This subroutine accepts as input an upper Hessenberg matrix\n*     H and performs an orthogonal similarity transformation\n*     designed to detect and deflate fully converged eigenvalues from\n*     a trailing principal submatrix.  On output H has been over-\n*     written by a new Hessenberg matrix that is a perturbation of\n*     an orthogonal similarity transformation of H.  It is to be\n*     hoped that the final version of H has many zero subdiagonal\n*     entries.\n*\n*     ******************************************************************\n\n*     WANTT   (input) LOGICAL\n*          If .TRUE., then the Hessenberg matrix H is fully updated\n*          so that the quasi-triangular Schur factor may be\n*          computed (in cooperation with the calling subroutine).\n*          If .FALSE., then only enough of H is updated to preserve\n*          the eigenvalues.\n*\n*     WANTZ   (input) LOGICAL\n*          If .TRUE., then the orthogonal matrix Z is updated so\n*          so that the orthogonal Schur factor may be computed\n*          (in cooperation with the calling subroutine).\n*          If .FALSE., then Z is not referenced.\n*\n*     N       (input) INTEGER\n*          The order of the matrix H and (if WANTZ is .TRUE.) the\n*          order of the orthogonal matrix Z.\n*\n*     KTOP    (input) INTEGER\n*          It is assumed that either KTOP = 1 or H(KTOP,KTOP-1)=0.\n*          KBOT and KTOP together determine an isolated block\n*          along the diagonal of the Hessenberg matrix.\n*\n*     KBOT    (input) INTEGER\n*          It is assumed without a check that either\n*          KBOT = N or H(KBOT+1,KBOT)=0.  KBOT and KTOP together\n*          determine an isolated block along the diagonal of the\n*          Hessenberg matrix.\n*\n*     NW      (input) INTEGER\n*          Deflation window size.  1 .LE. NW .LE. (KBOT-KTOP+1).\n*\n*     H       (input/output) DOUBLE PRECISION array, dimension (LDH,N)\n*          On input the initial N-by-N section of H stores the\n*          Hessenberg matrix undergoing aggressive early deflation.\n*          On output H has been transformed by an orthogonal\n*          similarity transformation, perturbed, and the returned\n*          to Hessenberg form that (it is to be hoped) has some\n*          zero subdiagonal entries.\n*\n*     LDH     (input) integer\n*          Leading dimension of H just as declared in the calling\n*          subroutine.  N .LE. LDH\n*\n*     ILOZ    (input) INTEGER\n*     IHIZ    (input) INTEGER\n*          Specify the rows of Z to which transformations must be\n*          applied if WANTZ is .TRUE.. 1 .LE. ILOZ .LE. IHIZ .LE. N.\n*\n*     Z       (input/output) DOUBLE PRECISION array, dimension (LDZ,N)\n*          IF WANTZ is .TRUE., then on output, the orthogonal\n*          similarity transformation mentioned above has been\n*          accumulated into Z(ILOZ:IHIZ,ILO:IHI) from the right.\n*          If WANTZ is .FALSE., then Z is unreferenced.\n*\n*     LDZ     (input) integer\n*          The leading dimension of Z just as declared in the\n*          calling subroutine.  1 .LE. LDZ.\n*\n*     NS      (output) integer\n*          The number of unconverged (ie approximate) eigenvalues\n*          returned in SR and SI that may be used as shifts by the\n*          calling subroutine.\n*\n*     ND      (output) integer\n*          The number of converged eigenvalues uncovered by this\n*          subroutine.\n*\n*     SR      (output) DOUBLE PRECISION array, dimension (KBOT)\n*     SI      (output) DOUBLE PRECISION array, dimension (KBOT)\n*          On output, the real and imaginary parts of approximate\n*          eigenvalues that may be used for shifts are stored in\n*          SR(KBOT-ND-NS+1) through SR(KBOT-ND) and\n*          SI(KBOT-ND-NS+1) through SI(KBOT-ND), respectively.\n*          The real and imaginary parts of converged eigenvalues\n*          are stored in SR(KBOT-ND+1) through SR(KBOT) and\n*          SI(KBOT-ND+1) through SI(KBOT), respectively.\n*\n*     V       (workspace) DOUBLE PRECISION array, dimension (LDV,NW)\n*          An NW-by-NW work array.\n*\n*     LDV     (input) integer scalar\n*          The leading dimension of V just as declared in the\n*          calling subroutine.  NW .LE. LDV\n*\n*     NH      (input) integer scalar\n*          The number of columns of T.  NH.GE.NW.\n*\n*     T       (workspace) DOUBLE PRECISION array, dimension (LDT,NW)\n*\n*     LDT     (input) integer\n*          The leading dimension of T just as declared in the\n*          calling subroutine.  NW .LE. LDT\n*\n*     NV      (input) integer\n*          The number of rows of work array WV available for\n*          workspace.  NV.GE.NW.\n*\n*     WV      (workspace) DOUBLE PRECISION array, dimension (LDWV,NW)\n*\n*     LDWV    (input) integer\n*          The leading dimension of W just as declared in the\n*          calling subroutine.  NW .LE. LDV\n*\n*     WORK    (workspace) DOUBLE PRECISION array, dimension (LWORK)\n*          On exit, WORK(1) is set to an estimate of the optimal value\n*          of LWORK for the given values of N, NW, KTOP and KBOT.\n*\n*     LWORK   (input) integer\n*          The dimension of the work array WORK.  LWORK = 2*NW\n*          suffices, but greater efficiency may result from larger\n*          values of LWORK.\n*\n*          If LWORK = -1, then a workspace query is assumed; DLAQR2\n*          only estimates the optimal workspace size for the given\n*          values of N, NW, KTOP and KBOT.  The estimate is returned\n*          in WORK(1).  No error message related to LWORK is issued\n*          by XERBLA.  Neither H nor Z are accessed.\n*\n\n*     ================================================================\n*     Based on contributions by\n*        Karen Braman and Ralph Byers, Department of Mathematics,\n*        University of Kansas, USA\n*\n*     ================================================================\n\n");
    return Qnil;
  }
  if (argc != 12)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 12)", argc);
  rb_wantt = argv[0];
  rb_wantz = argv[1];
  rb_ktop = argv[2];
  rb_kbot = argv[3];
  rb_nw = argv[4];
  rb_h = argv[5];
  rb_iloz = argv[6];
  rb_ihiz = argv[7];
  rb_z = argv[8];
  rb_nh = argv[9];
  rb_nv = argv[10];
  rb_lwork = argv[11];

  ktop = NUM2INT(rb_ktop);
  wantz = (rb_wantz == Qtrue);
  nh = NUM2INT(rb_nh);
  if (!NA_IsNArray(rb_z))
    rb_raise(rb_eArgError, "z (9th argument) must be NArray");
  if (NA_RANK(rb_z) != 2)
    rb_raise(rb_eArgError, "rank of z (9th argument) must be %d", 2);
  n = NA_SHAPE1(rb_z);
  ldz = NA_SHAPE0(rb_z);
  if (NA_TYPE(rb_z) != NA_DFLOAT)
    rb_z = na_change_type(rb_z, NA_DFLOAT);
  z = NA_PTR_TYPE(rb_z, doublereal*);
  kbot = NUM2INT(rb_kbot);
  lwork = NUM2INT(rb_lwork);
  iloz = NUM2INT(rb_iloz);
  nv = NUM2INT(rb_nv);
  wantt = (rb_wantt == Qtrue);
  ihiz = NUM2INT(rb_ihiz);
  if (!NA_IsNArray(rb_h))
    rb_raise(rb_eArgError, "h (6th argument) must be NArray");
  if (NA_RANK(rb_h) != 2)
    rb_raise(rb_eArgError, "rank of h (6th argument) must be %d", 2);
  if (NA_SHAPE1(rb_h) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of h must be the same as shape 1 of z");
  ldh = NA_SHAPE0(rb_h);
  if (NA_TYPE(rb_h) != NA_DFLOAT)
    rb_h = na_change_type(rb_h, NA_DFLOAT);
  h = NA_PTR_TYPE(rb_h, doublereal*);
  nw = NUM2INT(rb_nw);
  ldv = nw;
  ldwv = nw;
  ldt = nw;
  {
    int shape[1];
    shape[0] = MAX(1,kbot);
    rb_sr = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  sr = NA_PTR_TYPE(rb_sr, doublereal*);
  {
    int shape[1];
    shape[0] = MAX(1,kbot);
    rb_si = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  si = NA_PTR_TYPE(rb_si, doublereal*);
  {
    int shape[2];
    shape[0] = ldh;
    shape[1] = n;
    rb_h_out__ = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  h_out__ = NA_PTR_TYPE(rb_h_out__, doublereal*);
  MEMCPY(h_out__, h, doublereal, NA_TOTAL(rb_h));
  rb_h = rb_h_out__;
  h = h_out__;
  {
    int shape[2];
    shape[0] = ldz;
    shape[1] = n;
    rb_z_out__ = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  z_out__ = NA_PTR_TYPE(rb_z_out__, doublereal*);
  MEMCPY(z_out__, z, doublereal, NA_TOTAL(rb_z));
  rb_z = rb_z_out__;
  z = z_out__;
  v = ALLOC_N(doublereal, (ldv)*(MAX(1,nw)));
  t = ALLOC_N(doublereal, (ldt)*(MAX(1,nw)));
  wv = ALLOC_N(doublereal, (ldwv)*(MAX(1,nw)));
  work = ALLOC_N(doublereal, (MAX(1,lwork)));

  dlaqr2_(&wantt, &wantz, &n, &ktop, &kbot, &nw, h, &ldh, &iloz, &ihiz, z, &ldz, &ns, &nd, sr, si, v, &ldv, &nh, t, &ldt, &nv, wv, &ldwv, work, &lwork);

  free(v);
  free(t);
  free(wv);
  free(work);
  rb_ns = INT2NUM(ns);
  rb_nd = INT2NUM(nd);
  return rb_ary_new3(6, rb_ns, rb_nd, rb_sr, rb_si, rb_h, rb_z);
}

void
init_lapack_dlaqr2(VALUE mLapack){
  rb_define_module_function(mLapack, "dlaqr2", rb_dlaqr2, -1);
}

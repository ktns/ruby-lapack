#include "rb_lapack.h"

extern VOID slaqr5_(logical *wantt, logical *wantz, integer *kacc22, integer *n, integer *ktop, integer *kbot, integer *nshfts, real *sr, real *si, real *h, integer *ldh, integer *iloz, integer *ihiz, real *z, integer *ldz, real *v, integer *ldv, real *u, integer *ldu, integer *nv, real *wv, integer *ldwv, integer *nh, real *wh, integer *ldwh);

static VALUE sHelp, sUsage;

static VALUE
rblapack_slaqr5(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_wantt;
  logical wantt; 
  VALUE rblapack_wantz;
  logical wantz; 
  VALUE rblapack_kacc22;
  integer kacc22; 
  VALUE rblapack_ktop;
  integer ktop; 
  VALUE rblapack_kbot;
  integer kbot; 
  VALUE rblapack_sr;
  real *sr; 
  VALUE rblapack_si;
  real *si; 
  VALUE rblapack_h;
  real *h; 
  VALUE rblapack_iloz;
  integer iloz; 
  VALUE rblapack_ihiz;
  integer ihiz; 
  VALUE rblapack_z;
  real *z; 
  VALUE rblapack_ldz;
  integer ldz; 
  VALUE rblapack_nv;
  integer nv; 
  VALUE rblapack_nh;
  integer nh; 
  VALUE rblapack_sr_out__;
  real *sr_out__;
  VALUE rblapack_si_out__;
  real *si_out__;
  VALUE rblapack_h_out__;
  real *h_out__;
  VALUE rblapack_z_out__;
  real *z_out__;
  real *v;
  real *u;
  real *wv;
  real *wh;

  integer nshfts;
  integer ldh;
  integer n;
  integer ldv;
  integer ldu;
  integer ldwv;
  integer ldwh;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  sr, si, h, z = NumRu::Lapack.slaqr5( wantt, wantz, kacc22, ktop, kbot, sr, si, h, iloz, ihiz, z, ldz, nv, nh, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE SLAQR5( WANTT, WANTZ, KACC22, N, KTOP, KBOT, NSHFTS, SR, SI, H, LDH, ILOZ, IHIZ, Z, LDZ, V, LDV, U, LDU, NV, WV, LDWV, NH, WH, LDWH )\n\n*     This auxiliary subroutine called by SLAQR0 performs a\n*     single small-bulge multi-shift QR sweep.\n*\n\n*      WANTT  (input) logical scalar\n*             WANTT = .true. if the quasi-triangular Schur factor\n*             is being computed.  WANTT is set to .false. otherwise.\n*\n*      WANTZ  (input) logical scalar\n*             WANTZ = .true. if the orthogonal Schur factor is being\n*             computed.  WANTZ is set to .false. otherwise.\n*\n*      KACC22 (input) integer with value 0, 1, or 2.\n*             Specifies the computation mode of far-from-diagonal\n*             orthogonal updates.\n*        = 0: SLAQR5 does not accumulate reflections and does not\n*             use matrix-matrix multiply to update far-from-diagonal\n*             matrix entries.\n*        = 1: SLAQR5 accumulates reflections and uses matrix-matrix\n*             multiply to update the far-from-diagonal matrix entries.\n*        = 2: SLAQR5 accumulates reflections, uses matrix-matrix\n*             multiply to update the far-from-diagonal matrix entries,\n*             and takes advantage of 2-by-2 block structure during\n*             matrix multiplies.\n*\n*      N      (input) integer scalar\n*             N is the order of the Hessenberg matrix H upon which this\n*             subroutine operates.\n*\n*      KTOP   (input) integer scalar\n*      KBOT   (input) integer scalar\n*             These are the first and last rows and columns of an\n*             isolated diagonal block upon which the QR sweep is to be\n*             applied. It is assumed without a check that\n*                       either KTOP = 1  or   H(KTOP,KTOP-1) = 0\n*             and\n*                       either KBOT = N  or   H(KBOT+1,KBOT) = 0.\n*\n*      NSHFTS (input) integer scalar\n*             NSHFTS gives the number of simultaneous shifts.  NSHFTS\n*             must be positive and even.\n*\n*      SR     (input/output) REAL array of size (NSHFTS)\n*      SI     (input/output) REAL array of size (NSHFTS)\n*             SR contains the real parts and SI contains the imaginary\n*             parts of the NSHFTS shifts of origin that define the\n*             multi-shift QR sweep.  On output SR and SI may be\n*             reordered.\n*\n*      H      (input/output) REAL array of size (LDH,N)\n*             On input H contains a Hessenberg matrix.  On output a\n*             multi-shift QR sweep with shifts SR(J)+i*SI(J) is applied\n*             to the isolated diagonal block in rows and columns KTOP\n*             through KBOT.\n*\n*      LDH    (input) integer scalar\n*             LDH is the leading dimension of H just as declared in the\n*             calling procedure.  LDH.GE.MAX(1,N).\n*\n*      ILOZ   (input) INTEGER\n*      IHIZ   (input) INTEGER\n*             Specify the rows of Z to which transformations must be\n*             applied if WANTZ is .TRUE.. 1 .LE. ILOZ .LE. IHIZ .LE. N\n*\n*      Z      (input/output) REAL array of size (LDZ,IHI)\n*             If WANTZ = .TRUE., then the QR Sweep orthogonal\n*             similarity transformation is accumulated into\n*             Z(ILOZ:IHIZ,ILO:IHI) from the right.\n*             If WANTZ = .FALSE., then Z is unreferenced.\n*\n*      LDZ    (input) integer scalar\n*             LDA is the leading dimension of Z just as declared in\n*             the calling procedure. LDZ.GE.N.\n*\n*      V      (workspace) REAL array of size (LDV,NSHFTS/2)\n*\n*      LDV    (input) integer scalar\n*             LDV is the leading dimension of V as declared in the\n*             calling procedure.  LDV.GE.3.\n*\n*      U      (workspace) REAL array of size\n*             (LDU,3*NSHFTS-3)\n*\n*      LDU    (input) integer scalar\n*             LDU is the leading dimension of U just as declared in the\n*             in the calling subroutine.  LDU.GE.3*NSHFTS-3.\n*\n*      NH     (input) integer scalar\n*             NH is the number of columns in array WH available for\n*             workspace. NH.GE.1.\n*\n*      WH     (workspace) REAL array of size (LDWH,NH)\n*\n*      LDWH   (input) integer scalar\n*             Leading dimension of WH just as declared in the\n*             calling procedure.  LDWH.GE.3*NSHFTS-3.\n*\n*      NV     (input) integer scalar\n*             NV is the number of rows in WV agailable for workspace.\n*             NV.GE.1.\n*\n*      WV     (workspace) REAL array of size\n*             (LDWV,3*NSHFTS-3)\n*\n*      LDWV   (input) integer scalar\n*             LDWV is the leading dimension of WV as declared in the\n*             in the calling subroutine.  LDWV.GE.NV.\n*\n\n*     ================================================================\n*     Based on contributions by\n*        Karen Braman and Ralph Byers, Department of Mathematics,\n*        University of Kansas, USA\n*\n*     ================================================================\n*     Reference:\n*\n*     K. Braman, R. Byers and R. Mathias, The Multi-Shift QR\n*     Algorithm Part I: Maintaining Well Focused Shifts, and\n*     Level 3 Performance, SIAM Journal of Matrix Analysis,\n*     volume 23, pages 929--947, 2002.\n*\n*     ================================================================\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  sr, si, h, z = NumRu::Lapack.slaqr5( wantt, wantz, kacc22, ktop, kbot, sr, si, h, iloz, ihiz, z, ldz, nv, nh, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 14)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 14)", argc);
  rblapack_wantt = argv[0];
  rblapack_wantz = argv[1];
  rblapack_kacc22 = argv[2];
  rblapack_ktop = argv[3];
  rblapack_kbot = argv[4];
  rblapack_sr = argv[5];
  rblapack_si = argv[6];
  rblapack_h = argv[7];
  rblapack_iloz = argv[8];
  rblapack_ihiz = argv[9];
  rblapack_z = argv[10];
  rblapack_ldz = argv[11];
  rblapack_nv = argv[12];
  rblapack_nh = argv[13];
  if (rb_options != Qnil) {
  }

  if (!NA_IsNArray(rblapack_si))
    rb_raise(rb_eArgError, "si (7th argument) must be NArray");
  if (NA_RANK(rblapack_si) != 1)
    rb_raise(rb_eArgError, "rank of si (7th argument) must be %d", 1);
  nshfts = NA_SHAPE0(rblapack_si);
  if (NA_TYPE(rblapack_si) != NA_SFLOAT)
    rblapack_si = na_change_type(rblapack_si, NA_SFLOAT);
  si = NA_PTR_TYPE(rblapack_si, real*);
  kacc22 = NUM2INT(rblapack_kacc22);
  ktop = NUM2INT(rblapack_ktop);
  wantz = (rblapack_wantz == Qtrue);
  if (!NA_IsNArray(rblapack_sr))
    rb_raise(rb_eArgError, "sr (6th argument) must be NArray");
  if (NA_RANK(rblapack_sr) != 1)
    rb_raise(rb_eArgError, "rank of sr (6th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_sr) != nshfts)
    rb_raise(rb_eRuntimeError, "shape 0 of sr must be the same as shape 0 of si");
  if (NA_TYPE(rblapack_sr) != NA_SFLOAT)
    rblapack_sr = na_change_type(rblapack_sr, NA_SFLOAT);
  sr = NA_PTR_TYPE(rblapack_sr, real*);
  kbot = NUM2INT(rblapack_kbot);
  nh = NUM2INT(rblapack_nh);
  if (!NA_IsNArray(rblapack_h))
    rb_raise(rb_eArgError, "h (8th argument) must be NArray");
  if (NA_RANK(rblapack_h) != 2)
    rb_raise(rb_eArgError, "rank of h (8th argument) must be %d", 2);
  n = NA_SHAPE1(rblapack_h);
  ldh = NA_SHAPE0(rblapack_h);
  if (NA_TYPE(rblapack_h) != NA_SFLOAT)
    rblapack_h = na_change_type(rblapack_h, NA_SFLOAT);
  h = NA_PTR_TYPE(rblapack_h, real*);
  nv = NUM2INT(rblapack_nv);
  ihiz = NUM2INT(rblapack_ihiz);
  wantt = (rblapack_wantt == Qtrue);
  ldv = 3;
  iloz = NUM2INT(rblapack_iloz);
  ldz = n;
  if (!NA_IsNArray(rblapack_z))
    rb_raise(rb_eArgError, "z (11th argument) must be NArray");
  if (NA_RANK(rblapack_z) != 2)
    rb_raise(rb_eArgError, "rank of z (11th argument) must be %d", 2);
  if (NA_SHAPE1(rblapack_z) != (wantz ? ihiz : 0))
    rb_raise(rb_eRuntimeError, "shape 1 of z must be %d", wantz ? ihiz : 0);
  if (NA_SHAPE0(rblapack_z) != (wantz ? ldz : 0))
    rb_raise(rb_eRuntimeError, "shape 0 of z must be %d", wantz ? ldz : 0);
  if (NA_TYPE(rblapack_z) != NA_SFLOAT)
    rblapack_z = na_change_type(rblapack_z, NA_SFLOAT);
  z = NA_PTR_TYPE(rblapack_z, real*);
  ldwv = nv;
  ldu = 3*nshfts-3;
  ldwh = 3*nshfts-3;
  {
    int shape[1];
    shape[0] = nshfts;
    rblapack_sr_out__ = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  sr_out__ = NA_PTR_TYPE(rblapack_sr_out__, real*);
  MEMCPY(sr_out__, sr, real, NA_TOTAL(rblapack_sr));
  rblapack_sr = rblapack_sr_out__;
  sr = sr_out__;
  {
    int shape[1];
    shape[0] = nshfts;
    rblapack_si_out__ = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  si_out__ = NA_PTR_TYPE(rblapack_si_out__, real*);
  MEMCPY(si_out__, si, real, NA_TOTAL(rblapack_si));
  rblapack_si = rblapack_si_out__;
  si = si_out__;
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
    shape[0] = wantz ? ldz : 0;
    shape[1] = wantz ? ihiz : 0;
    rblapack_z_out__ = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  z_out__ = NA_PTR_TYPE(rblapack_z_out__, real*);
  MEMCPY(z_out__, z, real, NA_TOTAL(rblapack_z));
  rblapack_z = rblapack_z_out__;
  z = z_out__;
  v = ALLOC_N(real, (ldv)*(nshfts/2));
  u = ALLOC_N(real, (ldu)*(3*nshfts-3));
  wv = ALLOC_N(real, (ldwv)*(3*nshfts-3));
  wh = ALLOC_N(real, (ldwh)*(MAX(1,nh)));

  slaqr5_(&wantt, &wantz, &kacc22, &n, &ktop, &kbot, &nshfts, sr, si, h, &ldh, &iloz, &ihiz, z, &ldz, v, &ldv, u, &ldu, &nv, wv, &ldwv, &nh, wh, &ldwh);

  free(v);
  free(u);
  free(wv);
  free(wh);
  return rb_ary_new3(4, rblapack_sr, rblapack_si, rblapack_h, rblapack_z);
}

void
init_lapack_slaqr5(VALUE mLapack){
  rb_define_module_function(mLapack, "slaqr5", rblapack_slaqr5, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}

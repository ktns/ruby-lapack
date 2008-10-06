#include "rb_lapack.h"

static VALUE
rb_dlaqr5(int argc, VALUE *argv, VALUE self){
  VALUE rb_wantt;
  logical wantt; 
  VALUE rb_wantz;
  logical wantz; 
  VALUE rb_kacc22;
  integer kacc22; 
  VALUE rb_ktop;
  integer ktop; 
  VALUE rb_kbot;
  integer kbot; 
  VALUE rb_sr;
  doublereal *sr; 
  VALUE rb_si;
  doublereal *si; 
  VALUE rb_h;
  doublereal *h; 
  VALUE rb_iloz;
  integer iloz; 
  VALUE rb_ihiz;
  integer ihiz; 
  VALUE rb_z;
  doublereal *z; 
  VALUE rb_ldz;
  integer ldz; 
  VALUE rb_nv;
  integer nv; 
  VALUE rb_nh;
  integer nh; 
  VALUE rb_h_out__;
  doublereal *h_out__;
  VALUE rb_z_out__;
  doublereal *z_out__;
  doublereal *v;
  doublereal *u;
  doublereal *wv;
  doublereal *wh;

  integer nshfts;
  integer ldh;
  integer n;
  integer ldv;
  integer ldu;
  integer ldwv;
  integer ldwh;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  h, z = NumRu::Lapack.dlaqr5( wantt, wantz, kacc22, ktop, kbot, sr, si, h, iloz, ihiz, z, ldz, nv, nh)\n    or\n  NumRu::Lapack.dlaqr5  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE DLAQR5( WANTT, WANTZ, KACC22, N, KTOP, KBOT, NSHFTS, SR, SI, H, LDH, ILOZ, IHIZ, Z, LDZ, V, LDV, U, LDU, NV, WV, LDWV, NH, WH, LDWH )\n\n*     This auxiliary subroutine called by DLAQR0 performs a\n*     single small-bulge multi-shift QR sweep.\n*\n\n*      WANTT  (input) logical scalar\n*             WANTT = .true. if the quasi-triangular Schur factor\n*             is being computed.  WANTT is set to .false. otherwise.\n*\n*      WANTZ  (input) logical scalar\n*             WANTZ = .true. if the orthogonal Schur factor is being\n*             computed.  WANTZ is set to .false. otherwise.\n*\n*      KACC22 (input) integer with value 0, 1, or 2.\n*             Specifies the computation mode of far-from-diagonal\n*             orthogonal updates.\n*        = 0: DLAQR5 does not accumulate reflections and does not\n*             use matrix-matrix multiply to update far-from-diagonal\n*             matrix entries.\n*        = 1: DLAQR5 accumulates reflections and uses matrix-matrix\n*             multiply to update the far-from-diagonal matrix entries.\n*        = 2: DLAQR5 accumulates reflections, uses matrix-matrix\n*             multiply to update the far-from-diagonal matrix entries,\n*             and takes advantage of 2-by-2 block structure during\n*             matrix multiplies.\n*\n*      N      (input) integer scalar\n*             N is the order of the Hessenberg matrix H upon which this\n*             subroutine operates.\n*\n*      KTOP   (input) integer scalar\n*      KBOT   (input) integer scalar\n*             These are the first and last rows and columns of an\n*             isolated diagonal block upon which the QR sweep is to be\n*             applied. It is assumed without a check that\n*                       either KTOP = 1  or   H(KTOP,KTOP-1) = 0\n*             and\n*                       either KBOT = N  or   H(KBOT+1,KBOT) = 0.\n*\n*      NSHFTS (input) integer scalar\n*             NSHFTS gives the number of simultaneous shifts.  NSHFTS\n*             must be positive and even.\n*\n*      SR     (input) DOUBLE PRECISION array of size (NSHFTS)\n*      SI     (input) DOUBLE PRECISION array of size (NSHFTS)\n*             SR contains the real parts and SI contains the imaginary\n*             parts of the NSHFTS shifts of origin that define the\n*             multi-shift QR sweep.\n*\n*      H      (input/output) DOUBLE PRECISION array of size (LDH,N)\n*             On input H contains a Hessenberg matrix.  On output a\n*             multi-shift QR sweep with shifts SR(J)+i*SI(J) is applied\n*             to the isolated diagonal block in rows and columns KTOP\n*             through KBOT.\n*\n*      LDH    (input) integer scalar\n*             LDH is the leading dimension of H just as declared in the\n*             calling procedure.  LDH.GE.MAX(1,N).\n*\n*      ILOZ   (input) INTEGER\n*      IHIZ   (input) INTEGER\n*             Specify the rows of Z to which transformations must be\n*             applied if WANTZ is .TRUE.. 1 .LE. ILOZ .LE. IHIZ .LE. N\n*\n*      Z      (input/output) DOUBLE PRECISION array of size (LDZ,IHI)\n*             If WANTZ = .TRUE., then the QR Sweep orthogonal\n*             similarity transformation is accumulated into\n*             Z(ILOZ:IHIZ,ILO:IHI) from the right.\n*             If WANTZ = .FALSE., then Z is unreferenced.\n*\n*      LDZ    (input) integer scalar\n*             LDA is the leading dimension of Z just as declared in\n*             the calling procedure. LDZ.GE.N.\n*\n*      V      (workspace) DOUBLE PRECISION array of size (LDV,NSHFTS/2)\n*\n*      LDV    (input) integer scalar\n*             LDV is the leading dimension of V as declared in the\n*             calling procedure.  LDV.GE.3.\n*\n*      U      (workspace) DOUBLE PRECISION array of size\n*             (LDU,3*NSHFTS-3)\n*\n*      LDU    (input) integer scalar\n*             LDU is the leading dimension of U just as declared in the\n*             in the calling subroutine.  LDU.GE.3*NSHFTS-3.\n*\n*      NH     (input) integer scalar\n*             NH is the number of columns in array WH available for\n*             workspace. NH.GE.1.\n*\n*      WH     (workspace) DOUBLE PRECISION array of size (LDWH,NH)\n*\n*      LDWH   (input) integer scalar\n*             Leading dimension of WH just as declared in the\n*             calling procedure.  LDWH.GE.3*NSHFTS-3.\n*\n*      NV     (input) integer scalar\n*             NV is the number of rows in WV agailable for workspace.\n*             NV.GE.1.\n*\n*      WV     (workspace) DOUBLE PRECISION array of size\n*             (LDWV,3*NSHFTS-3)\n*\n*      LDWV   (input) integer scalar\n*             LDWV is the leading dimension of WV as declared in the\n*             in the calling subroutine.  LDWV.GE.NV.\n*\n*\n\n*     ================================================================\n*     Based on contributions by\n*        Karen Braman and Ralph Byers, Department of Mathematics,\n*        University of Kansas, USA\n*\n*     ============================================================\n*     Reference:\n*\n*     K. Braman, R. Byers and R. Mathias, The Multi-Shift QR\n*     Algorithm Part I: Maintaining Well Focused Shifts, and\n*     Level 3 Performance, SIAM Journal of Matrix Analysis,\n*     volume 23, pages 929--947, 2002.\n*\n*     ============================================================\n\n");
    return Qnil;
  }
  if (argc != 14)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 14)", argc);
  rb_wantt = argv[0];
  rb_wantz = argv[1];
  rb_kacc22 = argv[2];
  rb_ktop = argv[3];
  rb_kbot = argv[4];
  rb_sr = argv[5];
  rb_si = argv[6];
  rb_h = argv[7];
  rb_iloz = argv[8];
  rb_ihiz = argv[9];
  rb_z = argv[10];
  rb_ldz = argv[11];
  rb_nv = argv[12];
  rb_nh = argv[13];

  wantt = (rb_wantt == Qtrue);
  wantz = (rb_wantz == Qtrue);
  kacc22 = NUM2INT(rb_kacc22);
  ktop = NUM2INT(rb_ktop);
  kbot = NUM2INT(rb_kbot);
  iloz = NUM2INT(rb_iloz);
  ihiz = NUM2INT(rb_ihiz);
  ldz = NUM2INT(rb_ldz);
  nv = NUM2INT(rb_nv);
  nh = NUM2INT(rb_nh);
  if (!NA_IsNArray(rb_sr))
    rb_raise(rb_eArgError, "sr (6th argument) must be NArray");
  if (NA_RANK(rb_sr) != 1)
    rb_raise(rb_eArgError, "rank of sr (6th argument) must be %d", 1);
  nshfts = NA_SHAPE0(rb_sr);
  if (NA_TYPE(rb_sr) != NA_DFLOAT)
    rb_sr = na_change_type(rb_sr, NA_DFLOAT);
  sr = NA_PTR_TYPE(rb_sr, doublereal*);
  if (!NA_IsNArray(rb_si))
    rb_raise(rb_eArgError, "si (7th argument) must be NArray");
  if (NA_RANK(rb_si) != 1)
    rb_raise(rb_eArgError, "rank of si (7th argument) must be %d", 1);
  if (NA_SHAPE0(rb_si) != nshfts)
    rb_raise(rb_eRuntimeError, "shape 0 of si must be the same as shape 0 of sr");
  if (NA_TYPE(rb_si) != NA_DFLOAT)
    rb_si = na_change_type(rb_si, NA_DFLOAT);
  si = NA_PTR_TYPE(rb_si, doublereal*);
  if (!NA_IsNArray(rb_h))
    rb_raise(rb_eArgError, "h (8th argument) must be NArray");
  if (NA_RANK(rb_h) != 2)
    rb_raise(rb_eArgError, "rank of h (8th argument) must be %d", 2);
  ldh = NA_SHAPE0(rb_h);
  n = NA_SHAPE1(rb_h);
  if (NA_TYPE(rb_h) != NA_DFLOAT)
    rb_h = na_change_type(rb_h, NA_DFLOAT);
  h = NA_PTR_TYPE(rb_h, doublereal*);
  if (!NA_IsNArray(rb_z))
    rb_raise(rb_eArgError, "z (11th argument) must be NArray");
  if (NA_RANK(rb_z) != 2)
    rb_raise(rb_eArgError, "rank of z (11th argument) must be %d", 2);
  ldz = n;
  if (NA_SHAPE0(rb_z) != (wantz ? ldz : 0))
    rb_raise(rb_eRuntimeError, "shape 0 of z must be %d", wantz ? ldz : 0);
  if (NA_SHAPE1(rb_z) != (wantz ? ihiz : 0))
    rb_raise(rb_eRuntimeError, "shape 1 of z must be %d", wantz ? ihiz : 0);
  if (NA_TYPE(rb_z) != NA_DFLOAT)
    rb_z = na_change_type(rb_z, NA_DFLOAT);
  z = NA_PTR_TYPE(rb_z, doublereal*);
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
    shape[0] = wantz ? ldz : 0;
    shape[1] = wantz ? ihiz : 0;
    rb_z_out__ = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  z_out__ = NA_PTR_TYPE(rb_z_out__, doublereal*);
  MEMCPY(z_out__, z, doublereal, NA_TOTAL(rb_z));
  rb_z = rb_z_out__;
  z = z_out__;
  ldv = 3;
  v = ALLOC_N(doublereal, (ldv)*(nshfts/2));
  ldu = 3*nshfts-3;
  u = ALLOC_N(doublereal, (ldu)*(3*nshfts-3));
  ldwv = nv;
  wv = ALLOC_N(doublereal, (ldwv)*(3*nshfts-3));
  ldwh = 3*nshfts-3;
  wh = ALLOC_N(doublereal, (ldwh)*(MAX(1,nh)));

  dlaqr5_(&wantt, &wantz, &kacc22, &n, &ktop, &kbot, &nshfts, sr, si, h, &ldh, &iloz, &ihiz, z, &ldz, v, &ldv, u, &ldu, &nv, wv, &ldwv, &nh, wh, &ldwh);

  free(v);
  free(u);
  free(wv);
  free(wh);
  return rb_ary_new3(2, rb_h, rb_z);
}

void
init_lapack_dlaqr5(VALUE mLapack){
  rb_define_module_function(mLapack, "dlaqr5", rb_dlaqr5, -1);
}

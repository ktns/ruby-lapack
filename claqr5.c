#include "rb_lapack.h"

extern VOID claqr5_(logical *wantt, logical *wantz, integer *kacc22, integer *n, integer *ktop, integer *kbot, integer *nshfts, complex *s, complex *h, integer *ldh, integer *iloz, integer *ihiz, complex *z, integer *ldz, complex *v, integer *ldv, complex *u, integer *ldu, integer *nv, complex *wv, integer *ldwv, integer *nh, complex *wh, integer *ldwh);

static VALUE sHelp, sUsage;

static VALUE
rblapack_claqr5(int argc, VALUE *argv, VALUE self){
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
  VALUE rblapack_s;
  complex *s; 
  VALUE rblapack_h;
  complex *h; 
  VALUE rblapack_iloz;
  integer iloz; 
  VALUE rblapack_ihiz;
  integer ihiz; 
  VALUE rblapack_z;
  complex *z; 
  VALUE rblapack_ldz;
  integer ldz; 
  VALUE rblapack_nv;
  integer nv; 
  VALUE rblapack_nh;
  integer nh; 
  VALUE rblapack_s_out__;
  complex *s_out__;
  VALUE rblapack_h_out__;
  complex *h_out__;
  VALUE rblapack_z_out__;
  complex *z_out__;
  complex *v;
  complex *u;
  complex *wv;
  complex *wh;

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
      printf("%s\n", "USAGE:\n  s, h, z = NumRu::Lapack.claqr5( wantt, wantz, kacc22, ktop, kbot, s, h, iloz, ihiz, z, ldz, nv, nh, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE CLAQR5( WANTT, WANTZ, KACC22, N, KTOP, KBOT, NSHFTS, S, H, LDH, ILOZ, IHIZ, Z, LDZ, V, LDV, U, LDU, NV, WV, LDWV, NH, WH, LDWH )\n\n*     This auxiliary subroutine called by CLAQR0 performs a\n*     single small-bulge multi-shift QR sweep.\n*\n\n*      WANTT  (input) logical scalar\n*             WANTT = .true. if the triangular Schur factor\n*             is being computed.  WANTT is set to .false. otherwise.\n*\n*      WANTZ  (input) logical scalar\n*             WANTZ = .true. if the unitary Schur factor is being\n*             computed.  WANTZ is set to .false. otherwise.\n*\n*      KACC22 (input) integer with value 0, 1, or 2.\n*             Specifies the computation mode of far-from-diagonal\n*             orthogonal updates.\n*        = 0: CLAQR5 does not accumulate reflections and does not\n*             use matrix-matrix multiply to update far-from-diagonal\n*             matrix entries.\n*        = 1: CLAQR5 accumulates reflections and uses matrix-matrix\n*             multiply to update the far-from-diagonal matrix entries.\n*        = 2: CLAQR5 accumulates reflections, uses matrix-matrix\n*             multiply to update the far-from-diagonal matrix entries,\n*             and takes advantage of 2-by-2 block structure during\n*             matrix multiplies.\n*\n*      N      (input) integer scalar\n*             N is the order of the Hessenberg matrix H upon which this\n*             subroutine operates.\n*\n*      KTOP   (input) integer scalar\n*      KBOT   (input) integer scalar\n*             These are the first and last rows and columns of an\n*             isolated diagonal block upon which the QR sweep is to be\n*             applied. It is assumed without a check that\n*                       either KTOP = 1  or   H(KTOP,KTOP-1) = 0\n*             and\n*                       either KBOT = N  or   H(KBOT+1,KBOT) = 0.\n*\n*      NSHFTS (input) integer scalar\n*             NSHFTS gives the number of simultaneous shifts.  NSHFTS\n*             must be positive and even.\n*\n*      S      (input/output) COMPLEX array of size (NSHFTS)\n*             S contains the shifts of origin that define the multi-\n*             shift QR sweep.  On output S may be reordered.\n*\n*      H      (input/output) COMPLEX array of size (LDH,N)\n*             On input H contains a Hessenberg matrix.  On output a\n*             multi-shift QR sweep with shifts SR(J)+i*SI(J) is applied\n*             to the isolated diagonal block in rows and columns KTOP\n*             through KBOT.\n*\n*      LDH    (input) integer scalar\n*             LDH is the leading dimension of H just as declared in the\n*             calling procedure.  LDH.GE.MAX(1,N).\n*\n*      ILOZ   (input) INTEGER\n*      IHIZ   (input) INTEGER\n*             Specify the rows of Z to which transformations must be\n*             applied if WANTZ is .TRUE.. 1 .LE. ILOZ .LE. IHIZ .LE. N\n*\n*      Z      (input/output) COMPLEX array of size (LDZ,IHI)\n*             If WANTZ = .TRUE., then the QR Sweep unitary\n*             similarity transformation is accumulated into\n*             Z(ILOZ:IHIZ,ILO:IHI) from the right.\n*             If WANTZ = .FALSE., then Z is unreferenced.\n*\n*      LDZ    (input) integer scalar\n*             LDA is the leading dimension of Z just as declared in\n*             the calling procedure. LDZ.GE.N.\n*\n*      V      (workspace) COMPLEX array of size (LDV,NSHFTS/2)\n*\n*      LDV    (input) integer scalar\n*             LDV is the leading dimension of V as declared in the\n*             calling procedure.  LDV.GE.3.\n*\n*      U      (workspace) COMPLEX array of size\n*             (LDU,3*NSHFTS-3)\n*\n*      LDU    (input) integer scalar\n*             LDU is the leading dimension of U just as declared in the\n*             in the calling subroutine.  LDU.GE.3*NSHFTS-3.\n*\n*      NH     (input) integer scalar\n*             NH is the number of columns in array WH available for\n*             workspace. NH.GE.1.\n*\n*      WH     (workspace) COMPLEX array of size (LDWH,NH)\n*\n*      LDWH   (input) integer scalar\n*             Leading dimension of WH just as declared in the\n*             calling procedure.  LDWH.GE.3*NSHFTS-3.\n*\n*      NV     (input) integer scalar\n*             NV is the number of rows in WV agailable for workspace.\n*             NV.GE.1.\n*\n*      WV     (workspace) COMPLEX array of size\n*             (LDWV,3*NSHFTS-3)\n*\n*      LDWV   (input) integer scalar\n*             LDWV is the leading dimension of WV as declared in the\n*             in the calling subroutine.  LDWV.GE.NV.\n*\n\n*     ================================================================\n*     Based on contributions by\n*        Karen Braman and Ralph Byers, Department of Mathematics,\n*        University of Kansas, USA\n*\n*     ================================================================\n*     Reference:\n*\n*     K. Braman, R. Byers and R. Mathias, The Multi-Shift QR\n*     Algorithm Part I: Maintaining Well Focused Shifts, and\n*     Level 3 Performance, SIAM Journal of Matrix Analysis,\n*     volume 23, pages 929--947, 2002.\n*\n*     ================================================================\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  s, h, z = NumRu::Lapack.claqr5( wantt, wantz, kacc22, ktop, kbot, s, h, iloz, ihiz, z, ldz, nv, nh, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 13)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 13)", argc);
  rblapack_wantt = argv[0];
  rblapack_wantz = argv[1];
  rblapack_kacc22 = argv[2];
  rblapack_ktop = argv[3];
  rblapack_kbot = argv[4];
  rblapack_s = argv[5];
  rblapack_h = argv[6];
  rblapack_iloz = argv[7];
  rblapack_ihiz = argv[8];
  rblapack_z = argv[9];
  rblapack_ldz = argv[10];
  rblapack_nv = argv[11];
  rblapack_nh = argv[12];
  if (rb_options != Qnil) {
  }

  kacc22 = NUM2INT(rblapack_kacc22);
  ktop = NUM2INT(rblapack_ktop);
  wantz = (rblapack_wantz == Qtrue);
  kbot = NUM2INT(rblapack_kbot);
  ldz = NUM2INT(rblapack_ldz);
  nh = NUM2INT(rblapack_nh);
  if (!NA_IsNArray(rblapack_s))
    rb_raise(rb_eArgError, "s (6th argument) must be NArray");
  if (NA_RANK(rblapack_s) != 1)
    rb_raise(rb_eArgError, "rank of s (6th argument) must be %d", 1);
  nshfts = NA_SHAPE0(rblapack_s);
  if (NA_TYPE(rblapack_s) != NA_SCOMPLEX)
    rblapack_s = na_change_type(rblapack_s, NA_SCOMPLEX);
  s = NA_PTR_TYPE(rblapack_s, complex*);
  if (!NA_IsNArray(rblapack_h))
    rb_raise(rb_eArgError, "h (7th argument) must be NArray");
  if (NA_RANK(rblapack_h) != 2)
    rb_raise(rb_eArgError, "rank of h (7th argument) must be %d", 2);
  n = NA_SHAPE1(rblapack_h);
  ldh = NA_SHAPE0(rblapack_h);
  if (NA_TYPE(rblapack_h) != NA_SCOMPLEX)
    rblapack_h = na_change_type(rblapack_h, NA_SCOMPLEX);
  h = NA_PTR_TYPE(rblapack_h, complex*);
  ldv = 3;
  wantt = (rblapack_wantt == Qtrue);
  nv = NUM2INT(rblapack_nv);
  ihiz = NUM2INT(rblapack_ihiz);
  iloz = NUM2INT(rblapack_iloz);
  ldu = 3*nshfts-3;
  ldwv = nv;
  if (!NA_IsNArray(rblapack_z))
    rb_raise(rb_eArgError, "z (10th argument) must be NArray");
  if (NA_RANK(rblapack_z) != 2)
    rb_raise(rb_eArgError, "rank of z (10th argument) must be %d", 2);
  if (NA_SHAPE1(rblapack_z) != (wantz ? ihiz : 0))
    rb_raise(rb_eRuntimeError, "shape 1 of z must be %d", wantz ? ihiz : 0);
  if (NA_SHAPE0(rblapack_z) != (wantz ? ldz : 0))
    rb_raise(rb_eRuntimeError, "shape 0 of z must be %d", wantz ? ldz : 0);
  if (NA_TYPE(rblapack_z) != NA_SCOMPLEX)
    rblapack_z = na_change_type(rblapack_z, NA_SCOMPLEX);
  z = NA_PTR_TYPE(rblapack_z, complex*);
  ldwh = 3*nshfts-3;
  {
    int shape[1];
    shape[0] = nshfts;
    rblapack_s_out__ = na_make_object(NA_SCOMPLEX, 1, shape, cNArray);
  }
  s_out__ = NA_PTR_TYPE(rblapack_s_out__, complex*);
  MEMCPY(s_out__, s, complex, NA_TOTAL(rblapack_s));
  rblapack_s = rblapack_s_out__;
  s = s_out__;
  {
    int shape[2];
    shape[0] = ldh;
    shape[1] = n;
    rblapack_h_out__ = na_make_object(NA_SCOMPLEX, 2, shape, cNArray);
  }
  h_out__ = NA_PTR_TYPE(rblapack_h_out__, complex*);
  MEMCPY(h_out__, h, complex, NA_TOTAL(rblapack_h));
  rblapack_h = rblapack_h_out__;
  h = h_out__;
  {
    int shape[2];
    shape[0] = wantz ? ldz : 0;
    shape[1] = wantz ? ihiz : 0;
    rblapack_z_out__ = na_make_object(NA_SCOMPLEX, 2, shape, cNArray);
  }
  z_out__ = NA_PTR_TYPE(rblapack_z_out__, complex*);
  MEMCPY(z_out__, z, complex, NA_TOTAL(rblapack_z));
  rblapack_z = rblapack_z_out__;
  z = z_out__;
  v = ALLOC_N(complex, (ldv)*(nshfts/2));
  u = ALLOC_N(complex, (ldu)*(3*nshfts-3));
  wv = ALLOC_N(complex, (ldwv)*(3*nshfts-3));
  wh = ALLOC_N(complex, (ldwh)*(MAX(1,nh)));

  claqr5_(&wantt, &wantz, &kacc22, &n, &ktop, &kbot, &nshfts, s, h, &ldh, &iloz, &ihiz, z, &ldz, v, &ldv, u, &ldu, &nv, wv, &ldwv, &nh, wh, &ldwh);

  free(v);
  free(u);
  free(wv);
  free(wh);
  return rb_ary_new3(3, rblapack_s, rblapack_h, rblapack_z);
}

void
init_lapack_claqr5(VALUE mLapack){
  rb_define_module_function(mLapack, "claqr5", rblapack_claqr5, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}

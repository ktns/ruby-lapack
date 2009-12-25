#include "rb_lapack.h"

static VALUE
rb_claqr3(int argc, VALUE *argv, VALUE self){
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
  complex *h; 
  VALUE rb_iloz;
  integer iloz; 
  VALUE rb_ihiz;
  integer ihiz; 
  VALUE rb_z;
  complex *z; 
  VALUE rb_nh;
  integer nh; 
  VALUE rb_ldt;
  integer ldt; 
  VALUE rb_nv;
  integer nv; 
  VALUE rb_ldwv;
  integer ldwv; 
  VALUE rb_lwork;
  integer lwork; 
  VALUE rb_ns;
  integer ns; 
  VALUE rb_nd;
  integer nd; 
  VALUE rb_sh;
  complex *sh; 
  VALUE rb_h_out__;
  complex *h_out__;
  VALUE rb_z_out__;
  complex *z_out__;
  complex *v;
  complex *t;
  complex *wv;
  complex *work;

  integer ldh;
  integer n;
  integer ldz;
  integer ldv;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  ns, nd, sh, h, z = NumRu::Lapack.claqr3( wantt, wantz, ktop, kbot, nw, h, iloz, ihiz, z, nh, ldt, nv, ldwv, lwork)\n    or\n  NumRu::Lapack.claqr3  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE CLAQR3( WANTT, WANTZ, N, KTOP, KBOT, NW, H, LDH, ILOZ, IHIZ, Z, LDZ, NS, ND, SH, V, LDV, NH, T, LDT, NV, WV, LDWV, WORK, LWORK )\n\n*     Aggressive early deflation:\n*\n*     This subroutine accepts as input an upper Hessenberg matrix\n*     H and performs an unitary similarity transformation\n*     designed to detect and deflate fully converged eigenvalues from\n*     a trailing principal submatrix.  On output H has been over-\n*     written by a new Hessenberg matrix that is a perturbation of\n*     an unitary similarity transformation of H.  It is to be\n*     hoped that the final version of H has many zero subdiagonal\n*     entries.\n*\n*     ******************************************************************\n\n*     WANTT   (input) LOGICAL\n*          If .TRUE., then the Hessenberg matrix H is fully updated\n*          so that the triangular Schur factor may be\n*          computed (in cooperation with the calling subroutine).\n*          If .FALSE., then only enough of H is updated to preserve\n*          the eigenvalues.\n*\n*     WANTZ   (input) LOGICAL\n*          If .TRUE., then the unitary matrix Z is updated so\n*          so that the unitary Schur factor may be computed\n*          (in cooperation with the calling subroutine).\n*          If .FALSE., then Z is not referenced.\n*\n*     N       (input) INTEGER\n*          The order of the matrix H and (if WANTZ is .TRUE.) the\n*          order of the unitary matrix Z.\n*\n*     KTOP    (input) INTEGER\n*          It is assumed that either KTOP = 1 or H(KTOP,KTOP-1)=0.\n*          KBOT and KTOP together determine an isolated block\n*          along the diagonal of the Hessenberg matrix.\n*\n*     KBOT    (input) INTEGER\n*          It is assumed without a check that either\n*          KBOT = N or H(KBOT+1,KBOT)=0.  KBOT and KTOP together\n*          determine an isolated block along the diagonal of the\n*          Hessenberg matrix.\n*\n*     NW      (input) INTEGER\n*          Deflation window size.  1 .LE. NW .LE. (KBOT-KTOP+1).\n*\n*     H       (input/output) COMPLEX array, dimension (LDH,N)\n*          On input the initial N-by-N section of H stores the\n*          Hessenberg matrix undergoing aggressive early deflation.\n*          On output H has been transformed by a unitary\n*          similarity transformation, perturbed, and the returned\n*          to Hessenberg form that (it is to be hoped) has some\n*          zero subdiagonal entries.\n*\n*     LDH     (input) integer\n*          Leading dimension of H just as declared in the calling\n*          subroutine.  N .LE. LDH\n*\n*     ILOZ    (input) INTEGER\n*     IHIZ    (input) INTEGER\n*          Specify the rows of Z to which transformations must be\n*          applied if WANTZ is .TRUE.. 1 .LE. ILOZ .LE. IHIZ .LE. N.\n*\n*     Z       (input/output) COMPLEX array, dimension (LDZ,N)\n*          IF WANTZ is .TRUE., then on output, the unitary\n*          similarity transformation mentioned above has been\n*          accumulated into Z(ILOZ:IHIZ,ILO:IHI) from the right.\n*          If WANTZ is .FALSE., then Z is unreferenced.\n*\n*     LDZ     (input) integer\n*          The leading dimension of Z just as declared in the\n*          calling subroutine.  1 .LE. LDZ.\n*\n*     NS      (output) integer\n*          The number of unconverged (ie approximate) eigenvalues\n*          returned in SR and SI that may be used as shifts by the\n*          calling subroutine.\n*\n*     ND      (output) integer\n*          The number of converged eigenvalues uncovered by this\n*          subroutine.\n*\n*     SH      (output) COMPLEX array, dimension KBOT\n*          On output, approximate eigenvalues that may\n*          be used for shifts are stored in SH(KBOT-ND-NS+1)\n*          through SR(KBOT-ND).  Converged eigenvalues are\n*          stored in SH(KBOT-ND+1) through SH(KBOT).\n*\n*     V       (workspace) COMPLEX array, dimension (LDV,NW)\n*          An NW-by-NW work array.\n*\n*     LDV     (input) integer scalar\n*          The leading dimension of V just as declared in the\n*          calling subroutine.  NW .LE. LDV\n*\n*     NH      (input) integer scalar\n*          The number of columns of T.  NH.GE.NW.\n*\n*     T       (workspace) COMPLEX array, dimension (LDT,NW)\n*\n*     LDT     (input) integer\n*          The leading dimension of T just as declared in the\n*          calling subroutine.  NW .LE. LDT\n*\n*     NV      (input) integer\n*          The number of rows of work array WV available for\n*          workspace.  NV.GE.NW.\n*\n*     WV      (workspace) COMPLEX array, dimension (LDWV,NW)\n*\n*     LDWV    (input) integer\n*          The leading dimension of W just as declared in the\n*          calling subroutine.  NW .LE. LDV\n*\n*     WORK    (workspace) COMPLEX array, dimension LWORK.\n*          On exit, WORK(1) is set to an estimate of the optimal value\n*          of LWORK for the given values of N, NW, KTOP and KBOT.\n*\n*     LWORK   (input) integer\n*          The dimension of the work array WORK.  LWORK = 2*NW\n*          suffices, but greater efficiency may result from larger\n*          values of LWORK.\n*\n*          If LWORK = -1, then a workspace query is assumed; CLAQR3\n*          only estimates the optimal workspace size for the given\n*          values of N, NW, KTOP and KBOT.  The estimate is returned\n*          in WORK(1).  No error message related to LWORK is issued\n*          by XERBLA.  Neither H nor Z are accessed.\n*\n\n*     ================================================================\n*     Based on contributions by\n*        Karen Braman and Ralph Byers, Department of Mathematics,\n*        University of Kansas, USA\n*\n*     ================================================================\n\n");
    return Qnil;
  }
  if (argc != 14)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 14)", argc);
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
  rb_ldt = argv[10];
  rb_nv = argv[11];
  rb_ldwv = argv[12];
  rb_lwork = argv[13];

  wantt = (rb_wantt == Qtrue);
  wantz = (rb_wantz == Qtrue);
  ktop = NUM2INT(rb_ktop);
  kbot = NUM2INT(rb_kbot);
  nw = NUM2INT(rb_nw);
  iloz = NUM2INT(rb_iloz);
  ihiz = NUM2INT(rb_ihiz);
  nh = NUM2INT(rb_nh);
  ldt = NUM2INT(rb_ldt);
  nv = NUM2INT(rb_nv);
  ldwv = NUM2INT(rb_ldwv);
  lwork = NUM2INT(rb_lwork);
  if (!NA_IsNArray(rb_h))
    rb_raise(rb_eArgError, "h (6th argument) must be NArray");
  if (NA_RANK(rb_h) != 2)
    rb_raise(rb_eArgError, "rank of h (6th argument) must be %d", 2);
  ldh = NA_SHAPE0(rb_h);
  n = NA_SHAPE1(rb_h);
  if (NA_TYPE(rb_h) != NA_SCOMPLEX)
    rb_h = na_change_type(rb_h, NA_SCOMPLEX);
  h = NA_PTR_TYPE(rb_h, complex*);
  if (!NA_IsNArray(rb_z))
    rb_raise(rb_eArgError, "z (9th argument) must be NArray");
  if (NA_RANK(rb_z) != 2)
    rb_raise(rb_eArgError, "rank of z (9th argument) must be %d", 2);
  ldz = NA_SHAPE0(rb_z);
  if (NA_SHAPE1(rb_z) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of z must be the same as shape 1 of h");
  if (NA_TYPE(rb_z) != NA_SCOMPLEX)
    rb_z = na_change_type(rb_z, NA_SCOMPLEX);
  z = NA_PTR_TYPE(rb_z, complex*);
  {
    int shape[1];
    shape[0] = MAX(1,kbot);
    rb_sh = na_make_object(NA_SCOMPLEX, 1, shape, cNArray);
  }
  sh = NA_PTR_TYPE(rb_sh, complex*);
  {
    int shape[2];
    shape[0] = ldh;
    shape[1] = n;
    rb_h_out__ = na_make_object(NA_SCOMPLEX, 2, shape, cNArray);
  }
  h_out__ = NA_PTR_TYPE(rb_h_out__, complex*);
  MEMCPY(h_out__, h, complex, NA_TOTAL(rb_h));
  rb_h = rb_h_out__;
  h = h_out__;
  {
    int shape[2];
    shape[0] = ldz;
    shape[1] = n;
    rb_z_out__ = na_make_object(NA_SCOMPLEX, 2, shape, cNArray);
  }
  z_out__ = NA_PTR_TYPE(rb_z_out__, complex*);
  MEMCPY(z_out__, z, complex, NA_TOTAL(rb_z));
  rb_z = rb_z_out__;
  z = z_out__;
  ldv = nw;
  v = ALLOC_N(complex, (ldv)*(MAX(1,nw)));
  ldt = nw;
  t = ALLOC_N(complex, (ldv)*(MAX(1,nw)));
  ldwv = nw;
  wv = ALLOC_N(complex, (ldv)*(MAX(1,nw)));
  work = ALLOC_N(complex, (MAX(1,lwork)));

  claqr3_(&wantt, &wantz, &n, &ktop, &kbot, &nw, h, &ldh, &iloz, &ihiz, z, &ldz, &ns, &nd, sh, v, &ldv, &nh, t, &ldt, &nv, wv, &ldwv, work, &lwork);

  free(v);
  free(t);
  free(wv);
  free(work);
  rb_ns = INT2NUM(ns);
  rb_nd = INT2NUM(nd);
  return rb_ary_new3(5, rb_ns, rb_nd, rb_sh, rb_h, rb_z);
}

void
init_lapack_claqr3(VALUE mLapack){
  rb_define_module_function(mLapack, "claqr3", rb_claqr3, -1);
}

#include "rb_lapack.h"

extern VOID zlaqr4_(logical *wantt, logical *wantz, integer *n, integer *ilo, integer *ihi, doublecomplex *h, integer *ldh, doublecomplex *w, integer *iloz, integer *ihiz, doublecomplex *z, integer *ldz, doublecomplex *work, integer *lwork, integer *info);

static VALUE sHelp, sUsage;

static VALUE
rblapack_zlaqr4(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_wantt;
  logical wantt; 
  VALUE rblapack_wantz;
  logical wantz; 
  VALUE rblapack_ilo;
  integer ilo; 
  VALUE rblapack_h;
  doublecomplex *h; 
  VALUE rblapack_iloz;
  integer iloz; 
  VALUE rblapack_ihiz;
  integer ihiz; 
  VALUE rblapack_z;
  doublecomplex *z; 
  VALUE rblapack_lwork;
  integer lwork; 
  VALUE rblapack_w;
  doublecomplex *w; 
  VALUE rblapack_work;
  doublecomplex *work; 
  VALUE rblapack_info;
  integer info; 
  VALUE rblapack_h_out__;
  doublecomplex *h_out__;
  VALUE rblapack_z_out__;
  doublecomplex *z_out__;

  integer ldh;
  integer n;
  integer ldz;
  integer ihi;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  w, work, info, h, z = NumRu::Lapack.zlaqr4( wantt, wantz, ilo, h, iloz, ihiz, z, lwork, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE ZLAQR4( WANTT, WANTZ, N, ILO, IHI, H, LDH, W, ILOZ, IHIZ, Z, LDZ, WORK, LWORK, INFO )\n\n*     Purpose\n*     =======\n*\n*     ZLAQR4 computes the eigenvalues of a Hessenberg matrix H\n*     and, optionally, the matrices T and Z from the Schur decomposition\n*     H = Z T Z**H, where T is an upper triangular matrix (the\n*     Schur form), and Z is the unitary matrix of Schur vectors.\n*\n*     Optionally Z may be postmultiplied into an input unitary\n*     matrix Q so that this routine can give the Schur factorization\n*     of a matrix A which has been reduced to the Hessenberg form H\n*     by the unitary matrix Q:  A = Q*H*Q**H = (QZ)*H*(QZ)**H.\n*\n\n*     Arguments\n*     =========\n*\n*     WANTT   (input) LOGICAL\n*          = .TRUE. : the full Schur form T is required;\n*          = .FALSE.: only eigenvalues are required.\n*\n*     WANTZ   (input) LOGICAL\n*          = .TRUE. : the matrix of Schur vectors Z is required;\n*          = .FALSE.: Schur vectors are not required.\n*\n*     N     (input) INTEGER\n*           The order of the matrix H.  N .GE. 0.\n*\n*     ILO   (input) INTEGER\n*     IHI   (input) INTEGER\n*           It is assumed that H is already upper triangular in rows\n*           and columns 1:ILO-1 and IHI+1:N and, if ILO.GT.1,\n*           H(ILO,ILO-1) is zero. ILO and IHI are normally set by a\n*           previous call to ZGEBAL, and then passed to ZGEHRD when the\n*           matrix output by ZGEBAL is reduced to Hessenberg form.\n*           Otherwise, ILO and IHI should be set to 1 and N,\n*           respectively.  If N.GT.0, then 1.LE.ILO.LE.IHI.LE.N.\n*           If N = 0, then ILO = 1 and IHI = 0.\n*\n*     H     (input/output) COMPLEX*16 array, dimension (LDH,N)\n*           On entry, the upper Hessenberg matrix H.\n*           On exit, if INFO = 0 and WANTT is .TRUE., then H\n*           contains the upper triangular matrix T from the Schur\n*           decomposition (the Schur form). If INFO = 0 and WANT is\n*           .FALSE., then the contents of H are unspecified on exit.\n*           (The output value of H when INFO.GT.0 is given under the\n*           description of INFO below.)\n*\n*           This subroutine may explicitly set H(i,j) = 0 for i.GT.j and\n*           j = 1, 2, ... ILO-1 or j = IHI+1, IHI+2, ... N.\n*\n*     LDH   (input) INTEGER\n*           The leading dimension of the array H. LDH .GE. max(1,N).\n*\n*     W        (output) COMPLEX*16 array, dimension (N)\n*           The computed eigenvalues of H(ILO:IHI,ILO:IHI) are stored\n*           in W(ILO:IHI). If WANTT is .TRUE., then the eigenvalues are\n*           stored in the same order as on the diagonal of the Schur\n*           form returned in H, with W(i) = H(i,i).\n*\n*     Z     (input/output) COMPLEX*16 array, dimension (LDZ,IHI)\n*           If WANTZ is .FALSE., then Z is not referenced.\n*           If WANTZ is .TRUE., then Z(ILO:IHI,ILOZ:IHIZ) is\n*           replaced by Z(ILO:IHI,ILOZ:IHIZ)*U where U is the\n*           orthogonal Schur factor of H(ILO:IHI,ILO:IHI).\n*           (The output value of Z when INFO.GT.0 is given under\n*           the description of INFO below.)\n*\n*     LDZ   (input) INTEGER\n*           The leading dimension of the array Z.  if WANTZ is .TRUE.\n*           then LDZ.GE.MAX(1,IHIZ).  Otherwize, LDZ.GE.1.\n*\n*     WORK  (workspace/output) COMPLEX*16 array, dimension LWORK\n*           On exit, if LWORK = -1, WORK(1) returns an estimate of\n*           the optimal value for LWORK.\n*\n*     LWORK (input) INTEGER\n*           The dimension of the array WORK.  LWORK .GE. max(1,N)\n*           is sufficient, but LWORK typically as large as 6*N may\n*           be required for optimal performance.  A workspace query\n*           to determine the optimal workspace size is recommended.\n*\n*           If LWORK = -1, then ZLAQR4 does a workspace query.\n*           In this case, ZLAQR4 checks the input parameters and\n*           estimates the optimal workspace size for the given\n*           values of N, ILO and IHI.  The estimate is returned\n*           in WORK(1).  No error message related to LWORK is\n*           issued by XERBLA.  Neither H nor Z are accessed.\n*\n*\n*     INFO  (output) INTEGER\n*             =  0:  successful exit\n*           .GT. 0:  if INFO = i, ZLAQR4 failed to compute all of\n*                the eigenvalues.  Elements 1:ilo-1 and i+1:n of WR\n*                and WI contain those eigenvalues which have been\n*                successfully computed.  (Failures are rare.)\n*\n*                If INFO .GT. 0 and WANT is .FALSE., then on exit,\n*                the remaining unconverged eigenvalues are the eigen-\n*                values of the upper Hessenberg matrix rows and\n*                columns ILO through INFO of the final, output\n*                value of H.\n*\n*                If INFO .GT. 0 and WANTT is .TRUE., then on exit\n*\n*           (*)  (initial value of H)*U  = U*(final value of H)\n*\n*                where U is a unitary matrix.  The final\n*                value of  H is upper Hessenberg and triangular in\n*                rows and columns INFO+1 through IHI.\n*\n*                If INFO .GT. 0 and WANTZ is .TRUE., then on exit\n*\n*                  (final value of Z(ILO:IHI,ILOZ:IHIZ)\n*                   =  (initial value of Z(ILO:IHI,ILOZ:IHIZ)*U\n*\n*                where U is the unitary matrix in (*) (regard-\n*                less of the value of WANTT.)\n*\n*                If INFO .GT. 0 and WANTZ is .FALSE., then Z is not\n*                accessed.\n*\n\n*     ================================================================\n*     Based on contributions by\n*        Karen Braman and Ralph Byers, Department of Mathematics,\n*        University of Kansas, USA\n*\n*     ================================================================\n*     References:\n*       K. Braman, R. Byers and R. Mathias, The Multi-Shift QR\n*       Algorithm Part I: Maintaining Well Focused Shifts, and Level 3\n*       Performance, SIAM Journal of Matrix Analysis, volume 23, pages\n*       929--947, 2002.\n*\n*       K. Braman, R. Byers and R. Mathias, The Multi-Shift QR\n*       Algorithm Part II: Aggressive Early Deflation, SIAM Journal\n*       of Matrix Analysis, volume 23, pages 948--973, 2002.\n*\n*     ================================================================\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  w, work, info, h, z = NumRu::Lapack.zlaqr4( wantt, wantz, ilo, h, iloz, ihiz, z, lwork, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 8)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 8)", argc);
  rblapack_wantt = argv[0];
  rblapack_wantz = argv[1];
  rblapack_ilo = argv[2];
  rblapack_h = argv[3];
  rblapack_iloz = argv[4];
  rblapack_ihiz = argv[5];
  rblapack_z = argv[6];
  rblapack_lwork = argv[7];
  if (rb_options != Qnil) {
  }

  ilo = NUM2INT(rblapack_ilo);
  wantz = (rblapack_wantz == Qtrue);
  if (!NA_IsNArray(rblapack_h))
    rb_raise(rb_eArgError, "h (4th argument) must be NArray");
  if (NA_RANK(rblapack_h) != 2)
    rb_raise(rb_eArgError, "rank of h (4th argument) must be %d", 2);
  n = NA_SHAPE1(rblapack_h);
  ldh = NA_SHAPE0(rblapack_h);
  if (NA_TYPE(rblapack_h) != NA_DCOMPLEX)
    rblapack_h = na_change_type(rblapack_h, NA_DCOMPLEX);
  h = NA_PTR_TYPE(rblapack_h, doublecomplex*);
  if (!NA_IsNArray(rblapack_z))
    rb_raise(rb_eArgError, "z (7th argument) must be NArray");
  if (NA_RANK(rblapack_z) != 2)
    rb_raise(rb_eArgError, "rank of z (7th argument) must be %d", 2);
  ihi = NA_SHAPE1(rblapack_z);
  ldz = NA_SHAPE0(rblapack_z);
  if (ldz != (wantz ? MAX(1,ihiz) : 1))
    rb_raise(rb_eRuntimeError, "shape 0 of z must be %d", wantz ? MAX(1,ihiz) : 1);
  if (NA_TYPE(rblapack_z) != NA_DCOMPLEX)
    rblapack_z = na_change_type(rblapack_z, NA_DCOMPLEX);
  z = NA_PTR_TYPE(rblapack_z, doublecomplex*);
  ihiz = NUM2INT(rblapack_ihiz);
  iloz = NUM2INT(rblapack_iloz);
  wantt = (rblapack_wantt == Qtrue);
  lwork = NUM2INT(rblapack_lwork);
  ldz = wantz ? MAX(1,ihiz) : 1;
  {
    int shape[1];
    shape[0] = n;
    rblapack_w = na_make_object(NA_DCOMPLEX, 1, shape, cNArray);
  }
  w = NA_PTR_TYPE(rblapack_w, doublecomplex*);
  {
    int shape[1];
    shape[0] = MAX(1,lwork);
    rblapack_work = na_make_object(NA_DCOMPLEX, 1, shape, cNArray);
  }
  work = NA_PTR_TYPE(rblapack_work, doublecomplex*);
  {
    int shape[2];
    shape[0] = ldh;
    shape[1] = n;
    rblapack_h_out__ = na_make_object(NA_DCOMPLEX, 2, shape, cNArray);
  }
  h_out__ = NA_PTR_TYPE(rblapack_h_out__, doublecomplex*);
  MEMCPY(h_out__, h, doublecomplex, NA_TOTAL(rblapack_h));
  rblapack_h = rblapack_h_out__;
  h = h_out__;
  {
    int shape[2];
    shape[0] = ldz;
    shape[1] = ihi;
    rblapack_z_out__ = na_make_object(NA_DCOMPLEX, 2, shape, cNArray);
  }
  z_out__ = NA_PTR_TYPE(rblapack_z_out__, doublecomplex*);
  MEMCPY(z_out__, z, doublecomplex, NA_TOTAL(rblapack_z));
  rblapack_z = rblapack_z_out__;
  z = z_out__;

  zlaqr4_(&wantt, &wantz, &n, &ilo, &ihi, h, &ldh, w, &iloz, &ihiz, z, &ldz, work, &lwork, &info);

  rblapack_info = INT2NUM(info);
  return rb_ary_new3(5, rblapack_w, rblapack_work, rblapack_info, rblapack_h, rblapack_z);
}

void
init_lapack_zlaqr4(VALUE mLapack){
  rb_define_module_function(mLapack, "zlaqr4", rblapack_zlaqr4, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}

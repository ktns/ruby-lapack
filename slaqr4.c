#include "rb_lapack.h"

extern VOID slaqr4_(logical *wantt, logical *wantz, integer *n, integer *ilo, integer *ihi, real *h, integer *ldh, real *wr, real *wi, integer *iloz, integer *ihiz, real *z, integer *ldz, real *work, integer *lwork, integer *info);

static VALUE
rb_slaqr4(int argc, VALUE *argv, VALUE self){
  VALUE rb_wantt;
  logical wantt; 
  VALUE rb_wantz;
  logical wantz; 
  VALUE rb_ilo;
  integer ilo; 
  VALUE rb_h;
  real *h; 
  VALUE rb_iloz;
  integer iloz; 
  VALUE rb_ihiz;
  integer ihiz; 
  VALUE rb_z;
  real *z; 
  VALUE rb_lwork;
  integer lwork; 
  VALUE rb_wr;
  real *wr; 
  VALUE rb_wi;
  real *wi; 
  VALUE rb_work;
  real *work; 
  VALUE rb_info;
  integer info; 
  VALUE rb_h_out__;
  real *h_out__;
  VALUE rb_z_out__;
  real *z_out__;

  integer ldh;
  integer n;
  integer ldz;
  integer ihi;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  wr, wi, work, info, h, z = NumRu::Lapack.slaqr4( wantt, wantz, ilo, h, iloz, ihiz, z, lwork)\n    or\n  NumRu::Lapack.slaqr4  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE SLAQR4( WANTT, WANTZ, N, ILO, IHI, H, LDH, WR, WI, ILOZ, IHIZ, Z, LDZ, WORK, LWORK, INFO )\n\n*     Purpose\n*     =======\n*\n*     SLAQR4 computes the eigenvalues of a Hessenberg matrix H\n*     and, optionally, the matrices T and Z from the Schur decomposition\n*     H = Z T Z**T, where T is an upper quasi-triangular matrix (the\n*     Schur form), and Z is the orthogonal matrix of Schur vectors.\n*\n*     Optionally Z may be postmultiplied into an input orthogonal\n*     matrix Q so that this routine can give the Schur factorization\n*     of a matrix A which has been reduced to the Hessenberg form H\n*     by the orthogonal matrix Q:  A = Q*H*Q**T = (QZ)*T*(QZ)**T.\n*\n\n*     Arguments\n*     =========\n*\n*     WANTT   (input) LOGICAL\n*          = .TRUE. : the full Schur form T is required;\n*          = .FALSE.: only eigenvalues are required.\n*\n*     WANTZ   (input) LOGICAL\n*          = .TRUE. : the matrix of Schur vectors Z is required;\n*          = .FALSE.: Schur vectors are not required.\n*\n*     N     (input) INTEGER\n*           The order of the matrix H.  N .GE. 0.\n*\n*     ILO   (input) INTEGER\n*     IHI   (input) INTEGER\n*           It is assumed that H is already upper triangular in rows\n*           and columns 1:ILO-1 and IHI+1:N and, if ILO.GT.1,\n*           H(ILO,ILO-1) is zero. ILO and IHI are normally set by a\n*           previous call to SGEBAL, and then passed to SGEHRD when the\n*           matrix output by SGEBAL is reduced to Hessenberg form.\n*           Otherwise, ILO and IHI should be set to 1 and N,\n*           respectively.  If N.GT.0, then 1.LE.ILO.LE.IHI.LE.N.\n*           If N = 0, then ILO = 1 and IHI = 0.\n*\n*     H     (input/output) REAL array, dimension (LDH,N)\n*           On entry, the upper Hessenberg matrix H.\n*           On exit, if INFO = 0 and WANTT is .TRUE., then H contains\n*           the upper quasi-triangular matrix T from the Schur\n*           decomposition (the Schur form); 2-by-2 diagonal blocks\n*           (corresponding to complex conjugate pairs of eigenvalues)\n*           are returned in standard form, with H(i,i) = H(i+1,i+1)\n*           and H(i+1,i)*H(i,i+1).LT.0. If INFO = 0 and WANTT is\n*           .FALSE., then the contents of H are unspecified on exit.\n*           (The output value of H when INFO.GT.0 is given under the\n*           description of INFO below.)\n*\n*           This subroutine may explicitly set H(i,j) = 0 for i.GT.j and\n*           j = 1, 2, ... ILO-1 or j = IHI+1, IHI+2, ... N.\n*\n*     LDH   (input) INTEGER\n*           The leading dimension of the array H. LDH .GE. max(1,N).\n*\n*     WR    (output) REAL array, dimension (IHI)\n*     WI    (output) REAL array, dimension (IHI)\n*           The real and imaginary parts, respectively, of the computed\n*           eigenvalues of H(ILO:IHI,ILO:IHI) are stored in WR(ILO:IHI)\n*           and WI(ILO:IHI). If two eigenvalues are computed as a\n*           complex conjugate pair, they are stored in consecutive\n*           elements of WR and WI, say the i-th and (i+1)th, with\n*           WI(i) .GT. 0 and WI(i+1) .LT. 0. If WANTT is .TRUE., then\n*           the eigenvalues are stored in the same order as on the\n*           diagonal of the Schur form returned in H, with\n*           WR(i) = H(i,i) and, if H(i:i+1,i:i+1) is a 2-by-2 diagonal\n*           block, WI(i) = sqrt(-H(i+1,i)*H(i,i+1)) and\n*           WI(i+1) = -WI(i).\n*\n*     ILOZ     (input) INTEGER\n*     IHIZ     (input) INTEGER\n*           Specify the rows of Z to which transformations must be\n*           applied if WANTZ is .TRUE..\n*           1 .LE. ILOZ .LE. ILO; IHI .LE. IHIZ .LE. N.\n*\n*     Z     (input/output) REAL array, dimension (LDZ,IHI)\n*           If WANTZ is .FALSE., then Z is not referenced.\n*           If WANTZ is .TRUE., then Z(ILO:IHI,ILOZ:IHIZ) is\n*           replaced by Z(ILO:IHI,ILOZ:IHIZ)*U where U is the\n*           orthogonal Schur factor of H(ILO:IHI,ILO:IHI).\n*           (The output value of Z when INFO.GT.0 is given under\n*           the description of INFO below.)\n*\n*     LDZ   (input) INTEGER\n*           The leading dimension of the array Z.  if WANTZ is .TRUE.\n*           then LDZ.GE.MAX(1,IHIZ).  Otherwize, LDZ.GE.1.\n*\n*     WORK  (workspace/output) REAL array, dimension LWORK\n*           On exit, if LWORK = -1, WORK(1) returns an estimate of\n*           the optimal value for LWORK.\n*\n*     LWORK (input) INTEGER\n*           The dimension of the array WORK.  LWORK .GE. max(1,N)\n*           is sufficient, but LWORK typically as large as 6*N may\n*           be required for optimal performance.  A workspace query\n*           to determine the optimal workspace size is recommended.\n*\n*           If LWORK = -1, then SLAQR4 does a workspace query.\n*           In this case, SLAQR4 checks the input parameters and\n*           estimates the optimal workspace size for the given\n*           values of N, ILO and IHI.  The estimate is returned\n*           in WORK(1).  No error message related to LWORK is\n*           issued by XERBLA.  Neither H nor Z are accessed.\n*\n*\n*     INFO  (output) INTEGER\n*             =  0:  successful exit\n*           .GT. 0:  if INFO = i, SLAQR4 failed to compute all of\n*                the eigenvalues.  Elements 1:ilo-1 and i+1:n of WR\n*                and WI contain those eigenvalues which have been\n*                successfully computed.  (Failures are rare.)\n*\n*                If INFO .GT. 0 and WANT is .FALSE., then on exit,\n*                the remaining unconverged eigenvalues are the eigen-\n*                values of the upper Hessenberg matrix rows and\n*                columns ILO through INFO of the final, output\n*                value of H.\n*\n*                If INFO .GT. 0 and WANTT is .TRUE., then on exit\n*\n*           (*)  (initial value of H)*U  = U*(final value of H)\n*\n*                where U is an orthogonal matrix.  The final\n*                value of H is upper Hessenberg and quasi-triangular\n*                in rows and columns INFO+1 through IHI.\n*\n*                If INFO .GT. 0 and WANTZ is .TRUE., then on exit\n*\n*                  (final value of Z(ILO:IHI,ILOZ:IHIZ)\n*                   =  (initial value of Z(ILO:IHI,ILOZ:IHIZ)*U\n*\n*                where U is the orthogonal matrix in (*) (regard-\n*                less of the value of WANTT.)\n*\n*                If INFO .GT. 0 and WANTZ is .FALSE., then Z is not\n*                accessed.\n*\n\n*     ================================================================\n*     Based on contributions by\n*        Karen Braman and Ralph Byers, Department of Mathematics,\n*        University of Kansas, USA\n*\n*     ================================================================\n*     References:\n*       K. Braman, R. Byers and R. Mathias, The Multi-Shift QR\n*       Algorithm Part I: Maintaining Well Focused Shifts, and Level 3\n*       Performance, SIAM Journal of Matrix Analysis, volume 23, pages\n*       929--947, 2002.\n*\n*       K. Braman, R. Byers and R. Mathias, The Multi-Shift QR\n*       Algorithm Part II: Aggressive Early Deflation, SIAM Journal\n*       of Matrix Analysis, volume 23, pages 948--973, 2002.\n*\n*     ================================================================\n\n");
    return Qnil;
  }
  if (argc != 8)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 8)", argc);
  rb_wantt = argv[0];
  rb_wantz = argv[1];
  rb_ilo = argv[2];
  rb_h = argv[3];
  rb_iloz = argv[4];
  rb_ihiz = argv[5];
  rb_z = argv[6];
  rb_lwork = argv[7];

  ilo = NUM2INT(rb_ilo);
  wantz = (rb_wantz == Qtrue);
  if (!NA_IsNArray(rb_z))
    rb_raise(rb_eArgError, "z (7th argument) must be NArray");
  if (NA_RANK(rb_z) != 2)
    rb_raise(rb_eArgError, "rank of z (7th argument) must be %d", 2);
  ihi = NA_SHAPE1(rb_z);
  ldz = NA_SHAPE0(rb_z);
  if (NA_TYPE(rb_z) != NA_SFLOAT)
    rb_z = na_change_type(rb_z, NA_SFLOAT);
  z = NA_PTR_TYPE(rb_z, real*);
  lwork = NUM2INT(rb_lwork);
  iloz = NUM2INT(rb_iloz);
  wantt = (rb_wantt == Qtrue);
  ihiz = NUM2INT(rb_ihiz);
  if (!NA_IsNArray(rb_h))
    rb_raise(rb_eArgError, "h (4th argument) must be NArray");
  if (NA_RANK(rb_h) != 2)
    rb_raise(rb_eArgError, "rank of h (4th argument) must be %d", 2);
  n = NA_SHAPE1(rb_h);
  ldh = NA_SHAPE0(rb_h);
  if (NA_TYPE(rb_h) != NA_SFLOAT)
    rb_h = na_change_type(rb_h, NA_SFLOAT);
  h = NA_PTR_TYPE(rb_h, real*);
  {
    int shape[1];
    shape[0] = ihi;
    rb_wr = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  wr = NA_PTR_TYPE(rb_wr, real*);
  {
    int shape[1];
    shape[0] = ihi;
    rb_wi = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  wi = NA_PTR_TYPE(rb_wi, real*);
  {
    int shape[1];
    shape[0] = MAX(1,lwork);
    rb_work = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  work = NA_PTR_TYPE(rb_work, real*);
  {
    int shape[2];
    shape[0] = ldh;
    shape[1] = n;
    rb_h_out__ = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  h_out__ = NA_PTR_TYPE(rb_h_out__, real*);
  MEMCPY(h_out__, h, real, NA_TOTAL(rb_h));
  rb_h = rb_h_out__;
  h = h_out__;
  {
    int shape[2];
    shape[0] = ldz;
    shape[1] = ihi;
    rb_z_out__ = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  z_out__ = NA_PTR_TYPE(rb_z_out__, real*);
  MEMCPY(z_out__, z, real, NA_TOTAL(rb_z));
  rb_z = rb_z_out__;
  z = z_out__;

  slaqr4_(&wantt, &wantz, &n, &ilo, &ihi, h, &ldh, wr, wi, &iloz, &ihiz, z, &ldz, work, &lwork, &info);

  rb_info = INT2NUM(info);
  return rb_ary_new3(6, rb_wr, rb_wi, rb_work, rb_info, rb_h, rb_z);
}

void
init_lapack_slaqr4(VALUE mLapack){
  rb_define_module_function(mLapack, "slaqr4", rb_slaqr4, -1);
}

#include "rb_lapack.h"

extern VOID chseqr_(char *job, char *compz, integer *n, integer *ilo, integer *ihi, complex *h, integer *ldh, complex *w, complex *z, integer *ldz, complex *work, integer *lwork, integer *info);

static VALUE
rb_chseqr(int argc, VALUE *argv, VALUE self){
  VALUE rb_job;
  char job; 
  VALUE rb_compz;
  char compz; 
  VALUE rb_ilo;
  integer ilo; 
  VALUE rb_ihi;
  integer ihi; 
  VALUE rb_h;
  complex *h; 
  VALUE rb_z;
  complex *z; 
  VALUE rb_ldz;
  integer ldz; 
  VALUE rb_lwork;
  integer lwork; 
  VALUE rb_w;
  complex *w; 
  VALUE rb_work;
  complex *work; 
  VALUE rb_info;
  integer info; 
  VALUE rb_h_out__;
  complex *h_out__;
  VALUE rb_z_out__;
  complex *z_out__;

  integer ldh;
  integer n;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  w, work, info, h, z = NumRu::Lapack.chseqr( job, compz, ilo, ihi, h, z, ldz, lwork)\n    or\n  NumRu::Lapack.chseqr  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE CHSEQR( JOB, COMPZ, N, ILO, IHI, H, LDH, W, Z, LDZ, WORK, LWORK, INFO )\n\n*     Purpose\n*     =======\n*\n*     CHSEQR computes the eigenvalues of a Hessenberg matrix H\n*     and, optionally, the matrices T and Z from the Schur decomposition\n*     H = Z T Z**H, where T is an upper triangular matrix (the\n*     Schur form), and Z is the unitary matrix of Schur vectors.\n*\n*     Optionally Z may be postmultiplied into an input unitary\n*     matrix Q so that this routine can give the Schur factorization\n*     of a matrix A which has been reduced to the Hessenberg form H\n*     by the unitary matrix Q:  A = Q*H*Q**H = (QZ)*H*(QZ)**H.\n*\n\n*     Arguments\n*     =========\n*\n*     JOB   (input) CHARACTER*1\n*           = 'E':  compute eigenvalues only;\n*           = 'S':  compute eigenvalues and the Schur form T.\n*\n*     COMPZ (input) CHARACTER*1\n*           = 'N':  no Schur vectors are computed;\n*           = 'I':  Z is initialized to the unit matrix and the matrix Z\n*                   of Schur vectors of H is returned;\n*           = 'V':  Z must contain an unitary matrix Q on entry, and\n*                   the product Q*Z is returned.\n*\n*     N     (input) INTEGER\n*           The order of the matrix H.  N .GE. 0.\n*\n*     ILO   (input) INTEGER\n*     IHI   (input) INTEGER\n*           It is assumed that H is already upper triangular in rows\n*           and columns 1:ILO-1 and IHI+1:N. ILO and IHI are normally\n*           set by a previous call to CGEBAL, and then passed to CGEHRD\n*           when the matrix output by CGEBAL is reduced to Hessenberg\n*           form. Otherwise ILO and IHI should be set to 1 and N\n*           respectively.  If N.GT.0, then 1.LE.ILO.LE.IHI.LE.N.\n*           If N = 0, then ILO = 1 and IHI = 0.\n*\n*     H     (input/output) COMPLEX array, dimension (LDH,N)\n*           On entry, the upper Hessenberg matrix H.\n*           On exit, if INFO = 0 and JOB = 'S', H contains the upper\n*           triangular matrix T from the Schur decomposition (the\n*           Schur form). If INFO = 0 and JOB = 'E', the contents of\n*           H are unspecified on exit.  (The output value of H when\n*           INFO.GT.0 is given under the description of INFO below.)\n*\n*           Unlike earlier versions of CHSEQR, this subroutine may\n*           explicitly H(i,j) = 0 for i.GT.j and j = 1, 2, ... ILO-1\n*           or j = IHI+1, IHI+2, ... N.\n*\n*     LDH   (input) INTEGER\n*           The leading dimension of the array H. LDH .GE. max(1,N).\n*\n*     W        (output) COMPLEX array, dimension (N)\n*           The computed eigenvalues. If JOB = 'S', the eigenvalues are\n*           stored in the same order as on the diagonal of the Schur\n*           form returned in H, with W(i) = H(i,i).\n*\n*     Z     (input/output) COMPLEX array, dimension (LDZ,N)\n*           If COMPZ = 'N', Z is not referenced.\n*           If COMPZ = 'I', on entry Z need not be set and on exit,\n*           if INFO = 0, Z contains the unitary matrix Z of the Schur\n*           vectors of H.  If COMPZ = 'V', on entry Z must contain an\n*           N-by-N matrix Q, which is assumed to be equal to the unit\n*           matrix except for the submatrix Z(ILO:IHI,ILO:IHI). On exit,\n*           if INFO = 0, Z contains Q*Z.\n*           Normally Q is the unitary matrix generated by CUNGHR\n*           after the call to CGEHRD which formed the Hessenberg matrix\n*           H. (The output value of Z when INFO.GT.0 is given under\n*           the description of INFO below.)\n*\n*     LDZ   (input) INTEGER\n*           The leading dimension of the array Z.  if COMPZ = 'I' or\n*           COMPZ = 'V', then LDZ.GE.MAX(1,N).  Otherwize, LDZ.GE.1.\n*\n*     WORK  (workspace/output) COMPLEX array, dimension (LWORK)\n*           On exit, if INFO = 0, WORK(1) returns an estimate of\n*           the optimal value for LWORK.\n*\n*     LWORK (input) INTEGER\n*           The dimension of the array WORK.  LWORK .GE. max(1,N)\n*           is sufficient and delivers very good and sometimes\n*           optimal performance.  However, LWORK as large as 11*N\n*           may be required for optimal performance.  A workspace\n*           query is recommended to determine the optimal workspace\n*           size.\n*\n*           If LWORK = -1, then CHSEQR does a workspace query.\n*           In this case, CHSEQR checks the input parameters and\n*           estimates the optimal workspace size for the given\n*           values of N, ILO and IHI.  The estimate is returned\n*           in WORK(1).  No error message related to LWORK is\n*           issued by XERBLA.  Neither H nor Z are accessed.\n*\n*\n*     INFO  (output) INTEGER\n*             =  0:  successful exit\n*           .LT. 0:  if INFO = -i, the i-th argument had an illegal\n*                    value\n*           .GT. 0:  if INFO = i, CHSEQR failed to compute all of\n*                the eigenvalues.  Elements 1:ilo-1 and i+1:n of WR\n*                and WI contain those eigenvalues which have been\n*                successfully computed.  (Failures are rare.)\n*\n*                If INFO .GT. 0 and JOB = 'E', then on exit, the\n*                remaining unconverged eigenvalues are the eigen-\n*                values of the upper Hessenberg matrix rows and\n*                columns ILO through INFO of the final, output\n*                value of H.\n*\n*                If INFO .GT. 0 and JOB   = 'S', then on exit\n*\n*           (*)  (initial value of H)*U  = U*(final value of H)\n*\n*                where U is a unitary matrix.  The final\n*                value of  H is upper Hessenberg and triangular in\n*                rows and columns INFO+1 through IHI.\n*\n*                If INFO .GT. 0 and COMPZ = 'V', then on exit\n*\n*                  (final value of Z)  =  (initial value of Z)*U\n*\n*                where U is the unitary matrix in (*) (regard-\n*                less of the value of JOB.)\n*\n*                If INFO .GT. 0 and COMPZ = 'I', then on exit\n*                      (final value of Z)  = U\n*                where U is the unitary matrix in (*) (regard-\n*                less of the value of JOB.)\n*\n*                If INFO .GT. 0 and COMPZ = 'N', then Z is not\n*                accessed.\n*\n\n*     ================================================================\n*             Default values supplied by\n*             ILAENV(ISPEC,'CHSEQR',JOB(:1)//COMPZ(:1),N,ILO,IHI,LWORK).\n*             It is suggested that these defaults be adjusted in order\n*             to attain best performance in each particular\n*             computational environment.\n*\n*            ISPEC=12: The CLAHQR vs CLAQR0 crossover point.\n*                      Default: 75. (Must be at least 11.)\n*\n*            ISPEC=13: Recommended deflation window size.\n*                      This depends on ILO, IHI and NS.  NS is the\n*                      number of simultaneous shifts returned\n*                      by ILAENV(ISPEC=15).  (See ISPEC=15 below.)\n*                      The default for (IHI-ILO+1).LE.500 is NS.\n*                      The default for (IHI-ILO+1).GT.500 is 3*NS/2.\n*\n*            ISPEC=14: Nibble crossover point. (See IPARMQ for\n*                      details.)  Default: 14% of deflation window\n*                      size.\n*\n*            ISPEC=15: Number of simultaneous shifts in a multishift\n*                      QR iteration.\n*\n*                      If IHI-ILO+1 is ...\n*\n*                      greater than      ...but less    ... the\n*                      or equal to ...      than        default is\n*\n*                           1               30          NS =   2(+)\n*                          30               60          NS =   4(+)\n*                          60              150          NS =  10(+)\n*                         150              590          NS =  **\n*                         590             3000          NS =  64\n*                        3000             6000          NS = 128\n*                        6000             infinity      NS = 256\n*\n*                  (+)  By default some or all matrices of this order\n*                       are passed to the implicit double shift routine\n*                       CLAHQR and this parameter is ignored.  See\n*                       ISPEC=12 above and comments in IPARMQ for\n*                       details.\n*\n*                 (**)  The asterisks (**) indicate an ad-hoc\n*                       function of N increasing from 10 to 64.\n*\n*            ISPEC=16: Select structured matrix multiply.\n*                      If the number of simultaneous shifts (specified\n*                      by ISPEC=15) is less than 14, then the default\n*                      for ISPEC=16 is 0.  Otherwise the default for\n*                      ISPEC=16 is 2.\n*\n*     ================================================================\n*     Based on contributions by\n*        Karen Braman and Ralph Byers, Department of Mathematics,\n*        University of Kansas, USA\n*\n*     ================================================================\n*     References:\n*       K. Braman, R. Byers and R. Mathias, The Multi-Shift QR\n*       Algorithm Part I: Maintaining Well Focused Shifts, and Level 3\n*       Performance, SIAM Journal of Matrix Analysis, volume 23, pages\n*       929--947, 2002.\n*\n*       K. Braman, R. Byers and R. Mathias, The Multi-Shift QR\n*       Algorithm Part II: Aggressive Early Deflation, SIAM Journal\n*       of Matrix Analysis, volume 23, pages 948--973, 2002.\n*\n*     ================================================================\n\n");
    return Qnil;
  }
  if (argc != 8)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 8)", argc);
  rb_job = argv[0];
  rb_compz = argv[1];
  rb_ilo = argv[2];
  rb_ihi = argv[3];
  rb_h = argv[4];
  rb_z = argv[5];
  rb_ldz = argv[6];
  rb_lwork = argv[7];

  ilo = NUM2INT(rb_ilo);
  ldz = NUM2INT(rb_ldz);
  compz = StringValueCStr(rb_compz)[0];
  ihi = NUM2INT(rb_ihi);
  if (!NA_IsNArray(rb_h))
    rb_raise(rb_eArgError, "h (5th argument) must be NArray");
  if (NA_RANK(rb_h) != 2)
    rb_raise(rb_eArgError, "rank of h (5th argument) must be %d", 2);
  n = NA_SHAPE1(rb_h);
  ldh = NA_SHAPE0(rb_h);
  if (NA_TYPE(rb_h) != NA_SCOMPLEX)
    rb_h = na_change_type(rb_h, NA_SCOMPLEX);
  h = NA_PTR_TYPE(rb_h, complex*);
  job = StringValueCStr(rb_job)[0];
  lwork = NUM2INT(rb_lwork);
  if (!NA_IsNArray(rb_z))
    rb_raise(rb_eArgError, "z (6th argument) must be NArray");
  if (NA_RANK(rb_z) != 2)
    rb_raise(rb_eArgError, "rank of z (6th argument) must be %d", 2);
  if (NA_SHAPE1(rb_z) != (lsame_(&compz,"N") ? 0 : n))
    rb_raise(rb_eRuntimeError, "shape 1 of z must be %d", lsame_(&compz,"N") ? 0 : n);
  if (NA_SHAPE0(rb_z) != (lsame_(&compz,"N") ? 0 : ldz))
    rb_raise(rb_eRuntimeError, "shape 0 of z must be %d", lsame_(&compz,"N") ? 0 : ldz);
  if (NA_TYPE(rb_z) != NA_SCOMPLEX)
    rb_z = na_change_type(rb_z, NA_SCOMPLEX);
  z = NA_PTR_TYPE(rb_z, complex*);
  {
    int shape[1];
    shape[0] = n;
    rb_w = na_make_object(NA_SCOMPLEX, 1, shape, cNArray);
  }
  w = NA_PTR_TYPE(rb_w, complex*);
  {
    int shape[1];
    shape[0] = MAX(1,lwork);
    rb_work = na_make_object(NA_SCOMPLEX, 1, shape, cNArray);
  }
  work = NA_PTR_TYPE(rb_work, complex*);
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
    shape[0] = lsame_(&compz,"N") ? 0 : ldz;
    shape[1] = lsame_(&compz,"N") ? 0 : n;
    rb_z_out__ = na_make_object(NA_SCOMPLEX, 2, shape, cNArray);
  }
  z_out__ = NA_PTR_TYPE(rb_z_out__, complex*);
  MEMCPY(z_out__, z, complex, NA_TOTAL(rb_z));
  rb_z = rb_z_out__;
  z = z_out__;

  chseqr_(&job, &compz, &n, &ilo, &ihi, h, &ldh, w, z, &ldz, work, &lwork, &info);

  rb_info = INT2NUM(info);
  return rb_ary_new3(5, rb_w, rb_work, rb_info, rb_h, rb_z);
}

void
init_lapack_chseqr(VALUE mLapack){
  rb_define_module_function(mLapack, "chseqr", rb_chseqr, -1);
}

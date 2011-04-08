#include "rb_lapack.h"

extern integer iparmq_(integer *ispec, char *name, char *opts, integer *n, integer *ilo, integer *ihi, integer *lwork);

static VALUE sHelp, sUsage;

static VALUE
rblapack_iparmq(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_ispec;
  integer ispec; 
  VALUE rblapack_name;
  char name; 
  VALUE rblapack_opts;
  char opts; 
  VALUE rblapack_n;
  integer n; 
  VALUE rblapack_ilo;
  integer ilo; 
  VALUE rblapack_ihi;
  integer ihi; 
  VALUE rblapack_lwork;
  integer lwork; 
  VALUE rblapack___out__;
  integer __out__; 


  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  __out__ = NumRu::Lapack.iparmq( ispec, name, opts, n, ilo, ihi, lwork, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      INTEGER FUNCTION IPARMQ( ISPEC, NAME, OPTS, N, ILO, IHI, LWORK )\n\n*  Purpose\n*  =======\n*\n*       This program sets problem and machine dependent parameters\n*       useful for xHSEQR and its subroutines. It is called whenever \n*       ILAENV is called with 12 <= ISPEC <= 16\n*\n\n*  Arguments\n*  =========\n*\n*       ISPEC  (input) integer scalar\n*              ISPEC specifies which tunable parameter IPARMQ should\n*              return.\n*\n*              ISPEC=12: (INMIN)  Matrices of order nmin or less\n*                        are sent directly to xLAHQR, the implicit\n*                        double shift QR algorithm.  NMIN must be\n*                        at least 11.\n*\n*              ISPEC=13: (INWIN)  Size of the deflation window.\n*                        This is best set greater than or equal to\n*                        the number of simultaneous shifts NS.\n*                        Larger matrices benefit from larger deflation\n*                        windows.\n*\n*              ISPEC=14: (INIBL) Determines when to stop nibbling and\n*                        invest in an (expensive) multi-shift QR sweep.\n*                        If the aggressive early deflation subroutine\n*                        finds LD converged eigenvalues from an order\n*                        NW deflation window and LD.GT.(NW*NIBBLE)/100,\n*                        then the next QR sweep is skipped and early\n*                        deflation is applied immediately to the\n*                        remaining active diagonal block.  Setting\n*                        IPARMQ(ISPEC=14) = 0 causes TTQRE to skip a\n*                        multi-shift QR sweep whenever early deflation\n*                        finds a converged eigenvalue.  Setting\n*                        IPARMQ(ISPEC=14) greater than or equal to 100\n*                        prevents TTQRE from skipping a multi-shift\n*                        QR sweep.\n*\n*              ISPEC=15: (NSHFTS) The number of simultaneous shifts in\n*                        a multi-shift QR iteration.\n*\n*              ISPEC=16: (IACC22) IPARMQ is set to 0, 1 or 2 with the\n*                        following meanings.\n*                        0:  During the multi-shift QR sweep,\n*                            xLAQR5 does not accumulate reflections and\n*                            does not use matrix-matrix multiply to\n*                            update the far-from-diagonal matrix\n*                            entries.\n*                        1:  During the multi-shift QR sweep,\n*                            xLAQR5 and/or xLAQRaccumulates reflections and uses\n*                            matrix-matrix multiply to update the\n*                            far-from-diagonal matrix entries.\n*                        2:  During the multi-shift QR sweep.\n*                            xLAQR5 accumulates reflections and takes\n*                            advantage of 2-by-2 block structure during\n*                            matrix-matrix multiplies.\n*                        (If xTRMM is slower than xGEMM, then\n*                        IPARMQ(ISPEC=16)=1 may be more efficient than\n*                        IPARMQ(ISPEC=16)=2 despite the greater level of\n*                        arithmetic work implied by the latter choice.)\n*\n*       NAME    (input) character string\n*               Name of the calling subroutine\n*\n*       OPTS    (input) character string\n*               This is a concatenation of the string arguments to\n*               TTQRE.\n*\n*       N       (input) integer scalar\n*               N is the order of the Hessenberg matrix H.\n*\n*       ILO     (input) INTEGER\n*       IHI     (input) INTEGER\n*               It is assumed that H is already upper triangular\n*               in rows and columns 1:ILO-1 and IHI+1:N.\n*\n*       LWORK   (input) integer scalar\n*               The amount of workspace available.\n*\n\n*  Further Details\n*  ===============\n*\n*       Little is known about how best to choose these parameters.\n*       It is possible to use different values of the parameters\n*       for each of CHSEQR, DHSEQR, SHSEQR and ZHSEQR.\n*\n*       It is probably best to choose different parameters for\n*       different matrices and different parameters at different\n*       times during the iteration, but this has not been\n*       implemented --- yet.\n*\n*\n*       The best choices of most of the parameters depend\n*       in an ill-understood way on the relative execution\n*       rate of xLAQR3 and xLAQR5 and on the nature of each\n*       particular eigenvalue problem.  Experiment may be the\n*       only practical way to determine which choices are most\n*       effective.\n*\n*       Following is a list of default values supplied by IPARMQ.\n*       These defaults may be adjusted in order to attain better\n*       performance in any particular computational environment.\n*\n*       IPARMQ(ISPEC=12) The xLAHQR vs xLAQR0 crossover point.\n*                        Default: 75. (Must be at least 11.)\n*\n*       IPARMQ(ISPEC=13) Recommended deflation window size.\n*                        This depends on ILO, IHI and NS, the\n*                        number of simultaneous shifts returned\n*                        by IPARMQ(ISPEC=15).  The default for\n*                        (IHI-ILO+1).LE.500 is NS.  The default\n*                        for (IHI-ILO+1).GT.500 is 3*NS/2.\n*\n*       IPARMQ(ISPEC=14) Nibble crossover point.  Default: 14.\n*\n*       IPARMQ(ISPEC=15) Number of simultaneous shifts, NS.\n*                        a multi-shift QR iteration.\n*\n*                        If IHI-ILO+1 is ...\n*\n*                        greater than      ...but less    ... the\n*                        or equal to ...      than        default is\n*\n*                                0               30       NS =   2+\n*                               30               60       NS =   4+\n*                               60              150       NS =  10\n*                              150              590       NS =  **\n*                              590             3000       NS =  64\n*                             3000             6000       NS = 128\n*                             6000             infinity   NS = 256\n*\n*                    (+)  By default matrices of this order are\n*                         passed to the implicit double shift routine\n*                         xLAHQR.  See IPARMQ(ISPEC=12) above.   These\n*                         values of NS are used only in case of a rare\n*                         xLAHQR failure.\n*\n*                    (**) The asterisks (**) indicate an ad-hoc\n*                         function increasing from 10 to 64.\n*\n*       IPARMQ(ISPEC=16) Select structured matrix multiply.\n*                        (See ISPEC=16 above for details.)\n*                        Default: 3.\n*\n*     ================================================================\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  __out__ = NumRu::Lapack.iparmq( ispec, name, opts, n, ilo, ihi, lwork, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 7)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 7)", argc);
  rblapack_ispec = argv[0];
  rblapack_name = argv[1];
  rblapack_opts = argv[2];
  rblapack_n = argv[3];
  rblapack_ilo = argv[4];
  rblapack_ihi = argv[5];
  rblapack_lwork = argv[6];
  if (rb_options != Qnil) {
  }

  name = StringValueCStr(rblapack_name)[0];
  ilo = NUM2INT(rblapack_ilo);
  opts = StringValueCStr(rblapack_opts)[0];
  n = NUM2INT(rblapack_n);
  lwork = NUM2INT(rblapack_lwork);
  ispec = NUM2INT(rblapack_ispec);
  ihi = NUM2INT(rblapack_ihi);

  __out__ = iparmq_(&ispec, &name, &opts, &n, &ilo, &ihi, &lwork);

  rblapack___out__ = INT2NUM(__out__);
  return rblapack___out__;
}

void
init_lapack_iparmq(VALUE mLapack){
  rb_define_module_function(mLapack, "iparmq", rblapack_iparmq, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}

#include "rb_lapack.h"

extern VOID dlahqr_(logical *wantt, logical *wantz, integer *n, integer *ilo, integer *ihi, doublereal *h, integer *ldh, doublereal *wr, doublereal *wi, integer *iloz, integer *ihiz, doublereal *z, integer *ldz, integer *info);

static VALUE sHelp, sUsage;

static VALUE
rblapack_dlahqr(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_wantt;
  logical wantt; 
  VALUE rblapack_wantz;
  logical wantz; 
  VALUE rblapack_ilo;
  integer ilo; 
  VALUE rblapack_ihi;
  integer ihi; 
  VALUE rblapack_h;
  doublereal *h; 
  VALUE rblapack_iloz;
  integer iloz; 
  VALUE rblapack_ihiz;
  integer ihiz; 
  VALUE rblapack_z;
  doublereal *z; 
  VALUE rblapack_ldz;
  integer ldz; 
  VALUE rblapack_wr;
  doublereal *wr; 
  VALUE rblapack_wi;
  doublereal *wi; 
  VALUE rblapack_info;
  integer info; 
  VALUE rblapack_h_out__;
  doublereal *h_out__;
  VALUE rblapack_z_out__;
  doublereal *z_out__;

  integer ldh;
  integer n;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  wr, wi, info, h, z = NumRu::Lapack.dlahqr( wantt, wantz, ilo, ihi, h, iloz, ihiz, z, ldz, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE DLAHQR( WANTT, WANTZ, N, ILO, IHI, H, LDH, WR, WI, ILOZ, IHIZ, Z, LDZ, INFO )\n\n*     Purpose\n*     =======\n*\n*     DLAHQR is an auxiliary routine called by DHSEQR to update the\n*     eigenvalues and Schur decomposition already computed by DHSEQR, by\n*     dealing with the Hessenberg submatrix in rows and columns ILO to\n*     IHI.\n*\n\n*     Arguments\n*     =========\n*\n*     WANTT   (input) LOGICAL\n*          = .TRUE. : the full Schur form T is required;\n*          = .FALSE.: only eigenvalues are required.\n*\n*     WANTZ   (input) LOGICAL\n*          = .TRUE. : the matrix of Schur vectors Z is required;\n*          = .FALSE.: Schur vectors are not required.\n*\n*     N       (input) INTEGER\n*          The order of the matrix H.  N >= 0.\n*\n*     ILO     (input) INTEGER\n*     IHI     (input) INTEGER\n*          It is assumed that H is already upper quasi-triangular in\n*          rows and columns IHI+1:N, and that H(ILO,ILO-1) = 0 (unless\n*          ILO = 1). DLAHQR works primarily with the Hessenberg\n*          submatrix in rows and columns ILO to IHI, but applies\n*          transformations to all of H if WANTT is .TRUE..\n*          1 <= ILO <= max(1,IHI); IHI <= N.\n*\n*     H       (input/output) DOUBLE PRECISION array, dimension (LDH,N)\n*          On entry, the upper Hessenberg matrix H.\n*          On exit, if INFO is zero and if WANTT is .TRUE., H is upper\n*          quasi-triangular in rows and columns ILO:IHI, with any\n*          2-by-2 diagonal blocks in standard form. If INFO is zero\n*          and WANTT is .FALSE., the contents of H are unspecified on\n*          exit.  The output state of H if INFO is nonzero is given\n*          below under the description of INFO.\n*\n*     LDH     (input) INTEGER\n*          The leading dimension of the array H. LDH >= max(1,N).\n*\n*     WR      (output) DOUBLE PRECISION array, dimension (N)\n*     WI      (output) DOUBLE PRECISION array, dimension (N)\n*          The real and imaginary parts, respectively, of the computed\n*          eigenvalues ILO to IHI are stored in the corresponding\n*          elements of WR and WI. If two eigenvalues are computed as a\n*          complex conjugate pair, they are stored in consecutive\n*          elements of WR and WI, say the i-th and (i+1)th, with\n*          WI(i) > 0 and WI(i+1) < 0. If WANTT is .TRUE., the\n*          eigenvalues are stored in the same order as on the diagonal\n*          of the Schur form returned in H, with WR(i) = H(i,i), and, if\n*          H(i:i+1,i:i+1) is a 2-by-2 diagonal block,\n*          WI(i) = sqrt(H(i+1,i)*H(i,i+1)) and WI(i+1) = -WI(i).\n*\n*     ILOZ    (input) INTEGER\n*     IHIZ    (input) INTEGER\n*          Specify the rows of Z to which transformations must be\n*          applied if WANTZ is .TRUE..\n*          1 <= ILOZ <= ILO; IHI <= IHIZ <= N.\n*\n*     Z       (input/output) DOUBLE PRECISION array, dimension (LDZ,N)\n*          If WANTZ is .TRUE., on entry Z must contain the current\n*          matrix Z of transformations accumulated by DHSEQR, and on\n*          exit Z has been updated; transformations are applied only to\n*          the submatrix Z(ILOZ:IHIZ,ILO:IHI).\n*          If WANTZ is .FALSE., Z is not referenced.\n*\n*     LDZ     (input) INTEGER\n*          The leading dimension of the array Z. LDZ >= max(1,N).\n*\n*     INFO    (output) INTEGER\n*           =   0: successful exit\n*          .GT. 0: If INFO = i, DLAHQR failed to compute all the\n*                  eigenvalues ILO to IHI in a total of 30 iterations\n*                  per eigenvalue; elements i+1:ihi of WR and WI\n*                  contain those eigenvalues which have been\n*                  successfully computed.\n*\n*                  If INFO .GT. 0 and WANTT is .FALSE., then on exit,\n*                  the remaining unconverged eigenvalues are the\n*                  eigenvalues of the upper Hessenberg matrix rows\n*                  and columns ILO thorugh INFO of the final, output\n*                  value of H.\n*\n*                  If INFO .GT. 0 and WANTT is .TRUE., then on exit\n*          (*)       (initial value of H)*U  = U*(final value of H)\n*                  where U is an orthognal matrix.    The final\n*                  value of H is upper Hessenberg and triangular in\n*                  rows and columns INFO+1 through IHI.\n*\n*                  If INFO .GT. 0 and WANTZ is .TRUE., then on exit\n*                      (final value of Z)  = (initial value of Z)*U\n*                  where U is the orthogonal matrix in (*)\n*                  (regardless of the value of WANTT.)\n*\n\n*     Further Details\n*     ===============\n*\n*     02-96 Based on modifications by\n*     David Day, Sandia National Laboratory, USA\n*\n*     12-04 Further modifications by\n*     Ralph Byers, University of Kansas, USA\n*     This is a modified version of DLAHQR from LAPACK version 3.0.\n*     It is (1) more robust against overflow and underflow and\n*     (2) adopts the more conservative Ahues & Tisseur stopping\n*     criterion (LAWN 122, 1997).\n*\n*     =========================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  wr, wi, info, h, z = NumRu::Lapack.dlahqr( wantt, wantz, ilo, ihi, h, iloz, ihiz, z, ldz, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 9)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 9)", argc);
  rblapack_wantt = argv[0];
  rblapack_wantz = argv[1];
  rblapack_ilo = argv[2];
  rblapack_ihi = argv[3];
  rblapack_h = argv[4];
  rblapack_iloz = argv[5];
  rblapack_ihiz = argv[6];
  rblapack_z = argv[7];
  rblapack_ldz = argv[8];
  if (rb_options != Qnil) {
  }

  ilo = NUM2INT(rblapack_ilo);
  wantz = (rblapack_wantz == Qtrue);
  ldz = NUM2INT(rblapack_ldz);
  ihi = NUM2INT(rblapack_ihi);
  if (!NA_IsNArray(rblapack_h))
    rb_raise(rb_eArgError, "h (5th argument) must be NArray");
  if (NA_RANK(rblapack_h) != 2)
    rb_raise(rb_eArgError, "rank of h (5th argument) must be %d", 2);
  n = NA_SHAPE1(rblapack_h);
  ldh = NA_SHAPE0(rblapack_h);
  if (NA_TYPE(rblapack_h) != NA_DFLOAT)
    rblapack_h = na_change_type(rblapack_h, NA_DFLOAT);
  h = NA_PTR_TYPE(rblapack_h, doublereal*);
  wantt = (rblapack_wantt == Qtrue);
  ihiz = NUM2INT(rblapack_ihiz);
  iloz = NUM2INT(rblapack_iloz);
  if (!NA_IsNArray(rblapack_z))
    rb_raise(rb_eArgError, "z (8th argument) must be NArray");
  if (NA_RANK(rblapack_z) != 2)
    rb_raise(rb_eArgError, "rank of z (8th argument) must be %d", 2);
  if (NA_SHAPE1(rblapack_z) != (wantz ? n : 0))
    rb_raise(rb_eRuntimeError, "shape 1 of z must be %d", wantz ? n : 0);
  if (NA_SHAPE0(rblapack_z) != (wantz ? ldz : 0))
    rb_raise(rb_eRuntimeError, "shape 0 of z must be %d", wantz ? ldz : 0);
  if (NA_TYPE(rblapack_z) != NA_DFLOAT)
    rblapack_z = na_change_type(rblapack_z, NA_DFLOAT);
  z = NA_PTR_TYPE(rblapack_z, doublereal*);
  {
    int shape[1];
    shape[0] = n;
    rblapack_wr = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  wr = NA_PTR_TYPE(rblapack_wr, doublereal*);
  {
    int shape[1];
    shape[0] = n;
    rblapack_wi = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  wi = NA_PTR_TYPE(rblapack_wi, doublereal*);
  {
    int shape[2];
    shape[0] = ldh;
    shape[1] = n;
    rblapack_h_out__ = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  h_out__ = NA_PTR_TYPE(rblapack_h_out__, doublereal*);
  MEMCPY(h_out__, h, doublereal, NA_TOTAL(rblapack_h));
  rblapack_h = rblapack_h_out__;
  h = h_out__;
  {
    int shape[2];
    shape[0] = wantz ? ldz : 0;
    shape[1] = wantz ? n : 0;
    rblapack_z_out__ = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  z_out__ = NA_PTR_TYPE(rblapack_z_out__, doublereal*);
  MEMCPY(z_out__, z, doublereal, NA_TOTAL(rblapack_z));
  rblapack_z = rblapack_z_out__;
  z = z_out__;

  dlahqr_(&wantt, &wantz, &n, &ilo, &ihi, h, &ldh, wr, wi, &iloz, &ihiz, z, &ldz, &info);

  rblapack_info = INT2NUM(info);
  return rb_ary_new3(5, rblapack_wr, rblapack_wi, rblapack_info, rblapack_h, rblapack_z);
}

void
init_lapack_dlahqr(VALUE mLapack){
  rb_define_module_function(mLapack, "dlahqr", rblapack_dlahqr, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}

#include "rb_lapack.h"

extern VOID clahqr_(logical *wantt, logical *wantz, integer *n, integer *ilo, integer *ihi, complex *h, integer *ldh, complex *w, integer *iloz, integer *ihiz, complex *z, integer *ldz, integer *info);

static VALUE sHelp, sUsage;

static VALUE
rblapack_clahqr(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_wantt;
  logical wantt; 
  VALUE rblapack_wantz;
  logical wantz; 
  VALUE rblapack_ilo;
  integer ilo; 
  VALUE rblapack_ihi;
  integer ihi; 
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
  VALUE rblapack_w;
  complex *w; 
  VALUE rblapack_info;
  integer info; 
  VALUE rblapack_h_out__;
  complex *h_out__;
  VALUE rblapack_z_out__;
  complex *z_out__;

  integer ldh;
  integer n;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  w, info, h, z = NumRu::Lapack.clahqr( wantt, wantz, ilo, ihi, h, iloz, ihiz, z, ldz, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE CLAHQR( WANTT, WANTZ, N, ILO, IHI, H, LDH, W, ILOZ, IHIZ, Z, LDZ, INFO )\n\n*     Purpose\n*     =======\n*\n*     CLAHQR is an auxiliary routine called by CHSEQR to update the\n*     eigenvalues and Schur decomposition already computed by CHSEQR, by\n*     dealing with the Hessenberg submatrix in rows and columns ILO to\n*     IHI.\n*\n\n*     Arguments\n*     =========\n*\n*     WANTT   (input) LOGICAL\n*          = .TRUE. : the full Schur form T is required;\n*          = .FALSE.: only eigenvalues are required.\n*\n*     WANTZ   (input) LOGICAL\n*          = .TRUE. : the matrix of Schur vectors Z is required;\n*          = .FALSE.: Schur vectors are not required.\n*\n*     N       (input) INTEGER\n*          The order of the matrix H.  N >= 0.\n*\n*     ILO     (input) INTEGER\n*     IHI     (input) INTEGER\n*          It is assumed that H is already upper triangular in rows and\n*          columns IHI+1:N, and that H(ILO,ILO-1) = 0 (unless ILO = 1).\n*          CLAHQR works primarily with the Hessenberg submatrix in rows\n*          and columns ILO to IHI, but applies transformations to all of\n*          H if WANTT is .TRUE..\n*          1 <= ILO <= max(1,IHI); IHI <= N.\n*\n*     H       (input/output) COMPLEX array, dimension (LDH,N)\n*          On entry, the upper Hessenberg matrix H.\n*          On exit, if INFO is zero and if WANTT is .TRUE., then H\n*          is upper triangular in rows and columns ILO:IHI.  If INFO\n*          is zero and if WANTT is .FALSE., then the contents of H\n*          are unspecified on exit.  The output state of H in case\n*          INF is positive is below under the description of INFO.\n*\n*     LDH     (input) INTEGER\n*          The leading dimension of the array H. LDH >= max(1,N).\n*\n*     W       (output) COMPLEX array, dimension (N)\n*          The computed eigenvalues ILO to IHI are stored in the\n*          corresponding elements of W. If WANTT is .TRUE., the\n*          eigenvalues are stored in the same order as on the diagonal\n*          of the Schur form returned in H, with W(i) = H(i,i).\n*\n*     ILOZ    (input) INTEGER\n*     IHIZ    (input) INTEGER\n*          Specify the rows of Z to which transformations must be\n*          applied if WANTZ is .TRUE..\n*          1 <= ILOZ <= ILO; IHI <= IHIZ <= N.\n*\n*     Z       (input/output) COMPLEX array, dimension (LDZ,N)\n*          If WANTZ is .TRUE., on entry Z must contain the current\n*          matrix Z of transformations accumulated by CHSEQR, and on\n*          exit Z has been updated; transformations are applied only to\n*          the submatrix Z(ILOZ:IHIZ,ILO:IHI).\n*          If WANTZ is .FALSE., Z is not referenced.\n*\n*     LDZ     (input) INTEGER\n*          The leading dimension of the array Z. LDZ >= max(1,N).\n*\n*     INFO    (output) INTEGER\n*           =   0: successful exit\n*          .GT. 0: if INFO = i, CLAHQR failed to compute all the\n*                  eigenvalues ILO to IHI in a total of 30 iterations\n*                  per eigenvalue; elements i+1:ihi of W contain\n*                  those eigenvalues which have been successfully\n*                  computed.\n*\n*                  If INFO .GT. 0 and WANTT is .FALSE., then on exit,\n*                  the remaining unconverged eigenvalues are the\n*                  eigenvalues of the upper Hessenberg matrix\n*                  rows and columns ILO thorugh INFO of the final,\n*                  output value of H.\n*\n*                  If INFO .GT. 0 and WANTT is .TRUE., then on exit\n*          (*)       (initial value of H)*U  = U*(final value of H)\n*                  where U is an orthognal matrix.    The final\n*                  value of H is upper Hessenberg and triangular in\n*                  rows and columns INFO+1 through IHI.\n*\n*                  If INFO .GT. 0 and WANTZ is .TRUE., then on exit\n*                      (final value of Z)  = (initial value of Z)*U\n*                  where U is the orthogonal matrix in (*)\n*                  (regardless of the value of WANTT.)\n*\n\n*     Further Details\n*     ===============\n*\n*     02-96 Based on modifications by\n*     David Day, Sandia National Laboratory, USA\n*\n*     12-04 Further modifications by\n*     Ralph Byers, University of Kansas, USA\n*     This is a modified version of CLAHQR from LAPACK version 3.0.\n*     It is (1) more robust against overflow and underflow and\n*     (2) adopts the more conservative Ahues & Tisseur stopping\n*     criterion (LAWN 122, 1997).\n*\n*     =========================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  w, info, h, z = NumRu::Lapack.clahqr( wantt, wantz, ilo, ihi, h, iloz, ihiz, z, ldz, [:usage => usage, :help => help])\n");
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
  if (NA_TYPE(rblapack_h) != NA_SCOMPLEX)
    rblapack_h = na_change_type(rblapack_h, NA_SCOMPLEX);
  h = NA_PTR_TYPE(rblapack_h, complex*);
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
  if (NA_TYPE(rblapack_z) != NA_SCOMPLEX)
    rblapack_z = na_change_type(rblapack_z, NA_SCOMPLEX);
  z = NA_PTR_TYPE(rblapack_z, complex*);
  {
    int shape[1];
    shape[0] = n;
    rblapack_w = na_make_object(NA_SCOMPLEX, 1, shape, cNArray);
  }
  w = NA_PTR_TYPE(rblapack_w, complex*);
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
    shape[1] = wantz ? n : 0;
    rblapack_z_out__ = na_make_object(NA_SCOMPLEX, 2, shape, cNArray);
  }
  z_out__ = NA_PTR_TYPE(rblapack_z_out__, complex*);
  MEMCPY(z_out__, z, complex, NA_TOTAL(rblapack_z));
  rblapack_z = rblapack_z_out__;
  z = z_out__;

  clahqr_(&wantt, &wantz, &n, &ilo, &ihi, h, &ldh, w, &iloz, &ihiz, z, &ldz, &info);

  rblapack_info = INT2NUM(info);
  return rb_ary_new3(4, rblapack_w, rblapack_info, rblapack_h, rblapack_z);
}

void
init_lapack_clahqr(VALUE mLapack){
  rb_define_module_function(mLapack, "clahqr", rblapack_clahqr, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}

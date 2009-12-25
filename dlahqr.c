#include "rb_lapack.h"

static VALUE
rb_dlahqr(int argc, VALUE *argv, VALUE self){
  VALUE rb_wantt;
  logical wantt; 
  VALUE rb_wantz;
  logical wantz; 
  VALUE rb_ilo;
  integer ilo; 
  VALUE rb_ihi;
  integer ihi; 
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
  VALUE rb_wr;
  doublereal *wr; 
  VALUE rb_wi;
  doublereal *wi; 
  VALUE rb_info;
  integer info; 
  VALUE rb_h_out__;
  doublereal *h_out__;
  VALUE rb_z_out__;
  doublereal *z_out__;

  integer ldh;
  integer n;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  wr, wi, info, h, z = NumRu::Lapack.dlahqr( wantt, wantz, ilo, ihi, h, iloz, ihiz, z, ldz)\n    or\n  NumRu::Lapack.dlahqr  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE DLAHQR( WANTT, WANTZ, N, ILO, IHI, H, LDH, WR, WI, ILOZ, IHIZ, Z, LDZ, INFO )\n\n*     Purpose\n*     =======\n*\n*     DLAHQR is an auxiliary routine called by DHSEQR to update the\n*     eigenvalues and Schur decomposition already computed by DHSEQR, by\n*     dealing with the Hessenberg submatrix in rows and columns ILO to\n*     IHI.\n*\n\n*     Arguments\n*     =========\n*\n*     WANTT   (input) LOGICAL\n*          = .TRUE. : the full Schur form T is required;\n*          = .FALSE.: only eigenvalues are required.\n*\n*     WANTZ   (input) LOGICAL\n*          = .TRUE. : the matrix of Schur vectors Z is required;\n*          = .FALSE.: Schur vectors are not required.\n*\n*     N       (input) INTEGER\n*          The order of the matrix H.  N >= 0.\n*\n*     ILO     (input) INTEGER\n*     IHI     (input) INTEGER\n*          It is assumed that H is already upper quasi-triangular in\n*          rows and columns IHI+1:N, and that H(ILO,ILO-1) = 0 (unless\n*          ILO = 1). DLAHQR works primarily with the Hessenberg\n*          submatrix in rows and columns ILO to IHI, but applies\n*          transformations to all of H if WANTT is .TRUE..\n*          1 <= ILO <= max(1,IHI); IHI <= N.\n*\n*     H       (input/output) DOUBLE PRECISION array, dimension (LDH,N)\n*          On entry, the upper Hessenberg matrix H.\n*          On exit, if INFO is zero and if WANTT is .TRUE., H is upper\n*          quasi-triangular in rows and columns ILO:IHI, with any\n*          2-by-2 diagonal blocks in standard form. If INFO is zero\n*          and WANTT is .FALSE., the contents of H are unspecified on\n*          exit.  The output state of H if INFO is nonzero is given\n*          below under the description of INFO.\n*\n*     LDH     (input) INTEGER\n*          The leading dimension of the array H. LDH >= max(1,N).\n*\n*     WR      (output) DOUBLE PRECISION array, dimension (N)\n*     WI      (output) DOUBLE PRECISION array, dimension (N)\n*          The real and imaginary parts, respectively, of the computed\n*          eigenvalues ILO to IHI are stored in the corresponding\n*          elements of WR and WI. If two eigenvalues are computed as a\n*          complex conjugate pair, they are stored in consecutive\n*          elements of WR and WI, say the i-th and (i+1)th, with\n*          WI(i) > 0 and WI(i+1) < 0. If WANTT is .TRUE., the\n*          eigenvalues are stored in the same order as on the diagonal\n*          of the Schur form returned in H, with WR(i) = H(i,i), and, if\n*          H(i:i+1,i:i+1) is a 2-by-2 diagonal block,\n*          WI(i) = sqrt(H(i+1,i)*H(i,i+1)) and WI(i+1) = -WI(i).\n*\n*     ILOZ    (input) INTEGER\n*     IHIZ    (input) INTEGER\n*          Specify the rows of Z to which transformations must be\n*          applied if WANTZ is .TRUE..\n*          1 <= ILOZ <= ILO; IHI <= IHIZ <= N.\n*\n*     Z       (input/output) DOUBLE PRECISION array, dimension (LDZ,N)\n*          If WANTZ is .TRUE., on entry Z must contain the current\n*          matrix Z of transformations accumulated by DHSEQR, and on\n*          exit Z has been updated; transformations are applied only to\n*          the submatrix Z(ILOZ:IHIZ,ILO:IHI).\n*          If WANTZ is .FALSE., Z is not referenced.\n*\n*     LDZ     (input) INTEGER\n*          The leading dimension of the array Z. LDZ >= max(1,N).\n*\n*     INFO    (output) INTEGER\n*           =   0: successful exit\n*          .GT. 0: If INFO = i, DLAHQR failed to compute all the\n*                  eigenvalues ILO to IHI in a total of 30 iterations\n*                  per eigenvalue; elements i+1:ihi of WR and WI\n*                  contain those eigenvalues which have been\n*                  successfully computed.\n*\n*                  If INFO .GT. 0 and WANTT is .FALSE., then on exit,\n*                  the remaining unconverged eigenvalues are the\n*                  eigenvalues of the upper Hessenberg matrix rows\n*                  and columns ILO thorugh INFO of the final, output\n*                  value of H.\n*\n*                  If INFO .GT. 0 and WANTT is .TRUE., then on exit\n*          (*)       (initial value of H)*U  = U*(final value of H)\n*                  where U is an orthognal matrix.    The final\n*                  value of H is upper Hessenberg and triangular in\n*                  rows and columns INFO+1 through IHI.\n*\n*                  If INFO .GT. 0 and WANTZ is .TRUE., then on exit\n*                      (final value of Z)  = (initial value of Z)*U\n*                  where U is the orthogonal matrix in (*)\n*                  (regardless of the value of WANTT.)\n*\n\n*     Further Details\n*     ===============\n*\n*     02-96 Based on modifications by\n*     David Day, Sandia National Laboratory, USA\n*\n*     12-04 Further modifications by\n*     Ralph Byers, University of Kansas, USA\n*     This is a modified version of DLAHQR from LAPACK version 3.0.\n*     It is (1) more robust against overflow and underflow and\n*     (2) adopts the more conservative Ahues & Tisseur stopping\n*     criterion (LAWN 122, 1997).\n*\n*     =========================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 9)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 9)", argc);
  rb_wantt = argv[0];
  rb_wantz = argv[1];
  rb_ilo = argv[2];
  rb_ihi = argv[3];
  rb_h = argv[4];
  rb_iloz = argv[5];
  rb_ihiz = argv[6];
  rb_z = argv[7];
  rb_ldz = argv[8];

  wantt = (rb_wantt == Qtrue);
  wantz = (rb_wantz == Qtrue);
  ilo = NUM2INT(rb_ilo);
  ihi = NUM2INT(rb_ihi);
  iloz = NUM2INT(rb_iloz);
  ihiz = NUM2INT(rb_ihiz);
  ldz = NUM2INT(rb_ldz);
  if (!NA_IsNArray(rb_h))
    rb_raise(rb_eArgError, "h (5th argument) must be NArray");
  if (NA_RANK(rb_h) != 2)
    rb_raise(rb_eArgError, "rank of h (5th argument) must be %d", 2);
  ldh = NA_SHAPE0(rb_h);
  n = NA_SHAPE1(rb_h);
  if (NA_TYPE(rb_h) != NA_DFLOAT)
    rb_h = na_change_type(rb_h, NA_DFLOAT);
  h = NA_PTR_TYPE(rb_h, doublereal*);
  if (!NA_IsNArray(rb_z))
    rb_raise(rb_eArgError, "z (8th argument) must be NArray");
  if (NA_RANK(rb_z) != 2)
    rb_raise(rb_eArgError, "rank of z (8th argument) must be %d", 2);
  if (NA_SHAPE0(rb_z) != (wantz ? ldz : 0))
    rb_raise(rb_eRuntimeError, "shape 0 of z must be %d", wantz ? ldz : 0);
  if (NA_SHAPE1(rb_z) != (wantz ? n : 0))
    rb_raise(rb_eRuntimeError, "shape 1 of z must be %d", wantz ? n : 0);
  if (NA_TYPE(rb_z) != NA_DFLOAT)
    rb_z = na_change_type(rb_z, NA_DFLOAT);
  z = NA_PTR_TYPE(rb_z, doublereal*);
  {
    int shape[1];
    shape[0] = n;
    rb_wr = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  wr = NA_PTR_TYPE(rb_wr, doublereal*);
  {
    int shape[1];
    shape[0] = n;
    rb_wi = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  wi = NA_PTR_TYPE(rb_wi, doublereal*);
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
    shape[1] = wantz ? n : 0;
    rb_z_out__ = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  z_out__ = NA_PTR_TYPE(rb_z_out__, doublereal*);
  MEMCPY(z_out__, z, doublereal, NA_TOTAL(rb_z));
  rb_z = rb_z_out__;
  z = z_out__;

  dlahqr_(&wantt, &wantz, &n, &ilo, &ihi, h, &ldh, wr, wi, &iloz, &ihiz, z, &ldz, &info);

  rb_info = INT2NUM(info);
  return rb_ary_new3(5, rb_wr, rb_wi, rb_info, rb_h, rb_z);
}

void
init_lapack_dlahqr(VALUE mLapack){
  rb_define_module_function(mLapack, "dlahqr", rb_dlahqr, -1);
}

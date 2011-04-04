#include "rb_lapack.h"

extern VOID dtrevc_(char *side, char *howmny, logical *select, integer *n, doublereal *t, integer *ldt, doublereal *vl, integer *ldvl, doublereal *vr, integer *ldvr, integer *mm, integer *m, doublereal *work, integer *info);

static VALUE
rb_dtrevc(int argc, VALUE *argv, VALUE self){
  VALUE rb_side;
  char side; 
  VALUE rb_howmny;
  char howmny; 
  VALUE rb_select;
  logical *select; 
  VALUE rb_t;
  doublereal *t; 
  VALUE rb_vl;
  doublereal *vl; 
  VALUE rb_vr;
  doublereal *vr; 
  VALUE rb_m;
  integer m; 
  VALUE rb_info;
  integer info; 
  VALUE rb_select_out__;
  logical *select_out__;
  VALUE rb_vl_out__;
  doublereal *vl_out__;
  VALUE rb_vr_out__;
  doublereal *vr_out__;
  doublereal *work;

  integer n;
  integer ldt;
  integer ldvl;
  integer mm;
  integer ldvr;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  m, info, select, vl, vr = NumRu::Lapack.dtrevc( side, howmny, select, t, vl, vr)\n    or\n  NumRu::Lapack.dtrevc  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE DTREVC( SIDE, HOWMNY, SELECT, N, T, LDT, VL, LDVL, VR, LDVR, MM, M, WORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  DTREVC computes some or all of the right and/or left eigenvectors of\n*  a real upper quasi-triangular matrix T.\n*  Matrices of this type are produced by the Schur factorization of\n*  a real general matrix:  A = Q*T*Q**T, as computed by DHSEQR.\n*  \n*  The right eigenvector x and the left eigenvector y of T corresponding\n*  to an eigenvalue w are defined by:\n*  \n*     T*x = w*x,     (y**H)*T = w*(y**H)\n*  \n*  where y**H denotes the conjugate transpose of y.\n*  The eigenvalues are not input to this routine, but are read directly\n*  from the diagonal blocks of T.\n*  \n*  This routine returns the matrices X and/or Y of right and left\n*  eigenvectors of T, or the products Q*X and/or Q*Y, where Q is an\n*  input matrix.  If Q is the orthogonal factor that reduces a matrix\n*  A to Schur form T, then Q*X and Q*Y are the matrices of right and\n*  left eigenvectors of A.\n*\n\n*  Arguments\n*  =========\n*\n*  SIDE    (input) CHARACTER*1\n*          = 'R':  compute right eigenvectors only;\n*          = 'L':  compute left eigenvectors only;\n*          = 'B':  compute both right and left eigenvectors.\n*\n*  HOWMNY  (input) CHARACTER*1\n*          = 'A':  compute all right and/or left eigenvectors;\n*          = 'B':  compute all right and/or left eigenvectors,\n*                  backtransformed by the matrices in VR and/or VL;\n*          = 'S':  compute selected right and/or left eigenvectors,\n*                  as indicated by the logical array SELECT.\n*\n*  SELECT  (input/output) LOGICAL array, dimension (N)\n*          If HOWMNY = 'S', SELECT specifies the eigenvectors to be\n*          computed.\n*          If w(j) is a real eigenvalue, the corresponding real\n*          eigenvector is computed if SELECT(j) is .TRUE..\n*          If w(j) and w(j+1) are the real and imaginary parts of a\n*          complex eigenvalue, the corresponding complex eigenvector is\n*          computed if either SELECT(j) or SELECT(j+1) is .TRUE., and\n*          on exit SELECT(j) is set to .TRUE. and SELECT(j+1) is set to\n*          .FALSE..\n*          Not referenced if HOWMNY = 'A' or 'B'.\n*\n*  N       (input) INTEGER\n*          The order of the matrix T. N >= 0.\n*\n*  T       (input) DOUBLE PRECISION array, dimension (LDT,N)\n*          The upper quasi-triangular matrix T in Schur canonical form.\n*\n*  LDT     (input) INTEGER\n*          The leading dimension of the array T. LDT >= max(1,N).\n*\n*  VL      (input/output) DOUBLE PRECISION array, dimension (LDVL,MM)\n*          On entry, if SIDE = 'L' or 'B' and HOWMNY = 'B', VL must\n*          contain an N-by-N matrix Q (usually the orthogonal matrix Q\n*          of Schur vectors returned by DHSEQR).\n*          On exit, if SIDE = 'L' or 'B', VL contains:\n*          if HOWMNY = 'A', the matrix Y of left eigenvectors of T;\n*          if HOWMNY = 'B', the matrix Q*Y;\n*          if HOWMNY = 'S', the left eigenvectors of T specified by\n*                           SELECT, stored consecutively in the columns\n*                           of VL, in the same order as their\n*                           eigenvalues.\n*          A complex eigenvector corresponding to a complex eigenvalue\n*          is stored in two consecutive columns, the first holding the\n*          real part, and the second the imaginary part.\n*          Not referenced if SIDE = 'R'.\n*\n*  LDVL    (input) INTEGER\n*          The leading dimension of the array VL.  LDVL >= 1, and if\n*          SIDE = 'L' or 'B', LDVL >= N.\n*\n*  VR      (input/output) DOUBLE PRECISION array, dimension (LDVR,MM)\n*          On entry, if SIDE = 'R' or 'B' and HOWMNY = 'B', VR must\n*          contain an N-by-N matrix Q (usually the orthogonal matrix Q\n*          of Schur vectors returned by DHSEQR).\n*          On exit, if SIDE = 'R' or 'B', VR contains:\n*          if HOWMNY = 'A', the matrix X of right eigenvectors of T;\n*          if HOWMNY = 'B', the matrix Q*X;\n*          if HOWMNY = 'S', the right eigenvectors of T specified by\n*                           SELECT, stored consecutively in the columns\n*                           of VR, in the same order as their\n*                           eigenvalues.\n*          A complex eigenvector corresponding to a complex eigenvalue\n*          is stored in two consecutive columns, the first holding the\n*          real part and the second the imaginary part.\n*          Not referenced if SIDE = 'L'.\n*\n*  LDVR    (input) INTEGER\n*          The leading dimension of the array VR.  LDVR >= 1, and if\n*          SIDE = 'R' or 'B', LDVR >= N.\n*\n*  MM      (input) INTEGER\n*          The number of columns in the arrays VL and/or VR. MM >= M.\n*\n*  M       (output) INTEGER\n*          The number of columns in the arrays VL and/or VR actually\n*          used to store the eigenvectors.\n*          If HOWMNY = 'A' or 'B', M is set to N.\n*          Each selected real eigenvector occupies one column and each\n*          selected complex eigenvector occupies two columns.\n*\n*  WORK    (workspace) DOUBLE PRECISION array, dimension (3*N)\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit\n*          < 0:  if INFO = -i, the i-th argument had an illegal value\n*\n\n*  Further Details\n*  ===============\n*\n*  The algorithm used in this program is basically backward (forward)\n*  substitution, with scaling to make the the code robust against\n*  possible overflow.\n*\n*  Each eigenvector is normalized so that the element of largest\n*  magnitude has magnitude 1; here the magnitude of a complex number\n*  (x,y) is taken to be |x| + |y|.\n*\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 6)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 6)", argc);
  rb_side = argv[0];
  rb_howmny = argv[1];
  rb_select = argv[2];
  rb_t = argv[3];
  rb_vl = argv[4];
  rb_vr = argv[5];

  if (!NA_IsNArray(rb_vl))
    rb_raise(rb_eArgError, "vl (5th argument) must be NArray");
  if (NA_RANK(rb_vl) != 2)
    rb_raise(rb_eArgError, "rank of vl (5th argument) must be %d", 2);
  mm = NA_SHAPE1(rb_vl);
  ldvl = NA_SHAPE0(rb_vl);
  if (NA_TYPE(rb_vl) != NA_DFLOAT)
    rb_vl = na_change_type(rb_vl, NA_DFLOAT);
  vl = NA_PTR_TYPE(rb_vl, doublereal*);
  side = StringValueCStr(rb_side)[0];
  if (!NA_IsNArray(rb_vr))
    rb_raise(rb_eArgError, "vr (6th argument) must be NArray");
  if (NA_RANK(rb_vr) != 2)
    rb_raise(rb_eArgError, "rank of vr (6th argument) must be %d", 2);
  if (NA_SHAPE1(rb_vr) != mm)
    rb_raise(rb_eRuntimeError, "shape 1 of vr must be the same as shape 1 of vl");
  ldvr = NA_SHAPE0(rb_vr);
  if (NA_TYPE(rb_vr) != NA_DFLOAT)
    rb_vr = na_change_type(rb_vr, NA_DFLOAT);
  vr = NA_PTR_TYPE(rb_vr, doublereal*);
  if (!NA_IsNArray(rb_select))
    rb_raise(rb_eArgError, "select (3th argument) must be NArray");
  if (NA_RANK(rb_select) != 1)
    rb_raise(rb_eArgError, "rank of select (3th argument) must be %d", 1);
  n = NA_SHAPE0(rb_select);
  if (NA_TYPE(rb_select) != NA_LINT)
    rb_select = na_change_type(rb_select, NA_LINT);
  select = NA_PTR_TYPE(rb_select, logical*);
  if (!NA_IsNArray(rb_t))
    rb_raise(rb_eArgError, "t (4th argument) must be NArray");
  if (NA_RANK(rb_t) != 2)
    rb_raise(rb_eArgError, "rank of t (4th argument) must be %d", 2);
  if (NA_SHAPE1(rb_t) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of t must be the same as shape 0 of select");
  ldt = NA_SHAPE0(rb_t);
  if (NA_TYPE(rb_t) != NA_DFLOAT)
    rb_t = na_change_type(rb_t, NA_DFLOAT);
  t = NA_PTR_TYPE(rb_t, doublereal*);
  howmny = StringValueCStr(rb_howmny)[0];
  {
    int shape[1];
    shape[0] = n;
    rb_select_out__ = na_make_object(NA_LINT, 1, shape, cNArray);
  }
  select_out__ = NA_PTR_TYPE(rb_select_out__, logical*);
  MEMCPY(select_out__, select, logical, NA_TOTAL(rb_select));
  rb_select = rb_select_out__;
  select = select_out__;
  {
    int shape[2];
    shape[0] = ldvl;
    shape[1] = mm;
    rb_vl_out__ = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  vl_out__ = NA_PTR_TYPE(rb_vl_out__, doublereal*);
  MEMCPY(vl_out__, vl, doublereal, NA_TOTAL(rb_vl));
  rb_vl = rb_vl_out__;
  vl = vl_out__;
  {
    int shape[2];
    shape[0] = ldvr;
    shape[1] = mm;
    rb_vr_out__ = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  vr_out__ = NA_PTR_TYPE(rb_vr_out__, doublereal*);
  MEMCPY(vr_out__, vr, doublereal, NA_TOTAL(rb_vr));
  rb_vr = rb_vr_out__;
  vr = vr_out__;
  work = ALLOC_N(doublereal, (3*n));

  dtrevc_(&side, &howmny, select, &n, t, &ldt, vl, &ldvl, vr, &ldvr, &mm, &m, work, &info);

  free(work);
  rb_m = INT2NUM(m);
  rb_info = INT2NUM(info);
  return rb_ary_new3(5, rb_m, rb_info, rb_select, rb_vl, rb_vr);
}

void
init_lapack_dtrevc(VALUE mLapack){
  rb_define_module_function(mLapack, "dtrevc", rb_dtrevc, -1);
}

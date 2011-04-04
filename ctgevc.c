#include "rb_lapack.h"

extern VOID ctgevc_(char *side, char *howmny, logical *select, integer *n, complex *s, integer *lds, complex *p, integer *ldp, complex *vl, integer *ldvl, complex *vr, integer *ldvr, integer *mm, integer *m, complex *work, real *rwork, integer *info);

static VALUE
rb_ctgevc(int argc, VALUE *argv, VALUE self){
  VALUE rb_side;
  char side; 
  VALUE rb_howmny;
  char howmny; 
  VALUE rb_select;
  logical *select; 
  VALUE rb_s;
  complex *s; 
  VALUE rb_p;
  complex *p; 
  VALUE rb_vl;
  complex *vl; 
  VALUE rb_vr;
  complex *vr; 
  VALUE rb_m;
  integer m; 
  VALUE rb_info;
  integer info; 
  VALUE rb_vl_out__;
  complex *vl_out__;
  VALUE rb_vr_out__;
  complex *vr_out__;
  complex *work;
  real *rwork;

  integer n;
  integer lds;
  integer ldp;
  integer ldvl;
  integer mm;
  integer ldvr;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  m, info, vl, vr = NumRu::Lapack.ctgevc( side, howmny, select, s, p, vl, vr)\n    or\n  NumRu::Lapack.ctgevc  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE CTGEVC( SIDE, HOWMNY, SELECT, N, S, LDS, P, LDP, VL, LDVL, VR, LDVR, MM, M, WORK, RWORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  CTGEVC computes some or all of the right and/or left eigenvectors of\n*  a pair of complex matrices (S,P), where S and P are upper triangular.\n*  Matrix pairs of this type are produced by the generalized Schur\n*  factorization of a complex matrix pair (A,B):\n*  \n*     A = Q*S*Z**H,  B = Q*P*Z**H\n*  \n*  as computed by CGGHRD + CHGEQZ.\n*  \n*  The right eigenvector x and the left eigenvector y of (S,P)\n*  corresponding to an eigenvalue w are defined by:\n*  \n*     S*x = w*P*x,  (y**H)*S = w*(y**H)*P,\n*  \n*  where y**H denotes the conjugate tranpose of y.\n*  The eigenvalues are not input to this routine, but are computed\n*  directly from the diagonal elements of S and P.\n*  \n*  This routine returns the matrices X and/or Y of right and left\n*  eigenvectors of (S,P), or the products Z*X and/or Q*Y,\n*  where Z and Q are input matrices.\n*  If Q and Z are the unitary factors from the generalized Schur\n*  factorization of a matrix pair (A,B), then Z*X and Q*Y\n*  are the matrices of right and left eigenvectors of (A,B).\n*\n\n*  Arguments\n*  =========\n*\n*  SIDE    (input) CHARACTER*1\n*          = 'R': compute right eigenvectors only;\n*          = 'L': compute left eigenvectors only;\n*          = 'B': compute both right and left eigenvectors.\n*\n*  HOWMNY  (input) CHARACTER*1\n*          = 'A': compute all right and/or left eigenvectors;\n*          = 'B': compute all right and/or left eigenvectors,\n*                 backtransformed by the matrices in VR and/or VL;\n*          = 'S': compute selected right and/or left eigenvectors,\n*                 specified by the logical array SELECT.\n*\n*  SELECT  (input) LOGICAL array, dimension (N)\n*          If HOWMNY='S', SELECT specifies the eigenvectors to be\n*          computed.  The eigenvector corresponding to the j-th\n*          eigenvalue is computed if SELECT(j) = .TRUE..\n*          Not referenced if HOWMNY = 'A' or 'B'.\n*\n*  N       (input) INTEGER\n*          The order of the matrices S and P.  N >= 0.\n*\n*  S       (input) COMPLEX array, dimension (LDS,N)\n*          The upper triangular matrix S from a generalized Schur\n*          factorization, as computed by CHGEQZ.\n*\n*  LDS     (input) INTEGER\n*          The leading dimension of array S.  LDS >= max(1,N).\n*\n*  P       (input) COMPLEX array, dimension (LDP,N)\n*          The upper triangular matrix P from a generalized Schur\n*          factorization, as computed by CHGEQZ.  P must have real\n*          diagonal elements.\n*\n*  LDP     (input) INTEGER\n*          The leading dimension of array P.  LDP >= max(1,N).\n*\n*  VL      (input/output) COMPLEX array, dimension (LDVL,MM)\n*          On entry, if SIDE = 'L' or 'B' and HOWMNY = 'B', VL must\n*          contain an N-by-N matrix Q (usually the unitary matrix Q\n*          of left Schur vectors returned by CHGEQZ).\n*          On exit, if SIDE = 'L' or 'B', VL contains:\n*          if HOWMNY = 'A', the matrix Y of left eigenvectors of (S,P);\n*          if HOWMNY = 'B', the matrix Q*Y;\n*          if HOWMNY = 'S', the left eigenvectors of (S,P) specified by\n*                      SELECT, stored consecutively in the columns of\n*                      VL, in the same order as their eigenvalues.\n*          Not referenced if SIDE = 'R'.\n*\n*  LDVL    (input) INTEGER\n*          The leading dimension of array VL.  LDVL >= 1, and if\n*          SIDE = 'L' or 'l' or 'B' or 'b', LDVL >= N.\n*\n*  VR      (input/output) COMPLEX array, dimension (LDVR,MM)\n*          On entry, if SIDE = 'R' or 'B' and HOWMNY = 'B', VR must\n*          contain an N-by-N matrix Q (usually the unitary matrix Z\n*          of right Schur vectors returned by CHGEQZ).\n*          On exit, if SIDE = 'R' or 'B', VR contains:\n*          if HOWMNY = 'A', the matrix X of right eigenvectors of (S,P);\n*          if HOWMNY = 'B', the matrix Z*X;\n*          if HOWMNY = 'S', the right eigenvectors of (S,P) specified by\n*                      SELECT, stored consecutively in the columns of\n*                      VR, in the same order as their eigenvalues.\n*          Not referenced if SIDE = 'L'.\n*\n*  LDVR    (input) INTEGER\n*          The leading dimension of the array VR.  LDVR >= 1, and if\n*          SIDE = 'R' or 'B', LDVR >= N.\n*\n*  MM      (input) INTEGER\n*          The number of columns in the arrays VL and/or VR. MM >= M.\n*\n*  M       (output) INTEGER\n*          The number of columns in the arrays VL and/or VR actually\n*          used to store the eigenvectors.  If HOWMNY = 'A' or 'B', M\n*          is set to N.  Each selected eigenvector occupies one column.\n*\n*  WORK    (workspace) COMPLEX array, dimension (2*N)\n*\n*  RWORK   (workspace) REAL array, dimension (2*N)\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit.\n*          < 0:  if INFO = -i, the i-th argument had an illegal value.\n*\n\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 7)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 7)", argc);
  rb_side = argv[0];
  rb_howmny = argv[1];
  rb_select = argv[2];
  rb_s = argv[3];
  rb_p = argv[4];
  rb_vl = argv[5];
  rb_vr = argv[6];

  if (!NA_IsNArray(rb_vl))
    rb_raise(rb_eArgError, "vl (6th argument) must be NArray");
  if (NA_RANK(rb_vl) != 2)
    rb_raise(rb_eArgError, "rank of vl (6th argument) must be %d", 2);
  mm = NA_SHAPE1(rb_vl);
  ldvl = NA_SHAPE0(rb_vl);
  if (NA_TYPE(rb_vl) != NA_SCOMPLEX)
    rb_vl = na_change_type(rb_vl, NA_SCOMPLEX);
  vl = NA_PTR_TYPE(rb_vl, complex*);
  side = StringValueCStr(rb_side)[0];
  if (!NA_IsNArray(rb_p))
    rb_raise(rb_eArgError, "p (5th argument) must be NArray");
  if (NA_RANK(rb_p) != 2)
    rb_raise(rb_eArgError, "rank of p (5th argument) must be %d", 2);
  n = NA_SHAPE1(rb_p);
  ldp = NA_SHAPE0(rb_p);
  if (NA_TYPE(rb_p) != NA_SCOMPLEX)
    rb_p = na_change_type(rb_p, NA_SCOMPLEX);
  p = NA_PTR_TYPE(rb_p, complex*);
  if (!NA_IsNArray(rb_vr))
    rb_raise(rb_eArgError, "vr (7th argument) must be NArray");
  if (NA_RANK(rb_vr) != 2)
    rb_raise(rb_eArgError, "rank of vr (7th argument) must be %d", 2);
  if (NA_SHAPE1(rb_vr) != mm)
    rb_raise(rb_eRuntimeError, "shape 1 of vr must be the same as shape 1 of vl");
  ldvr = NA_SHAPE0(rb_vr);
  if (NA_TYPE(rb_vr) != NA_SCOMPLEX)
    rb_vr = na_change_type(rb_vr, NA_SCOMPLEX);
  vr = NA_PTR_TYPE(rb_vr, complex*);
  if (!NA_IsNArray(rb_s))
    rb_raise(rb_eArgError, "s (4th argument) must be NArray");
  if (NA_RANK(rb_s) != 2)
    rb_raise(rb_eArgError, "rank of s (4th argument) must be %d", 2);
  if (NA_SHAPE1(rb_s) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of s must be the same as shape 1 of p");
  lds = NA_SHAPE0(rb_s);
  if (NA_TYPE(rb_s) != NA_SCOMPLEX)
    rb_s = na_change_type(rb_s, NA_SCOMPLEX);
  s = NA_PTR_TYPE(rb_s, complex*);
  if (!NA_IsNArray(rb_select))
    rb_raise(rb_eArgError, "select (3th argument) must be NArray");
  if (NA_RANK(rb_select) != 1)
    rb_raise(rb_eArgError, "rank of select (3th argument) must be %d", 1);
  if (NA_SHAPE0(rb_select) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of select must be the same as shape 1 of p");
  if (NA_TYPE(rb_select) != NA_LINT)
    rb_select = na_change_type(rb_select, NA_LINT);
  select = NA_PTR_TYPE(rb_select, logical*);
  howmny = StringValueCStr(rb_howmny)[0];
  {
    int shape[2];
    shape[0] = ldvl;
    shape[1] = mm;
    rb_vl_out__ = na_make_object(NA_SCOMPLEX, 2, shape, cNArray);
  }
  vl_out__ = NA_PTR_TYPE(rb_vl_out__, complex*);
  MEMCPY(vl_out__, vl, complex, NA_TOTAL(rb_vl));
  rb_vl = rb_vl_out__;
  vl = vl_out__;
  {
    int shape[2];
    shape[0] = ldvr;
    shape[1] = mm;
    rb_vr_out__ = na_make_object(NA_SCOMPLEX, 2, shape, cNArray);
  }
  vr_out__ = NA_PTR_TYPE(rb_vr_out__, complex*);
  MEMCPY(vr_out__, vr, complex, NA_TOTAL(rb_vr));
  rb_vr = rb_vr_out__;
  vr = vr_out__;
  work = ALLOC_N(complex, (2*n));
  rwork = ALLOC_N(real, (2*n));

  ctgevc_(&side, &howmny, select, &n, s, &lds, p, &ldp, vl, &ldvl, vr, &ldvr, &mm, &m, work, rwork, &info);

  free(work);
  free(rwork);
  rb_m = INT2NUM(m);
  rb_info = INT2NUM(info);
  return rb_ary_new3(4, rb_m, rb_info, rb_vl, rb_vr);
}

void
init_lapack_ctgevc(VALUE mLapack){
  rb_define_module_function(mLapack, "ctgevc", rb_ctgevc, -1);
}

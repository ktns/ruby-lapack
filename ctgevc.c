#include "rb_lapack.h"

extern VOID ctgevc_(char *side, char *howmny, logical *select, integer *n, complex *s, integer *lds, complex *p, integer *ldp, complex *vl, integer *ldvl, complex *vr, integer *ldvr, integer *mm, integer *m, complex *work, real *rwork, integer *info);

static VALUE sHelp, sUsage;

static VALUE
rblapack_ctgevc(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_side;
  char side; 
  VALUE rblapack_howmny;
  char howmny; 
  VALUE rblapack_select;
  logical *select; 
  VALUE rblapack_s;
  complex *s; 
  VALUE rblapack_p;
  complex *p; 
  VALUE rblapack_vl;
  complex *vl; 
  VALUE rblapack_vr;
  complex *vr; 
  VALUE rblapack_m;
  integer m; 
  VALUE rblapack_info;
  integer info; 
  VALUE rblapack_vl_out__;
  complex *vl_out__;
  VALUE rblapack_vr_out__;
  complex *vr_out__;
  complex *work;
  real *rwork;

  integer n;
  integer lds;
  integer ldp;
  integer ldvl;
  integer mm;
  integer ldvr;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  m, info, vl, vr = NumRu::Lapack.ctgevc( side, howmny, select, s, p, vl, vr, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE CTGEVC( SIDE, HOWMNY, SELECT, N, S, LDS, P, LDP, VL, LDVL, VR, LDVR, MM, M, WORK, RWORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  CTGEVC computes some or all of the right and/or left eigenvectors of\n*  a pair of complex matrices (S,P), where S and P are upper triangular.\n*  Matrix pairs of this type are produced by the generalized Schur\n*  factorization of a complex matrix pair (A,B):\n*  \n*     A = Q*S*Z**H,  B = Q*P*Z**H\n*  \n*  as computed by CGGHRD + CHGEQZ.\n*  \n*  The right eigenvector x and the left eigenvector y of (S,P)\n*  corresponding to an eigenvalue w are defined by:\n*  \n*     S*x = w*P*x,  (y**H)*S = w*(y**H)*P,\n*  \n*  where y**H denotes the conjugate tranpose of y.\n*  The eigenvalues are not input to this routine, but are computed\n*  directly from the diagonal elements of S and P.\n*  \n*  This routine returns the matrices X and/or Y of right and left\n*  eigenvectors of (S,P), or the products Z*X and/or Q*Y,\n*  where Z and Q are input matrices.\n*  If Q and Z are the unitary factors from the generalized Schur\n*  factorization of a matrix pair (A,B), then Z*X and Q*Y\n*  are the matrices of right and left eigenvectors of (A,B).\n*\n\n*  Arguments\n*  =========\n*\n*  SIDE    (input) CHARACTER*1\n*          = 'R': compute right eigenvectors only;\n*          = 'L': compute left eigenvectors only;\n*          = 'B': compute both right and left eigenvectors.\n*\n*  HOWMNY  (input) CHARACTER*1\n*          = 'A': compute all right and/or left eigenvectors;\n*          = 'B': compute all right and/or left eigenvectors,\n*                 backtransformed by the matrices in VR and/or VL;\n*          = 'S': compute selected right and/or left eigenvectors,\n*                 specified by the logical array SELECT.\n*\n*  SELECT  (input) LOGICAL array, dimension (N)\n*          If HOWMNY='S', SELECT specifies the eigenvectors to be\n*          computed.  The eigenvector corresponding to the j-th\n*          eigenvalue is computed if SELECT(j) = .TRUE..\n*          Not referenced if HOWMNY = 'A' or 'B'.\n*\n*  N       (input) INTEGER\n*          The order of the matrices S and P.  N >= 0.\n*\n*  S       (input) COMPLEX array, dimension (LDS,N)\n*          The upper triangular matrix S from a generalized Schur\n*          factorization, as computed by CHGEQZ.\n*\n*  LDS     (input) INTEGER\n*          The leading dimension of array S.  LDS >= max(1,N).\n*\n*  P       (input) COMPLEX array, dimension (LDP,N)\n*          The upper triangular matrix P from a generalized Schur\n*          factorization, as computed by CHGEQZ.  P must have real\n*          diagonal elements.\n*\n*  LDP     (input) INTEGER\n*          The leading dimension of array P.  LDP >= max(1,N).\n*\n*  VL      (input/output) COMPLEX array, dimension (LDVL,MM)\n*          On entry, if SIDE = 'L' or 'B' and HOWMNY = 'B', VL must\n*          contain an N-by-N matrix Q (usually the unitary matrix Q\n*          of left Schur vectors returned by CHGEQZ).\n*          On exit, if SIDE = 'L' or 'B', VL contains:\n*          if HOWMNY = 'A', the matrix Y of left eigenvectors of (S,P);\n*          if HOWMNY = 'B', the matrix Q*Y;\n*          if HOWMNY = 'S', the left eigenvectors of (S,P) specified by\n*                      SELECT, stored consecutively in the columns of\n*                      VL, in the same order as their eigenvalues.\n*          Not referenced if SIDE = 'R'.\n*\n*  LDVL    (input) INTEGER\n*          The leading dimension of array VL.  LDVL >= 1, and if\n*          SIDE = 'L' or 'l' or 'B' or 'b', LDVL >= N.\n*\n*  VR      (input/output) COMPLEX array, dimension (LDVR,MM)\n*          On entry, if SIDE = 'R' or 'B' and HOWMNY = 'B', VR must\n*          contain an N-by-N matrix Q (usually the unitary matrix Z\n*          of right Schur vectors returned by CHGEQZ).\n*          On exit, if SIDE = 'R' or 'B', VR contains:\n*          if HOWMNY = 'A', the matrix X of right eigenvectors of (S,P);\n*          if HOWMNY = 'B', the matrix Z*X;\n*          if HOWMNY = 'S', the right eigenvectors of (S,P) specified by\n*                      SELECT, stored consecutively in the columns of\n*                      VR, in the same order as their eigenvalues.\n*          Not referenced if SIDE = 'L'.\n*\n*  LDVR    (input) INTEGER\n*          The leading dimension of the array VR.  LDVR >= 1, and if\n*          SIDE = 'R' or 'B', LDVR >= N.\n*\n*  MM      (input) INTEGER\n*          The number of columns in the arrays VL and/or VR. MM >= M.\n*\n*  M       (output) INTEGER\n*          The number of columns in the arrays VL and/or VR actually\n*          used to store the eigenvectors.  If HOWMNY = 'A' or 'B', M\n*          is set to N.  Each selected eigenvector occupies one column.\n*\n*  WORK    (workspace) COMPLEX array, dimension (2*N)\n*\n*  RWORK   (workspace) REAL array, dimension (2*N)\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit.\n*          < 0:  if INFO = -i, the i-th argument had an illegal value.\n*\n\n*  =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  m, info, vl, vr = NumRu::Lapack.ctgevc( side, howmny, select, s, p, vl, vr, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 7)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 7)", argc);
  rblapack_side = argv[0];
  rblapack_howmny = argv[1];
  rblapack_select = argv[2];
  rblapack_s = argv[3];
  rblapack_p = argv[4];
  rblapack_vl = argv[5];
  rblapack_vr = argv[6];
  if (rb_options != Qnil) {
  }

  if (!NA_IsNArray(rblapack_vl))
    rb_raise(rb_eArgError, "vl (6th argument) must be NArray");
  if (NA_RANK(rblapack_vl) != 2)
    rb_raise(rb_eArgError, "rank of vl (6th argument) must be %d", 2);
  mm = NA_SHAPE1(rblapack_vl);
  ldvl = NA_SHAPE0(rblapack_vl);
  if (NA_TYPE(rblapack_vl) != NA_SCOMPLEX)
    rblapack_vl = na_change_type(rblapack_vl, NA_SCOMPLEX);
  vl = NA_PTR_TYPE(rblapack_vl, complex*);
  side = StringValueCStr(rblapack_side)[0];
  if (!NA_IsNArray(rblapack_p))
    rb_raise(rb_eArgError, "p (5th argument) must be NArray");
  if (NA_RANK(rblapack_p) != 2)
    rb_raise(rb_eArgError, "rank of p (5th argument) must be %d", 2);
  n = NA_SHAPE1(rblapack_p);
  ldp = NA_SHAPE0(rblapack_p);
  if (NA_TYPE(rblapack_p) != NA_SCOMPLEX)
    rblapack_p = na_change_type(rblapack_p, NA_SCOMPLEX);
  p = NA_PTR_TYPE(rblapack_p, complex*);
  if (!NA_IsNArray(rblapack_vr))
    rb_raise(rb_eArgError, "vr (7th argument) must be NArray");
  if (NA_RANK(rblapack_vr) != 2)
    rb_raise(rb_eArgError, "rank of vr (7th argument) must be %d", 2);
  if (NA_SHAPE1(rblapack_vr) != mm)
    rb_raise(rb_eRuntimeError, "shape 1 of vr must be the same as shape 1 of vl");
  ldvr = NA_SHAPE0(rblapack_vr);
  if (NA_TYPE(rblapack_vr) != NA_SCOMPLEX)
    rblapack_vr = na_change_type(rblapack_vr, NA_SCOMPLEX);
  vr = NA_PTR_TYPE(rblapack_vr, complex*);
  if (!NA_IsNArray(rblapack_s))
    rb_raise(rb_eArgError, "s (4th argument) must be NArray");
  if (NA_RANK(rblapack_s) != 2)
    rb_raise(rb_eArgError, "rank of s (4th argument) must be %d", 2);
  if (NA_SHAPE1(rblapack_s) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of s must be the same as shape 1 of p");
  lds = NA_SHAPE0(rblapack_s);
  if (NA_TYPE(rblapack_s) != NA_SCOMPLEX)
    rblapack_s = na_change_type(rblapack_s, NA_SCOMPLEX);
  s = NA_PTR_TYPE(rblapack_s, complex*);
  if (!NA_IsNArray(rblapack_select))
    rb_raise(rb_eArgError, "select (3th argument) must be NArray");
  if (NA_RANK(rblapack_select) != 1)
    rb_raise(rb_eArgError, "rank of select (3th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_select) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of select must be the same as shape 1 of p");
  if (NA_TYPE(rblapack_select) != NA_LINT)
    rblapack_select = na_change_type(rblapack_select, NA_LINT);
  select = NA_PTR_TYPE(rblapack_select, logical*);
  howmny = StringValueCStr(rblapack_howmny)[0];
  {
    int shape[2];
    shape[0] = ldvl;
    shape[1] = mm;
    rblapack_vl_out__ = na_make_object(NA_SCOMPLEX, 2, shape, cNArray);
  }
  vl_out__ = NA_PTR_TYPE(rblapack_vl_out__, complex*);
  MEMCPY(vl_out__, vl, complex, NA_TOTAL(rblapack_vl));
  rblapack_vl = rblapack_vl_out__;
  vl = vl_out__;
  {
    int shape[2];
    shape[0] = ldvr;
    shape[1] = mm;
    rblapack_vr_out__ = na_make_object(NA_SCOMPLEX, 2, shape, cNArray);
  }
  vr_out__ = NA_PTR_TYPE(rblapack_vr_out__, complex*);
  MEMCPY(vr_out__, vr, complex, NA_TOTAL(rblapack_vr));
  rblapack_vr = rblapack_vr_out__;
  vr = vr_out__;
  work = ALLOC_N(complex, (2*n));
  rwork = ALLOC_N(real, (2*n));

  ctgevc_(&side, &howmny, select, &n, s, &lds, p, &ldp, vl, &ldvl, vr, &ldvr, &mm, &m, work, rwork, &info);

  free(work);
  free(rwork);
  rblapack_m = INT2NUM(m);
  rblapack_info = INT2NUM(info);
  return rb_ary_new3(4, rblapack_m, rblapack_info, rblapack_vl, rblapack_vr);
}

void
init_lapack_ctgevc(VALUE mLapack){
  rb_define_module_function(mLapack, "ctgevc", rblapack_ctgevc, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}

#include "rb_lapack.h"

extern VOID stgevc_(char *side, char *howmny, logical *select, integer *n, real *s, integer *lds, real *p, integer *ldp, real *vl, integer *ldvl, real *vr, integer *ldvr, integer *mm, integer *m, real *work, integer *info);

static VALUE sHelp, sUsage;

static VALUE
rblapack_stgevc(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_side;
  char side; 
  VALUE rblapack_howmny;
  char howmny; 
  VALUE rblapack_select;
  logical *select; 
  VALUE rblapack_s;
  real *s; 
  VALUE rblapack_p;
  real *p; 
  VALUE rblapack_vl;
  real *vl; 
  VALUE rblapack_vr;
  real *vr; 
  VALUE rblapack_m;
  integer m; 
  VALUE rblapack_info;
  integer info; 
  VALUE rblapack_vl_out__;
  real *vl_out__;
  VALUE rblapack_vr_out__;
  real *vr_out__;
  real *work;

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
      printf("%s\n", "USAGE:\n  m, info, vl, vr = NumRu::Lapack.stgevc( side, howmny, select, s, p, vl, vr, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE STGEVC( SIDE, HOWMNY, SELECT, N, S, LDS, P, LDP, VL, LDVL, VR, LDVR, MM, M, WORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  STGEVC computes some or all of the right and/or left eigenvectors of\n*  a pair of real matrices (S,P), where S is a quasi-triangular matrix\n*  and P is upper triangular.  Matrix pairs of this type are produced by\n*  the generalized Schur factorization of a matrix pair (A,B):\n*\n*     A = Q*S*Z**T,  B = Q*P*Z**T\n*\n*  as computed by SGGHRD + SHGEQZ.\n*\n*  The right eigenvector x and the left eigenvector y of (S,P)\n*  corresponding to an eigenvalue w are defined by:\n*  \n*     S*x = w*P*x,  (y**H)*S = w*(y**H)*P,\n*  \n*  where y**H denotes the conjugate tranpose of y.\n*  The eigenvalues are not input to this routine, but are computed\n*  directly from the diagonal blocks of S and P.\n*  \n*  This routine returns the matrices X and/or Y of right and left\n*  eigenvectors of (S,P), or the products Z*X and/or Q*Y,\n*  where Z and Q are input matrices.\n*  If Q and Z are the orthogonal factors from the generalized Schur\n*  factorization of a matrix pair (A,B), then Z*X and Q*Y\n*  are the matrices of right and left eigenvectors of (A,B).\n* \n\n*  Arguments\n*  =========\n*\n*  SIDE    (input) CHARACTER*1\n*          = 'R': compute right eigenvectors only;\n*          = 'L': compute left eigenvectors only;\n*          = 'B': compute both right and left eigenvectors.\n*\n*  HOWMNY  (input) CHARACTER*1\n*          = 'A': compute all right and/or left eigenvectors;\n*          = 'B': compute all right and/or left eigenvectors,\n*                 backtransformed by the matrices in VR and/or VL;\n*          = 'S': compute selected right and/or left eigenvectors,\n*                 specified by the logical array SELECT.\n*\n*  SELECT  (input) LOGICAL array, dimension (N)\n*          If HOWMNY='S', SELECT specifies the eigenvectors to be\n*          computed.  If w(j) is a real eigenvalue, the corresponding\n*          real eigenvector is computed if SELECT(j) is .TRUE..\n*          If w(j) and w(j+1) are the real and imaginary parts of a\n*          complex eigenvalue, the corresponding complex eigenvector\n*          is computed if either SELECT(j) or SELECT(j+1) is .TRUE.,\n*          and on exit SELECT(j) is set to .TRUE. and SELECT(j+1) is\n*          set to .FALSE..\n*          Not referenced if HOWMNY = 'A' or 'B'.\n*\n*  N       (input) INTEGER\n*          The order of the matrices S and P.  N >= 0.\n*\n*  S       (input) REAL array, dimension (LDS,N)\n*          The upper quasi-triangular matrix S from a generalized Schur\n*          factorization, as computed by SHGEQZ.\n*\n*  LDS     (input) INTEGER\n*          The leading dimension of array S.  LDS >= max(1,N).\n*\n*  P       (input) REAL array, dimension (LDP,N)\n*          The upper triangular matrix P from a generalized Schur\n*          factorization, as computed by SHGEQZ.\n*          2-by-2 diagonal blocks of P corresponding to 2-by-2 blocks\n*          of S must be in positive diagonal form.\n*\n*  LDP     (input) INTEGER\n*          The leading dimension of array P.  LDP >= max(1,N).\n*\n*  VL      (input/output) REAL array, dimension (LDVL,MM)\n*          On entry, if SIDE = 'L' or 'B' and HOWMNY = 'B', VL must\n*          contain an N-by-N matrix Q (usually the orthogonal matrix Q\n*          of left Schur vectors returned by SHGEQZ).\n*          On exit, if SIDE = 'L' or 'B', VL contains:\n*          if HOWMNY = 'A', the matrix Y of left eigenvectors of (S,P);\n*          if HOWMNY = 'B', the matrix Q*Y;\n*          if HOWMNY = 'S', the left eigenvectors of (S,P) specified by\n*                      SELECT, stored consecutively in the columns of\n*                      VL, in the same order as their eigenvalues.\n*\n*          A complex eigenvector corresponding to a complex eigenvalue\n*          is stored in two consecutive columns, the first holding the\n*          real part, and the second the imaginary part.\n*\n*          Not referenced if SIDE = 'R'.\n*\n*  LDVL    (input) INTEGER\n*          The leading dimension of array VL.  LDVL >= 1, and if\n*          SIDE = 'L' or 'B', LDVL >= N.\n*\n*  VR      (input/output) REAL array, dimension (LDVR,MM)\n*          On entry, if SIDE = 'R' or 'B' and HOWMNY = 'B', VR must\n*          contain an N-by-N matrix Z (usually the orthogonal matrix Z\n*          of right Schur vectors returned by SHGEQZ).\n*\n*          On exit, if SIDE = 'R' or 'B', VR contains:\n*          if HOWMNY = 'A', the matrix X of right eigenvectors of (S,P);\n*          if HOWMNY = 'B' or 'b', the matrix Z*X;\n*          if HOWMNY = 'S' or 's', the right eigenvectors of (S,P)\n*                      specified by SELECT, stored consecutively in the\n*                      columns of VR, in the same order as their\n*                      eigenvalues.\n*\n*          A complex eigenvector corresponding to a complex eigenvalue\n*          is stored in two consecutive columns, the first holding the\n*          real part and the second the imaginary part.\n*          \n*          Not referenced if SIDE = 'L'.\n*\n*  LDVR    (input) INTEGER\n*          The leading dimension of the array VR.  LDVR >= 1, and if\n*          SIDE = 'R' or 'B', LDVR >= N.\n*\n*  MM      (input) INTEGER\n*          The number of columns in the arrays VL and/or VR. MM >= M.\n*\n*  M       (output) INTEGER\n*          The number of columns in the arrays VL and/or VR actually\n*          used to store the eigenvectors.  If HOWMNY = 'A' or 'B', M\n*          is set to N.  Each selected real eigenvector occupies one\n*          column and each selected complex eigenvector occupies two\n*          columns.\n*\n*  WORK    (workspace) REAL array, dimension (6*N)\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit.\n*          < 0:  if INFO = -i, the i-th argument had an illegal value.\n*          > 0:  the 2-by-2 block (INFO:INFO+1) does not have a complex\n*                eigenvalue.\n*\n\n*  Further Details\n*  ===============\n*\n*  Allocation of workspace:\n*  ---------- -- ---------\n*\n*     WORK( j ) = 1-norm of j-th column of A, above the diagonal\n*     WORK( N+j ) = 1-norm of j-th column of B, above the diagonal\n*     WORK( 2*N+1:3*N ) = real part of eigenvector\n*     WORK( 3*N+1:4*N ) = imaginary part of eigenvector\n*     WORK( 4*N+1:5*N ) = real part of back-transformed eigenvector\n*     WORK( 5*N+1:6*N ) = imaginary part of back-transformed eigenvector\n*\n*  Rowwise vs. columnwise solution methods:\n*  ------- --  ---------- -------- -------\n*\n*  Finding a generalized eigenvector consists basically of solving the\n*  singular triangular system\n*\n*   (A - w B) x = 0     (for right) or:   (A - w B)**H y = 0  (for left)\n*\n*  Consider finding the i-th right eigenvector (assume all eigenvalues\n*  are real). The equation to be solved is:\n*       n                   i\n*  0 = sum  C(j,k) v(k)  = sum  C(j,k) v(k)     for j = i,. . .,1\n*      k=j                 k=j\n*\n*  where  C = (A - w B)  (The components v(i+1:n) are 0.)\n*\n*  The \"rowwise\" method is:\n*\n*  (1)  v(i) := 1\n*  for j = i-1,. . .,1:\n*                          i\n*      (2) compute  s = - sum C(j,k) v(k)   and\n*                        k=j+1\n*\n*      (3) v(j) := s / C(j,j)\n*\n*  Step 2 is sometimes called the \"dot product\" step, since it is an\n*  inner product between the j-th row and the portion of the eigenvector\n*  that has been computed so far.\n*\n*  The \"columnwise\" method consists basically in doing the sums\n*  for all the rows in parallel.  As each v(j) is computed, the\n*  contribution of v(j) times the j-th column of C is added to the\n*  partial sums.  Since FORTRAN arrays are stored columnwise, this has\n*  the advantage that at each step, the elements of C that are accessed\n*  are adjacent to one another, whereas with the rowwise method, the\n*  elements accessed at a step are spaced LDS (and LDP) words apart.\n*\n*  When finding left eigenvectors, the matrix in question is the\n*  transpose of the one in storage, so the rowwise method then\n*  actually accesses columns of A and B at each step, and so is the\n*  preferred method.\n*\n*  =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  m, info, vl, vr = NumRu::Lapack.stgevc( side, howmny, select, s, p, vl, vr, [:usage => usage, :help => help])\n");
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
  if (NA_TYPE(rblapack_vl) != NA_SFLOAT)
    rblapack_vl = na_change_type(rblapack_vl, NA_SFLOAT);
  vl = NA_PTR_TYPE(rblapack_vl, real*);
  side = StringValueCStr(rblapack_side)[0];
  if (!NA_IsNArray(rblapack_p))
    rb_raise(rb_eArgError, "p (5th argument) must be NArray");
  if (NA_RANK(rblapack_p) != 2)
    rb_raise(rb_eArgError, "rank of p (5th argument) must be %d", 2);
  n = NA_SHAPE1(rblapack_p);
  ldp = NA_SHAPE0(rblapack_p);
  if (NA_TYPE(rblapack_p) != NA_SFLOAT)
    rblapack_p = na_change_type(rblapack_p, NA_SFLOAT);
  p = NA_PTR_TYPE(rblapack_p, real*);
  if (!NA_IsNArray(rblapack_vr))
    rb_raise(rb_eArgError, "vr (7th argument) must be NArray");
  if (NA_RANK(rblapack_vr) != 2)
    rb_raise(rb_eArgError, "rank of vr (7th argument) must be %d", 2);
  if (NA_SHAPE1(rblapack_vr) != mm)
    rb_raise(rb_eRuntimeError, "shape 1 of vr must be the same as shape 1 of vl");
  ldvr = NA_SHAPE0(rblapack_vr);
  if (NA_TYPE(rblapack_vr) != NA_SFLOAT)
    rblapack_vr = na_change_type(rblapack_vr, NA_SFLOAT);
  vr = NA_PTR_TYPE(rblapack_vr, real*);
  if (!NA_IsNArray(rblapack_s))
    rb_raise(rb_eArgError, "s (4th argument) must be NArray");
  if (NA_RANK(rblapack_s) != 2)
    rb_raise(rb_eArgError, "rank of s (4th argument) must be %d", 2);
  if (NA_SHAPE1(rblapack_s) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of s must be the same as shape 1 of p");
  lds = NA_SHAPE0(rblapack_s);
  if (NA_TYPE(rblapack_s) != NA_SFLOAT)
    rblapack_s = na_change_type(rblapack_s, NA_SFLOAT);
  s = NA_PTR_TYPE(rblapack_s, real*);
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
    rblapack_vl_out__ = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  vl_out__ = NA_PTR_TYPE(rblapack_vl_out__, real*);
  MEMCPY(vl_out__, vl, real, NA_TOTAL(rblapack_vl));
  rblapack_vl = rblapack_vl_out__;
  vl = vl_out__;
  {
    int shape[2];
    shape[0] = ldvr;
    shape[1] = mm;
    rblapack_vr_out__ = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  vr_out__ = NA_PTR_TYPE(rblapack_vr_out__, real*);
  MEMCPY(vr_out__, vr, real, NA_TOTAL(rblapack_vr));
  rblapack_vr = rblapack_vr_out__;
  vr = vr_out__;
  work = ALLOC_N(real, (6*n));

  stgevc_(&side, &howmny, select, &n, s, &lds, p, &ldp, vl, &ldvl, vr, &ldvr, &mm, &m, work, &info);

  free(work);
  rblapack_m = INT2NUM(m);
  rblapack_info = INT2NUM(info);
  return rb_ary_new3(4, rblapack_m, rblapack_info, rblapack_vl, rblapack_vr);
}

void
init_lapack_stgevc(VALUE mLapack){
  rb_define_module_function(mLapack, "stgevc", rblapack_stgevc, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}

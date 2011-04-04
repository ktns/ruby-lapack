#include "rb_lapack.h"

extern VOID dtgevc_(char *side, char *howmny, logical *select, integer *n, doublereal *s, integer *lds, doublereal *p, integer *ldp, doublereal *vl, integer *ldvl, doublereal *vr, integer *ldvr, integer *mm, integer *m, doublereal *work, integer *info);

static VALUE
rb_dtgevc(int argc, VALUE *argv, VALUE self){
  VALUE rb_side;
  char side; 
  VALUE rb_howmny;
  char howmny; 
  VALUE rb_select;
  logical *select; 
  VALUE rb_s;
  doublereal *s; 
  VALUE rb_p;
  doublereal *p; 
  VALUE rb_vl;
  doublereal *vl; 
  VALUE rb_vr;
  doublereal *vr; 
  VALUE rb_m;
  integer m; 
  VALUE rb_info;
  integer info; 
  VALUE rb_vl_out__;
  doublereal *vl_out__;
  VALUE rb_vr_out__;
  doublereal *vr_out__;
  doublereal *work;

  integer n;
  integer lds;
  integer ldp;
  integer ldvl;
  integer mm;
  integer ldvr;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  m, info, vl, vr = NumRu::Lapack.dtgevc( side, howmny, select, s, p, vl, vr)\n    or\n  NumRu::Lapack.dtgevc  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE DTGEVC( SIDE, HOWMNY, SELECT, N, S, LDS, P, LDP, VL, LDVL, VR, LDVR, MM, M, WORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  DTGEVC computes some or all of the right and/or left eigenvectors of\n*  a pair of real matrices (S,P), where S is a quasi-triangular matrix\n*  and P is upper triangular.  Matrix pairs of this type are produced by\n*  the generalized Schur factorization of a matrix pair (A,B):\n*\n*     A = Q*S*Z**T,  B = Q*P*Z**T\n*\n*  as computed by DGGHRD + DHGEQZ.\n*\n*  The right eigenvector x and the left eigenvector y of (S,P)\n*  corresponding to an eigenvalue w are defined by:\n*  \n*     S*x = w*P*x,  (y**H)*S = w*(y**H)*P,\n*  \n*  where y**H denotes the conjugate tranpose of y.\n*  The eigenvalues are not input to this routine, but are computed\n*  directly from the diagonal blocks of S and P.\n*  \n*  This routine returns the matrices X and/or Y of right and left\n*  eigenvectors of (S,P), or the products Z*X and/or Q*Y,\n*  where Z and Q are input matrices.\n*  If Q and Z are the orthogonal factors from the generalized Schur\n*  factorization of a matrix pair (A,B), then Z*X and Q*Y\n*  are the matrices of right and left eigenvectors of (A,B).\n* \n\n*  Arguments\n*  =========\n*\n*  SIDE    (input) CHARACTER*1\n*          = 'R': compute right eigenvectors only;\n*          = 'L': compute left eigenvectors only;\n*          = 'B': compute both right and left eigenvectors.\n*\n*  HOWMNY  (input) CHARACTER*1\n*          = 'A': compute all right and/or left eigenvectors;\n*          = 'B': compute all right and/or left eigenvectors,\n*                 backtransformed by the matrices in VR and/or VL;\n*          = 'S': compute selected right and/or left eigenvectors,\n*                 specified by the logical array SELECT.\n*\n*  SELECT  (input) LOGICAL array, dimension (N)\n*          If HOWMNY='S', SELECT specifies the eigenvectors to be\n*          computed.  If w(j) is a real eigenvalue, the corresponding\n*          real eigenvector is computed if SELECT(j) is .TRUE..\n*          If w(j) and w(j+1) are the real and imaginary parts of a\n*          complex eigenvalue, the corresponding complex eigenvector\n*          is computed if either SELECT(j) or SELECT(j+1) is .TRUE.,\n*          and on exit SELECT(j) is set to .TRUE. and SELECT(j+1) is\n*          set to .FALSE..\n*          Not referenced if HOWMNY = 'A' or 'B'.\n*\n*  N       (input) INTEGER\n*          The order of the matrices S and P.  N >= 0.\n*\n*  S       (input) DOUBLE PRECISION array, dimension (LDS,N)\n*          The upper quasi-triangular matrix S from a generalized Schur\n*          factorization, as computed by DHGEQZ.\n*\n*  LDS     (input) INTEGER\n*          The leading dimension of array S.  LDS >= max(1,N).\n*\n*  P       (input) DOUBLE PRECISION array, dimension (LDP,N)\n*          The upper triangular matrix P from a generalized Schur\n*          factorization, as computed by DHGEQZ.\n*          2-by-2 diagonal blocks of P corresponding to 2-by-2 blocks\n*          of S must be in positive diagonal form.\n*\n*  LDP     (input) INTEGER\n*          The leading dimension of array P.  LDP >= max(1,N).\n*\n*  VL      (input/output) DOUBLE PRECISION array, dimension (LDVL,MM)\n*          On entry, if SIDE = 'L' or 'B' and HOWMNY = 'B', VL must\n*          contain an N-by-N matrix Q (usually the orthogonal matrix Q\n*          of left Schur vectors returned by DHGEQZ).\n*          On exit, if SIDE = 'L' or 'B', VL contains:\n*          if HOWMNY = 'A', the matrix Y of left eigenvectors of (S,P);\n*          if HOWMNY = 'B', the matrix Q*Y;\n*          if HOWMNY = 'S', the left eigenvectors of (S,P) specified by\n*                      SELECT, stored consecutively in the columns of\n*                      VL, in the same order as their eigenvalues.\n*\n*          A complex eigenvector corresponding to a complex eigenvalue\n*          is stored in two consecutive columns, the first holding the\n*          real part, and the second the imaginary part.\n*\n*          Not referenced if SIDE = 'R'.\n*\n*  LDVL    (input) INTEGER\n*          The leading dimension of array VL.  LDVL >= 1, and if\n*          SIDE = 'L' or 'B', LDVL >= N.\n*\n*  VR      (input/output) DOUBLE PRECISION array, dimension (LDVR,MM)\n*          On entry, if SIDE = 'R' or 'B' and HOWMNY = 'B', VR must\n*          contain an N-by-N matrix Z (usually the orthogonal matrix Z\n*          of right Schur vectors returned by DHGEQZ).\n*\n*          On exit, if SIDE = 'R' or 'B', VR contains:\n*          if HOWMNY = 'A', the matrix X of right eigenvectors of (S,P);\n*          if HOWMNY = 'B' or 'b', the matrix Z*X;\n*          if HOWMNY = 'S' or 's', the right eigenvectors of (S,P)\n*                      specified by SELECT, stored consecutively in the\n*                      columns of VR, in the same order as their\n*                      eigenvalues.\n*\n*          A complex eigenvector corresponding to a complex eigenvalue\n*          is stored in two consecutive columns, the first holding the\n*          real part and the second the imaginary part.\n*          \n*          Not referenced if SIDE = 'L'.\n*\n*  LDVR    (input) INTEGER\n*          The leading dimension of the array VR.  LDVR >= 1, and if\n*          SIDE = 'R' or 'B', LDVR >= N.\n*\n*  MM      (input) INTEGER\n*          The number of columns in the arrays VL and/or VR. MM >= M.\n*\n*  M       (output) INTEGER\n*          The number of columns in the arrays VL and/or VR actually\n*          used to store the eigenvectors.  If HOWMNY = 'A' or 'B', M\n*          is set to N.  Each selected real eigenvector occupies one\n*          column and each selected complex eigenvector occupies two\n*          columns.\n*\n*  WORK    (workspace) DOUBLE PRECISION array, dimension (6*N)\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit.\n*          < 0:  if INFO = -i, the i-th argument had an illegal value.\n*          > 0:  the 2-by-2 block (INFO:INFO+1) does not have a complex\n*                eigenvalue.\n*\n\n*  Further Details\n*  ===============\n*\n*  Allocation of workspace:\n*  ---------- -- ---------\n*\n*     WORK( j ) = 1-norm of j-th column of A, above the diagonal\n*     WORK( N+j ) = 1-norm of j-th column of B, above the diagonal\n*     WORK( 2*N+1:3*N ) = real part of eigenvector\n*     WORK( 3*N+1:4*N ) = imaginary part of eigenvector\n*     WORK( 4*N+1:5*N ) = real part of back-transformed eigenvector\n*     WORK( 5*N+1:6*N ) = imaginary part of back-transformed eigenvector\n*\n*  Rowwise vs. columnwise solution methods:\n*  ------- --  ---------- -------- -------\n*\n*  Finding a generalized eigenvector consists basically of solving the\n*  singular triangular system\n*\n*   (A - w B) x = 0     (for right) or:   (A - w B)**H y = 0  (for left)\n*\n*  Consider finding the i-th right eigenvector (assume all eigenvalues\n*  are real). The equation to be solved is:\n*       n                   i\n*  0 = sum  C(j,k) v(k)  = sum  C(j,k) v(k)     for j = i,. . .,1\n*      k=j                 k=j\n*\n*  where  C = (A - w B)  (The components v(i+1:n) are 0.)\n*\n*  The \"rowwise\" method is:\n*\n*  (1)  v(i) := 1\n*  for j = i-1,. . .,1:\n*                          i\n*      (2) compute  s = - sum C(j,k) v(k)   and\n*                        k=j+1\n*\n*      (3) v(j) := s / C(j,j)\n*\n*  Step 2 is sometimes called the \"dot product\" step, since it is an\n*  inner product between the j-th row and the portion of the eigenvector\n*  that has been computed so far.\n*\n*  The \"columnwise\" method consists basically in doing the sums\n*  for all the rows in parallel.  As each v(j) is computed, the\n*  contribution of v(j) times the j-th column of C is added to the\n*  partial sums.  Since FORTRAN arrays are stored columnwise, this has\n*  the advantage that at each step, the elements of C that are accessed\n*  are adjacent to one another, whereas with the rowwise method, the\n*  elements accessed at a step are spaced LDS (and LDP) words apart.\n*\n*  When finding left eigenvectors, the matrix in question is the\n*  transpose of the one in storage, so the rowwise method then\n*  actually accesses columns of A and B at each step, and so is the\n*  preferred method.\n*\n*  =====================================================================\n*\n\n");
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
  if (NA_TYPE(rb_vl) != NA_DFLOAT)
    rb_vl = na_change_type(rb_vl, NA_DFLOAT);
  vl = NA_PTR_TYPE(rb_vl, doublereal*);
  side = StringValueCStr(rb_side)[0];
  if (!NA_IsNArray(rb_p))
    rb_raise(rb_eArgError, "p (5th argument) must be NArray");
  if (NA_RANK(rb_p) != 2)
    rb_raise(rb_eArgError, "rank of p (5th argument) must be %d", 2);
  n = NA_SHAPE1(rb_p);
  ldp = NA_SHAPE0(rb_p);
  if (NA_TYPE(rb_p) != NA_DFLOAT)
    rb_p = na_change_type(rb_p, NA_DFLOAT);
  p = NA_PTR_TYPE(rb_p, doublereal*);
  if (!NA_IsNArray(rb_vr))
    rb_raise(rb_eArgError, "vr (7th argument) must be NArray");
  if (NA_RANK(rb_vr) != 2)
    rb_raise(rb_eArgError, "rank of vr (7th argument) must be %d", 2);
  if (NA_SHAPE1(rb_vr) != mm)
    rb_raise(rb_eRuntimeError, "shape 1 of vr must be the same as shape 1 of vl");
  ldvr = NA_SHAPE0(rb_vr);
  if (NA_TYPE(rb_vr) != NA_DFLOAT)
    rb_vr = na_change_type(rb_vr, NA_DFLOAT);
  vr = NA_PTR_TYPE(rb_vr, doublereal*);
  if (!NA_IsNArray(rb_s))
    rb_raise(rb_eArgError, "s (4th argument) must be NArray");
  if (NA_RANK(rb_s) != 2)
    rb_raise(rb_eArgError, "rank of s (4th argument) must be %d", 2);
  if (NA_SHAPE1(rb_s) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of s must be the same as shape 1 of p");
  lds = NA_SHAPE0(rb_s);
  if (NA_TYPE(rb_s) != NA_DFLOAT)
    rb_s = na_change_type(rb_s, NA_DFLOAT);
  s = NA_PTR_TYPE(rb_s, doublereal*);
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
  work = ALLOC_N(doublereal, (6*n));

  dtgevc_(&side, &howmny, select, &n, s, &lds, p, &ldp, vl, &ldvl, vr, &ldvr, &mm, &m, work, &info);

  free(work);
  rb_m = INT2NUM(m);
  rb_info = INT2NUM(info);
  return rb_ary_new3(4, rb_m, rb_info, rb_vl, rb_vr);
}

void
init_lapack_dtgevc(VALUE mLapack){
  rb_define_module_function(mLapack, "dtgevc", rb_dtgevc, -1);
}

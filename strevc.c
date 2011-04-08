#include "rb_lapack.h"

extern VOID strevc_(char *side, char *howmny, logical *select, integer *n, real *t, integer *ldt, real *vl, integer *ldvl, real *vr, integer *ldvr, integer *mm, integer *m, real *work, integer *info);

static VALUE sHelp, sUsage;

static VALUE
rblapack_strevc(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_side;
  char side; 
  VALUE rblapack_howmny;
  char howmny; 
  VALUE rblapack_select;
  logical *select; 
  VALUE rblapack_t;
  real *t; 
  VALUE rblapack_vl;
  real *vl; 
  VALUE rblapack_vr;
  real *vr; 
  VALUE rblapack_m;
  integer m; 
  VALUE rblapack_info;
  integer info; 
  VALUE rblapack_select_out__;
  logical *select_out__;
  VALUE rblapack_vl_out__;
  real *vl_out__;
  VALUE rblapack_vr_out__;
  real *vr_out__;
  real *work;

  integer n;
  integer ldt;
  integer ldvl;
  integer mm;
  integer ldvr;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  m, info, select, vl, vr = NumRu::Lapack.strevc( side, howmny, select, t, vl, vr, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE STREVC( SIDE, HOWMNY, SELECT, N, T, LDT, VL, LDVL, VR, LDVR, MM, M, WORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  STREVC computes some or all of the right and/or left eigenvectors of\n*  a real upper quasi-triangular matrix T.\n*  Matrices of this type are produced by the Schur factorization of\n*  a real general matrix:  A = Q*T*Q**T, as computed by SHSEQR.\n*  \n*  The right eigenvector x and the left eigenvector y of T corresponding\n*  to an eigenvalue w are defined by:\n*  \n*     T*x = w*x,     (y**H)*T = w*(y**H)\n*  \n*  where y**H denotes the conjugate transpose of y.\n*  The eigenvalues are not input to this routine, but are read directly\n*  from the diagonal blocks of T.\n*  \n*  This routine returns the matrices X and/or Y of right and left\n*  eigenvectors of T, or the products Q*X and/or Q*Y, where Q is an\n*  input matrix.  If Q is the orthogonal factor that reduces a matrix\n*  A to Schur form T, then Q*X and Q*Y are the matrices of right and\n*  left eigenvectors of A.\n*\n\n*  Arguments\n*  =========\n*\n*  SIDE    (input) CHARACTER*1\n*          = 'R':  compute right eigenvectors only;\n*          = 'L':  compute left eigenvectors only;\n*          = 'B':  compute both right and left eigenvectors.\n*\n*  HOWMNY  (input) CHARACTER*1\n*          = 'A':  compute all right and/or left eigenvectors;\n*          = 'B':  compute all right and/or left eigenvectors,\n*                  backtransformed by the matrices in VR and/or VL;\n*          = 'S':  compute selected right and/or left eigenvectors,\n*                  as indicated by the logical array SELECT.\n*\n*  SELECT  (input/output) LOGICAL array, dimension (N)\n*          If HOWMNY = 'S', SELECT specifies the eigenvectors to be\n*          computed.\n*          If w(j) is a real eigenvalue, the corresponding real\n*          eigenvector is computed if SELECT(j) is .TRUE..\n*          If w(j) and w(j+1) are the real and imaginary parts of a\n*          complex eigenvalue, the corresponding complex eigenvector is\n*          computed if either SELECT(j) or SELECT(j+1) is .TRUE., and\n*          on exit SELECT(j) is set to .TRUE. and SELECT(j+1) is set to\n*          .FALSE..\n*          Not referenced if HOWMNY = 'A' or 'B'.\n*\n*  N       (input) INTEGER\n*          The order of the matrix T. N >= 0.\n*\n*  T       (input) REAL array, dimension (LDT,N)\n*          The upper quasi-triangular matrix T in Schur canonical form.\n*\n*  LDT     (input) INTEGER\n*          The leading dimension of the array T. LDT >= max(1,N).\n*\n*  VL      (input/output) REAL array, dimension (LDVL,MM)\n*          On entry, if SIDE = 'L' or 'B' and HOWMNY = 'B', VL must\n*          contain an N-by-N matrix Q (usually the orthogonal matrix Q\n*          of Schur vectors returned by SHSEQR).\n*          On exit, if SIDE = 'L' or 'B', VL contains:\n*          if HOWMNY = 'A', the matrix Y of left eigenvectors of T;\n*          if HOWMNY = 'B', the matrix Q*Y;\n*          if HOWMNY = 'S', the left eigenvectors of T specified by\n*                           SELECT, stored consecutively in the columns\n*                           of VL, in the same order as their\n*                           eigenvalues.\n*          A complex eigenvector corresponding to a complex eigenvalue\n*          is stored in two consecutive columns, the first holding the\n*          real part, and the second the imaginary part.\n*          Not referenced if SIDE = 'R'.\n*\n*  LDVL    (input) INTEGER\n*          The leading dimension of the array VL.  LDVL >= 1, and if\n*          SIDE = 'L' or 'B', LDVL >= N.\n*\n*  VR      (input/output) REAL array, dimension (LDVR,MM)\n*          On entry, if SIDE = 'R' or 'B' and HOWMNY = 'B', VR must\n*          contain an N-by-N matrix Q (usually the orthogonal matrix Q\n*          of Schur vectors returned by SHSEQR).\n*          On exit, if SIDE = 'R' or 'B', VR contains:\n*          if HOWMNY = 'A', the matrix X of right eigenvectors of T;\n*          if HOWMNY = 'B', the matrix Q*X;\n*          if HOWMNY = 'S', the right eigenvectors of T specified by\n*                           SELECT, stored consecutively in the columns\n*                           of VR, in the same order as their\n*                           eigenvalues.\n*          A complex eigenvector corresponding to a complex eigenvalue\n*          is stored in two consecutive columns, the first holding the\n*          real part and the second the imaginary part.\n*          Not referenced if SIDE = 'L'.\n*\n*  LDVR    (input) INTEGER\n*          The leading dimension of the array VR.  LDVR >= 1, and if\n*          SIDE = 'R' or 'B', LDVR >= N.\n*\n*  MM      (input) INTEGER\n*          The number of columns in the arrays VL and/or VR. MM >= M.\n*\n*  M       (output) INTEGER\n*          The number of columns in the arrays VL and/or VR actually\n*          used to store the eigenvectors.\n*          If HOWMNY = 'A' or 'B', M is set to N.\n*          Each selected real eigenvector occupies one column and each\n*          selected complex eigenvector occupies two columns.\n*\n*  WORK    (workspace) REAL array, dimension (3*N)\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit\n*          < 0:  if INFO = -i, the i-th argument had an illegal value\n*\n\n*  Further Details\n*  ===============\n*\n*  The algorithm used in this program is basically backward (forward)\n*  substitution, with scaling to make the the code robust against\n*  possible overflow.\n*\n*  Each eigenvector is normalized so that the element of largest\n*  magnitude has magnitude 1; here the magnitude of a complex number\n*  (x,y) is taken to be |x| + |y|.\n*\n*  =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  m, info, select, vl, vr = NumRu::Lapack.strevc( side, howmny, select, t, vl, vr, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 6)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 6)", argc);
  rblapack_side = argv[0];
  rblapack_howmny = argv[1];
  rblapack_select = argv[2];
  rblapack_t = argv[3];
  rblapack_vl = argv[4];
  rblapack_vr = argv[5];
  if (rb_options != Qnil) {
  }

  if (!NA_IsNArray(rblapack_vl))
    rb_raise(rb_eArgError, "vl (5th argument) must be NArray");
  if (NA_RANK(rblapack_vl) != 2)
    rb_raise(rb_eArgError, "rank of vl (5th argument) must be %d", 2);
  mm = NA_SHAPE1(rblapack_vl);
  ldvl = NA_SHAPE0(rblapack_vl);
  if (NA_TYPE(rblapack_vl) != NA_SFLOAT)
    rblapack_vl = na_change_type(rblapack_vl, NA_SFLOAT);
  vl = NA_PTR_TYPE(rblapack_vl, real*);
  side = StringValueCStr(rblapack_side)[0];
  if (!NA_IsNArray(rblapack_vr))
    rb_raise(rb_eArgError, "vr (6th argument) must be NArray");
  if (NA_RANK(rblapack_vr) != 2)
    rb_raise(rb_eArgError, "rank of vr (6th argument) must be %d", 2);
  if (NA_SHAPE1(rblapack_vr) != mm)
    rb_raise(rb_eRuntimeError, "shape 1 of vr must be the same as shape 1 of vl");
  ldvr = NA_SHAPE0(rblapack_vr);
  if (NA_TYPE(rblapack_vr) != NA_SFLOAT)
    rblapack_vr = na_change_type(rblapack_vr, NA_SFLOAT);
  vr = NA_PTR_TYPE(rblapack_vr, real*);
  if (!NA_IsNArray(rblapack_select))
    rb_raise(rb_eArgError, "select (3th argument) must be NArray");
  if (NA_RANK(rblapack_select) != 1)
    rb_raise(rb_eArgError, "rank of select (3th argument) must be %d", 1);
  n = NA_SHAPE0(rblapack_select);
  if (NA_TYPE(rblapack_select) != NA_LINT)
    rblapack_select = na_change_type(rblapack_select, NA_LINT);
  select = NA_PTR_TYPE(rblapack_select, logical*);
  if (!NA_IsNArray(rblapack_t))
    rb_raise(rb_eArgError, "t (4th argument) must be NArray");
  if (NA_RANK(rblapack_t) != 2)
    rb_raise(rb_eArgError, "rank of t (4th argument) must be %d", 2);
  if (NA_SHAPE1(rblapack_t) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of t must be the same as shape 0 of select");
  ldt = NA_SHAPE0(rblapack_t);
  if (NA_TYPE(rblapack_t) != NA_SFLOAT)
    rblapack_t = na_change_type(rblapack_t, NA_SFLOAT);
  t = NA_PTR_TYPE(rblapack_t, real*);
  howmny = StringValueCStr(rblapack_howmny)[0];
  {
    int shape[1];
    shape[0] = n;
    rblapack_select_out__ = na_make_object(NA_LINT, 1, shape, cNArray);
  }
  select_out__ = NA_PTR_TYPE(rblapack_select_out__, logical*);
  MEMCPY(select_out__, select, logical, NA_TOTAL(rblapack_select));
  rblapack_select = rblapack_select_out__;
  select = select_out__;
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
  work = ALLOC_N(real, (3*n));

  strevc_(&side, &howmny, select, &n, t, &ldt, vl, &ldvl, vr, &ldvr, &mm, &m, work, &info);

  free(work);
  rblapack_m = INT2NUM(m);
  rblapack_info = INT2NUM(info);
  return rb_ary_new3(5, rblapack_m, rblapack_info, rblapack_select, rblapack_vl, rblapack_vr);
}

void
init_lapack_strevc(VALUE mLapack){
  rb_define_module_function(mLapack, "strevc", rblapack_strevc, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}

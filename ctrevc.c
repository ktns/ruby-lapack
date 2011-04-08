#include "rb_lapack.h"

extern VOID ctrevc_(char *side, char *howmny, logical *select, integer *n, complex *t, integer *ldt, complex *vl, integer *ldvl, complex *vr, integer *ldvr, integer *mm, integer *m, complex *work, real *rwork, integer *info);

static VALUE sHelp, sUsage;

static VALUE
rblapack_ctrevc(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_side;
  char side; 
  VALUE rblapack_howmny;
  char howmny; 
  VALUE rblapack_select;
  logical *select; 
  VALUE rblapack_t;
  complex *t; 
  VALUE rblapack_vl;
  complex *vl; 
  VALUE rblapack_vr;
  complex *vr; 
  VALUE rblapack_m;
  integer m; 
  VALUE rblapack_info;
  integer info; 
  VALUE rblapack_t_out__;
  complex *t_out__;
  VALUE rblapack_vl_out__;
  complex *vl_out__;
  VALUE rblapack_vr_out__;
  complex *vr_out__;
  complex *work;
  real *rwork;

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
      printf("%s\n", "USAGE:\n  m, info, t, vl, vr = NumRu::Lapack.ctrevc( side, howmny, select, t, vl, vr, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE CTREVC( SIDE, HOWMNY, SELECT, N, T, LDT, VL, LDVL, VR, LDVR, MM, M, WORK, RWORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  CTREVC computes some or all of the right and/or left eigenvectors of\n*  a complex upper triangular matrix T.\n*  Matrices of this type are produced by the Schur factorization of\n*  a complex general matrix:  A = Q*T*Q**H, as computed by CHSEQR.\n*  \n*  The right eigenvector x and the left eigenvector y of T corresponding\n*  to an eigenvalue w are defined by:\n*  \n*               T*x = w*x,     (y**H)*T = w*(y**H)\n*  \n*  where y**H denotes the conjugate transpose of the vector y.\n*  The eigenvalues are not input to this routine, but are read directly\n*  from the diagonal of T.\n*  \n*  This routine returns the matrices X and/or Y of right and left\n*  eigenvectors of T, or the products Q*X and/or Q*Y, where Q is an\n*  input matrix.  If Q is the unitary factor that reduces a matrix A to\n*  Schur form T, then Q*X and Q*Y are the matrices of right and left\n*  eigenvectors of A.\n*\n\n*  Arguments\n*  =========\n*\n*  SIDE    (input) CHARACTER*1\n*          = 'R':  compute right eigenvectors only;\n*          = 'L':  compute left eigenvectors only;\n*          = 'B':  compute both right and left eigenvectors.\n*\n*  HOWMNY  (input) CHARACTER*1\n*          = 'A':  compute all right and/or left eigenvectors;\n*          = 'B':  compute all right and/or left eigenvectors,\n*                  backtransformed using the matrices supplied in\n*                  VR and/or VL;\n*          = 'S':  compute selected right and/or left eigenvectors,\n*                  as indicated by the logical array SELECT.\n*\n*  SELECT  (input) LOGICAL array, dimension (N)\n*          If HOWMNY = 'S', SELECT specifies the eigenvectors to be\n*          computed.\n*          The eigenvector corresponding to the j-th eigenvalue is\n*          computed if SELECT(j) = .TRUE..\n*          Not referenced if HOWMNY = 'A' or 'B'.\n*\n*  N       (input) INTEGER\n*          The order of the matrix T. N >= 0.\n*\n*  T       (input/output) COMPLEX array, dimension (LDT,N)\n*          The upper triangular matrix T.  T is modified, but restored\n*          on exit.\n*\n*  LDT     (input) INTEGER\n*          The leading dimension of the array T. LDT >= max(1,N).\n*\n*  VL      (input/output) COMPLEX array, dimension (LDVL,MM)\n*          On entry, if SIDE = 'L' or 'B' and HOWMNY = 'B', VL must\n*          contain an N-by-N matrix Q (usually the unitary matrix Q of\n*          Schur vectors returned by CHSEQR).\n*          On exit, if SIDE = 'L' or 'B', VL contains:\n*          if HOWMNY = 'A', the matrix Y of left eigenvectors of T;\n*          if HOWMNY = 'B', the matrix Q*Y;\n*          if HOWMNY = 'S', the left eigenvectors of T specified by\n*                           SELECT, stored consecutively in the columns\n*                           of VL, in the same order as their\n*                           eigenvalues.\n*          Not referenced if SIDE = 'R'.\n*\n*  LDVL    (input) INTEGER\n*          The leading dimension of the array VL.  LDVL >= 1, and if\n*          SIDE = 'L' or 'B', LDVL >= N.\n*\n*  VR      (input/output) COMPLEX array, dimension (LDVR,MM)\n*          On entry, if SIDE = 'R' or 'B' and HOWMNY = 'B', VR must\n*          contain an N-by-N matrix Q (usually the unitary matrix Q of\n*          Schur vectors returned by CHSEQR).\n*          On exit, if SIDE = 'R' or 'B', VR contains:\n*          if HOWMNY = 'A', the matrix X of right eigenvectors of T;\n*          if HOWMNY = 'B', the matrix Q*X;\n*          if HOWMNY = 'S', the right eigenvectors of T specified by\n*                           SELECT, stored consecutively in the columns\n*                           of VR, in the same order as their\n*                           eigenvalues.\n*          Not referenced if SIDE = 'L'.\n*\n*  LDVR    (input) INTEGER\n*          The leading dimension of the array VR.  LDVR >= 1, and if\n*          SIDE = 'R' or 'B'; LDVR >= N.\n*\n*  MM      (input) INTEGER\n*          The number of columns in the arrays VL and/or VR. MM >= M.\n*\n*  M       (output) INTEGER\n*          The number of columns in the arrays VL and/or VR actually\n*          used to store the eigenvectors.  If HOWMNY = 'A' or 'B', M\n*          is set to N.  Each selected eigenvector occupies one\n*          column.\n*\n*  WORK    (workspace) COMPLEX array, dimension (2*N)\n*\n*  RWORK   (workspace) REAL array, dimension (N)\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit\n*          < 0:  if INFO = -i, the i-th argument had an illegal value\n*\n\n*  Further Details\n*  ===============\n*\n*  The algorithm used in this program is basically backward (forward)\n*  substitution, with scaling to make the the code robust against\n*  possible overflow.\n*\n*  Each eigenvector is normalized so that the element of largest\n*  magnitude has magnitude 1; here the magnitude of a complex number\n*  (x,y) is taken to be |x| + |y|.\n*\n*  =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  m, info, t, vl, vr = NumRu::Lapack.ctrevc( side, howmny, select, t, vl, vr, [:usage => usage, :help => help])\n");
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
  if (NA_TYPE(rblapack_vl) != NA_SCOMPLEX)
    rblapack_vl = na_change_type(rblapack_vl, NA_SCOMPLEX);
  vl = NA_PTR_TYPE(rblapack_vl, complex*);
  side = StringValueCStr(rblapack_side)[0];
  if (!NA_IsNArray(rblapack_vr))
    rb_raise(rb_eArgError, "vr (6th argument) must be NArray");
  if (NA_RANK(rblapack_vr) != 2)
    rb_raise(rb_eArgError, "rank of vr (6th argument) must be %d", 2);
  if (NA_SHAPE1(rblapack_vr) != mm)
    rb_raise(rb_eRuntimeError, "shape 1 of vr must be the same as shape 1 of vl");
  ldvr = NA_SHAPE0(rblapack_vr);
  if (NA_TYPE(rblapack_vr) != NA_SCOMPLEX)
    rblapack_vr = na_change_type(rblapack_vr, NA_SCOMPLEX);
  vr = NA_PTR_TYPE(rblapack_vr, complex*);
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
  if (NA_TYPE(rblapack_t) != NA_SCOMPLEX)
    rblapack_t = na_change_type(rblapack_t, NA_SCOMPLEX);
  t = NA_PTR_TYPE(rblapack_t, complex*);
  howmny = StringValueCStr(rblapack_howmny)[0];
  {
    int shape[2];
    shape[0] = ldt;
    shape[1] = n;
    rblapack_t_out__ = na_make_object(NA_SCOMPLEX, 2, shape, cNArray);
  }
  t_out__ = NA_PTR_TYPE(rblapack_t_out__, complex*);
  MEMCPY(t_out__, t, complex, NA_TOTAL(rblapack_t));
  rblapack_t = rblapack_t_out__;
  t = t_out__;
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
  rwork = ALLOC_N(real, (n));

  ctrevc_(&side, &howmny, select, &n, t, &ldt, vl, &ldvl, vr, &ldvr, &mm, &m, work, rwork, &info);

  free(work);
  free(rwork);
  rblapack_m = INT2NUM(m);
  rblapack_info = INT2NUM(info);
  return rb_ary_new3(5, rblapack_m, rblapack_info, rblapack_t, rblapack_vl, rblapack_vr);
}

void
init_lapack_ctrevc(VALUE mLapack){
  rb_define_module_function(mLapack, "ctrevc", rblapack_ctrevc, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}

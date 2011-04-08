#include "rb_lapack.h"

extern VOID strexc_(char *compq, integer *n, real *t, integer *ldt, real *q, integer *ldq, integer *ifst, integer *ilst, real *work, integer *info);

static VALUE sHelp, sUsage;

static VALUE
rblapack_strexc(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_compq;
  char compq; 
  VALUE rblapack_t;
  real *t; 
  VALUE rblapack_q;
  real *q; 
  VALUE rblapack_ifst;
  integer ifst; 
  VALUE rblapack_ilst;
  integer ilst; 
  VALUE rblapack_info;
  integer info; 
  VALUE rblapack_t_out__;
  real *t_out__;
  VALUE rblapack_q_out__;
  real *q_out__;
  real *work;

  integer ldt;
  integer n;
  integer ldq;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  info, t, q, ifst, ilst = NumRu::Lapack.strexc( compq, t, q, ifst, ilst, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE STREXC( COMPQ, N, T, LDT, Q, LDQ, IFST, ILST, WORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  STREXC reorders the real Schur factorization of a real matrix\n*  A = Q*T*Q**T, so that the diagonal block of T with row index IFST is\n*  moved to row ILST.\n*\n*  The real Schur form T is reordered by an orthogonal similarity\n*  transformation Z**T*T*Z, and optionally the matrix Q of Schur vectors\n*  is updated by postmultiplying it with Z.\n*\n*  T must be in Schur canonical form (as returned by SHSEQR), that is,\n*  block upper triangular with 1-by-1 and 2-by-2 diagonal blocks; each\n*  2-by-2 diagonal block has its diagonal elements equal and its\n*  off-diagonal elements of opposite sign.\n*\n\n*  Arguments\n*  =========\n*\n*  COMPQ   (input) CHARACTER*1\n*          = 'V':  update the matrix Q of Schur vectors;\n*          = 'N':  do not update Q.\n*\n*  N       (input) INTEGER\n*          The order of the matrix T. N >= 0.\n*\n*  T       (input/output) REAL array, dimension (LDT,N)\n*          On entry, the upper quasi-triangular matrix T, in Schur\n*          Schur canonical form.\n*          On exit, the reordered upper quasi-triangular matrix, again\n*          in Schur canonical form.\n*\n*  LDT     (input) INTEGER\n*          The leading dimension of the array T. LDT >= max(1,N).\n*\n*  Q       (input/output) REAL array, dimension (LDQ,N)\n*          On entry, if COMPQ = 'V', the matrix Q of Schur vectors.\n*          On exit, if COMPQ = 'V', Q has been postmultiplied by the\n*          orthogonal transformation matrix Z which reorders T.\n*          If COMPQ = 'N', Q is not referenced.\n*\n*  LDQ     (input) INTEGER\n*          The leading dimension of the array Q.  LDQ >= max(1,N).\n*\n*  IFST    (input/output) INTEGER\n*  ILST    (input/output) INTEGER\n*          Specify the reordering of the diagonal blocks of T.\n*          The block with row index IFST is moved to row ILST, by a\n*          sequence of transpositions between adjacent blocks.\n*          On exit, if IFST pointed on entry to the second row of a\n*          2-by-2 block, it is changed to point to the first row; ILST\n*          always points to the first row of the block in its final\n*          position (which may differ from its input value by +1 or -1).\n*          1 <= IFST <= N; 1 <= ILST <= N.\n*\n*  WORK    (workspace) REAL array, dimension (N)\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit\n*          < 0:  if INFO = -i, the i-th argument had an illegal value\n*          = 1:  two adjacent blocks were too close to swap (the problem\n*                is very ill-conditioned); T may have been partially\n*                reordered, and ILST points to the first row of the\n*                current position of the block being moved.\n*\n\n*  =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  info, t, q, ifst, ilst = NumRu::Lapack.strexc( compq, t, q, ifst, ilst, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 5)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 5)", argc);
  rblapack_compq = argv[0];
  rblapack_t = argv[1];
  rblapack_q = argv[2];
  rblapack_ifst = argv[3];
  rblapack_ilst = argv[4];
  if (rb_options != Qnil) {
  }

  compq = StringValueCStr(rblapack_compq)[0];
  ilst = NUM2INT(rblapack_ilst);
  if (!NA_IsNArray(rblapack_q))
    rb_raise(rb_eArgError, "q (3th argument) must be NArray");
  if (NA_RANK(rblapack_q) != 2)
    rb_raise(rb_eArgError, "rank of q (3th argument) must be %d", 2);
  n = NA_SHAPE1(rblapack_q);
  ldq = NA_SHAPE0(rblapack_q);
  if (NA_TYPE(rblapack_q) != NA_SFLOAT)
    rblapack_q = na_change_type(rblapack_q, NA_SFLOAT);
  q = NA_PTR_TYPE(rblapack_q, real*);
  if (!NA_IsNArray(rblapack_t))
    rb_raise(rb_eArgError, "t (2th argument) must be NArray");
  if (NA_RANK(rblapack_t) != 2)
    rb_raise(rb_eArgError, "rank of t (2th argument) must be %d", 2);
  if (NA_SHAPE1(rblapack_t) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of t must be the same as shape 1 of q");
  ldt = NA_SHAPE0(rblapack_t);
  if (NA_TYPE(rblapack_t) != NA_SFLOAT)
    rblapack_t = na_change_type(rblapack_t, NA_SFLOAT);
  t = NA_PTR_TYPE(rblapack_t, real*);
  ifst = NUM2INT(rblapack_ifst);
  {
    int shape[2];
    shape[0] = ldt;
    shape[1] = n;
    rblapack_t_out__ = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  t_out__ = NA_PTR_TYPE(rblapack_t_out__, real*);
  MEMCPY(t_out__, t, real, NA_TOTAL(rblapack_t));
  rblapack_t = rblapack_t_out__;
  t = t_out__;
  {
    int shape[2];
    shape[0] = ldq;
    shape[1] = n;
    rblapack_q_out__ = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  q_out__ = NA_PTR_TYPE(rblapack_q_out__, real*);
  MEMCPY(q_out__, q, real, NA_TOTAL(rblapack_q));
  rblapack_q = rblapack_q_out__;
  q = q_out__;
  work = ALLOC_N(real, (n));

  strexc_(&compq, &n, t, &ldt, q, &ldq, &ifst, &ilst, work, &info);

  free(work);
  rblapack_info = INT2NUM(info);
  rblapack_ifst = INT2NUM(ifst);
  rblapack_ilst = INT2NUM(ilst);
  return rb_ary_new3(5, rblapack_info, rblapack_t, rblapack_q, rblapack_ifst, rblapack_ilst);
}

void
init_lapack_strexc(VALUE mLapack){
  rb_define_module_function(mLapack, "strexc", rblapack_strexc, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}

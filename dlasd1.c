#include "rb_lapack.h"

extern VOID dlasd1_(integer *nl, integer *nr, integer *sqre, doublereal *d, doublereal *alpha, doublereal *beta, doublereal *u, integer *ldu, doublereal *vt, integer *ldvt, integer *idxq, integer *iwork, doublereal *work, integer *info);

static VALUE sHelp, sUsage;

static VALUE
rblapack_dlasd1(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_nl;
  integer nl; 
  VALUE rblapack_nr;
  integer nr; 
  VALUE rblapack_sqre;
  integer sqre; 
  VALUE rblapack_d;
  doublereal *d; 
  VALUE rblapack_alpha;
  doublereal alpha; 
  VALUE rblapack_beta;
  doublereal beta; 
  VALUE rblapack_u;
  doublereal *u; 
  VALUE rblapack_vt;
  doublereal *vt; 
  VALUE rblapack_idxq;
  integer *idxq; 
  VALUE rblapack_info;
  integer info; 
  VALUE rblapack_d_out__;
  doublereal *d_out__;
  VALUE rblapack_u_out__;
  doublereal *u_out__;
  VALUE rblapack_vt_out__;
  doublereal *vt_out__;
  integer *iwork;
  doublereal *work;

  integer n;
  integer ldu;
  integer ldvt;
  integer m;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  idxq, info, d, alpha, beta, u, vt = NumRu::Lapack.dlasd1( nl, nr, sqre, d, alpha, beta, u, vt, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE DLASD1( NL, NR, SQRE, D, ALPHA, BETA, U, LDU, VT, LDVT, IDXQ, IWORK, WORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  DLASD1 computes the SVD of an upper bidiagonal N-by-M matrix B,\n*  where N = NL + NR + 1 and M = N + SQRE. DLASD1 is called from DLASD0.\n*\n*  A related subroutine DLASD7 handles the case in which the singular\n*  values (and the singular vectors in factored form) are desired.\n*\n*  DLASD1 computes the SVD as follows:\n*\n*                ( D1(in)  0    0     0 )\n*    B = U(in) * (   Z1'   a   Z2'    b ) * VT(in)\n*                (   0     0   D2(in) 0 )\n*\n*      = U(out) * ( D(out) 0) * VT(out)\n*\n*  where Z' = (Z1' a Z2' b) = u' VT', and u is a vector of dimension M\n*  with ALPHA and BETA in the NL+1 and NL+2 th entries and zeros\n*  elsewhere; and the entry b is empty if SQRE = 0.\n*\n*  The left singular vectors of the original matrix are stored in U, and\n*  the transpose of the right singular vectors are stored in VT, and the\n*  singular values are in D.  The algorithm consists of three stages:\n*\n*     The first stage consists of deflating the size of the problem\n*     when there are multiple singular values or when there are zeros in\n*     the Z vector.  For each such occurence the dimension of the\n*     secular equation problem is reduced by one.  This stage is\n*     performed by the routine DLASD2.\n*\n*     The second stage consists of calculating the updated\n*     singular values. This is done by finding the square roots of the\n*     roots of the secular equation via the routine DLASD4 (as called\n*     by DLASD3). This routine also calculates the singular vectors of\n*     the current problem.\n*\n*     The final stage consists of computing the updated singular vectors\n*     directly using the updated singular values.  The singular vectors\n*     for the current problem are multiplied with the singular vectors\n*     from the overall problem.\n*\n\n*  Arguments\n*  =========\n*\n*  NL     (input) INTEGER\n*         The row dimension of the upper block.  NL >= 1.\n*\n*  NR     (input) INTEGER\n*         The row dimension of the lower block.  NR >= 1.\n*\n*  SQRE   (input) INTEGER\n*         = 0: the lower block is an NR-by-NR square matrix.\n*         = 1: the lower block is an NR-by-(NR+1) rectangular matrix.\n*\n*         The bidiagonal matrix has row dimension N = NL + NR + 1,\n*         and column dimension M = N + SQRE.\n*\n*  D      (input/output) DOUBLE PRECISION array,\n*                        dimension (N = NL+NR+1).\n*         On entry D(1:NL,1:NL) contains the singular values of the\n*         upper block; and D(NL+2:N) contains the singular values of\n*         the lower block. On exit D(1:N) contains the singular values\n*         of the modified matrix.\n*\n*  ALPHA  (input/output) DOUBLE PRECISION\n*         Contains the diagonal element associated with the added row.\n*\n*  BETA   (input/output) DOUBLE PRECISION\n*         Contains the off-diagonal element associated with the added\n*         row.\n*\n*  U      (input/output) DOUBLE PRECISION array, dimension(LDU,N)\n*         On entry U(1:NL, 1:NL) contains the left singular vectors of\n*         the upper block; U(NL+2:N, NL+2:N) contains the left singular\n*         vectors of the lower block. On exit U contains the left\n*         singular vectors of the bidiagonal matrix.\n*\n*  LDU    (input) INTEGER\n*         The leading dimension of the array U.  LDU >= max( 1, N ).\n*\n*  VT     (input/output) DOUBLE PRECISION array, dimension(LDVT,M)\n*         where M = N + SQRE.\n*         On entry VT(1:NL+1, 1:NL+1)' contains the right singular\n*         vectors of the upper block; VT(NL+2:M, NL+2:M)' contains\n*         the right singular vectors of the lower block. On exit\n*         VT' contains the right singular vectors of the\n*         bidiagonal matrix.\n*\n*  LDVT   (input) INTEGER\n*         The leading dimension of the array VT.  LDVT >= max( 1, M ).\n*\n*  IDXQ  (output) INTEGER array, dimension(N)\n*         This contains the permutation which will reintegrate the\n*         subproblem just solved back into sorted order, i.e.\n*         D( IDXQ( I = 1, N ) ) will be in ascending order.\n*\n*  IWORK  (workspace) INTEGER array, dimension( 4 * N )\n*\n*  WORK   (workspace) DOUBLE PRECISION array, dimension( 3*M**2 + 2*M )\n*\n*  INFO   (output) INTEGER\n*          = 0:  successful exit.\n*          < 0:  if INFO = -i, the i-th argument had an illegal value.\n*          > 0:  if INFO = 1, a singular value did not converge\n*\n\n*  Further Details\n*  ===============\n*\n*  Based on contributions by\n*     Ming Gu and Huan Ren, Computer Science Division, University of\n*     California at Berkeley, USA\n*\n*  =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  idxq, info, d, alpha, beta, u, vt = NumRu::Lapack.dlasd1( nl, nr, sqre, d, alpha, beta, u, vt, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 8)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 8)", argc);
  rblapack_nl = argv[0];
  rblapack_nr = argv[1];
  rblapack_sqre = argv[2];
  rblapack_d = argv[3];
  rblapack_alpha = argv[4];
  rblapack_beta = argv[5];
  rblapack_u = argv[6];
  rblapack_vt = argv[7];
  if (rb_options != Qnil) {
  }

  alpha = NUM2DBL(rblapack_alpha);
  if (!NA_IsNArray(rblapack_u))
    rb_raise(rb_eArgError, "u (7th argument) must be NArray");
  if (NA_RANK(rblapack_u) != 2)
    rb_raise(rb_eArgError, "rank of u (7th argument) must be %d", 2);
  n = NA_SHAPE1(rblapack_u);
  if (n != (nl+nr+1))
    rb_raise(rb_eRuntimeError, "shape 1 of u must be %d", nl+nr+1);
  ldu = NA_SHAPE0(rblapack_u);
  if (NA_TYPE(rblapack_u) != NA_DFLOAT)
    rblapack_u = na_change_type(rblapack_u, NA_DFLOAT);
  u = NA_PTR_TYPE(rblapack_u, doublereal*);
  if (!NA_IsNArray(rblapack_d))
    rb_raise(rb_eArgError, "d (4th argument) must be NArray");
  if (NA_RANK(rblapack_d) != 1)
    rb_raise(rb_eArgError, "rank of d (4th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_d) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of d must be the same as shape 1 of u");
  if (NA_TYPE(rblapack_d) != NA_DFLOAT)
    rblapack_d = na_change_type(rblapack_d, NA_DFLOAT);
  d = NA_PTR_TYPE(rblapack_d, doublereal*);
  nr = NUM2INT(rblapack_nr);
  beta = NUM2DBL(rblapack_beta);
  nl = NUM2INT(rblapack_nl);
  if (!NA_IsNArray(rblapack_vt))
    rb_raise(rb_eArgError, "vt (8th argument) must be NArray");
  if (NA_RANK(rblapack_vt) != 2)
    rb_raise(rb_eArgError, "rank of vt (8th argument) must be %d", 2);
  m = NA_SHAPE1(rblapack_vt);
  if (m != (n + sqre))
    rb_raise(rb_eRuntimeError, "shape 1 of vt must be %d", n + sqre);
  ldvt = NA_SHAPE0(rblapack_vt);
  if (NA_TYPE(rblapack_vt) != NA_DFLOAT)
    rblapack_vt = na_change_type(rblapack_vt, NA_DFLOAT);
  vt = NA_PTR_TYPE(rblapack_vt, doublereal*);
  sqre = NUM2INT(rblapack_sqre);
  n = nl+nr+1;
  m = n + sqre;
  {
    int shape[1];
    shape[0] = n;
    rblapack_idxq = na_make_object(NA_LINT, 1, shape, cNArray);
  }
  idxq = NA_PTR_TYPE(rblapack_idxq, integer*);
  {
    int shape[1];
    shape[0] = n;
    rblapack_d_out__ = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  d_out__ = NA_PTR_TYPE(rblapack_d_out__, doublereal*);
  MEMCPY(d_out__, d, doublereal, NA_TOTAL(rblapack_d));
  rblapack_d = rblapack_d_out__;
  d = d_out__;
  {
    int shape[2];
    shape[0] = ldu;
    shape[1] = n;
    rblapack_u_out__ = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  u_out__ = NA_PTR_TYPE(rblapack_u_out__, doublereal*);
  MEMCPY(u_out__, u, doublereal, NA_TOTAL(rblapack_u));
  rblapack_u = rblapack_u_out__;
  u = u_out__;
  {
    int shape[2];
    shape[0] = ldvt;
    shape[1] = m;
    rblapack_vt_out__ = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  vt_out__ = NA_PTR_TYPE(rblapack_vt_out__, doublereal*);
  MEMCPY(vt_out__, vt, doublereal, NA_TOTAL(rblapack_vt));
  rblapack_vt = rblapack_vt_out__;
  vt = vt_out__;
  iwork = ALLOC_N(integer, (4 * n));
  work = ALLOC_N(doublereal, (3*pow(m,2) + 2*m));

  dlasd1_(&nl, &nr, &sqre, d, &alpha, &beta, u, &ldu, vt, &ldvt, idxq, iwork, work, &info);

  free(iwork);
  free(work);
  rblapack_info = INT2NUM(info);
  rblapack_alpha = rb_float_new((double)alpha);
  rblapack_beta = rb_float_new((double)beta);
  return rb_ary_new3(7, rblapack_idxq, rblapack_info, rblapack_d, rblapack_alpha, rblapack_beta, rblapack_u, rblapack_vt);
}

void
init_lapack_dlasd1(VALUE mLapack){
  rb_define_module_function(mLapack, "dlasd1", rblapack_dlasd1, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}

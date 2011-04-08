#include "rb_lapack.h"

extern VOID dgesvj_(char *joba, char *jobu, char *jobv, integer *m, integer *n, doublereal *a, integer *lda, doublereal *sva, integer *mv, doublereal *v, integer *ldv, doublereal *work, integer *lwork, integer *info);

static VALUE sHelp, sUsage;

static VALUE
rblapack_dgesvj(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_joba;
  char joba; 
  VALUE rblapack_jobu;
  char jobu; 
  VALUE rblapack_jobv;
  char jobv; 
  VALUE rblapack_m;
  integer m; 
  VALUE rblapack_a;
  doublereal *a; 
  VALUE rblapack_mv;
  integer mv; 
  VALUE rblapack_v;
  doublereal *v; 
  VALUE rblapack_work;
  doublereal *work; 
  VALUE rblapack_sva;
  doublereal *sva; 
  VALUE rblapack_info;
  integer info; 
  VALUE rblapack_a_out__;
  doublereal *a_out__;
  VALUE rblapack_v_out__;
  doublereal *v_out__;
  VALUE rblapack_work_out__;
  doublereal *work_out__;

  integer lda;
  integer n;
  integer ldv;
  integer lwork;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  sva, info, a, v, work = NumRu::Lapack.dgesvj( joba, jobu, jobv, m, a, mv, v, work, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE DGESVJ( JOBA, JOBU, JOBV, M, N, A, LDA, SVA, MV, V, LDV, WORK, LWORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  DGESVJ computes the singular value decomposition (SVD) of a real\n*  M-by-N matrix A, where M >= N. The SVD of A is written as\n*                                     [++]   [xx]   [x0]   [xx]\n*               A = U * SIGMA * V^t,  [++] = [xx] * [ox] * [xx]\n*                                     [++]   [xx]\n*  where SIGMA is an N-by-N diagonal matrix, U is an M-by-N orthonormal\n*  matrix, and V is an N-by-N orthogonal matrix. The diagonal elements\n*  of SIGMA are the singular values of A. The columns of U and V are the\n*  left and the right singular vectors of A, respectively.\n*\n*  Further Details\n*  ~~~~~~~~~~~~~~~\n*  The orthogonal N-by-N matrix V is obtained as a product of Jacobi plane\n*  rotations. The rotations are implemented as fast scaled rotations of\n*  Anda and Park [1]. In the case of underflow of the Jacobi angle, a\n*  modified Jacobi transformation of Drmac [4] is used. Pivot strategy uses\n*  column interchanges of de Rijk [2]. The relative accuracy of the computed\n*  singular values and the accuracy of the computed singular vectors (in\n*  angle metric) is as guaranteed by the theory of Demmel and Veselic [3].\n*  The condition number that determines the accuracy in the full rank case\n*  is essentially min_{D=diag} kappa(A*D), where kappa(.) is the\n*  spectral condition number. The best performance of this Jacobi SVD\n*  procedure is achieved if used in an  accelerated version of Drmac and\n*  Veselic [5,6], and it is the kernel routine in the SIGMA library [7].\n*  Some tunning parameters (marked with [TP]) are available for the\n*  implementer.\n*  The computational range for the nonzero singular values is the  machine\n*  number interval ( UNDERFLOW , OVERFLOW ). In extreme cases, even\n*  denormalized singular values can be computed with the corresponding\n*  gradual loss of accurate digits.\n*\n*  Contributors\n*  ~~~~~~~~~~~~\n*  Zlatko Drmac (Zagreb, Croatia) and Kresimir Veselic (Hagen, Germany)\n*\n*  References\n*  ~~~~~~~~~~\n* [1] A. A. Anda and H. Park: Fast plane rotations with dynamic scaling.\n*     SIAM J. matrix Anal. Appl., Vol. 15 (1994), pp. 162-174.\n* [2] P. P. M. De Rijk: A one-sided Jacobi algorithm for computing the\n*     singular value decomposition on a vector computer.\n*     SIAM J. Sci. Stat. Comp., Vol. 10 (1998), pp. 359-371.\n* [3] J. Demmel and K. Veselic: Jacobi method is more accurate than QR.\n* [4] Z. Drmac: Implementation of Jacobi rotations for accurate singular\n*     value computation in floating point arithmetic.\n*     SIAM J. Sci. Comp., Vol. 18 (1997), pp. 1200-1222.\n* [5] Z. Drmac and K. Veselic: New fast and accurate Jacobi SVD algorithm I.\n*     SIAM J. Matrix Anal. Appl. Vol. 35, No. 2 (2008), pp. 1322-1342.\n*     LAPACK Working note 169.\n* [6] Z. Drmac and K. Veselic: New fast and accurate Jacobi SVD algorithm II.\n*     SIAM J. Matrix Anal. Appl. Vol. 35, No. 2 (2008), pp. 1343-1362.\n*     LAPACK Working note 170.\n* [7] Z. Drmac: SIGMA - mathematical software library for accurate SVD, PSV,\n*     QSVD, (H,K)-SVD computations.\n*     Department of Mathematics, University of Zagreb, 2008.\n*\n*  Bugs, Examples and Comments\n*  ~~~~~~~~~~~~~~~~~~~~~~~~~~~\n*  Please report all bugs and send interesting test examples and comments to\n*  drmac@math.hr. Thank you.\n*\n\n*  Arguments\n*  =========\n*\n*  JOBA    (input) CHARACTER* 1\n*          Specifies the structure of A.\n*          = 'L': The input matrix A is lower triangular;\n*          = 'U': The input matrix A is upper triangular;\n*          = 'G': The input matrix A is general M-by-N matrix, M >= N.\n*\n*  JOBU    (input) CHARACTER*1\n*          Specifies whether to compute the left singular vectors\n*          (columns of U):\n*          = 'U': The left singular vectors corresponding to the nonzero\n*                 singular values are computed and returned in the leading\n*                 columns of A. See more details in the description of A.\n*                 The default numerical orthogonality threshold is set to\n*                 approximately TOL=CTOL*EPS, CTOL=DSQRT(M), EPS=DLAMCH('E').\n*          = 'C': Analogous to JOBU='U', except that user can control the\n*                 level of numerical orthogonality of the computed left\n*                 singular vectors. TOL can be set to TOL = CTOL*EPS, where\n*                 CTOL is given on input in the array WORK.\n*                 No CTOL smaller than ONE is allowed. CTOL greater\n*                 than 1 / EPS is meaningless. The option 'C'\n*                 can be used if M*EPS is satisfactory orthogonality\n*                 of the computed left singular vectors, so CTOL=M could\n*                 save few sweeps of Jacobi rotations.\n*                 See the descriptions of A and WORK(1).\n*          = 'N': The matrix U is not computed. However, see the\n*                 description of A.\n*\n*  JOBV    (input) CHARACTER*1\n*          Specifies whether to compute the right singular vectors, that\n*          is, the matrix V:\n*          = 'V' : the matrix V is computed and returned in the array V\n*          = 'A' : the Jacobi rotations are applied to the MV-by-N\n*                  array V. In other words, the right singular vector\n*                  matrix V is not computed explicitly, instead it is\n*                  applied to an MV-by-N matrix initially stored in the\n*                  first MV rows of V.\n*          = 'N' : the matrix V is not computed and the array V is not\n*                  referenced\n*\n*  M       (input) INTEGER\n*          The number of rows of the input matrix A.  M >= 0.\n*\n*  N       (input) INTEGER\n*          The number of columns of the input matrix A.\n*          M >= N >= 0.\n*\n*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)\n*          On entry, the M-by-N matrix A.\n*          On exit :\n*          If JOBU .EQ. 'U' .OR. JOBU .EQ. 'C' :\n*                 If INFO .EQ. 0 :\n*                 RANKA orthonormal columns of U are returned in the\n*                 leading RANKA columns of the array A. Here RANKA <= N\n*                 is the number of computed singular values of A that are\n*                 above the underflow threshold DLAMCH('S'). The singular\n*                 vectors corresponding to underflowed or zero singular\n*                 values are not computed. The value of RANKA is returned\n*                 in the array WORK as RANKA=NINT(WORK(2)). Also see the\n*                 descriptions of SVA and WORK. The computed columns of U\n*                 are mutually numerically orthogonal up to approximately\n*                 TOL=DSQRT(M)*EPS (default); or TOL=CTOL*EPS (JOBU.EQ.'C'),\n*                 see the description of JOBU.\n*                 If INFO .GT. 0 :\n*                 the procedure DGESVJ did not converge in the given number\n*                 of iterations (sweeps). In that case, the computed\n*                 columns of U may not be orthogonal up to TOL. The output\n*                 U (stored in A), SIGMA (given by the computed singular\n*                 values in SVA(1:N)) and V is still a decomposition of the\n*                 input matrix A in the sense that the residual\n*                 ||A-SCALE*U*SIGMA*V^T||_2 / ||A||_2 is small.\n*\n*          If JOBU .EQ. 'N' :\n*                 If INFO .EQ. 0 :\n*                 Note that the left singular vectors are 'for free' in the\n*                 one-sided Jacobi SVD algorithm. However, if only the\n*                 singular values are needed, the level of numerical\n*                 orthogonality of U is not an issue and iterations are\n*                 stopped when the columns of the iterated matrix are\n*                 numerically orthogonal up to approximately M*EPS. Thus,\n*                 on exit, A contains the columns of U scaled with the\n*                 corresponding singular values.\n*                 If INFO .GT. 0 :\n*                 the procedure DGESVJ did not converge in the given number\n*                 of iterations (sweeps).\n*\n*  LDA     (input) INTEGER\n*          The leading dimension of the array A.  LDA >= max(1,M).\n*\n*  SVA     (workspace/output) DOUBLE PRECISION array, dimension (N)\n*          On exit :\n*          If INFO .EQ. 0 :\n*          depending on the value SCALE = WORK(1), we have:\n*                 If SCALE .EQ. ONE :\n*                 SVA(1:N) contains the computed singular values of A.\n*                 During the computation SVA contains the Euclidean column\n*                 norms of the iterated matrices in the array A.\n*                 If SCALE .NE. ONE :\n*                 The singular values of A are SCALE*SVA(1:N), and this\n*                 factored representation is due to the fact that some of the\n*                 singular values of A might underflow or overflow.\n*          If INFO .GT. 0 :\n*          the procedure DGESVJ did not converge in the given number of\n*          iterations (sweeps) and SCALE*SVA(1:N) may not be accurate.\n*\n*  MV      (input) INTEGER\n*          If JOBV .EQ. 'A', then the product of Jacobi rotations in DGESVJ\n*          is applied to the first MV rows of V. See the description of JOBV.\n*\n*  V       (input/output) DOUBLE PRECISION array, dimension (LDV,N)\n*          If JOBV = 'V', then V contains on exit the N-by-N matrix of\n*                         the right singular vectors;\n*          If JOBV = 'A', then V contains the product of the computed right\n*                         singular vector matrix and the initial matrix in\n*                         the array V.\n*          If JOBV = 'N', then V is not referenced.\n*\n*  LDV     (input) INTEGER\n*          The leading dimension of the array V, LDV .GE. 1.\n*          If JOBV .EQ. 'V', then LDV .GE. max(1,N).\n*          If JOBV .EQ. 'A', then LDV .GE. max(1,MV) .\n*\n*  WORK    (input/workspace/output) DOUBLE PRECISION array, dimension max(4,M+N).\n*          On entry :\n*          If JOBU .EQ. 'C' :\n*          WORK(1) = CTOL, where CTOL defines the threshold for convergence.\n*                    The process stops if all columns of A are mutually\n*                    orthogonal up to CTOL*EPS, EPS=DLAMCH('E').\n*                    It is required that CTOL >= ONE, i.e. it is not\n*                    allowed to force the routine to obtain orthogonality\n*                    below EPS.\n*          On exit :\n*          WORK(1) = SCALE is the scaling factor such that SCALE*SVA(1:N)\n*                    are the computed singular values of A.\n*                    (See description of SVA().)\n*          WORK(2) = NINT(WORK(2)) is the number of the computed nonzero\n*                    singular values.\n*          WORK(3) = NINT(WORK(3)) is the number of the computed singular\n*                    values that are larger than the underflow threshold.\n*          WORK(4) = NINT(WORK(4)) is the number of sweeps of Jacobi\n*                    rotations needed for numerical convergence.\n*          WORK(5) = max_{i.NE.j} |COS(A(:,i),A(:,j))| in the last sweep.\n*                    This is useful information in cases when DGESVJ did\n*                    not converge, as it can be used to estimate whether\n*                    the output is stil useful and for post festum analysis.\n*          WORK(6) = the largest absolute value over all sines of the\n*                    Jacobi rotation angles in the last sweep. It can be\n*                    useful for a post festum analysis.\n*\n*  LWORK   (input) INTEGER\n*          length of WORK, WORK >= MAX(6,M+N)\n*\n*  INFO    (output) INTEGER\n*          = 0 : successful exit.\n*          < 0 : if INFO = -i, then the i-th argument had an illegal value\n*          > 0 : DGESVJ did not converge in the maximal allowed number (30)\n*                of sweeps. The output may still be useful. See the\n*                description of WORK.\n*\n\n*  =====================================================================\n*\n*     .. Local Parameters ..\n      DOUBLE PRECISION   ZERO, HALF, ONE, TWO\n      PARAMETER          ( ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0,\n     +                   TWO = 2.0D0 )\n      INTEGER            NSWEEP\n      PARAMETER          ( NSWEEP = 30 )\n*     ..\n*     .. Local Scalars ..\n      DOUBLE PRECISION   AAPP, AAPP0, AAPQ, AAQQ, APOAQ, AQOAP, BIG,\n     +                   BIGTHETA, CS, CTOL, EPSLN, LARGE, MXAAPQ,\n     +                   MXSINJ, ROOTBIG, ROOTEPS, ROOTSFMIN, ROOTTOL,\n     +                   SKL, SFMIN, SMALL, SN, T, TEMP1, THETA,\n     +                   THSIGN, TOL\n      INTEGER            BLSKIP, EMPTSW, i, ibr, IERR, igl, IJBLSK, ir1,\n     +                   ISWROT, jbc, jgl, KBL, LKAHEAD, MVL, N2, N34,\n     +                   N4, NBL, NOTROT, p, PSKIPPED, q, ROWSKIP,\n     +                   SWBAND\n      LOGICAL            APPLV, GOSCALE, LOWER, LSVEC, NOSCALE, ROTOK,\n     +                   RSVEC, UCTOL, UPPER\n*     ..\n*     .. Local Arrays ..\n      DOUBLE PRECISION   FASTR( 5 )\n*     ..\n*     .. Intrinsic Functions ..\n      INTRINSIC          DABS, DMAX1, DMIN1, DBLE, MIN0, DSIGN, DSQRT\n*     ..\n*     .. External Functions ..\n*     ..\n*     from BLAS\n      DOUBLE PRECISION   DDOT, DNRM2\n      EXTERNAL           DDOT, DNRM2\n      INTEGER            IDAMAX\n      EXTERNAL           IDAMAX\n*     from LAPACK\n      DOUBLE PRECISION   DLAMCH\n      EXTERNAL           DLAMCH\n      LOGICAL            LSAME\n      EXTERNAL           LSAME\n*     ..\n*     .. External Subroutines ..\n*     ..\n*     from BLAS\n      EXTERNAL           DAXPY, DCOPY, DROTM, DSCAL, DSWAP\n*     from LAPACK\n      EXTERNAL           DLASCL, DLASET, DLASSQ, XERBLA\n*\n      EXTERNAL           DGSVJ0, DGSVJ1\n*     ..\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  sva, info, a, v, work = NumRu::Lapack.dgesvj( joba, jobu, jobv, m, a, mv, v, work, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 8)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 8)", argc);
  rblapack_joba = argv[0];
  rblapack_jobu = argv[1];
  rblapack_jobv = argv[2];
  rblapack_m = argv[3];
  rblapack_a = argv[4];
  rblapack_mv = argv[5];
  rblapack_v = argv[6];
  rblapack_work = argv[7];
  if (rb_options != Qnil) {
  }

  if (!NA_IsNArray(rblapack_v))
    rb_raise(rb_eArgError, "v (7th argument) must be NArray");
  if (NA_RANK(rblapack_v) != 2)
    rb_raise(rb_eArgError, "rank of v (7th argument) must be %d", 2);
  n = NA_SHAPE1(rblapack_v);
  ldv = NA_SHAPE0(rblapack_v);
  if (NA_TYPE(rblapack_v) != NA_DFLOAT)
    rblapack_v = na_change_type(rblapack_v, NA_DFLOAT);
  v = NA_PTR_TYPE(rblapack_v, doublereal*);
  if (!NA_IsNArray(rblapack_a))
    rb_raise(rb_eArgError, "a (5th argument) must be NArray");
  if (NA_RANK(rblapack_a) != 2)
    rb_raise(rb_eArgError, "rank of a (5th argument) must be %d", 2);
  if (NA_SHAPE1(rblapack_a) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of a must be the same as shape 1 of v");
  lda = NA_SHAPE0(rblapack_a);
  if (NA_TYPE(rblapack_a) != NA_DFLOAT)
    rblapack_a = na_change_type(rblapack_a, NA_DFLOAT);
  a = NA_PTR_TYPE(rblapack_a, doublereal*);
  m = NUM2INT(rblapack_m);
  jobu = StringValueCStr(rblapack_jobu)[0];
  mv = NUM2INT(rblapack_mv);
  jobv = StringValueCStr(rblapack_jobv)[0];
  joba = StringValueCStr(rblapack_joba)[0];
  if (!NA_IsNArray(rblapack_work))
    rb_raise(rb_eArgError, "work (8th argument) must be NArray");
  if (NA_RANK(rblapack_work) != 1)
    rb_raise(rb_eArgError, "rank of work (8th argument) must be %d", 1);
  lwork = NA_SHAPE0(rblapack_work);
  if (lwork != (MAX(4,m+n)))
    rb_raise(rb_eRuntimeError, "shape 0 of work must be %d", MAX(4,m+n));
  if (NA_TYPE(rblapack_work) != NA_DFLOAT)
    rblapack_work = na_change_type(rblapack_work, NA_DFLOAT);
  work = NA_PTR_TYPE(rblapack_work, doublereal*);
  lwork = MAX(4,m+n);
  {
    int shape[1];
    shape[0] = n;
    rblapack_sva = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  sva = NA_PTR_TYPE(rblapack_sva, doublereal*);
  {
    int shape[2];
    shape[0] = lda;
    shape[1] = n;
    rblapack_a_out__ = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  a_out__ = NA_PTR_TYPE(rblapack_a_out__, doublereal*);
  MEMCPY(a_out__, a, doublereal, NA_TOTAL(rblapack_a));
  rblapack_a = rblapack_a_out__;
  a = a_out__;
  {
    int shape[2];
    shape[0] = ldv;
    shape[1] = n;
    rblapack_v_out__ = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  v_out__ = NA_PTR_TYPE(rblapack_v_out__, doublereal*);
  MEMCPY(v_out__, v, doublereal, NA_TOTAL(rblapack_v));
  rblapack_v = rblapack_v_out__;
  v = v_out__;
  {
    int shape[1];
    shape[0] = lwork;
    rblapack_work_out__ = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  work_out__ = NA_PTR_TYPE(rblapack_work_out__, doublereal*);
  MEMCPY(work_out__, work, doublereal, NA_TOTAL(rblapack_work));
  rblapack_work = rblapack_work_out__;
  work = work_out__;

  dgesvj_(&joba, &jobu, &jobv, &m, &n, a, &lda, sva, &mv, v, &ldv, work, &lwork, &info);

  rblapack_info = INT2NUM(info);
  return rb_ary_new3(5, rblapack_sva, rblapack_info, rblapack_a, rblapack_v, rblapack_work);
}

void
init_lapack_dgesvj(VALUE mLapack){
  rb_define_module_function(mLapack, "dgesvj", rblapack_dgesvj, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}

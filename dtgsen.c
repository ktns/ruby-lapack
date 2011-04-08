#include "rb_lapack.h"

extern VOID dtgsen_(integer *ijob, logical *wantq, logical *wantz, logical *select, integer *n, doublereal *a, integer *lda, doublereal *b, integer *ldb, doublereal *alphar, doublereal *alphai, doublereal *beta, doublereal *q, integer *ldq, doublereal *z, integer *ldz, integer *m, doublereal *pl, doublereal *pr, doublereal *dif, doublereal *work, integer *lwork, integer *iwork, integer *liwork, integer *info);

static VALUE sHelp, sUsage;

static VALUE
rblapack_dtgsen(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_ijob;
  integer ijob; 
  VALUE rblapack_wantq;
  logical wantq; 
  VALUE rblapack_wantz;
  logical wantz; 
  VALUE rblapack_select;
  logical *select; 
  VALUE rblapack_a;
  doublereal *a; 
  VALUE rblapack_b;
  doublereal *b; 
  VALUE rblapack_q;
  doublereal *q; 
  VALUE rblapack_z;
  doublereal *z; 
  VALUE rblapack_lwork;
  integer lwork; 
  VALUE rblapack_liwork;
  integer liwork; 
  VALUE rblapack_alphar;
  doublereal *alphar; 
  VALUE rblapack_alphai;
  doublereal *alphai; 
  VALUE rblapack_beta;
  doublereal *beta; 
  VALUE rblapack_m;
  integer m; 
  VALUE rblapack_pl;
  doublereal pl; 
  VALUE rblapack_pr;
  doublereal pr; 
  VALUE rblapack_dif;
  doublereal *dif; 
  VALUE rblapack_work;
  doublereal *work; 
  VALUE rblapack_iwork;
  integer *iwork; 
  VALUE rblapack_info;
  integer info; 
  VALUE rblapack_a_out__;
  doublereal *a_out__;
  VALUE rblapack_b_out__;
  doublereal *b_out__;
  VALUE rblapack_q_out__;
  doublereal *q_out__;
  VALUE rblapack_z_out__;
  doublereal *z_out__;

  integer n;
  integer lda;
  integer ldb;
  integer ldq;
  integer ldz;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  alphar, alphai, beta, m, pl, pr, dif, work, iwork, info, a, b, q, z = NumRu::Lapack.dtgsen( ijob, wantq, wantz, select, a, b, q, z, lwork, liwork, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE DTGSEN( IJOB, WANTQ, WANTZ, SELECT, N, A, LDA, B, LDB, ALPHAR, ALPHAI, BETA, Q, LDQ, Z, LDZ, M, PL, PR, DIF, WORK, LWORK, IWORK, LIWORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  DTGSEN reorders the generalized real Schur decomposition of a real\n*  matrix pair (A, B) (in terms of an orthonormal equivalence trans-\n*  formation Q' * (A, B) * Z), so that a selected cluster of eigenvalues\n*  appears in the leading diagonal blocks of the upper quasi-triangular\n*  matrix A and the upper triangular B. The leading columns of Q and\n*  Z form orthonormal bases of the corresponding left and right eigen-\n*  spaces (deflating subspaces). (A, B) must be in generalized real\n*  Schur canonical form (as returned by DGGES), i.e. A is block upper\n*  triangular with 1-by-1 and 2-by-2 diagonal blocks. B is upper\n*  triangular.\n*\n*  DTGSEN also computes the generalized eigenvalues\n*\n*              w(j) = (ALPHAR(j) + i*ALPHAI(j))/BETA(j)\n*\n*  of the reordered matrix pair (A, B).\n*\n*  Optionally, DTGSEN computes the estimates of reciprocal condition\n*  numbers for eigenvalues and eigenspaces. These are Difu[(A11,B11),\n*  (A22,B22)] and Difl[(A11,B11), (A22,B22)], i.e. the separation(s)\n*  between the matrix pairs (A11, B11) and (A22,B22) that correspond to\n*  the selected cluster and the eigenvalues outside the cluster, resp.,\n*  and norms of \"projections\" onto left and right eigenspaces w.r.t.\n*  the selected cluster in the (1,1)-block.\n*\n\n*  Arguments\n*  =========\n*\n*  IJOB    (input) INTEGER\n*          Specifies whether condition numbers are required for the\n*          cluster of eigenvalues (PL and PR) or the deflating subspaces\n*          (Difu and Difl):\n*           =0: Only reorder w.r.t. SELECT. No extras.\n*           =1: Reciprocal of norms of \"projections\" onto left and right\n*               eigenspaces w.r.t. the selected cluster (PL and PR).\n*           =2: Upper bounds on Difu and Difl. F-norm-based estimate\n*               (DIF(1:2)).\n*           =3: Estimate of Difu and Difl. 1-norm-based estimate\n*               (DIF(1:2)).\n*               About 5 times as expensive as IJOB = 2.\n*           =4: Compute PL, PR and DIF (i.e. 0, 1 and 2 above): Economic\n*               version to get it all.\n*           =5: Compute PL, PR and DIF (i.e. 0, 1 and 3 above)\n*\n*  WANTQ   (input) LOGICAL\n*          .TRUE. : update the left transformation matrix Q;\n*          .FALSE.: do not update Q.\n*\n*  WANTZ   (input) LOGICAL\n*          .TRUE. : update the right transformation matrix Z;\n*          .FALSE.: do not update Z.\n*\n*  SELECT  (input) LOGICAL array, dimension (N)\n*          SELECT specifies the eigenvalues in the selected cluster.\n*          To select a real eigenvalue w(j), SELECT(j) must be set to\n*          .TRUE.. To select a complex conjugate pair of eigenvalues\n*          w(j) and w(j+1), corresponding to a 2-by-2 diagonal block,\n*          either SELECT(j) or SELECT(j+1) or both must be set to\n*          .TRUE.; a complex conjugate pair of eigenvalues must be\n*          either both included in the cluster or both excluded.\n*\n*  N       (input) INTEGER\n*          The order of the matrices A and B. N >= 0.\n*\n*  A       (input/output) DOUBLE PRECISION array, dimension(LDA,N)\n*          On entry, the upper quasi-triangular matrix A, with (A, B) in\n*          generalized real Schur canonical form.\n*          On exit, A is overwritten by the reordered matrix A.\n*\n*  LDA     (input) INTEGER\n*          The leading dimension of the array A. LDA >= max(1,N).\n*\n*  B       (input/output) DOUBLE PRECISION array, dimension(LDB,N)\n*          On entry, the upper triangular matrix B, with (A, B) in\n*          generalized real Schur canonical form.\n*          On exit, B is overwritten by the reordered matrix B.\n*\n*  LDB     (input) INTEGER\n*          The leading dimension of the array B. LDB >= max(1,N).\n*\n*  ALPHAR  (output) DOUBLE PRECISION array, dimension (N)\n*  ALPHAI  (output) DOUBLE PRECISION array, dimension (N)\n*  BETA    (output) DOUBLE PRECISION array, dimension (N)\n*          On exit, (ALPHAR(j) + ALPHAI(j)*i)/BETA(j), j=1,...,N, will\n*          be the generalized eigenvalues.  ALPHAR(j) + ALPHAI(j)*i\n*          and BETA(j),j=1,...,N  are the diagonals of the complex Schur\n*          form (S,T) that would result if the 2-by-2 diagonal blocks of\n*          the real generalized Schur form of (A,B) were further reduced\n*          to triangular form using complex unitary transformations.\n*          If ALPHAI(j) is zero, then the j-th eigenvalue is real; if\n*          positive, then the j-th and (j+1)-st eigenvalues are a\n*          complex conjugate pair, with ALPHAI(j+1) negative.\n*\n*  Q       (input/output) DOUBLE PRECISION array, dimension (LDQ,N)\n*          On entry, if WANTQ = .TRUE., Q is an N-by-N matrix.\n*          On exit, Q has been postmultiplied by the left orthogonal\n*          transformation matrix which reorder (A, B); The leading M\n*          columns of Q form orthonormal bases for the specified pair of\n*          left eigenspaces (deflating subspaces).\n*          If WANTQ = .FALSE., Q is not referenced.\n*\n*  LDQ     (input) INTEGER\n*          The leading dimension of the array Q.  LDQ >= 1;\n*          and if WANTQ = .TRUE., LDQ >= N.\n*\n*  Z       (input/output) DOUBLE PRECISION array, dimension (LDZ,N)\n*          On entry, if WANTZ = .TRUE., Z is an N-by-N matrix.\n*          On exit, Z has been postmultiplied by the left orthogonal\n*          transformation matrix which reorder (A, B); The leading M\n*          columns of Z form orthonormal bases for the specified pair of\n*          left eigenspaces (deflating subspaces).\n*          If WANTZ = .FALSE., Z is not referenced.\n*\n*  LDZ     (input) INTEGER\n*          The leading dimension of the array Z. LDZ >= 1;\n*          If WANTZ = .TRUE., LDZ >= N.\n*\n*  M       (output) INTEGER\n*          The dimension of the specified pair of left and right eigen-\n*          spaces (deflating subspaces). 0 <= M <= N.\n*\n*  PL      (output) DOUBLE PRECISION\n*  PR      (output) DOUBLE PRECISION\n*          If IJOB = 1, 4 or 5, PL, PR are lower bounds on the\n*          reciprocal of the norm of \"projections\" onto left and right\n*          eigenspaces with respect to the selected cluster.\n*          0 < PL, PR <= 1.\n*          If M = 0 or M = N, PL = PR  = 1.\n*          If IJOB = 0, 2 or 3, PL and PR are not referenced.\n*\n*  DIF     (output) DOUBLE PRECISION array, dimension (2).\n*          If IJOB >= 2, DIF(1:2) store the estimates of Difu and Difl.\n*          If IJOB = 2 or 4, DIF(1:2) are F-norm-based upper bounds on\n*          Difu and Difl. If IJOB = 3 or 5, DIF(1:2) are 1-norm-based\n*          estimates of Difu and Difl.\n*          If M = 0 or N, DIF(1:2) = F-norm([A, B]).\n*          If IJOB = 0 or 1, DIF is not referenced.\n*\n*  WORK    (workspace/output) DOUBLE PRECISION array,\n*          dimension (MAX(1,LWORK)) \n*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.\n*\n*  LWORK   (input) INTEGER\n*          The dimension of the array WORK. LWORK >=  4*N+16.\n*          If IJOB = 1, 2 or 4, LWORK >= MAX(4*N+16, 2*M*(N-M)).\n*          If IJOB = 3 or 5, LWORK >= MAX(4*N+16, 4*M*(N-M)).\n*\n*          If LWORK = -1, then a workspace query is assumed; the routine\n*          only calculates the optimal size of the WORK array, returns\n*          this value as the first entry of the WORK array, and no error\n*          message related to LWORK is issued by XERBLA.\n*\n*  IWORK   (workspace/output) INTEGER array, dimension (MAX(1,LIWORK))\n*          On exit, if INFO = 0, IWORK(1) returns the optimal LIWORK.\n*\n*  LIWORK  (input) INTEGER\n*          The dimension of the array IWORK. LIWORK >= 1.\n*          If IJOB = 1, 2 or 4, LIWORK >=  N+6.\n*          If IJOB = 3 or 5, LIWORK >= MAX(2*M*(N-M), N+6).\n*\n*          If LIWORK = -1, then a workspace query is assumed; the\n*          routine only calculates the optimal size of the IWORK array,\n*          returns this value as the first entry of the IWORK array, and\n*          no error message related to LIWORK is issued by XERBLA.\n*\n*  INFO    (output) INTEGER\n*            =0: Successful exit.\n*            <0: If INFO = -i, the i-th argument had an illegal value.\n*            =1: Reordering of (A, B) failed because the transformed\n*                matrix pair (A, B) would be too far from generalized\n*                Schur form; the problem is very ill-conditioned.\n*                (A, B) may have been partially reordered.\n*                If requested, 0 is returned in DIF(*), PL and PR.\n*\n\n*  Further Details\n*  ===============\n*\n*  DTGSEN first collects the selected eigenvalues by computing\n*  orthogonal U and W that move them to the top left corner of (A, B).\n*  In other words, the selected eigenvalues are the eigenvalues of\n*  (A11, B11) in:\n*\n*                U'*(A, B)*W = (A11 A12) (B11 B12) n1\n*                              ( 0  A22),( 0  B22) n2\n*                                n1  n2    n1  n2\n*\n*  where N = n1+n2 and U' means the transpose of U. The first n1 columns\n*  of U and W span the specified pair of left and right eigenspaces\n*  (deflating subspaces) of (A, B).\n*\n*  If (A, B) has been obtained from the generalized real Schur\n*  decomposition of a matrix pair (C, D) = Q*(A, B)*Z', then the\n*  reordered generalized real Schur form of (C, D) is given by\n*\n*           (C, D) = (Q*U)*(U'*(A, B)*W)*(Z*W)',\n*\n*  and the first n1 columns of Q*U and Z*W span the corresponding\n*  deflating subspaces of (C, D) (Q and Z store Q*U and Z*W, resp.).\n*\n*  Note that if the selected eigenvalue is sufficiently ill-conditioned,\n*  then its value may differ significantly from its value before\n*  reordering.\n*\n*  The reciprocal condition numbers of the left and right eigenspaces\n*  spanned by the first n1 columns of U and W (or Q*U and Z*W) may\n*  be returned in DIF(1:2), corresponding to Difu and Difl, resp.\n*\n*  The Difu and Difl are defined as:\n*\n*       Difu[(A11, B11), (A22, B22)] = sigma-min( Zu )\n*  and\n*       Difl[(A11, B11), (A22, B22)] = Difu[(A22, B22), (A11, B11)],\n*\n*  where sigma-min(Zu) is the smallest singular value of the\n*  (2*n1*n2)-by-(2*n1*n2) matrix\n*\n*       Zu = [ kron(In2, A11)  -kron(A22', In1) ]\n*            [ kron(In2, B11)  -kron(B22', In1) ].\n*\n*  Here, Inx is the identity matrix of size nx and A22' is the\n*  transpose of A22. kron(X, Y) is the Kronecker product between\n*  the matrices X and Y.\n*\n*  When DIF(2) is small, small changes in (A, B) can cause large changes\n*  in the deflating subspace. An approximate (asymptotic) bound on the\n*  maximum angular error in the computed deflating subspaces is\n*\n*       EPS * norm((A, B)) / DIF(2),\n*\n*  where EPS is the machine precision.\n*\n*  The reciprocal norm of the projectors on the left and right\n*  eigenspaces associated with (A11, B11) may be returned in PL and PR.\n*  They are computed as follows. First we compute L and R so that\n*  P*(A, B)*Q is block diagonal, where\n*\n*       P = ( I -L ) n1           Q = ( I R ) n1\n*           ( 0  I ) n2    and        ( 0 I ) n2\n*             n1 n2                    n1 n2\n*\n*  and (L, R) is the solution to the generalized Sylvester equation\n*\n*       A11*R - L*A22 = -A12\n*       B11*R - L*B22 = -B12\n*\n*  Then PL = (F-norm(L)**2+1)**(-1/2) and PR = (F-norm(R)**2+1)**(-1/2).\n*  An approximate (asymptotic) bound on the average absolute error of\n*  the selected eigenvalues is\n*\n*       EPS * norm((A, B)) / PL.\n*\n*  There are also global error bounds which valid for perturbations up\n*  to a certain restriction:  A lower bound (x) on the smallest\n*  F-norm(E,F) for which an eigenvalue of (A11, B11) may move and\n*  coalesce with an eigenvalue of (A22, B22) under perturbation (E,F),\n*  (i.e. (A + E, B + F), is\n*\n*   x = min(Difu,Difl)/((1/(PL*PL)+1/(PR*PR))**(1/2)+2*max(1/PL,1/PR)).\n*\n*  An approximate bound on x can be computed from DIF(1:2), PL and PR.\n*\n*  If y = ( F-norm(E,F) / x) <= 1, the angles between the perturbed\n*  (L', R') and unperturbed (L, R) left and right deflating subspaces\n*  associated with the selected cluster in the (1,1)-blocks can be\n*  bounded as\n*\n*   max-angle(L, L') <= arctan( y * PL / (1 - y * (1 - PL * PL)**(1/2))\n*   max-angle(R, R') <= arctan( y * PR / (1 - y * (1 - PR * PR)**(1/2))\n*\n*  See LAPACK User's Guide section 4.11 or the following references\n*  for more information.\n*\n*  Note that if the default method for computing the Frobenius-norm-\n*  based estimate DIF is not wanted (see DLATDF), then the parameter\n*  IDIFJB (see below) should be changed from 3 to 4 (routine DLATDF\n*  (IJOB = 2 will be used)). See DTGSYL for more details.\n*\n*  Based on contributions by\n*     Bo Kagstrom and Peter Poromaa, Department of Computing Science,\n*     Umea University, S-901 87 Umea, Sweden.\n*\n*  References\n*  ==========\n*\n*  [1] B. Kagstrom; A Direct Method for Reordering Eigenvalues in the\n*      Generalized Real Schur Form of a Regular Matrix Pair (A, B), in\n*      M.S. Moonen et al (eds), Linear Algebra for Large Scale and\n*      Real-Time Applications, Kluwer Academic Publ. 1993, pp 195-218.\n*\n*  [2] B. Kagstrom and P. Poromaa; Computing Eigenspaces with Specified\n*      Eigenvalues of a Regular Matrix Pair (A, B) and Condition\n*      Estimation: Theory, Algorithms and Software,\n*      Report UMINF - 94.04, Department of Computing Science, Umea\n*      University, S-901 87 Umea, Sweden, 1994. Also as LAPACK Working\n*      Note 87. To appear in Numerical Algorithms, 1996.\n*\n*  [3] B. Kagstrom and P. Poromaa, LAPACK-Style Algorithms and Software\n*      for Solving the Generalized Sylvester Equation and Estimating the\n*      Separation between Regular Matrix Pairs, Report UMINF - 93.23,\n*      Department of Computing Science, Umea University, S-901 87 Umea,\n*      Sweden, December 1993, Revised April 1994, Also as LAPACK Working\n*      Note 75. To appear in ACM Trans. on Math. Software, Vol 22, No 1,\n*      1996.\n*\n*  =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  alphar, alphai, beta, m, pl, pr, dif, work, iwork, info, a, b, q, z = NumRu::Lapack.dtgsen( ijob, wantq, wantz, select, a, b, q, z, lwork, liwork, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 10)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 10)", argc);
  rblapack_ijob = argv[0];
  rblapack_wantq = argv[1];
  rblapack_wantz = argv[2];
  rblapack_select = argv[3];
  rblapack_a = argv[4];
  rblapack_b = argv[5];
  rblapack_q = argv[6];
  rblapack_z = argv[7];
  rblapack_lwork = argv[8];
  rblapack_liwork = argv[9];
  if (rb_options != Qnil) {
  }

  ijob = NUM2INT(rblapack_ijob);
  if (!NA_IsNArray(rblapack_a))
    rb_raise(rb_eArgError, "a (5th argument) must be NArray");
  if (NA_RANK(rblapack_a) != 2)
    rb_raise(rb_eArgError, "rank of a (5th argument) must be %d", 2);
  n = NA_SHAPE1(rblapack_a);
  lda = NA_SHAPE0(rblapack_a);
  if (NA_TYPE(rblapack_a) != NA_DFLOAT)
    rblapack_a = na_change_type(rblapack_a, NA_DFLOAT);
  a = NA_PTR_TYPE(rblapack_a, doublereal*);
  wantz = (rblapack_wantz == Qtrue);
  if (!NA_IsNArray(rblapack_b))
    rb_raise(rb_eArgError, "b (6th argument) must be NArray");
  if (NA_RANK(rblapack_b) != 2)
    rb_raise(rb_eArgError, "rank of b (6th argument) must be %d", 2);
  if (NA_SHAPE1(rblapack_b) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of b must be the same as shape 1 of a");
  ldb = NA_SHAPE0(rblapack_b);
  if (NA_TYPE(rblapack_b) != NA_DFLOAT)
    rblapack_b = na_change_type(rblapack_b, NA_DFLOAT);
  b = NA_PTR_TYPE(rblapack_b, doublereal*);
  liwork = NUM2INT(rblapack_liwork);
  wantq = (rblapack_wantq == Qtrue);
  if (!NA_IsNArray(rblapack_z))
    rb_raise(rb_eArgError, "z (8th argument) must be NArray");
  if (NA_RANK(rblapack_z) != 2)
    rb_raise(rb_eArgError, "rank of z (8th argument) must be %d", 2);
  if (NA_SHAPE1(rblapack_z) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of z must be the same as shape 1 of a");
  ldz = NA_SHAPE0(rblapack_z);
  if (NA_TYPE(rblapack_z) != NA_DFLOAT)
    rblapack_z = na_change_type(rblapack_z, NA_DFLOAT);
  z = NA_PTR_TYPE(rblapack_z, doublereal*);
  lwork = NUM2INT(rblapack_lwork);
  if (!NA_IsNArray(rblapack_q))
    rb_raise(rb_eArgError, "q (7th argument) must be NArray");
  if (NA_RANK(rblapack_q) != 2)
    rb_raise(rb_eArgError, "rank of q (7th argument) must be %d", 2);
  if (NA_SHAPE1(rblapack_q) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of q must be the same as shape 1 of a");
  ldq = NA_SHAPE0(rblapack_q);
  if (NA_TYPE(rblapack_q) != NA_DFLOAT)
    rblapack_q = na_change_type(rblapack_q, NA_DFLOAT);
  q = NA_PTR_TYPE(rblapack_q, doublereal*);
  if (!NA_IsNArray(rblapack_select))
    rb_raise(rb_eArgError, "select (4th argument) must be NArray");
  if (NA_RANK(rblapack_select) != 1)
    rb_raise(rb_eArgError, "rank of select (4th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_select) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of select must be the same as shape 1 of a");
  if (NA_TYPE(rblapack_select) != NA_LINT)
    rblapack_select = na_change_type(rblapack_select, NA_LINT);
  select = NA_PTR_TYPE(rblapack_select, logical*);
  {
    int shape[1];
    shape[0] = n;
    rblapack_alphar = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  alphar = NA_PTR_TYPE(rblapack_alphar, doublereal*);
  {
    int shape[1];
    shape[0] = n;
    rblapack_alphai = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  alphai = NA_PTR_TYPE(rblapack_alphai, doublereal*);
  {
    int shape[1];
    shape[0] = n;
    rblapack_beta = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  beta = NA_PTR_TYPE(rblapack_beta, doublereal*);
  {
    int shape[1];
    shape[0] = 2;
    rblapack_dif = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  dif = NA_PTR_TYPE(rblapack_dif, doublereal*);
  {
    int shape[1];
    shape[0] = MAX(1,lwork);
    rblapack_work = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  work = NA_PTR_TYPE(rblapack_work, doublereal*);
  {
    int shape[1];
    shape[0] = ijob==0 ? 0 : MAX(1,liwork);
    rblapack_iwork = na_make_object(NA_LINT, 1, shape, cNArray);
  }
  iwork = NA_PTR_TYPE(rblapack_iwork, integer*);
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
    shape[0] = ldb;
    shape[1] = n;
    rblapack_b_out__ = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  b_out__ = NA_PTR_TYPE(rblapack_b_out__, doublereal*);
  MEMCPY(b_out__, b, doublereal, NA_TOTAL(rblapack_b));
  rblapack_b = rblapack_b_out__;
  b = b_out__;
  {
    int shape[2];
    shape[0] = ldq;
    shape[1] = n;
    rblapack_q_out__ = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  q_out__ = NA_PTR_TYPE(rblapack_q_out__, doublereal*);
  MEMCPY(q_out__, q, doublereal, NA_TOTAL(rblapack_q));
  rblapack_q = rblapack_q_out__;
  q = q_out__;
  {
    int shape[2];
    shape[0] = ldz;
    shape[1] = n;
    rblapack_z_out__ = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  z_out__ = NA_PTR_TYPE(rblapack_z_out__, doublereal*);
  MEMCPY(z_out__, z, doublereal, NA_TOTAL(rblapack_z));
  rblapack_z = rblapack_z_out__;
  z = z_out__;

  dtgsen_(&ijob, &wantq, &wantz, select, &n, a, &lda, b, &ldb, alphar, alphai, beta, q, &ldq, z, &ldz, &m, &pl, &pr, dif, work, &lwork, iwork, &liwork, &info);

  rblapack_m = INT2NUM(m);
  rblapack_pl = rb_float_new((double)pl);
  rblapack_pr = rb_float_new((double)pr);
  rblapack_info = INT2NUM(info);
  return rb_ary_new3(14, rblapack_alphar, rblapack_alphai, rblapack_beta, rblapack_m, rblapack_pl, rblapack_pr, rblapack_dif, rblapack_work, rblapack_iwork, rblapack_info, rblapack_a, rblapack_b, rblapack_q, rblapack_z);
}

void
init_lapack_dtgsen(VALUE mLapack){
  rb_define_module_function(mLapack, "dtgsen", rblapack_dtgsen, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}

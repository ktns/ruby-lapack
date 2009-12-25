#include "rb_lapack.h"

static VALUE
rb_ilaenv(int argc, VALUE *argv, VALUE self){
  VALUE rb_ispec;
  integer ispec; 
  VALUE rb_name;
  char *name; 
  VALUE rb_opts;
  char *opts; 
  VALUE rb_n1;
  integer n1; 
  VALUE rb_n2;
  integer n2; 
  VALUE rb_n3;
  integer n3; 
  VALUE rb_n4;
  integer n4; 
  VALUE rb___out__;
  integer __out__; 


  if (argc == 0) {
    printf("%s\n", "USAGE:\n  __out__ = NumRu::Lapack.ilaenv( ispec, name, opts, n1, n2, n3, n4)\n    or\n  NumRu::Lapack.ilaenv  # print help\n\n\nFORTRAN MANUAL\n      INTEGER FUNCTION ILAENV( ISPEC, NAME, OPTS, N1, N2, N3, N4 )\n\n*  Purpose\n*  =======\n*\n*  ILAENV is called from the LAPACK routines to choose problem-dependent\n*  parameters for the local environment.  See ISPEC for a description of\n*  the parameters.\n*\n*  ILAENV returns an INTEGER\n*  if ILAENV >= 0: ILAENV returns the value of the parameter specified by ISPEC\n*  if ILAENV < 0:  if ILAENV = -k, the k-th argument had an illegal value.\n*\n*  This version provides a set of parameters which should give good,\n*  but not optimal, performance on many of the currently available\n*  computers.  Users are encouraged to modify this subroutine to set\n*  the tuning parameters for their particular machine using the option\n*  and problem size information in the arguments.\n*\n*  This routine will not function correctly if it is converted to all\n*  lower case.  Converting it to all upper case is allowed.\n*\n\n*  Arguments\n*  =========\n*\n*  ISPEC   (input) INTEGER\n*          Specifies the parameter to be returned as the value of\n*          ILAENV.\n*          = 1: the optimal blocksize; if this value is 1, an unblocked\n*               algorithm will give the best performance.\n*          = 2: the minimum block size for which the block routine\n*               should be used; if the usable block size is less than\n*               this value, an unblocked routine should be used.\n*          = 3: the crossover point (in a block routine, for N less\n*               than this value, an unblocked routine should be used)\n*          = 4: the number of shifts, used in the nonsymmetric\n*               eigenvalue routines (DEPRECATED)\n*          = 5: the minimum column dimension for blocking to be used;\n*               rectangular blocks must have dimension at least k by m,\n*               where k is given by ILAENV(2,...) and m by ILAENV(5,...)\n*          = 6: the crossover point for the SVD (when reducing an m by n\n*               matrix to bidiagonal form, if max(m,n)/min(m,n) exceeds\n*               this value, a QR factorization is used first to reduce\n*               the matrix to a triangular form.)\n*          = 7: the number of processors\n*          = 8: the crossover point for the multishift QR method\n*               for nonsymmetric eigenvalue problems (DEPRECATED)\n*          = 9: maximum size of the subproblems at the bottom of the\n*               computation tree in the divide-and-conquer algorithm\n*               (used by xGELSD and xGESDD)\n*          =10: ieee NaN arithmetic can be trusted not to trap\n*          =11: infinity arithmetic can be trusted not to trap\n*          12 <= ISPEC <= 16:\n*               xHSEQR or one of its subroutines,\n*               see IPARMQ for detailed explanation\n*\n*  NAME    (input) CHARACTER*(*)\n*          The name of the calling subroutine, in either upper case or\n*          lower case.\n*\n*  OPTS    (input) CHARACTER*(*)\n*          The character options to the subroutine NAME, concatenated\n*          into a single character string.  For example, UPLO = 'U',\n*          TRANS = 'T', and DIAG = 'N' for a triangular routine would\n*          be specified as OPTS = 'UTN'.\n*\n*  N1      (input) INTEGER\n*  N2      (input) INTEGER\n*  N3      (input) INTEGER\n*  N4      (input) INTEGER\n*          Problem dimensions for the subroutine NAME; these may not all\n*          be required.\n*\n\n*  Further Details\n*  ===============\n*\n*  The following conventions have been used when calling ILAENV from the\n*  LAPACK routines:\n*  1)  OPTS is a concatenation of all of the character options to\n*      subroutine NAME, in the same order that they appear in the\n*      argument list for NAME, even if they are not used in determining\n*      the value of the parameter specified by ISPEC.\n*  2)  The problem dimensions N1, N2, N3, N4 are specified in the order\n*      that they appear in the argument list for NAME.  N1 is used\n*      first, N2 second, and so on, and unused problem dimensions are\n*      passed a value of -1.\n*  3)  The parameter value returned by ILAENV is checked for validity in\n*      the calling subroutine.  For example, ILAENV is used to retrieve\n*      the optimal blocksize for STRTRI as follows:\n*\n*      NB = ILAENV( 1, 'STRTRI', UPLO // DIAG, N, -1, -1, -1 )\n*      IF( NB.LE.1 ) NB = MAX( 1, N )\n*\n*  =====================================================================\n*\n*     .. Local Scalars ..\n      INTEGER            I, IC, IZ, NB, NBMIN, NX\n      LOGICAL            CNAME, SNAME\n      CHARACTER          C1*1, C2*2, C4*2, C3*3, SUBNAM*6\n*     ..\n*     .. Intrinsic Functions ..\n      INTRINSIC          CHAR, ICHAR, INT, MIN, REAL\n*     ..\n*     .. External Functions ..\n      INTEGER            IEEECK, IPARMQ\n      EXTERNAL           IEEECK, IPARMQ\n*     ..\n\n");
    return Qnil;
  }
  if (argc != 7)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 7)", argc);
  rb_ispec = argv[0];
  rb_name = argv[1];
  rb_opts = argv[2];
  rb_n1 = argv[3];
  rb_n2 = argv[4];
  rb_n3 = argv[5];
  rb_n4 = argv[6];

  ispec = NUM2INT(rb_ispec);
  n1 = NUM2INT(rb_n1);
  n2 = NUM2INT(rb_n2);
  n3 = NUM2INT(rb_n3);
  n4 = NUM2INT(rb_n4);
  name = StringValueCStr(rb_name);
  opts = StringValueCStr(rb_opts);

  __out__ = ilaenv_(&ispec, name, opts, &n1, &n2, &n3, &n4);

  rb___out__ = INT2NUM(__out__);
  return rb___out__;
}

void
init_lapack_ilaenv(VALUE mLapack){
  rb_define_module_function(mLapack, "ilaenv", rb_ilaenv, -1);
}

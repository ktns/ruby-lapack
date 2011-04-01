$:.unshift File.dirname(__FILE__)

require "yaml"
require "digest/md5"
require "pp"
require "common"

CTYPES = {
  "INTEGER" => "integer",
  "CHARACTER" => "char",
  "REAL" => "real",
  "DOUBLE PRECISION" => "doublereal",
  "COMPLEX" => "complex",
  "COMPLEX*16" => "doublecomplex",
  "DOUBLE COMPLEX" => "doublecomplex",
  "LOGICAL" => "logical"
}

ARGS = {
  "csyequb" => {
    "uplo" => {:intent => "input", :type => "char"},
    "work" => {:intent => "input", :type => "real"}
  },
  "zsyequb" => {
    "uplo" => {:intent => "input", :type => "char"},
    "work" => {:intent => "input", :type => "doublereal"}
  },
  "zhfrk" => {
    "uplo" => {:intent => "input", :type => "char"},
    "trans" => {:intent => "input", :type => "char"},
    "n" => {:intent => "input", :type => "integer"},
    "k" => {:intent => "input", :type => "integer"},
    "alpha" => {:intent => "input", :type => "doublereal"},
    "beta" => {:intent => "input", :type => "doublereal"},
    "a" => {:intent => "input", :type => "doublereal"},
    "lda" => {:intent => "input", :type => "integer"},
    "c" => {:intent => "input", :type => "doublereal"},
  },
  "zpstf2" => {
    "work" => {:intent => "input", :type => "doublereal"}
  },
  "ssfrk" => {
    "uplo" => {:intent => "input", :type => "char"},
    "trans" => {:intent => "input", :type => "char"},
    "n" => {:intent => "input", :type => "integer"},
    "k" => {:intent => "input", :type => "integer"},
    "alpha" => {:intent => "input", :type => "real"},
    "beta" => {:intent => "input", :type => "real"},
    "a" => {:intent => "input", :type => "real"},
    "lda" => {:intent => "input", :type => "integer"},
    "c" => {:intent => "input", :type => "real"},
  },
  "dsfrk" => {
    "uplo" => {:intent => "input", :type => "char"},
    "trans" => {:intent => "input", :type => "char"},
    "n" => {:intent => "input", :type => "integer"},
    "k" => {:intent => "input", :type => "integer"},
    "alpha" => {:intent => "input", :type => "doublereal"},
    "beta" => {:intent => "input", :type => "doublereal"},
    "a" => {:intent => "input", :type => "doublereal"},
    "lda" => {:intent => "input", :type => "integer"},
    "c" => {:intent => "input", :type => "doublereal"},
  },
  "spstf2" => {
    "work" => {:intent => "input", :type => "real"}
  },
  "slarrf" => {
    "clgapl" => {:intent => "input", :type => "real"},
    "clgapr" => {:intent => "input", :type => "real"},
  },
  "cheequb" => {
    "uplo" => {:intent => "input", :type => "char"},
    "work" => {:intent => "input", :type => "complex"}
  },
  "zheequb" => {
    "uplo" => {:intent => "input", :type => "char"},
    "work" => {:intent => "input", :type => "doublecomplex"}
  },
  "dsyequb" => {
    "uplo" => {:intent => "input", :type => "char"},
    "work" => {:intent => "input", :type => "doublereal"}
  },
  "chesvxx" => {
    "uplo" => {:intent => "input", :type => "char"}
  },
  "zhesvxx" => {
    "uplo" => {:intent => "input", :type => "char"}
  },
  "dgesvj" => {
    "lwork" => {:intent => "input", :type => "integer"}
  },
  "dpstf2" => {
    "work" => {:intent => "input", :type => "doublereal"}
  },
  "cpstf2" => {
    "work" => {:intent => "input", :type => "real"}
  },
  "spstrf" => {
    "work" => {:intent => "input", :type => "real"}
  },
  "dpstrf" => {
    "work" => {:intent => "input", :type => "doublereal"}
  },
  "zpstrf" => {
    "work" => {:intent => "input", :type => "doublereal"}
  },
  "cpstrf" => {
    "work" => {:intent => "input", :type => "doublereal"}
  },
  "sgesvj" => {
    "lwork" => {:intent => "input", :type => "integer"}
  },
  "sla_lin_berr" => {
    "berr" => {:intent => "output", :type => "real"}
  },
  "dla_lin_berr" => {
    "berr" => {:intent => "output", :type => "doublereal"}
  },
  "ssyequb" => {
    "uplo" => {:intent => "input", :type => "char"},
    "work" => {:intent => "input", :type => "real"},
  },
  "slasq4" => {
    "n0in" => {:intent => "input", :type => "integer"},
  },
  "dlasq4" => {
    "n0in" => {:intent => "input", :type => "integer"},
  },
  "slaqr1" => {
    "si1" => {:intent => "input", :type => "real"},
    "sr2" => {:intent => "input", :type => "real"},
    "si2" => {:intent => "input", :type => "real"},
  },
  "dlaqr1" => {
    "si1" => {:intent => "input", :type => "doublereal"},
    "sr2" => {:intent => "input", :type => "doublereal"},
    "si2" => {:intent => "input", :type => "doublereal"},
  },
  "claqr1" => {
    "s2" => {:intent => "input", :type => "complex"},
  },
  "zlaqr1" => {
    "s2" => {:intent => "input", :type => "doublecomplex"},
  },
  "slarrf" => {
    "clgapl" => {:intent => "input", :type => "real"},
    "clgapr" => {:intent => "input", :type => "real"},
    "info" => {:intent => "output", :type => "integer"},
  },
  "dlarrf" => {
    "clgapl" => {:intent => "input", :type => "doublereal"},
    "clgapr" => {:intent => "input", :type => "doublereal"},
    "info" => {:intent => "output", :type => "integer"},
  }
}

DIMS = {
  "strttp" => {
    "ap" => ["n*(n+1)/2"]
  },
  "dtrttp" => {
    "ap" => ["n*(n+1)/2"]
  },
  "stftri" => {
    "a" => ["n*(n+1)/2"]
  },
  "dtftri" => {
    "a" => ["n*(n+1)/2"]
  },
  "zpstrf" => {
    "work" => ["2*n"]
  },
  "spftrf" => {
    "a" => ["n*(n+1)/2"]
  },
  "cpftrf" => {
    "a" => ["n*(n+1)/2"]
  },
  "cpftrs" => {
    "a" => ["n*(n+1)/2"]
  },
  "zpftrs" => {
    "a" => ["n*(n+1)/2"]
  },
  "strttf" => {
    "arf" => ["n*(n+1)/2"]
  },
  "dtrttf" => {
    "arf" => ["n*(n+1)/2"]
  },
  "sla_gbrfsx_extended" => {
    "ab" => ["ldab", "n"]
  },
  "stgex2" => {
    "work" => ["lwork"]
  },
  "dtgex2" => {
    "work" => ["lwork"]
  },
  "dlasd1" => {
    "d" => ["n"]
  },
  "slaruv" => {
    "x" => ["MAX(1,n)"]
  },
  "dlaruv" => {
    "x" => ["MAX(1,n)"]
  },
  "slasyf" => {
    "w" => ["ldw","MAX(1,nb)"]
  },
  "dlasyf" => {
    "w" => ["ldw","MAX(1,nb)"]
  },
  "clasyf" => {
    "w" => ["ldw","MAX(1,nb)"]
  },
  "zlasyf" => {
    "w" => ["ldw","MAX(1,nb)"]
  },
  "slaeda" => {
    "qptr" => ["ldqptr"]
  },
  "dlaeda" => {
    "qptr" => ["ldqptr"]
  },
  "slasdt" => {
    "inode" => ["MAX(1,n)"],
    "ndiml" => ["MAX(1,n)"],
    "ndimr" => ["MAX(1,n)"],
  },
  "dlasdt" => {
    "inode" => ["MAX(1,n)"],
    "ndiml" => ["MAX(1,n)"],
    "ndimr" => ["MAX(1,n)"],
  },
  "sgbequ" => {
    "r" => ["MAX(1,m)"]
  },
  "dgbequ" => {
    "r" => ["MAX(1,m)"]
  },
  "cgbequ" => {
    "r" => ["MAX(1,m)"]
  },
  "zgbequ" => {
    "r" => ["MAX(1,m)"]
  },
  "slaed9" => {
    "d" => ["MAX(1,n)"],
    "q" => ["ldq", "MAX(1,n)"]
  },
  "dlaed9" => {
    "d" => ["MAX(1,n)"],
    "q" => ["ldq", "MAX(1,n)"]
  },
  "slarnv" => {
    "x" => ["MAX(1,n)"],
  },
  "dlarnv" => {
    "x" => ["MAX(1,n)"],
  },
  "clarnv" => {
    "x" => ["MAX(1,n)"],
  },
  "zlarnv" => {
    "x" => ["MAX(1,n)"]
  },
  "slabrd" => {
    "d" => ["MAX(1,nb)"],
    "e" => ["MAX(1,nb)"],
    "tauq" => ["MAX(1,nb)"],
    "taup" => ["MAX(1,nb)"],
    "x" => ["ldx", "MAX(1,nb)"],
    "y" => ["ldy", "MAX(1,nb)"],
  },
  "dlabrd" => {
    "d" => ["MAX(1,nb)"],
    "e" => ["MAX(1,nb)"],
    "tauq" => ["MAX(1,nb)"],
    "taup" => ["MAX(1,nb)"],
    "x" => ["ldx", "MAX(1,nb)"],
    "y" => ["ldy", "MAX(1,nb)"],
  },
  "clabrd" => {
    "d" => ["MAX(1,nb)"],
    "e" => ["MAX(1,nb)"],
    "tauq" => ["MAX(1,nb)"],
    "taup" => ["MAX(1,nb)"],
    "x" => ["ldx", "MAX(1,nb)"],
    "y" => ["ldy", "MAX(1,nb)"],
  },
  "zlabrd" => {
    "d" => ["MAX(1,nb)"],
    "e" => ["MAX(1,nb)"],
    "tauq" => ["MAX(1,nb)"],
    "taup" => ["MAX(1,nb)"],
    "x" => ["ldx", "MAX(1,nb)"],
    "y" => ["ldy", "MAX(1,nb)"],
  },
  "slahr2" => {
    "tau" => ["MAX(1,nb)"],
    "t" => ["ldt","MAX(1,nb)"],
    "y" => ["ldy","MAX(1,nb)"],
  },
  "dlahr2" => {
    "tau" => ["MAX(1,nb)"],
    "t" => ["ldt","MAX(1,nb)"],
    "y" => ["ldy","MAX(1,nb)"],
  },
  "clahr2" => {
    "tau" => ["MAX(1,nb)"],
    "t" => ["ldt","MAX(1,nb)"],
    "y" => ["ldy","MAX(1,nb)"],
  },
  "zlahr2" => {
    "tau" => ["MAX(1,nb)"],
    "t" => ["ldt","MAX(1,nb)"],
    "y" => ["ldy","MAX(1,nb)"],
  },
  "slahrd" => {
    "tau" => ["MAX(1,nb)"],
    "t" => ["ldt","MAX(1,nb)"],
    "y" => ["ldy","MAX(1,nb)"],
  },
  "dlahrd" => {
    "tau" => ["MAX(1,nb)"],
    "t" => ["ldt","MAX(1,nb)"],
    "y" => ["ldy","MAX(1,nb)"],
  },
  "clahrd" => {
    "tau" => ["MAX(1,nb)"],
    "t" => ["ldt","MAX(1,nb)"],
    "y" => ["ldy","MAX(1,nb)"],
  },
  "zlahrd" => {
    "tau" => ["MAX(1,nb)"],
    "t" => ["ldt","MAX(1,nb)"],
    "y" => ["ldy","MAX(1,nb)"],
  },
  "slaqr2" => {
    "sr" => ["MAX(1,kbot)"],
    "si" => ["MAX(1,kbot)"],
    "v" => ["ldv", "MAX(1,nw)"],
    "t" => ["ldt", "MAX(1,nw)"],
    "wv" => ["ldwv", "MAX(1,nw)"],
  },
  "dlaqr2" => {
    "sr" => ["MAX(1,kbot)"],
    "si" => ["MAX(1,kbot)"],
    "v" => ["ldv", "MAX(1,nw)"],
    "t" => ["ldt", "MAX(1,nw)"],
    "wv" => ["ldwv", "MAX(1,nw)"],
  },
  "claqr2" => {
    "sh" => ["MAX(1,kbot)"],
    "v" => ["ldv", "MAX(1,nw)"],
    "t" => ["ldv", "MAX(1,nw)"],
    "wv" => ["ldv", "MAX(1,nw)"],
  },
  "zlaqr2" => {
    "sh" => ["MAX(1,kbot)"],
    "v" => ["ldv", "MAX(1,nw)"],
    "t" => ["ldv", "MAX(1,nw)"],
    "wv" => ["ldv", "MAX(1,nw)"],
  },
  "slaqr3" => {
    "sr" => ["MAX(1,kbot)"],
    "si" => ["MAX(1,kbot)"],
    "v" => ["ldv", "MAX(1,nw)"],
    "t" => ["ldt", "MAX(1,nw)"],
    "wv" => ["ldwv", "MAX(1,nw)"],
  },
  "dlaqr3" => {
    "sr" => ["MAX(1,kbot)"],
    "si" => ["MAX(1,kbot)"],
    "v" => ["ldv", "MAX(1,nw)"],
    "t" => ["ldt", "MAX(1,nw)"],
    "wv" => ["ldwv", "MAX(1,nw)"],
  },
  "claqr3" => {
    "sh" => ["MAX(1,kbot)"],
    "v" => ["ldv", "MAX(1,nw)"],
    "t" => ["ldv", "MAX(1,nw)"],
    "wv" => ["ldv", "MAX(1,nw)"],
  },
  "zlaqr3" => {
    "sh" => ["MAX(1,kbot)"],
    "v" => ["ldv", "MAX(1,nw)"],
    "t" => ["ldv", "MAX(1,nw)"],
    "wv" => ["ldv", "MAX(1,nw)"],
  },
  "slatrd" => {
    "w" => ["ldw", "MAX(n,nb)"]
  },
  "dlatrd" => {
    "w" => ["ldw", "MAX(n,nb)"]
  },
  "clatrd" => {
    "w" => ["ldw", "MAX(n,nb)"]
  },
  "zlatrd" => {
    "w" => ["ldw", "MAX(n,nb)"]
  },
  "clahef" => {
    "w" => ["ldw", "MAX(n,nb)"]
  },
  "zlahef" => {
    "w" => ["ldw", "MAX(n,nb)"]
  },
  "sopgtr" => {
    "ap" => ["ldap"],
    "tau" => ["ldtau"],
  },
  "dopgtr" => {
    "ap" => ["ldap"],
    "tau" => ["ldtau"],
  },
  "cupgtr" => {
    "ap" => ["ldap"],
    "tau" => ["ldtau"],
  },
  "zupgtr" => {
    "ap" => ["ldap"],
    "tau" => ["ldtau"],
  },
  "ssptrd" => {
    "ap" => ["ldap"]
  },
  "dsptrd" => {
    "ap" => ["ldap"]
  },
  "chptrd" => {
    "ap" => ["ldap"]
  },
  "zhptrd" => {
    "ap" => ["ldap"]
  },
  "sspgv" => {
    "ap" => ["ldap"]
  },
  "chpev" => {
    "ap" => ["ldap"],
  },
  "zhpev" => {
    "ap" => ["ldap"],
  },
  "chpevx" => {
    "ap" => ["ldap"],
  },
  "zhpevx" => {
    "ap" => ["ldap"],
  },
  "dspgv" => {
    "ap" => ["ldap"]
  },
  "sppequ" => {
    "ap" => ["ldap"]
  },
  "dppequ" => {
    "ap" => ["ldap"]
  },
  "cppequ" => {
    "ap" => ["ldap"]
  },
  "zppequ" => {
    "ap" => ["ldap"]
  },
  "sspgvd" => {
    "ap" => ["ldap"]
  },
  "dspgvd" => {
    "ap" => ["ldap"]
  },
  "stpcon" => {
    "ap" => ["ldap"]
  },
  "dtpcon" => {
    "ap" => ["ldap"]
  },
  "ctpcon" => {
    "ap" => ["ldap"]
  },
  "ztpcon" => {
    "ap" => ["ldap"]
  },
  "sspgvx" => {
    "ap" => ["ldap"],
  },
  "dspgvx" => {
    "ap" => ["ldap"],
  },
  "chpevd" => {
    "ap" => ["ldap"],
  },
  "zhpevd" => {
    "ap" => ["ldap"],
  },
  "sspev" => {
    "ap" => ["ldap"],
  },
  "dspev" => {
    "ap" => ["ldap"],
  },
  "sspevd" => {
    "ap" => ["ldap"],
  },
  "dspevd" => {
    "ap" => ["ldap"],
  },
  "sspevx" => {
    "ap" => ["ldap"],
  },
  "dspevx" => {
    "ap" => ["ldap"],
  },
  "sppcon" => {
    "ap" => ["ldap"]
  },
  "dppcon" => {
    "ap" => ["ldap"]
  },
  "cppcon" => {
    "ap" => ["ldap"]
  },
  "zppcon" => {
    "ap" => ["ldap"]
  },
  "chpgv" => {
    "ap" => ["ldap"]
  },
  "zhpgv" => {
    "ap" => ["ldap"]
  },
  "chpgvd" => {
    "ap" => ["ldap"]
  },
  "zhpgvd" => {
    "ap" => ["ldap"]
  },
  "chpgvx" => {
    "ap" => ["ldap"]
  },
  "zhpgvx" => {
    "ap" => ["ldap"]
  },
  "chptrf" => {
    "ap" => ["ldap"]
  },
  "zhptrf" => {
    "ap" => ["ldap"]
  },
  "ssptrf" => {
    "ap" => ["ldap"]
  },
  "dsptrf" => {
    "ap" => ["ldap"]
  },
  "csptrf" => {
    "ap" => ["ldap"]
  },
  "zsptrf" => {
    "ap" => ["ldap"]
  },
  "slaqr0" => {
    "z" => ["ldz","ihi"]
  },
  "dlaqr0" => {
    "z" => ["ldz","ihi"]
  },
  "claqr0" => {
    "z" => ["ldz","ihi"],
    "work" => ["MAX(1,lwork)"]
  },
  "zlaqr0" => {
    "z" => ["wantz ? ldz : 0","wantz ? ihi : 0"],
    "work" => ["MAX(1,lwork)"]
  },
  "slaqr4" => {
    "z" => ["ldz","ihi"]
  },
  "dlaqr4" => {
    "z" => ["ldz","ihi"]
  },
  "claqr4" => {
    "z" => ["ldz","ihi"]
  },
  "zlaqr4" => {
    "z" => ["ldz","ihi"]
  },
  "slaqr5" => {
    "z" => ["wantz ? ldz : 0","wantz ? ihiz : 0"],
    "wh" => ["ldwh", "MAX(1,nh)"]
  },
  "dlaqr5" => {
    "z" => ["wantz ? ldz : 0","wantz ? ihiz : 0"],
    "wh" => ["ldwh", "MAX(1,nh)"]
  },
  "claqr5" => {
    "z" => ["wantz ? ldz : 0","wantz ? ihiz : 0"],
    "wh" => ["ldwh", "MAX(1,nh)"]
  },
  "zlaqr5" => {
    "z" => ["wantz ? ldz : 0","wantz ? ihiz : 0"],
    "wh" => ["ldwh", "MAX(1,nh)"]
  },
  "slasd0" => {
    "u" => ["ldu", "n"]
  },
  "dlasd0" => {
    "u" => ["ldu", "n"]
  },
  "slasd4" => {
    "delta" => ["n"],
    "work" => ["n"]
  },
  "dlasd4" => {
    "delta" => ["n"],
    "work" => ["n"]
  },
  "slasdq" => {
    "e" => ['sqre==0 ? n-1 : sqre==1 ? n : 0'],
    "work" => ["4*n"]
  },
  "dlasdq" => {
    "e" => ['sqre==0 ? n-1 : sqre==1 ? n : 0'],
    "work" => ["4*n"]
  },
  "slaed3" => {
    "q2" => ["n","n"]
  },
  "dlaed3" => {
    "q2" => ["n","n"]
  },
  "slaed4" => {
    "delta" => ["n"]
  },
  "dlaed4" => {
    "delta" => ["n"]
  },
  "slaed8" => {
    "q" => ['icompq==0 ? 0 : ldq', 'icompq==0 ? 0 : n'],
    "q2" => ['icompq==0 ? 0 : ldq2', 'icompq==0 ? 0 : n']
  },
  "dlaed8" => {
    "q" => ['icompq==0 ? 0 : ldq', 'icompq==0 ? 0 : n'],
    "q2" => ['icompq==0 ? 0 : ldq2', 'icompq==0 ? 0 : n']
  },
  "claed8" => {
    "q2" => ['ldq2', 'n']
  },
  "zlaed8" => {
    "q2" => ['ldq2', 'n']
  },
  "slasq3" => {
    "z" => ['4*n0']
  },
  "dlasq3" => {
    "z" => ['4*n0']
  },
  "slasq4" => {
    "z" => ['4*n0']
  },
  "dlasq4" => {
    "z" => ['4*n0']
  },
  "slasq5" => {
    "z" => ['4*n0']
  },
  "dlasq5" => {
    "z" => ['4*n0']
  },
  "slasq6" => {
    "z" => ['4*n0']
  },
  "dlasq6" => {
    "z" => ['4*n0']
  },
  "slazq3" => {
    "z" => ['4*n0']
  },
  "dlazq3" => {
    "z" => ['4*n0']
  },
  "slazq4" => {
    "z" => ['4*n0']
  },
  "dlazq4" => {
    "z" => ['4*n0']
  },
  "slahqr" => {
    "z" => ['wantz ? ldz : 0', 'wantz ? n : 0']
  },
  "dlahqr" => {
    "z" => ['wantz ? ldz : 0', 'wantz ? n : 0']
  },
  "clahqr" => {
    "z" => ['wantz ? ldz : 0', 'wantz ? n : 0']
  },
  "zlahqr" => {
    "z" => ['wantz ? ldz : 0', 'wantz ? n : 0']
  },
  "cpbsvx" => {
    "afb" => ["ldafb","n"]
  },
  "dgesdd" => {
    "vt" => ["ldvt","n"]
  },
  "slasda" => {
    "u" => ["ldu", "MAX(1,smlsiz)"],
    "s" => ['icompq==1 ? n : icompq==0 ? 1 : 0']
  },
  "dlasda" => {
    "u" => ["ldu", "MAX(1,smlsiz)"],
    "s" => ['icompq==1 ? n : icompq==0 ? 1 : 0']
  },
  "clalsa" => {
    "rwork" => ['MAX(n,(smlsiz+1)*nrhs*3)']
  },
  "zlalsa" => {
    "rwork" => ['MAX(n,(smlsiz+1)*nrhs*3)']
  },
  "stgsen" => {
    "iwork" => ['ijob==0 ? 0 : MAX(1,liwork)']
  },
  "dtgsen" => {
    "iwork" => ['ijob==0 ? 0 : MAX(1,liwork)']
  },
  "ztgsen" => {
    "work" => ['ijob==0 ? 0 : MAX(1,lwork)'],
    "iwork" => ['ijob==0 ? 0 : MAX(1,liwork)']
  },
  "ctgsen" => {
    "work" => ['ijob==0 ? 0 : MAX(1,lwork)'],
    "iwork" => ['ijob==0 ? 0 : MAX(1,liwork)']
  },
  "slarfx" => {
    "work" => ['lsame_(&side,"L") ? n : lsame_(&side,"R") ? m : 0']
  },
  "dlarfx" => {
    "work" => ['lsame_(&side,"L") ? n : lsame_(&side,"R") ? m : 0']
  },
  "clarfx" => {
    "work" => ['lsame_(&side,"L") ? n : lsame_(&side,"R") ? m : 0']
  },
  "zlarfx" => {
    "work" => ['lsame_(&side,"L") ? n : lsame_(&side,"R") ? m : 0']
  },
  "slaebz" => {
    "nval" => ['(ijob==1||ijob==2) ? 0 : ijob==3 ? minp : 0'],
    "c" => ['ijob==1 ? 0 : (ijob==2||ijob==3) ? mmax : 0'],
    "nab" => ['mmax', '2']
  },
  "dlaebz" => {
    "nval" => ['(ijob==1||ijob==2) ? 0 : ijob==3 ? minp : 0'],
    "c" => ['ijob==1 ? 0 : (ijob==2||ijob==3) ? mmax : 0'],
    "nab" => ['mmax', '2']
  },
  "sbdsdc" => {
    "u" => ['lsame_(&compq,"I") ? ldu : 0','lsame_(&compq,"I") ? n : 0'],
    "vt" => ['lsame_(&compq,"I") ? ldvt : 0','lsame_(&compq,"I") ? n : 0'],
    "q" => ['lsame_(&compq,"I") ? ldq : 0'],
    "iq" => ['lsame_(&compq,"I") ? ldiq : 0'],
    "work" => ['MAX(1,lwork)']
  },
  "dbdsdc" => {
    "u" => ['lsame_(&compq,"I") ? ldu : 0','lsame_(&compq,"I") ? n : 0'],
    "vt" => ['lsame_(&compq,"I") ? ldvt : 0','lsame_(&compq,"I") ? n : 0'],
    "q" => ['lsame_(&compq,"I") ? ldq : 0'],
    "iq" => ['lsame_(&compq,"I") ? ldiq : 0'],
    "work" => ['MAX(1,lwork)']
  },
  "sgelsx" => {
    "work" => ['MAX((MIN(m,n))+3*n,2*(MIN(m,n))*nrhs)']
  },
  "dgelsx" => {
    "work" => ['MAX((MIN(m,n))+3*n,2*(MIN(m,n))*nrhs)']
  },
  "cgelsx" => {
    "work" => ['MIN(m,n) + MAX(n,2*(MIN(m,n))+nrhs)']
  },
  "zgelsx" => {
    "work" => ['MIN(m,n) + MAX(n,2*(MIN(m,n))+nrhs)']
  },
  "sgeevx" => {
    "iwork" => ['(lsame_(&sense,"N")||lsame_(&sense,"E")) ? 0 : 2*n-2']
  },
  "dgeevx" => {
    "iwork" => ['(lsame_(&sense,"N")||lsame_(&sense,"E")) ? 0 : 2*n-2']
  },
  "ctgex2" => {
    "q" => ['wantq ? ldq : 0', 'wantq ? n : 0'],
    "z" => ['wantq ? ldz : 0', 'wantq ? n : 0']
  },
  "ztgex2" => {
    "q" => ['wantq ? ldq : 0', 'wantq ? n : 0'],
    "z" => ['wantq ? ldz : 0', 'wantq ? n : 0']
  },
  "clalsd" => {
    "rwork" => ['9*n+2*n*smlsiz+8*n*nlvl+3*smlsiz*nrhs+(smlsiz+1)*(smlsiz+1)']
  },
  "zlalsd" => {
    "rwork" => ['9*n+2*n*smlsiz+8*n*nlvl+3*smlsiz*nrhs+(smlsiz+1)*(smlsiz+1)']
  },
  "ssbgvx" => {
    "work" => ['7*n'],
    "iwork" => ['5*n'],
  }
}


TYPES = {
  "slarrc" => {
    "vl" => "real",
    "vu" => "real",
    "d" => "real",
    "e" => "real",
    "pivmin" => "real"
  },
  "slarrb" => {
    "pivmin" => "real",
    "spdiam" => "real"
  },
  "slarre" => {
    "pivmin" => "real"
  },
  "slarrf" => {
    "pivmin" => "real",
    "spdiam" => "real"
  },
  "dlarrf" => {
    "spdiam" => "doublereal"
  },
  "slarrj" => {
    "pivmin" => "real",
    "spdiam" => "real"
  },
  "slarrv" => {
    "pivmin" => "real"
  },
  "clarrv" =>{
    "pivmin" => "real"
  },
  "zggbal" => {
    "work" => "doublereal"
  },
  "clarcm" => {
    "b" => "complex"
  },
  "zlarcm" => {
    "b" => "doublecomplex"
  },
  "clatdf" => {
    "z" => "complex",
    "rhs" => "complex",
  },
  "zlatdf" => {
    "z" => "doublecomplex",
    "rhs" => "doublecomplex",
  },
  "dggbal" => {
    "work" => "doublereal"
  },
  "zggevx" => {
    "rwork" => "doublereal"
  },
  "dlazq3" => {
    "dmin1" => "doublereal",
    "dmin2" => "doublereal",
    "dn" => "doublereal",
    "dn1" => "doublereal",
    "dn2" => "doublereal",
    "tau" => "doublereal",
  },
  "zlag2c" => {
    "a" => "doublecomplex",
    "sa" => "complex"
  },
  "clag2z" => {
    "a" => "doublecomplex",
    "sa" => "complex"
  },
  "cptts2" => {
    "b" => "complex"
  },
  "zptts2" => {
    "b" => "doublecomplex"
  },
  "cpttrs" => {
    "b" => "complex"
  },
  "zpttrs" => {
    "b" => "doublecomplex"
  },
  "zpbrfs" => {
    "ab" => "doublecomplex"
  }
}



SUBSTS = Hash.new
SUBSTS["dlasd1"] = {"n" => "nl+nr+1"}
%w(slaeda dlaeda).each{|n| SUBSTS[n] = {"n" => "ldqptr-2"}}
%w(sgbbrd dgbbrd cgbbrd zgbbrd).each{|n| SUBSTS[n] = {"m" => "ldab"}}
%w(sggsvp dggsvp cggsvp zggsvp sgeequ dgeequ cgeequ zgeequ cungr2 zungr2 cungl2 zungl2 sorgr2 dorgr2 sorgl2 dorgl2 stzrqf dtzrqf ctzrqf ztzrqf stzrzf dtzrzf ctzrzf ztzrzf sgerq2 dgerq2 cgerq2 zgerq2 sgelq2 dgelq2 cgelq2 zgelq2 slatrz dlatrz clatrz zlatrz).each{|n| SUBSTS[n] = {"m" => "lda"}}
%w(slaqps dlaqps claqps zlaqps).each{|n| SUBSTS[n] = {"kb" => "nb"}}
%w(sstevx dstevx sstein dstein cstein zstein).each{|n| SUBSTS[n] = {"m" => "n"}}
%w(sopgtr dopgtr cupgtr zupgtr).each{|n| SUBSTS[n] = {"n" => "ldtau+1"}}
%w(stgsna dtgsna ctgsna ztgsna strsna dtrsna ctrsna ztrsna).each{|n| SUBSTS[n] = {"mm" => "m"}}
%w(sggsvd dggsvd cggsvd zggsvd).each{|n| SUBSTS[n] = {"m" => "lda", "p" => "ldb"}}
SUBSTS["dlansp"] = {"lwork" => '(lsame_(&norm,"I") || lsame_(&norm,"1") || lsame_(&norm,"0")) ? n : 0'}
%w(ssptrd dsptrd chptrd zhptrd sppequ dppequ cppequ zppequ chpevd zhpevd sspev dspev sspevd dspevd chpev zhpev sspgv chpgvx zhpgvx dspgv sspgvd dspgvd chpgv zhpgv chpgvd zhpgvd chptrf zhptrf ssptrf dsptrf csptrf zsptrf stpcon dtpcon ctpcon ztpcon sppcon dppcon cppcon zppcon).each{|n| SUBSTS[n]  = {"n" => '(int)(sqrt((double)8*ldap+1)-1)/2'}}
%w(sspgvx dspgvx sspevx dspevx chpevx zhpevx).each{|n| SUBSTS[n] = {"n" => '(int)(sqrt((double)8*ldap+1)-1)/2', "m" => 'lsame_(&range,"A") ? n : lsame_(&range,"I") ? iu-il+1 : 0'}}
%w(sstevr dstevr ssyevr dsyevr dsyevx ssyevx).each{|n| SUBSTS[n] = {"m" => 'lsame_(&range,"I") ? iu-il+1 : n'}}
%w(ssygvx dsygvx cheevr zheevr chbevx zhbevx sstemr dstemr cstemr zstemr sstegr dstegr cstegr zstegr ssbevx dsbevx cheevx zheevx chegvx zhegvx ssbgvx dsbgvx).each{|n| SUBSTS[n] = {"m" => 'lsame_(&range,"A") ? n : lsame_(&range,"I") ? iu-il+1 : 0'}}
%w(stgex2 dtgex2).each{|n| SUBSTS[n] = {"lwork" => 'MAX(1,(MAX(n*(n2+n1),(n2+n1)*(n2+n1)*2)))'}}
%w(zlange zlanhs).each{|n| SUBSTS[n] = {"lwork" => 'lsame_(&norm,"I") ? n : 0'}}
%w(spprfs dpprfs cpprfs zpprfs chpsv zhpsv sspsv dspsv cspsv zspsv stprfs dtprfs ctprfs ztprfs).each{|n| SUBSTS[n] = {"n" => 'ldb'}}
SUBSTS["slasd0"] = {"m" => "sqre == 0 ? n : sqre == 1 ? n+1 : 0", "ldu" => "n", "ldvt" => "m"}
SUBSTS["dlasd0"] = {"m" => "sqre == 0 ? n : sqre == 1 ? n+1 : 0", "ldu" => "n", "ldu" => "n", "ldvt" => "n"}
SUBSTS["slasd3"] = {"ldu2" => "n", "n" => "nl + nr + 1", "ldvt2" => "n", "m" => "n+sqre"}
SUBSTS["dlasd3"] = {"ldu2" => "n", "n" => "nl + nr + 1", "ldvt2" => "n", "m" => "n + sqre"}
%w(slasda dlasda).each{|n| SUBSTS[n] = {"m" => "sqre == 0 ? n : sqre == 1 ? n+1 : 0", "ldu" => "n"}}
%w(slals0 dlals0 clals0 zlals0 slalsa dlalsa clalsa zlalsa).each{|n| SUBSTS[n] = {"ldbx" => "n"}}
%w(slasd2 dlasd2).each{|n| SUBSTS[n] = {"ldu2" => "n", "ldvt2" => "m"}}
%w(slasd6 dlasd6).each{|n| SUBSTS[n] = {"m" => "n + sqre", "n" => "nl + nr + 1"}}
SUBSTS["dgesdd"] = {"ldvt" => '(lsame_(&jobz,"A")||(lsame_(&jobz,"O")&&(m>=n))) ? n : lsame_(&jobz,"S") ? MIN(m,n) : 0'}
SUBSTS["zlaqr4"] = {"ldz" => 'wantz ? MAX(1,ihiz) : 1'}
%w(slaqr5 dlaqr5).each{|n| SUBSTS[n] = {"ldz" => 'n'}}
%w(sgelsd dgelsd).each{|n| SUBSTS[n] = {"c__0" => "0", "c__9" => "9", "smlsiz" => 'ilaenv_(&c__9,"'+n.upcase+'"," ",&c__0,&c__0,&c__0,&c__0)', "nlvl" => 'MAX(0,((int)(log(((double)(MIN(m,n)))/(smlsiz+1))/log(2.0))+1))', "liwork" => '3*(MIN(m,n))*nlvl+11*(MIN(m,n))'}}
%w(cgelsd zgelsd).each{|n| SUBSTS[n] = {"c__9" => "9", "c__0" => "0", "smlsiz" => 'ilaenv_(&c__9,"'+n.upcase+'"," ",&c__0,&c__0,&c__0,&c__0)', "nlvl" => 'MAX(0,(int)(log(1.0*MIN(m,n)/(smlsiz+1))/log(2.0)))', "lrwork" => 'm>=n ? 10*n+2*n*smlsiz+8*n*nlvl+3*smlsiz*nrhs+(smlsiz+1)*(smlsiz+1) : 10*m+2*m*smlsiz+8*m*nlvl+2*smlsiz*nrhs+(smlsiz+1)*(smlsiz+1)', "liwork" => 'MAX(1,3*(MIN(m,n))*nlvl+11*(MIN(m,n)))'}}
%w(sormbr dormbr cunmbr zunmbr).each{|n| SUBSTS[n] = {"nq" => 'lsame_(&side,"L") ? m : lsame_(&side,"R") ? n : 0'}}
%w(sbdsdc dbdsdc).each{|n| SUBSTS[n] = {"c__0" => "0", "c__9" => "9", "smlsiz" => 'ilaenv_(&c__9, "'+n.upcase+'", " ", &c__0, &c__0, &c__0, &c__0)', "ldq" => 'lsame_(&compq,"P") ? n*(11+2*smlsiz+8*(int)(log(((double)n)/(smlsiz+1))/log(2.0))) : 0', "ldiq" => 'lsame_(&compq,"P") ? n*(3+3*(int)(log(((double)n)/(smlsiz+1))/log(2.0))) : 0', "lwork" => 'lsame_(&compq,"N") ? 4*n : lsame_(&compq,"P") ? 6*n : lsame_(&compq,"I") ? 3*n*n+4*n : 0'}}
%w(clalsd zlalsd).each{|n| SUBSTS[n] = {"nlvl" => '( (int)( log(((double)n)/(smlsiz+1))/log(2.0) ) ) + 1'}}
%w(claed8 zlaed8).each{|n| SUBSTS[n] = {"ldq2" => "n"}}





def pow(str)
  str = str.gsub(/([A-Z\d]+?)\*\*([A-Z\d]+)/, 'pow(\\1,\\2)')
  str = str.gsub(/\(([^\)]+)\)\*\*([A-Z\d]+)/, 'pow(\\1,\\2)')
  str.gsub!(/lg (\w+)/, 'LG(\\1)')
  str.gsub!(/(\w) (LG\()/, '\\1*\\2')
  return str
end


def get_vname(str)
  v = pow(str).strip.downcase.gsub(/lg\(/,"LG(").gsub(/([^a-z])max/,'\\1MAX').gsub(/^max/,'MAX').gsub(/min/,"MIN")
  if /if / =~ v || (/\,/ =~ v && /(MIN|MAX|pow)\(/ !~ v) || / the / =~ v
    raise "vname is invalid #{v}"
  end
  v
end


def get_dims(str)
  if /^\((.*)\)\.?$/ =~ str
    str = $1.strip
  end
  if /(max|min)/i =~ str
    dims = Array.new
    while (str && str != "")
      if /^((((MAX|MIN)\s*\([^\,]+\,)|[^\,])+)\,?(.*)$/i =~ str
        dns = $1.strip
        str = $5.strip
        dims.push dns
        next
      else
        /^([^\,]+)(\,.*)$/ =~ str
        dims.push $1
        str = $2.strip
      end
      if /^\,(.*)/ =~ str
        str = $1.strip
      end
    end
  else
    dims = str.split(",").collect{|dim| dim.sub(/\.$/,"")}
  end
  dims.collect do |dim|
    dim.sub!(/;\Z/,"")
    if /\A\((.+)\)\Z/ =~ dim
      dim = $1
    end
    get_vname(dim)
  end
end

AO = {"and"=>"&&","or"=>"||"}
def get_cond(cond,v=nil)
  cond.sub!(/\.$/,"")
  if /^(.+?)\s*(and|or)\s*(.+)$/ =~ cond
    c0 = $1.strip
    ao = $2
    c1 = $3.strip
    if /^([A-Z\d]+)\s*=/ =~ c0
      v = get_vname($1)
    else
      v = nil
    end
    cond = "((#{get_cond(c0,v)}) #{AO[ao]} (#{get_cond(c1,v)}))"
  elsif /^(.+?)\s*=\s*\'(.+)'$/ =~ cond
    cond = "lsame_(&#{$1.downcase},\"#{$2}\")"
  else
    if /=.*=/ =~ cond
      conds = cond.split("=").collect{|c| c.strip.downcase}
      cond = Array.new
      (conds.length-1).times{|i|
        cond.push "(#{conds[i]}==#{conds[i+1]})"
      }
      cond = cond.join("&&")
    elsif v && (/^\'(.+)\'$/ =~ cond)
      cond = "lsame_(&#{v},\"#{$1}\")"
    elsif v && (/^\d+$/ =~ cond)
      cond = "#{v} == #{cond}"
    else
      cond = cond.gsub(/=/,"==").downcase
    end
  end
  cond
end



def read_file(fname)
  flag_sub = false
  subr = nil
  sub_type = nil
  flag_pur = false
  purpose = nil
  flag_arg = false
  args = nil
  flag_fd = false
  fd = nil
  File.foreach(fname){|line|

    if /This routine is not for general use\./ =~ line
      return nil
    end
    if /^\*\s+\.\. Parameters \.\./ =~ line || /^\*\s+\.\. Executable Statements \.\./ =~ line
      break
    end

    if flag_sub
      if (/^     \$\s* (.+)$/ =~ line) || (/^     \+\s* (.+)$/ =~ line) || (/^     \&\s* (.+)$/ =~ line)
        subr << " " << $1.chomp
      else
        flag_sub = false
      end
      next
    elsif /^      SUBROUTINE/ =~ line
      subr = line.chomp
      flag_sub = true
      next
    elsif /^      RECURSIVE SUBROUTINE/ =~ line
      subr = line.chomp
      flag_sub = true
      next
    elsif /^      ([A-Z\s\d\*]+)\s+FUNCTION/ =~ line
      subr = line.chomp
      flag_sub = true
    elsif /^      FUNCTION/ =~ line && /^[sd]laneg/ =~ File.basename(fname)
      subr = line.chomp
      flag_sub = true
    end

    if flag_pur
      if /^\*\s+Arguments$/ =~ line
        flag_pur = false
        args = line
        flag_arg = true
        next
      else
        case File.basename(fname)
        when /^[cz]la_lin_berr/
          if /^\*     N       \(input\) INTEGER$/ =~ line
            flag_pur = false
            args = line
            flag_arg = true
            next
          end
        when /^[sdcz]laqr1/
          if /^\*       N      \(input\) integer$/ =~ line
            flag_pur = false
            args = line
            flag_arg = true
            next
          end
        when /^[sdcz]laqr2/, /^[sdcz]laqr3/
          if /^\*     WANTT   \(input\) LOGICAL$/ =~ line
            flag_pur = false
            args = line
            flag_arg = true
            next
          end
        when /^[sdcz]laqr5/
          if /^\*      WANTT  \(input\) logical scalar$/ =~ line
            flag_pur = false
            args = line
            flag_arg = true
            next
          end
        when /^[sd]lasq4/, /^[sd]lazq4/
          if /^\*  I0    \(input\) INTEGER$/ =~ line
            flag_pur = false
            args = line
            flag_arg = true
            next
          end
        end
        purpose << line
        next
      end
    elsif /^\*\s+Purpose$/ =~ line
      purpose = line
      flag_pur = true
      next
    else
      case File.basename(fname)
      when /^[sdcz]laqr1/
        if /^\*       Given a 2-by-2 or 3-by-3 matrix H, [SDCZ]LAQR1 sets v to a$/ =~ line
          purpose = line
          flag_pur = true
        end
      when /^[sdcz]laqr2/
        if /^\*     This subroutine is identical to [SDCZ]LAQR3 except that it avoids$/ =~ line
          purpose = line
          flag_pur = true
        end
      when /^[sdcz]laqr3/
        if /^\*     Aggressive early deflation:$/ =~ line
          purpose = line
          flag_pur = true
        end
      when /^[sdcz]laqr5/
        if /^\*     This auxiliary subroutine called by [SDCZ]LAQR0 performs a$/ =~ line
          purpose = line
          flag_pur = true
        end
      end
    end

    if flag_arg
      if /^\*\s+Further Details$/ =~ line || /^\*\s+={40,}$/ =~ line
        flag_arg = false
        fd = line
        flag_fd = true
        next
      else
        args << line
        next
      end
    end
    
    if flag_fd
      fd << line
      next
    end
  }

  unless subr
    raise "subr not found #{fname}"
  end
  unless purpose
    raise "purpose not found #{fname}"
  end
  unless args
    raise "args not found #{fname}"
  end

  help = subr + "\n\n" + purpose + "\n" + args
  help << "\n" + fd if fd

  return {:subr => subr, :purpose => purpose, :args => args, :help => help}
end

def parse_file(fname)
  hash = read_file(fname)
  subr = hash[:subr]
  purpose = hash[:purpose]
  args = hash[:args]
  help = hash[:help]

  if /^      (?:RECURSIVE )?SUBROUTINE\s+([A-Z\d_]+)\s*\(([^\)]+)\)/ =~ subr
    sub_name = $1.downcase
    arg_names = $2
    sub_type = :subroutine
  elsif /^      ([A-Z\s\*\d]+[A-Z\d])\s+FUNCTION\s+([A-Z\d_]+)\s*\(([^\)]+)\)/ =~ subr
    f_type = $1.strip
    sub_name = $2.downcase
    arg_names = $3
    sub_type = :function
    if f_type == "CHARACTER*1"
      f_type = "CHARACTER"
    end
    func_type = CTYPES[f_type]
    unless func_type
      raise "func_type #{f_type} is not defined"
    end
  elsif /^      FUNCTION\s+([A-Z\d]+)\(([^\)]+)\)/ =~ subr
    sub_name = $1.downcase
    arg_names = $2
    sub_type = :function
    case File.basename(fname)
    when /^[sd]laneg/
      func_type = "integer"
    else
      raise "function name is invalid #{subr}"
    end
  else
    raise "subroutine or function name is invalid #{subr}"
  end
  arg_names = arg_names.split(",").collect{|arg| arg.strip.downcase}
  unless arg_names
    raise "arg_names is nil #{fname}"
  end

  if @@debug
    p sub_name
    p arg_names
  end

  flag = false
  flag1 = false
  ary = Array.new
  args.each{|line|
    case line
    when /^\*\s*$/, /^\*  Arguments$/, /^\*  =========$/
      flag = false
      flag1 = false
      next
    end

    if /^\*\s+([A-Z_\d,\s]+)\s+\((input|output|workspace|in|external procedure)[^\)]*\)\s+([A-Za-z]+)/ =~ line
      name = $1
      intent = $2
      type = $3
      name.strip!
      ary.push line.chomp
      if /array/i =~ line
        flag = true
      elsif (/^LD/=~name || /L.*WORK/=~name) && /input/=~intent && /INTEGER/i=~type
        flag1 = true
      end
    elsif flag
      line.sub!(/^\*\s*/,"")
      if /[Ii]f / =~ line || /where / =~ line || /^[a-z]/ =~ line || /[a-z]$/ =~ ary[-1] || /is not referenced/ =~ line || /dimension (of [A-Z\d]+ must be )?at least/ =~ ary[-1] || /^Otherwise, / =~ line
        if /^[^\s]+ is the/ =~ line || /^If [^\.]+ must contain/ =~ line || /where in / =~ line || /At each / =~ line || /^[^\.]+\, [A-Z\d]+ contains / =~ line || /[A-Z]+\(\d\) and [A-Z]+\(\d\) contain the/ =~ line
          flag = false
          next
        end
        if  /^\( lg\(\s*\w?\s*\) = smallest integer \w$/ =~ line || /^such that 2\^\w \>= \w \)$/ =~ line
          next
        end
        while /^(.*)On exit/ =~ line || /^(.*?)[^\.]+ on exit/ =~ line || /^(.*)On entry/ =~ line || /^(.*)The/ =~ line || /^(.*?)[^\.]+ the [^d]/=~ line || /^(.*?)[^\.]+ is an /=~ line
          line = $1.strip
          flag = false
        end
        if line != ""
          if /^[A-Z]/ =~ line && /[\.\,;]$/ !~ ary[-1] && /if$/ !~ ary[-1]
            ary[-1] << ";"
          end
          if /\)$/ =~ ary[-1] && /^\(/ =~ line
            ary[-1] << ";"
          end
          ary[-1] << " " << line.chomp
        end
      elsif /^(?:The d|D)imension must be at least (.*)$/ =~ line
        dim = $1
        ary[-1].sub!(/\.$/," ")
        ary[-1] << "(#{dim})"
        flag = false
      else
        flag = false
      end
    elsif flag1
      line.sub!(/^\*\s*/,"")
      ary[-1] << " " << line.chomp
    else
      if /^\*\s+If L.?WORK = -1, .*a workspace query/ =~ line
        line.sub!(/^\*\s*/,"")
        ary[-1] << " " << line.chomp
        flag1 = true
      end
    end
  }

  if @debug
    pp ary
  end

  args = Hash.new
  subst = SUBSTS[sub_name] || Hash.new
  ary.each{|line|
    line.strip!
    /^\*\s+([A-Z\d_,\s]+)\s+\(([^\)]+)\)\s*(.*)$/ =~ line
    name = $1.downcase.strip
    intent = $2
    type = $3.sub(/\.$/,"")
    hash =  Hash.new

    intent = "input" if intent == "in"
    intent = "input or input/output" if intent == "input or input/ouptut"
    case sub_name
    when /^[cz]larcm$/
      case name
      when "c"
        intent = "output"
      end
    when /^[cz]lacrm$/
      case name
      when "c"
        intent = "output"
      end
    end
    hash[:intent] =  intent
    if /^l.*work/ =~ name && (/The (dimension|length)? of (the )?(array|work|WORK)/ =~ type || /The amount of workspace available/ =~ type)
      if (/^(.*?) The (?:dimension|length)? of (?:the )?array\s+(?:WORK\.\s*)?(.+)/ =~ type) || (/^(.*?) The (?:dimension|length)? of (?:the )?(?:array|work|WORK)\.?\s+(.+)/ =~ type)
        type = $1.strip
        str = $2
      elsif /^(.+?) The amount of workspace available/ =~ type
        type = $1.strip
        str = nil
      else
        raise "invalid:  #{type}, #{name}, #{sub_name}"
      end
      aname = name.sub(/^l/,"")
      fff = false
      if (/If #{name.upcase} = -1, (then )?.*a workspace query/ =~ str) || (/If #{name.upcase} = -1 .+? returns the optimal/ =~ str)
        fff = true
      elsif /^([cz]tgsna|[cz]hetrf)$/=~sub_name && name=="lwork"
        fff = true
      elsif str
        unless /^[sd]tgex2$/ =~ sub_name && name == "lwork"
          warn "invalid #{str} #{name} #{sub_name}"
        end
      end
      if fff
        unless ar = args[aname]
          raise "array not found #{aname} #{name} #{sub_name}"
        end
        d = ar[:dims]
        unless d.length == 1
          raise "array length is not 1 #{aname} #{name} #{sub_name}"
        end
        if /^[a-z\d_]+$/ =~ d[0]
          d[0] = "MAX(1,#{d[0]})"
        end
      end
    elsif /^ld/ =~ name && (/(leading|first) dimension of/ =~ type || /(Leading|First) dimension of/ =~ type)
      if /^([^,]+) T[Hh]e (?:leading|first) dimension of\s+(?:(?:the )?vector)?(?:the array(?:s)?)?\s*(.+)$/ =~ type || /^(.+?) LD[A-Z]+ is the leading dimension of\s+(.+)$/ =~ type || /^(.+?) On entry, (?:LD[A-Z]+\s+specifies the )?(?:leading|first) dimension of (.+)$/ =~ type || /^([^,]+) (?:[Ll]eading|First) dimension of\s+(.+)$/ =~ type
        type = $1.strip
        str = $2.strip
      elsif /^([^,]+)\,\s+(.+?)\s+The (?:leading|first) dimension of\s+(.+)$/ =~ type
        type = $1.strip
        str0 = $2.strip
        str1 = $3.strip
        str = str1 + ", " + str0
      else
        raise "invalid #{type}, #{name} #{sub_name}"
      end
      if /^(?:the )?arrays (.+?),? and ([A-Z\d_]+)[\,\.]?\s*(.*)$/ =~ str
        anames = $1.strip
        aname1 = $2
        str = $3.strip
        anames = anames.split(",")
        anames.push aname1
        anames.collect!{|an| an.downcase}
      elsif /^([A-Z\d_]+) and ([A-Z\d_]+)\, must be at least (.*)$/ =~ str
        anames = [$1.downcase, $2.downcase]
        str = $3.strip
      elsif /^(?:the )?(?:array |matrix )?([A-Z\d_]+)\.?\,?\s*(.*)$/ =~ str
        anames = [$1.strip.downcase]
        str = $2.strip
        case sub_name
        when /^[cz]laqr[23]$/, /^[sd]laqr3$/, /^[sd]laqr2$/
          case name
          when "ldwv"
            anames = ["wv"]
          end
        when /^[sdcz]pbsvx$/, /^[sdcz]gbbrd$/, /^[sdcz]pbequ$/
          case name
          when "ldab"
            anames = ["ab"]
          end
        end
      else
        raise "error #{str} #{name} #{sub_name}"
      end
      ff0 = false
      ff1 = true
      anames.each{|an|
        if (ar = args[an])
          ff0 = true
          if /input/ =~ ar[:intent]
            ff1 = false
            break
          end
        end
      }
      unless ff0
        warn "arg not found [#{anames.join(",")}], #{name}, #{sub_name}"
      end
      if ff1 && str != ""
        if anames.length == 1
          aname = anames[0]
        elsif anames.include?(name.sub(/^ld/,""))
          aname = name.sub(/^ld/,"")
        else
          case sub_name
          when /^[sd]lasda$/
            case name
            when "ldgcol"
              aname = "givcol"
            end
          when /^[sd]lasd6$/
            case name
            when "ldgnum"
              aname = "givnum"
            end
          end
        end
        unless aname
          raise "cannot select anames [#{anames.join(",")}], #{name}, #{sub_name}"
        end
        un = name.upcase
        str.sub!(/\.$/,"")
        str.gsub!(/\.GE\./,"=")
        str.gsub!(/>=/,"=")
        str.gsub!(/=\s*>/,"=")
        str.gsub!(/(It )?must be at least/, "")
        str.gsub!(/#{un}\s*=/," ")
        str.sub!(/(just )?as declared in the (in the )?calling (subroutine|procedure)\./,"")
        str.strip
        if  /^([^;]+); if ([^,]+), ([^;]+); if ([^,]+), (.+)$/ =~ str
          if @@debug
            p 100
          end
          v = get_vname($1)
          cond0 = get_cond($2)
          v0 = get_vname($3)
          cond1 = get_cond($4)
          v1 = get_vname($5)
          str = "#{cond0} ? #{v0} : #{cond1} ? #{v1} : #{v}"
        elsif /^If ([^\,]+)\,\s+([^;]+); if ([^\,]+)\,\s+(.*)$/ =~ str || /^If ([^\,]+)\,\s+([^\.]+)\. If ([^\,]+)\,\s+(.*)$/ =~ str
          if @@debug
            p 110
          end
          v0 = get_vname($1)
          cond0 = get_cond($2)
          v1 = get_vname($3)
          cond1 = get_cond($4)
          str = "#{cond0} ? #{v0} : #{cond1} ? #{v1} : 0"
        elsif /^([^\,]+)\, and if ([^\,]+)\, (?:then )?(.*)$/ =~ str || /^([^;]+); if ([^\,]+)\, (?:then )?(.*)$/ =~ str
          if @@debug
            p 120
          end
          v0 = get_vname($1)
          cond1 = get_cond($2)
          v1 = get_vname($3)
          str = "#{cond1} ? #{v1} : #{v0}"
        elsif /^(.+?) if ([^;]+); (.+?) otherwise$/ =~ str
          if @@debug
            p 130
          end
          v0 = get_vname($1)
          cond0 = get_cond($2)
          v = get_vname($3)
          str = "#{cond0} ? #{v0} : #{v}"
        elsif /^If ([^,]+), then\s+(.+?)\.\s+In any case, (.+)$/ =~ str
          if @@debug
            p 140
          end
          cond0 = get_cond($1)
          v0 = get_vname($2)
          v = get_vname($3)
          str = "#{cond0} ? #{v0} : #{v}"
        elsif  /^([^;]+); and if ([^,]+), (.+)$/ =~ str
          if @@debug
            p 150
          end
          v = get_vname($1)
          cond0 = get_cond($2)
          v0 = get_vname($3)
          str = "#{cond0} ? #{v0} : #{v}"
        elsif /^(.+?) \.LE\. #{un}$/ =~ str
          if @@debug
            p 160
          end
          str = get_vname($1)
        elsif /^If ([^,]+), then\s+(.+)$/ =~ str
          if @@debug
            p 170
          end
          cond0 = get_cond($1)
          v0 = get_vname($2)
          str = "#{cond0} ? #{v0} : 0"
        elsif (/^[sd]bdsdc$/ =~ sub_name && /^ld(u|vt)$/ =~ name)
          if @@debug
            p 180
          end
          str = 'lsame_(&compq,"I") ? MAX(1,n) : 0'
        else
          if @@debug
            p 190
          end
          if /^[sdcz]laqr[23]$/ =~ sub_name && name == "ldwv"
            str = "nw"
          end
          begin
            str = get_vname(str)
          rescue
            warn "error #{str}, #{name}, #{sub_name}"
          end
        end
        if /^[sdcz]larrv$/ =~ sub_name && name == "ldz"
          str = "n"
        end
        subst[name] = str
      end
    elsif /^(.*?) array of size (.*)$/ =~ type || /^(.*?) arrays?,?(.*)$/i =~ type
      type = $1.strip
      str = $2.strip
      if /^([A-Z\s]+) work$/ =~ type
        type = $1.strip
      end
      if "CHARACTER(1)" == type
        type = "CHARACTER"
      end
      if /\Alength (.*)\z/ =~ str
        str = $1
      end
      d = DIMS[sub_name]
      dims = d[name] if d

      unless dims
        str.gsub!(/(#{name}) must be at least\s+/,'\\1 = ')
        str.gsub!(/must be at least\s+/,"")
        str.gsub!(/at least\s*/,"")
        str.gsub!(/dimension ([^\(]+) if/, '(\\1) if')
        str.gsub!(/dimension \>= /,"")
        str.gsub!(/the dimension of #{name.upcase}/,"")
        str.gsub!(/.?\s+dimension is/,"")
        str.gsub!(/(of )?dimensions?/,"")
        str.gsub!(/\.GE\.?/,">=")
        str.gsub!(/\.GT\.?/,">")
        str.gsub!(/\.LE\.?/,"<=")
        str.gsub!(/\.LT\.?/,"<")
        str.gsub!(/\.TRUE\.?/,"TRUE")
        str.gsub!(/\.FALSE\.?/,"FALSE")
        str.gsub!(/log_2\s*\(/,"1.0/log(2.0)*log((double)")
        str.gsub!(/INT\s*\(/,"(int)\(")

        str.strip!
        if intent == "input"
          str.sub!(/\s*if .*$/i, "")
        end
        if /(...rfsx|...svxx)/ =~ sub_name && name == "params"
          str = ""
        end
        if /[Ii]f/ =~ str || /where/ =~ str || /when/ =~ str
          dims = Array.new
          flag = true
          str.sub!(/;$/,"")
          if /^\([^;]+\); (\(.+?\) if [A-Z\d]+ = .+? or \(.+?\) if [A-Z\d]+ = .+)$/ =~ str
            str = $1
          end
          if @@debug
            p name
            p str
          end
          while (str && str!="")
            if (/^(?:and )?\((.*?)\)\s+if\s+([^\,\;]+)[\,\;]\s+(.*)$/ =~ str) && (dim = $1.strip) && (cond = $2.strip) && (str_tmp = $3.strip) && (/ if .* if / !~ cond) && (/\([^\)]+$/ !~ cond) && (/=/ !~ dim)
              if @@debug
                p 1
              end
              str = str_tmp
            elsif (/^\(([^;]+?)\)\s+if\s+(.+?)\s+(?:or|and)\s+(\(.*\) if.*)$/ =~ str) || (/^\(([^;]+?)\)\s+if\s+(.+?)\s+(\(.*\)\s+if.*)$/ =~ str)
              if @@debug
                p 2
              end
              dim = $1
              cond = $2.strip
              str = $3.strip
            elsif /^If ([^\,]+)\,\s+([^;\.]+)[;\.]?(.*)$/ =~ str
              if @@debug
                p 2.5
              end
              cond = $1.strip
              dim = $2
              str = $3.strip
            elsif /^\((.*?)\)\s+if\s+([^\,]+)\, (\(.*\) otherwise)$/ =~ str
              if @@debug
                p 3
              end
              dim = $1.strip
              cond = $2.strip
              str = $3.strip
            elsif (/^(?:and )?\((.+?)\)\s+if\s+(.*)$/ =~ str) && (dim = $1.strip) && (cond = $2.strip) && (/if/ !~ dim) && (/ if /i !~ cond) && (/^[^\(]+\)/ !~ dim)
              if @@debug
                p 4
              end
              str = nil
            elsif /\A\((.*)\) when (.*) and/ =~ str
              dim = $1
              cond = $2
              str = $'
            elsif /^(.*?)\s+otherwise$/ =~ str && (dim = $1.strip) && (/ if / !~ dim)
              if @@debug
                p 5
              end
              if cond == "not referenced"
                dims = dims.collect{|d| d + "0"}
              else
                get_dims(dim).each_with_index{|d,i|
                  dims[i] = dims[i] + d
                }
              end
              flag = false
              break
            elsif /^\((.+?)\);?\s+If ([A-Z_\d]+ = [^\,]+)\, (then )?#{name.upcase} is not referenced.?$/ =~ str
              if @@debug
                p 6
              end
              dims = get_dims($1)
              cond = get_cond($2)
              dims = dims.collect{|dim|
                "#{cond} ? 0 : #{dim}"
              }
              flag = false
              break
            elsif /^\((.+?)\); Not referenced if (.*)$/ =~ str
              if @@debug
                p 6.5
              end
              v0 = get_vname($1)
              cond0 = get_cond($2)
              dims = ["#{cond0} ? 0 : #{v0}"]
              flag = false
              break
            elsif /\((.+?)\)\, where ([A-Z\d]+) >= ([A-Z\d]+) when ([^;]+); otherwise, #{name.upcase} is not referenced\.?$/ =~ str
              if @@debug
                p 7
              end
              dims = get_dims($1)
              c0 = $2
              c1 = $3
              cond = get_cond($4)
              subst[get_vname(c0)] = "#{cond} ? #{get_vname(c1)} : 0"
              flag = false
              break
            elsif /\((.+)\)[\.\,] where ([A-Z\d]+) = (.+)\.?$/ =~ str
              if @@debug
                p 8
              end
              dims = get_dims($1)
              c0 = $2
              c1 = $3
              subst[get_vname(c0)] = get_vname(c1)
              flag = false
              break
            elsif /^(?:and )?(not referenced) if (.+).?$/ =~ str
              if @@debug
                p 9
              end
              dim = $1
              cond = $2
              str = nil
            elsif (/^\((.+?)\) where ([A-Z\d_]+) = (.*)$/ =~ str) && (dim = $1) && (c0 = $2) && (c1 = $3) && /=/ !~ c1
              dims = get_dims(dim)
              subst[get_vname(c0)] = get_vname(c1)
              flag = false
              break
            else
              if /^\((.+?)\) (.+?) = (.+?) (?:when|if) ([^,]+), and (.+?) when (.+)$/ =~ str
                if @@debug
                  p 10
                end
                dim = $1
                dim1 = $2
                unless dim == dim1
                  raise "error #{name} #{sub_name}"
                end
                c0 = $3
                cond0 = get_cond($4)
                c1 = $5
                cond1 = get_cond($6)
                dims = get_dims($1)
                subst[get_vname(dim)] = "#{cond0} ? #{get_vname(c0)} : #{cond1} ? #{get_vname(c1)} : 0"
                flag = false
                break
              elsif /^\((.+?)\) (.+?) = (.+?) (?:when|if) ([^,]+), and (.+?) otherwise.?$/ =~ str
                if @@debug
                  p 11
                end
                dim = $1
                dim1 = $2
                unless dim == dim1
                  raise "error #{name} #{sub_name}"
                end
                c0 = $3
                cond0 = get_cond($4)
                c1 = $5
                dims = get_dims(dim)
                subst[get_vname(dim)] = "#{cond0} ? #{get_vname(c0)} : #{get_vname(c1)}"
                flag = false
                break
              elsif /^\((.+?)\);\s+[Ii]f ([A-Z_\d]+ = [^\,]+)\, ([A-Z_\d]+) \>?= ([^\.]+)\.\s+Otherwise\, ([A-Z_\d]+) \>?= (.*)$/ =~ str
                if @@debug
                  p 14
                end
                dim = $1
                cond0 = $2
                c00 = $3.downcase
                c01 = $4.downcase
                c10 = $5.downcase
                c11 = $6.downcase
                unless c00==c10
                  raise "error #{name} #{sub_name}"
                end
                dims = get_dims(dim)
                cond0 = get_cond(cond0)
                subst[get_vname(c00)] = "#{cond0} ? #{get_vname(c01)} : #{get_vname(c11)}"
                flag = false
                break
              end
              if /^\((.+?)\)\.\s+[Ii]f ([A-Z_\d]+ = [^\,]+)\, ([A-Z_\d]+) = ([A-Z_\d]+); if ([A-Z_\d]+ = [^\,]+), ([A-Z_\d]+) = ([A-Z_\d]+)$/ =~ str
                if @@debug
                  p 13
                end
                dim = $1
                cond0 = $2
                c00 = $3.downcase
                c01 = $4.downcase
                cond1 = $5
                c10 = $6.downcase
                c11 = $7.downcase
              elsif /\((.+?)\);? ([A-Z_\d]+) = (.+?) if ([^;]+); ([A-Z_\d]+) = (.+?) if (.+)$/ =~ str
                if @@debug
                  p 13
                end
                dim = $1
                c00 = $2.downcase
                c01 = $3.downcase
                cond0 = $4
                c10 = $5.downcase
                c11 = $6.downcase
                cond1 = $7
              else
#                raise "error '#{str}' in #{name} (#{fname})"
                 warn "'#{str}' in #{name} (#{fname})"
                 dim = "dummy_" + str
              end
              if /dummy_/ =~ dim
                dims = dim
              else
                dims = get_dims(dim)
                cond0 = get_cond(cond0)
                cond1 = get_cond(cond1)
                if c00 == c10
                  subst[get_vname(c00)] = "#{cond0} ? #{get_vname(c01)} : #{cond1} ? #{get_vname(c11)} : 0"
                else
                  raise "error #{name} #{sub_name}"
                end
              end
              flag = false
              break
            end
            if flag
              cond = get_cond(cond)
              if dim == "not referenced"
                dims = dims.collect{|d| d + cond + " ? 0 : "}
              else
                ds = get_dims(dim)
                ds.each_with_index{|d,i|
                  dims[i] ||= ""
                  dims[i] << cond + " ? " + d + " : "
                }
              end
            end
          end
          dims = dims.collect{|d| d + "0"} if flag
        else
          dims = get_dims(str)
        end

        dims.each_with_index{|dim,i|
          if /\?/ =~ dim
            str = dim
            aa = Array.new
            while( /^[^\?]+ \? ([^:]+) : (.+)$/ =~ str )
              aa.push $1
              str = $2
            end
            aa.push str
            if aa.length > 2
              aa = aa.uniq
              if aa.length == 1 || (aa.length == 2 && aa[-1] == "0")
                dims[i] = aa[0]
              end
            elsif aa.length == 2
              if aa.uniq.length == 1
                dims[i] = aa[0]
              end
            else
              raise "error [#{aa.join(",")}] #{dim} #{name} #{sub_name}"
            end
          end
        }
      end
      if /^[sd]lasda$/ =~ sub_name  && name == "difl"
        subst["nlvl"].sub!(/\)$/,"")
      end
      hash[:dims] = dims
    elsif (/^(.+) work array$/ =~ type) || (/^(.+) array$/ =~ type)
      type = $1.strip
      dims = DIMS[sub_name] && DIMS[sub_name][name]
      unless dims
        raise "dimension is not defined (#{name} in #{fname})"
      end
    elsif /^CHARACTER\*\s*(\d+)$/ =~ type || /^CHARACTER\*?\((\d+)\)$/ =~ type || /^character string/ =~ type
      type = "CHARACTER"
      if $1 && $1 != "1"
        hash[:dims] = [$1]
      end
    elsif /^CHARACTER\*\(\*\)$/ =~ type
      type = "CHARACTER"
      hash[:dims] = ["*"]
    elsif /^(.*)\,\s*#{name}\s*=\s*>/i =~ type
      type = $1.strip
    elsif /^(.+?)\s+scalar$/ =~ type
      type = $1
    elsif /^(.+?)\s+with value/ =~ type
      type = $1
    end
    if (t = TYPES[sub_name]) && (t = t[name])
      hash[:type] = t
    else
      type.sub!(/ scalar/,"")
      if /^(.+?) FUNCTION of ([^\s]+) (.+?) arguments?$/ =~ type
        hash[:block_type] = CTYPES[$1] || raise("error block type is invalid #{$1} #{name} #{sub_name}")
        hash[:block_arg_num] = {"one"=>1,"two"=>2,"three"=>3}[$2] || raise("error #{$2}")
        hash[:block_arg_type] = CTYPES[$3] || raise("error block arg type is invalid #{$3} #{name} #{sub_name}")
      else
        hash[:type] = CTYPES[type.upcase]
        unless hash[:type]
          warn("type (#{type}) is not defined in #{name} (#{fname})")
        end
      end
    end
    if /,/ =~ name
      name.split(",").each{|n1|
        args[n1.strip] = hash.dup
      }
    else
      args[name] = hash
    end
  }
  case sub_name
  when /^[cz]laqr[04]$/
    args["iloz"] = {:type => "integer", :intent => "input"}
    args["ihiz"] = {:type => "integer", :intent => "input"}
  end
  if @@debug
    pp args
  end

  return  {:name => sub_name, :category => sub_type, :type => func_type, :argument_names => arg_names, :arguments => args, :fortran_help => help, :md5sum => get_md5sum(help), :substitutions => subst}


end


def create_hash(fname)
  hash = parse_file(fname)
  if hash.nil?
    warn "skip #{fname}"
    return nil
  end
  sub_name = hash[:name]
  arg_names = hash.delete(:argument_names)
  args = hash.delete(:arguments)

  args.each do |k,v|
    case v[:intent]
    when "input or output", "input or input/output", "input / output", "input/workspace/output", "input/workspace"
      v[:intent] = "input/output"
    when "workspace/output"
      v[:intent] = "output"
    end
  end


  unless sub_name
    warn "this has no subroutine (#{fname})"
    return nil
  end
  unless arg_names
    raise "no arg_names (#{fname})"
  end
  unless args && !args.empty?
    raise "no args (#{fname})"
  end

  hash[:arguments] = arg_names.map{|name| {name => args[name]} }

  return hash
end

def generate_code(fnames)
  nfnames = fnames.length
  sub_names = Array.new
  fnames.each_with_index{|fname,i|
    print "#{i+1}/#{nfnames}\n" if (i+1)%100==0
    hash = read_file(fname)
    next if hash.nil?
    help = hash[:help]
    basename = File.basename(fname, ".f")
    def_dir = File.join(File.dirname(__FILE__), "defs")
    def_fname = File.join(def_dir, basename)
    if File.exists?(def_fname)
      hash = nil
      File.open(def_fname){|file| hash = YAML.load(file.read) }
      md5sum = hash[:md5sum]
      next if get_md5sum(help) == md5sum && !@@force
    end

    hash = create_hash(fname)
    def hash.each
      [:name, :md5sum, :category].each do |k|
        yield(k, self[k])
      end
      if self[:category] == :function
        yield(:type, self[k])
      end
      [:arguments, :substitutions, :fortran_help].each do |k|
        yield(k, self[k])
      end
    end
    p "write #{hash[:name]}" if @@debug
    File.open(def_fname, "w") do |file|
      file.write hash.to_yaml
    end
  }
end

def get_md5sum(src)
  Digest::MD5.hexdigest(src.sub(/LAPACK routine \(version \d.\d\) --/,""))
end





@@debug = ARGV.delete("--debug")
@@force = ARGV.delete("--force")

dname = ARGV[0] || raise("Usage: ruby #$0 path_to_lapack_src")
if File.directory?(dname)
  fnames = Dir[ File.join(dname,"*.f") ]
elsif File.file?(dname)
  fnames = [dname]
#  @@debug = true
end

generate_code(fnames)



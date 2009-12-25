DataTypes = [
             ["S", "REAL"],
             ["D", "DOUBLE PRECISION"],
             ["C", "COMPLEX"],
             ["Z", "COMPLEX*16 or DOUBLE COMPLEX"],
             ["DS", "Data type in double but solving problem using single precision"],
             ["ZC", "Data type in complex*16 but solving problem using complex precision"]
]

MatrixTypes = [
              ["BD", "bidiagonal"],
              ["DI", "diagonal"],
              ["GB", "general band"],
              ["GE", "general (i.e., unsymmetric, in some cases rectangular)"],
              ["GG", "general matrices, generalized problem (i.e., a pair of general matrices)"],
              ["GT", "general tridiagonal"],
              ["HB", "(complex) Hermitian band"],
              ["HE", "(complex) Hermitian"],
              ["HG", "upper Hessenberg matrix, generalized problem (i.e a Hessenberg and a triangular matrix)"],
              ["HP", "(complex) Hermitian, packed storage"],
              ["HS", "upper Hessenberg"],
              ["OP", "(real) orthogonal, packed storage"],
              ["OR", "(real) orthogonal"],
              ["PB", "symmetric or Hermitian positive definite band"],
              ["PO", "symmetric or Hermitian positive definite"],
              ["PP", "symmetric or Hermitian positive definite, packed storage"],
              ["PT", "symmetric or Hermitian positive definite tridiagonal"],
              ["SB", "(real) symmetric band"],
              ["SP", "symmetric, packed storage"],
              ["ST", "(real) symmetric tridiagonal"],
              ["SY", "symmetric"],
              ["TB", "triangular band"],
              ["TG", "triangular matrices, generalized problem (i.e., a pair of triangular matrices)"],
              ["TP", "triangular, packed storage"],
              ["TR", "triangular (or in some cases quasi-triangular)"],
              ["TZ", "trapezoidal"],
              ["UN", "(complex) unitary"],
              ["UP", "(complex) unitary, packed storageBDbidiagonal"]
]


def parse_html(fname)
  hash = Hash.new
  name = nil
  File.foreach(fname){|line|
    if /^file <a href=".+">([a-z_\d]+)\.f<\/a>/ =~ line
      name = $1
    elsif name
      if /^for\s+(.*)$/ =~ line
        hash[name] = $1
      elsif /^,\s+(.*)$/ =~ line
        hash[name] ||= ""
        hash[name] << $1
      elsif /^gams/ =~ line 
        name = nil
      end
    end
  }
  return hash
end

require "numru/lapack"
include NumRu

prefix = "../doc"

path = ARGV[0] || raise("Usage: ruby #$0 path_to_document_html")
desc = Hash.new
%w(s d c z ds zc).each{|tn|
  fname =  File.join(path, tn+".html")
  desc.update parse_html(fname)
}

methods = Lapack.methods
DataTypes.each{|cdt, dt|
  cdt = cdt.downcase
  dmethods = Array.new
  methods.each{|m|
    dmethods.push m if /^#{cdt}/ =~ m
  }
  mts = Array.new
  MatrixTypes.each{|cmt, mt|
    cmt = cmt.downcase
    reg = /^#{cdt}#{cmt}/
    ms = Array.new
    dmethods.each{|m|
      next unless reg =~ m
      ms.push m
    }
    ms.sort!
    unless ms.empty?
      mts.push [cmt,mt]
      File.open(File.join(prefix,"#{cdt}#{cmt}.html"),"w"){|file|
        file.print <<"EOF"
<HTML>
  <HEAD>
    <TITLE>#{dt} routines for #{mt} matrix</TITLE>
  </HEAD>
  <BODY>
    <A NAME="top"></A>
    <H1>#{dt} routines for #{mt} matrix</H1>
    <UL>
EOF
        ms.each{|m|
          file.print <<"EOF"
      <LI><A HREF=\"##{m}\">#{m}</A> : #{desc[m]}</LI>
EOF
        }
        file.print <<"EOF"
    </UL>

EOF
        ms.each{|m|
          file.print <<"EOF"
  <A NAME="#{m}"></A>
  <H2>#{m}</H2>
  #{desc[m]}
  <PRE>
EOF
          stdout_org = STDOUT.dup
          STDOUT.flush
          STDOUT.reopen(file)
          Lapack.send(m)
          STDOUT.flush
          STDOUT.reopen(stdout_org)
          file.print <<"EOF"
    </PRE>
    <A HREF="#top">go to the page top</A>

EOF
        }
        file.print <<"EOF"
    <HR />
    <A HREF="#{cdt}.html">back to matrix types</A>
  </BODY>
</HTML>
EOF
      }
    end
  }

  unless mts.empty?
    File.open(File.join(prefix,"#{cdt}.html"),"w"){|file|
      file.print <<"EOF"
<HTML>
  <HEAD>
    <TITLE>#{dt} routines</TITLE>
  </HEAD>
  <BODY>
    <H1>#{dt} routines</H1>
    <UL>
EOF
      mts.each{|cmt,mt|
        file.print "      <LI><A HREF=\"#{cdt}#{cmt}.html\">#{cmt.upcase}: #{mt}</A></LI>\n"
      }
      file.print <<"EOF"
    </UL>
    <HR />
    <A HREF="index.html">back to data types</A>
  </BODY>
</HTML>
EOF
    }
  end
}

File.open(File.join(prefix,"index.html"),"w"){|file|
  file.print <<"EOF"
<HTML>
  <HEAD>
    <TITLE>LAPACK routines</TITLE>
  </HEAD>
  <BODY>
    <H1>Data types</H1>
    <UL>
EOF
  DataTypes.each{|cdt,dt|
    file.print "      <LI><A HREF=\"#{cdt.downcase}.html\">#{cdt}: #{dt}</A></LI>\n"
  }
  file.print <<"EOF"
    </UL>
  </BODY>
</HTML>
EOF
}

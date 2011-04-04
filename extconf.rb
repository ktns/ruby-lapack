require "mkmf"


def header_not_found(name)
  warn <<EOF
 #{name}.h was not found.
 If you have #{name}.h, try the following:
   % ruby extconf.rb --with-#{name}-include=path
EOF
  exit 1
end

def library_not_found(lname, fname=nil)
  if fname
    warn <<EOF
  #{fname} was not found.
  If you have #{lname} library, try the following:
    % ruby extconf.rb --with-#{lname}-lib=path --with-#{lname}-name=name
  e.g.
    If you have /usr/local/#{lname}/#{fname},
     % ruby extconf.rb --with-#{lname}-lib=/usr/local/#{lname} --with-#{lname}-name=#{fname}
EOF
    exit 1
  else
    warn <<EOF
  lib#{lname}.{a|so} was not found.
  If you have lib#{lname}.{a|so}, try the following:
    % ruby extconf.rb --with-#{lname}-lib=path
EOF
    exit 1
  end
end

def try_func(func, libs, headers = nil, &b)
  headers = cpp_include(headers)
    try_link(<<"SRC", libs, &b) or try_link(<<"SRC", libs, &b)
#{COMMON_HEADERS}
#{headers}
/*top*/
int main() { return 0; }
int MAIN__() { return main(); }
int t() { void ((*volatile p)()); p = (void ((*)()))#{func}; return 0; }
SRC
#{headers}
/*top*/
int main() { return 0; }
int MAIN__() { return main(); }
int t() { #{func}(); return 0; }
SRC
end
    

def find_library(lib, func=nil, name=nil)
  func = "main" if !func or func.empty?
  ldir = with_config(lib+'-lib')
  ldirs = ldir ? Array === ldir ? ldir : ldir.split(File::PATH_SEPARATOR) : []
  $LIBPATH = ldirs | $LIBPATH
  if /\.(a|so)$/ =~ name
    libs = $libs
    $LIBPATH.each{|path|
      f = File.join(path,name)
      if File.exist?(f)
        libs = f + " " + $libs
        break
      end
    }
  else
    name = LIBARG%lib
    libs = append_library($libs, lib)
  end
  paths = {}
  checking_for "#{func}() in #{name}" do
    libpath = $LIBPATH
    begin
      until r = try_func(func, libs) or paths.empty?
        $LIBPATH = libpath | [paths.shift]
      end
      if r
        $libs = libs
        libpath = nil
      end
    ensure
      $LIBPATH = libpath if libpath
    end
    r
  end
end




dir_config("lapack")
unless find_library("lapack")
  library_not_found("lapack",nil)

  warn "LAPACK will be tried to find"

  name = with_config("blas-name","blas_LINUX.a")
  unless have_library(name)
    lib_path = with_config("blas-lib","/usr/local/lib")
    _libarg = LIBARG
    LIBARG.replace "#{lib_path}/%s"
    unless have_library(name)
      library_not_found("blas",name)
    end
    LIBARG.replace _libarg
  end
  name = with_config("lapack-name","lapack_LINUX.a")
  unless have_library(name)
    lib_path = with_config("lapack-lib","/usr/local/lib")
    _libarg = LIBARG
    LIBARG.replace "#{lib_path}/%s"
    unless have_library(name)
      library_not_found("lapack",name)
    end
    LIBARG.replace _libarg
  end
end

sitearchdir = Config::CONFIG["sitearchdir"]
dir_config("narray", sitearchdir, sitearchdir)
unless find_header("narray.h") && have_header("narray_config.h")
  header_not_found("narray")
end
unless find_library("narray", nil, "narray.so")
  library_not_found("narray","narray.so")
end

create_makefile("numru/lapack")

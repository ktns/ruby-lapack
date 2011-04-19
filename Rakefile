require "rubygems"
require "rake/clean"
require "rake/gempackagetask"

version = 1.5
target_prefix = "numru"

# get destdir
if i = ARGV.index{|arg| /\ADESTDIR=/ =~ arg}
  destdir = ARGV[i].sub(/\ADESTDIR=/,"")
else
  destdir = ""
end

# get sitelibdir
if i = ARGV.index{|arg| /\ASITELIBDIR=/ =~ arg}
  libdir = ARGV[i].sub(/\ASITELIBDIR=/,"")
  unless File.exist?(libdir) && File.directory?(libdir)
    raise "SITELIBDIR is invalid: #{sitelibdir}"
  end
else
  libdir = Config::CONFIG["sitelibdir"]
end
# get sitearchdir
if i = ARGV.index{|arg| /\ASITEARCHDIR=/ =~ arg}
  archdir = ARGV[i].sub(/\ASITEARCHLIBDIR=/,"")
  unless File.exist?(archdir) && File.directory?(archdir)
    raise "SITEARCHDIR is invalid: #{sitearchdir}"
  end
else
  archdir = Config::CONFIG["sitearchdir"]
end



NAME = "lapack"
LIBS = FileList["lib/#{target_prefix}/*rb"]
DLLIB = "ext/#{NAME}.so"
so_file = File.join("lib", target_prefix, "#{NAME}.so")


task :default => so_file

desc "Building extensions"
file so_file => DLLIB do
  mkdir_p File.dirname(so_file)
  rm_f so_file
  ln_s DLLIB, so_file
end
file DLLIB => "ext/Makefile" do
  system("cd ext; make")
end
file "ext/Makefile" => "ext/rb_lapack.h" do
  system("cd ext; ruby extconf.rb")
end
file "ext/rb_lapack.h" => "dev/make_csrc.rb" do
  system("ruby dev/make_csrc.rb")
end

desc "Install files to system"
task :install => [:install_so, :install_rb]

task :install_so => DLLIB do
  dst = File.join(destdir, archdir, target_prefix)
  mkdir_p dst
  install DLLIB, dst, :mode => 0755
end

task :install_rb => LIBS do
  dst = File.join(destdir, libdir, target_prefix)
  mkdir_p dst
  LIBS.each do |lib|
    install lib, dst, :mode => 644
  end
end

desc "execute tests"
task :tests => so_file do
  system("testrb tests")
end


CLEAN.include("ext/*.o")
CLOBBER.include("ext/lapack.so")


PKG_FILES = FileList["lib/#{target_prefix}/*rb"]
PKG_FILES.include("ext/rb_lapack.h")
PKG_FILES.include("ext/f2c_minimal.h")
PKG_FILES.include("ext/*.c")
PKG_FILES.include("Rakefile")
PKG_FILES.include("COPYING", "GPL", "README.rdoc")
PKG_FILES.include("doc/*.html", "samples/**/*rb")
PKG_FILES.include("dev/*.rb", "dev/defs/*")
TEST_FILES = FileList["tests/**/*.rb"]

spec = Gem::Specification.new do |s|
  s.name = "ruby-lapack"
  s.version = version
  s.summary = "A Ruby wrapper of Lapack"
  s.description = <<EOL
Ruby-LAPACK is a Ruby wrapper of Lapack, which is a linear algebra package (http://www.netlib.org/lapack/).
EOL
  s.author = "Seiya Nishizawa"
  s.email = "seiya@gfd-dennou.org"
  s.homepage = "http://ruby.gfd-dennou.org/products/ruby-lapack/"
  s.has_rdoc = false
  s.files = PKG_FILES
  s.test_files = TEST_FILES
  s.add_dependency('narray')
  s.extensions = %w(ext/extconf.rb)
end


Rake::GemPackageTask.new(spec) do |pkg|
  pkg.need_tar_gz = true
  pkg.need_tar_bz2 = true
end



binary_pkg = "pkg/#{spec.name}-#{spec.version}-#{Config::CONFIG["arch"]}.gem"
desc "Build binary package"
task :binary_package => binary_pkg

file binary_pkg => so_file do
  files = PKG_FILES.dup
  files.include so_file
  spec.platform = Gem::Platform::CURRENT
  spec.files = files
  spec.extensions = []
  Gem::Builder.new(spec).build
  mv File.basename(binary_pkg), binary_pkg
end

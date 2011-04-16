require "rubygems"
require "rake/clean"
require "rake/gempackagetask"

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


task :default => DLLIB

desc "building extensions"
file DLLIB do
  system("cd ext; make")
end


desc "install files to system"
task :install => [:install_so, :install_rb]

task :install_so => DLLIB do
  install DLLIB, File.join(destdir, archdir, target_prefix), 0755
end

task :install_rb => LIBS do
  LIB.each do |lib|
    install lib, File.join(destdir, libdir, target_prefix), 644
  end
end


CLEAN.include("ext/*.o")
CLOBBER.include("ext/lapack.so")


PKG_FILES = FileList["lib/numru/*rb"]
PKG_FILES.include("ext/*.c", "ext/*h")
PKG_FILES.include("Rakefile")
PKG_FILES.include("COPYING", "GPL", "README.rdoc")
PKG_FILES.include("doc/*.html", "samples/**/*rb")
PKG_FILES.include("dev/*.rb", "dev/defs/*")
TEST_FILES = FileList["tests/**/*.rb"]

spec = Gem::Specification.new do |s|
  s.name = "ruby-lapack"
  s.version = "1.3"
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

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

task :install_so => DLLLIB do
  install DLLIB, File.join(destdir, archdir, target_prefix), 0755
end

task :install_rb => LIBS do
  LIB.each do |lib|
    install lib, File.join(destdir, libdir, target_prefix), 644
  end
end

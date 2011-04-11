$:.unshift File.dirname(__FILE__)

require "yaml"
require "pp"
require "common"

RBPREFIX = "rblapack_"

NATYPES = {
  "integer" => "NA_LINT",
  "real" => "NA_SFLOAT",
  "doublereal" => "NA_DFLOAT",
  "complex" => "NA_SCOMPLEX",
  "doublecomplex" => "NA_DCOMPLEX",
  "logical" => "NA_LINT",
}



def get_cobj(name, type, sub_name)
  case type
  when "integer"
    return "  #{name} = NUM2INT(#{RBPREFIX}#{name});\n"
  when "real"
    return "  #{name} = (real)NUM2DBL(#{RBPREFIX}#{name});\n"
  when "doublereal"
    return "  #{name} = NUM2DBL(#{RBPREFIX}#{name});\n"
  when "complex"
    code =<<"EOF"
  #{name}.r = (real)NUM2DBL(rb_funcall(#{RBPREFIX}#{name}, rb_intern("real"), 0));
  #{name}.i = (real)NUM2DBL(rb_funcall(#{RBPREFIX}#{name}, rb_intern("imag"), 0));
EOF
    return code
  when "doublecomplex"
    code =<<"EOF"
  #{name}.r = NUM2DBL(rb_funcall(#{RBPREFIX}#{name}, rb_intern("real"), 0));
  #{name}.i = NUM2DBL(rb_funcall(#{RBPREFIX}#{name}, rb_intern("imag"), 0));
EOF
    return code
  when "char"
    return "  #{name} = StringValueCStr(#{RBPREFIX}#{name})[0];\n"
  when "logical"
    return "  #{name} = (#{RBPREFIX}#{name} == Qtrue);\n"
  else
    raise "type (#{type}) is not defined in #{name} (#{sub_name})"
  end
end

def get_robj(name, type, flag=false)
  case type
  when "integer"
    cname =  flag ? "(*#{name})" : name
    return "  #{RBPREFIX}#{name} = INT2NUM(#{cname});\n"
  when "real", "doublereal"
    cname =  flag ? "(*#{name})" : name
    return "  #{RBPREFIX}#{name} = rb_float_new((double)#{cname});\n"
  when "complex", "doublecomplex"
    if flag
      r = "(#{name}->r)"
      i = "(#{name}->i)"
    else
      r = "(#{name}.r)"
      i = "(#{name}.i)"
    end
    return "  #{RBPREFIX}#{name} = rb_funcall(rb_gv_get(\"Complex\"), rb_intern(\"new\"), 2, rb_float_new((double)#{r}), rb_float_new((double)#{i}));\n"
  when "char"
    return "  #{RBPREFIX}#{name} = rb_str_new(&#{name},1);\n"
  when "logical"
    return "  #{RBPREFIX}#{name} = #{name} ? Qtrue : Qfalse;\n"
  else
    raise "type (#{type}) is not defined in #{name}"
  end
end


def get_input(name, type, dims, i, varset, sub_name, subst)
  if dims.nil?
    return get_cobj(name, type, sub_name)
  else
    if type == "char"
      return "  #{name} = StringValueCStr(#{RBPREFIX}#{name});\n"
    end
    if i.kind_of?(Integer)
      arg = "#{i+1}th argument"
    else
      arg = "option"
    end
    code =<<"EOF"
  if (!NA_IsNArray(#{RBPREFIX}#{name}))
    rb_raise(rb_eArgError, "#{name} (#{arg}) must be NArray");
  if (NA_RANK(#{RBPREFIX}#{name}) != #{dims.length})
    rb_raise(rb_eArgError, "rank of #{name} (#{arg}) must be %d", #{dims.length});
EOF
    ndim = dims.length
    ndim.times do |jj|
      j = ndim - jj - 1
      dim = dims[j]
      raise "bug: NA_SHAPE? cannot use {#{dim} in #{name}: #{sub_name}" if j>2
      if varset.include?(dim)
        code << <<"EOF"
  if (NA_SHAPE#{j}(#{RBPREFIX}#{name}) != #{dim})
    rb_raise(rb_eRuntimeError, "shape #{j} of #{name} must be #{dim}");
EOF
      elsif (shape = @shape[dim])
        code << <<"EOF"
  if (NA_SHAPE#{j}(#{RBPREFIX}#{name}) != #{dim})
    rb_raise(rb_eRuntimeError, "shape #{j} of #{name} must be the same as shape #{shape[:index]} of #{shape[:name]}");
EOF
      elsif /^[a-z][a-z_\d]*$/ !~ dim
        get_vars(dim).each{|d|
          unless varset.include?(d) || @shape.include?(d)
            raise "undefined #{d}  #{name} #{sub_name}"
          end
        }
        code << <<"EOF"
  if (NA_SHAPE#{j}(#{RBPREFIX}#{name}) != (#{dim}))
    rb_raise(rb_eRuntimeError, "shape #{j} of #{name} must be %d", #{dim});
EOF
      else
        code << "  #{dim} = NA_SHAPE#{j}(#{RBPREFIX}#{name});\n"
        @shape[dim] = {:name => name, :index => j}
        if s = subst[dim]
          code << <<EOF
  if (#{dim} != (#{s}))
    rb_raise(rb_eRuntimeError, "shape #{j} of #{name} must be %d", #{s});
EOF
          if /^[a-z][a-z_\d]*$/ =~ s
            code << "  #{s} = #{dim};\n"
            varset.push s
          end
        end
      end
    end
    natype =  NATYPES[type] || raise("na type is not deifned (#{type})")
    code << <<"EOF"
  if (NA_TYPE(#{RBPREFIX}#{name}) != #{natype})
    #{RBPREFIX}#{name} = na_change_type(#{RBPREFIX}#{name}, #{natype});
  #{name} = NA_PTR_TYPE(#{RBPREFIX}#{name}, #{type}*);
EOF
  end
end



def create_code(name, flag)
  def_fname = File.join(File.dirname(__FILE__), "defs", name)
  hash = nil
  begin
    File.open(def_fname) do |file|
      hash = YAML.load(file.read)
    end
  rescue
    p name
    raise $!
  end

  sub_name = hash[:name]
  sub_type = hash[:category]
  func_type = hash[:type] if sub_type == :function
  args = hash[:arguments]
  help = hash[:fortran_help]
  subst = hash[:substitutions]
  extras = hash[:extras]

  arg_names = Array.new
  args_new = Hash.new
  args.each do |arg|
    arg_names.push arg.keys[0]
    args_new.update(arg)
  end
  args = args_new

  unless sub_name
    warn "this has no subroutine (#{name})"
    return nil
  end
  unless arg_names
    raise "no arg_names (#{name})"
  end
  unless args && !args.empty?
    raise "no args (#{name})"
  end

  unless flag
    return true, sub_name
  end


  inputs = Array.new
  outputs = Array.new
  inouts = Array.new
  workspaces = Array.new
  options = Array.new
  block = nil
  arg_names.each{|aname|
    arg = args[aname]
    unless arg
      raise "arg #{aname} is not defined (#{sub_name} in #{aname})"
    end
    case arg[:intent]
    when "input"
      if arg[:option]
        options.push aname
      else
        inputs.push aname
      end
    when "output"
      outputs.push aname
    when "input/output"
      if arg[:option]
        options.push aname
      else
        inputs.push aname
      end
      inouts.push aname
    when "workspace"
      workspaces.push aname
    when "external procedure"
      if block
        raise "only one block is supported"
      end
      block = aname
    else
      raise "intent (#{arg[:intent]}) is invalid (#{aname}) #{sub_name}"
    end
    unless arg[:type]
      if arg[:intent] == "external procedure"
        arg[:type] = "L_fp"
      else
        raise "type is null (#{aname}) #{sub_name}"
      end
    end
    if d = arg[:dims]
      if d.length == 0
        raise "length of dim is 0 #{aname} #{sub_name}"
      end
    end
  }
  args.each{|key,arg|
    dims = arg[:dims]
    next unless dims
    dims.each{|dim|
      inputs.delete(dim)
    }
  }
  subst.keys.each do |k|
    inputs.delete(k)
  end

  if @@debug
    p "inputs"
    p inputs
    p "outputs"
    p outputs
    p "inouts"
    p inouts
    p "workspaces"
    p workspaces
    p "options"
    p options
    p "block"
    p block
  end

  code = ""

  cargs = arg_names.map do |an|
    arg = args[an]
    t = arg[:type]
    t + (t=="L_fp" ? " " : "* ") + an
  end.join(", ")
  case sub_type
  when :subroutine
    code << "extern VOID #{sub_name}_(#{cargs});\n\n"
  when:function
    outputs.push "__out__"
    args["__out__"] = {:type => func_type}
    if /complex/ =~ func_type || func_type == "char"
      code << "extern VOID #{sub_name}_(#{func_type} *__out__, #{cargs});\n\n"
    else
      code << "extern #{func_type} #{sub_name}_(#{cargs});\n\n"
    end
  else
    raise "category is invalid: #{sub_type} (#{sub_name})"
  end

  code << <<"EOF"

static VALUE
#{RBPREFIX}#{sub_name}(int argc, VALUE *argv, VALUE self){
EOF

  dimdefs = Array.new
  (inputs+options+outputs).each{|aname|
    arg = args[aname]
    code << <<"EOF"
  VALUE #{RBPREFIX}#{aname};
  #{arg[:type]} #{arg[:dims] ? "*" : ""}#{aname}; 
EOF
    dimdefs.push aname
  }
  inouts.each{|aname|
    arg = args[aname]
    if arg[:dims]
      code << <<"EOF"
  VALUE #{RBPREFIX}#{aname}_out__;
  #{arg[:type]} *#{aname}_out__;
EOF
    end
  }
  if extras
    extras.each do |k,v|
      code << "  #{v} #{k};\n"
      dimdefs.push k
    end
  end

  workspaces.each{|aname|
    arg = args[aname]
    code << "  #{arg[:type]} #{arg[:dims] ? "*" : ""}#{aname};\n"
  }
  code << "\n"
  (inputs+options+outputs+workspaces).each{|aname|
    arg = args[aname]
    if dims = arg[:dims]
      dims.each{|dim|
        if dim.kind_of?(Hash)
          p name
        end
        if /^[a-z][a-z_\d]*$/ !~ dim # if dim is a Numeric
          next
        end
        unless dimdefs.include?(dim) # untill it's defined
          code << "  integer #{dim};\n"
          dimdefs.push dim
        end
      }
    end
  }
  if ss = subst
    ss.each{|k,v|
      unless dimdefs.include?(k)
        code << "  integer #{k};\n"
        dimdefs.push k
      end
    }
  end



  code << "\n"

  if block
    block_help = "{|" + ['a','b','c'][0...args[block][:block_arg_num]].join(",") + "| ... }"
  else
    block_help = ""
  end

  usage_code = <<"EOF"
USAGE:
  #{(outputs+inouts).join(", ")} = NumRu::Lapack.#{sub_name}( #{inputs.join(", ")}, [#{(options+["usage","help"]).map{|on| ":"+on+" => "+on}.join(", ")}])#{block_help}
EOF

  help_code = <<"EOF"
#{usage_code}

FORTRAN MANUAL
#{help}
EOF
  ilen = inputs.length
  code << <<"EOF"
  VALUE #{RBPREFIX}options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    #{RBPREFIX}options = argv[argc];
    if (rb_hash_aref(#{RBPREFIX}options, sHelp) == Qtrue) {
      printf("%s\\n", "#{help_code.gsub(/\\/,'\\\\\\').gsub(/\n/,'\n').gsub(/"/,'\"')}");
      rb_exit(0);
    }
    if (rb_hash_aref(#{RBPREFIX}options, sUsage) == Qtrue) {
      printf("%s\\n", "#{usage_code.gsub(/\\/,'\\\\\\').gsub(/\n/,'\n').gsub(/"/,'\"')}");
      rb_exit(0);
    } 
  } else
    #{RBPREFIX}options = Qnil;
  if (argc != #{ilen})
    rb_raise(rb_eArgError,"wrong number of arguments (%d for #{ilen})", argc);
EOF
  inputs.each_with_index{|arg,i|
    code << "  #{RBPREFIX}#{arg} = argv[#{i}];\n"
  }
  code << "  if (#{RBPREFIX}options != Qnil) {\n"
  options.each do |opt|
    code << "    #{RBPREFIX}#{opt} = rb_hash_aref(#{RBPREFIX}options, ID2SYM(rb_intern(\"#{opt}\")));\n"
  end
  code << "  } else {\n"
  options.each do |opt|
    code << "    #{RBPREFIX}#{opt} = Qnil;\n"
  end
  code << "  }\n"

  code << "\n"

  order = Hash.new
  (inputs+options).each_with_index do |arg,i|
    aryd = Array.new
    aryp = Array.new
    if dim = args[arg][:dims]
      dim.each do |d|
        vs = get_vars(d)
        if vs.length==1 && vs[0] == d && !subst.keys.include?(d)
          aryp.push d
        else
          aryd.push vs
        end
      end
    end
    if vs = args[arg][:default]
      get_vars(vs).each do |v|
        aryd.push v
      end
    end
    aryd.flatten!
    aryp.uniq!
    aryd.uniq!
    order[arg] = {:depends => aryd, :type => :input, :order => i, :provides => aryp}
  end
  subst.each do |k,v|
    order[k] = {:depends => get_vars(v).uniq, :type => :subst, :value => v}
  end
  order = order.sort do |a0,a1|
    k0, v0 = a0
    k1, v1 = a1
    d0 = v0[:depends]
    d1 = v1[:depends]
    p0 = v0[:provides]
    p1 = v1[:provides]
    if d0.empty? && d1.empty?
      0
    elsif d0.empty?
      -1
    elsif d1.empty?
      1
    else
      flag0 = d0.include?(k1)
      flag1 = d1.include?(k0)
      if p0
        p0.each do |p|
          if d1.include?(p)
            flag1 = true
            break
          end
        end
      end
      if p1
        p1.each do |p|
          if d0.include?(p)
            flag0 = true
            break
          end
        end
      end
      if flag0 && flag1
        pp order
        p [k0, k1]
        pp v0
        pp v1
        raise "depends each other #{name}"
      end
      flag0 ? 1 : flag1 ? -1 : 0
    end
  end

 if @@debug
   p "order"
   pp order
 end

  varset = Array.new
  @shape = Hash.new
  order.each do |name, v|
    if v[:type] == :input
      arg = args[name]
      if arg[:option]
        code << <<EOF
  if (#{RBPREFIX}options == Qnil || #{RBPREFIX}#{name} == Qnil)
    #{name} = #{arg[:default] || "NULL"};
  else
  #{get_input(name, arg[:type], arg[:dims], :option, varset, sub_name, subst)}
EOF
      else
        code << get_input(name, arg[:type], arg[:dims], v[:order], varset, sub_name, subst)
      end
    else
      unless varset.include?(name)
        code << "  #{name} = #{v[:value]};\n"
      end
    end
    varset.push name
  end

  outputs.each{|name|
    arg = args[name]
    type = arg[:type]
    if dims = arg[:dims]
      code << <<"EOF"
  {
    int shape[#{dims.length}];
EOF
      dims.each_with_index{|dim,k|
        get_vars(dim).each{|d|
          unless varset.include?(d) || @shape[d]
            raise "undefined #{d}  #{name} #{sub_name}"
          end
        }
        code << "    shape[#{k}] = #{dim};\n"
      }
      code << <<"EOF"
    #{RBPREFIX}#{name} = na_make_object(#{NATYPES[type]}, #{dims.length}, shape, cNArray);
  }
  #{name} = NA_PTR_TYPE(#{RBPREFIX}#{name}, #{type}*);
EOF
    end
  }

  inouts.each{|name|
    arg = args[name]
    type = arg[:type]
    if dims = arg[:dims]
      if outdims = arg[:outdims]
        if outdims.length != dims.length
          raise "dimensions for input and output are different: #{dims.join(",")} and #{outdims.join(",")}"
        end
      end
      code << <<"EOF"
  {
    int shape[#{dims.length}];
EOF
      dims.each_with_index{|dim,k|
        ds = get_vars(dim)
        if outdims
          od = outdims[k]
          ds += get_vars(od)
        end
        ds.each{|d|
          unless varset.include?(d) || @shape[d]
            raise "undefined #{d}  #{name} #{sub_name}"
          end
        }
        d = outdims && dim != od ? "MAX(#{dim}, #{od})" : dim
        code << "    shape[#{k}] = #{d};\n"
      }
      code << <<"EOF"
    #{RBPREFIX}#{name}_out__ = na_make_object(#{NATYPES[type]}, #{dims.length}, shape, cNArray);
  }
  #{name}_out__ = NA_PTR_TYPE(#{RBPREFIX}#{name}_out__, #{type}*);
EOF
      if outdims
        sh = Array.new
        ndims = dims.length
        code << <<EOF
  {
    VALUE __shape__[#{ndims+1}];
EOF
        ndims.times do |n|
          d = dims[n]
          od = outdims[n]
          if d == od
            code << "    __shape__[#{n}] = Qtrue;\n"
          else
            code << "    __shape__[#{n}] = #{d} < #{od} ? rb_range_new(#{RBPREFIX}ZERO, INT2NUM(#{d}), Qtrue) : Qtrue;\n"
          end # if d == od
        end # ndims.times do
        code << <<"EOF"
    __shape__[#{ndims}] = #{RBPREFIX}#{name};
    na_aset(#{ndims+1}, __shape__, #{RBPREFIX}#{name}_out__);
  }
EOF
      else
        code << <<"EOF"
  MEMCPY(#{name}_out__, #{name}, #{type}, NA_TOTAL(#{RBPREFIX}#{name}));
EOF
      end
      code << <<"EOF"
  #{RBPREFIX}#{name} = #{RBPREFIX}#{name}_out__;
  #{name} = #{name}_out__;
EOF
    end
  }

  workspaces.each{|name|
    arg = args[name]
    if dims = arg[:dims]
      type = arg[:type]
      dims.each{|dim|
        get_vars(dim).each{|d|
          unless varset.include?(d) || @shape[d]
            raise "undefined #{d}  #{name} #{sub_name}"
          end
        }
      }
      len = dims.collect{|dim| "(#{dim})"}.join("*")
      code << "  #{name} = ALLOC_N(#{type}, #{len});\n"
    end
  }
  code << "\n"


  cargs = arg_names.collect{|name|
    block== name ? RBPREFIX+name : args[name][:dims] ? name : "&"+name
  }
  if sub_type == :function
    if /complex/ =~ func_type || func_type == "char"
      code << "  #{sub_name}_(&__out__, #{cargs.join(", ")});\n\n"
    else
      code << "  __out__ = #{sub_name}_(#{cargs.join(", ")});\n\n"
    end
  else
    code << "  #{sub_name}_(#{cargs.join(", ")});\n\n"
  end

  workspaces.each{|name|
    arg = args[name]
    if dims = arg[:dims]
      type = arg[:type]
      len = dims.collect{|dim| "(#{dim})"}.join("*")
      code << "  free(#{name});\n"
    end
  }

  out = outputs + inouts

  out.each{|name|
    arg = args[name]
    if dims = arg[:dims]
      if arg[:type] == "char"
        code << "  #{RBPREFIX}#{name} = rb_str_new2(&#{name});\n"
      elsif outdims = arg[:outdims]
        ndims = dims.length
        code << <<EOF
  {
    VALUE __shape__[#{ndims}];
EOF
        ndims.times do |n|
          d = dims[n]
          od = outdims[n]
          if d == od
            code << "    __shape__[#{n}] = Qtrue;\n"
          else
            code << "    __shape__[#{n}] = #{d} < #{od} ? Qtrue : rb_range_new(#{RBPREFIX}ZERO, INT2NUM(#{od}), Qtrue);\n"
          end # if d == od
        end # ndims.times do
        code << <<"EOF"
    #{RBPREFIX}#{name} = na_aref(#{ndims}, __shape__, #{RBPREFIX}#{name});
  }
EOF
      end
    else
      if name == "__out__" && arg[:type].nil?
        p sub_name
        p hash[:fortran_help].split("\n")[0]
      end
      code << get_robj(name, arg[:type])
    end
  }

  case out.length
  when 0
    result = "Qnil"
  when 1
    result = RBPREFIX+out[0];
  else
    result = "rb_ary_new3(#{out.length}, #{out.collect{|op| RBPREFIX+op}.join(", ")})"
  end

  code << <<"EOF"
  return #{result};
}

EOF

  code << <<"EOF"
void
init_lapack_#{sub_name}(VALUE mLapack){
  rb_define_module_function(mLapack, \"#{sub_name}\", #{RBPREFIX}#{sub_name}, -1);
}
EOF

  code_all = "#include \"rb_lapack.h\"\n\n"

  if block
    arg = args[block]
    type = arg[:block_type]
    atype = arg[:block_arg_type]
    anum = arg[:block_arg_num]
    cas = Array.new
    ras = Array.new
    anum.times{|n|
      cas.push "#{atype} *arg#{n}"
      ras.push "#{RBPREFIX}arg#{n}"
    }
    code_all << <<EOF
static #{type}
#{RBPREFIX}#{block}(#{cas.join(", ")}){
  VALUE #{ras.join(", ")};

  VALUE #{RBPREFIX}ret;
  #{type} ret;

EOF
    anum.times{|n|
      code_all << get_robj("arg#{n}", atype, true)
    }
    code_all << <<EOF

  #{RBPREFIX}ret = rb_yield_values(#{anum}, #{ras.join(", ")});

EOF
    code_all << get_cobj("ret", type, sub_name)
    code_all << <<EOF
  return ret;
}

EOF
  end

  code_all << code

  return [code_all, sub_name]
end

def generate_code(fnames, names)
  nfnames = fnames.length
  sub_names = Array.new
  fnames.each_with_index{|fname,i|
    print "#{i+1}/#{nfnames}\n" if (i+1)%100==0
    name = File.basename(fname)
    flag = names.nil? || names.include?(name)
    code, sub_name = create_code(name, flag)
    if code
      sub_names.push sub_name
      if flag
        File.open(sub_name+".c","w"){|file|
          file.print code
        }
      end
    end
  }

  File.open("rb_lapack.h","w"){|file|
    file.print <<"EOF"
#include <string.h>
#include <math.h>
#include "ruby.h"
#include "narray.h"
#include "f2c_minimal.h"

#define MAX(a,b) (a > b ? a : b)
#define MIN(a,b) (a < b ? a : b)
#define LG(n) ((int)ceil(log((double)n)/log(2.0)))

extern logical lsame_(char *ca, char *cb);
extern integer ilatrans_(char* trans);
extern integer ilaenv_(integer* ispec, char* name, char* opts, integer* n1, integer* n2, integer* n3, integer* n4);


VALUE sHelp, sUsage;
VALUE #{RBPREFIX}ZERO;

EOF
  }


  File.open("rb_lapack.c","w"){|file|
    file.print <<"EOF"
#include "ruby.h"
#include "rb_lapack.h"

EOF

    sub_names.each{|sname|
      file.print "extern void init_lapack_#{sname}(VALUE mLapack);\n"
    }

    file.print <<"EOF"

void Init_lapack(){
  VALUE mNumRu;
  VALUE mLapack;

  rb_require("narray");

  mNumRu = rb_define_module("NumRu");
  mLapack = rb_define_module_under(mNumRu, "Lapack");

  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));

  #{RBPREFIX}ZERO = INT2NUM(0);

EOF
    sub_names.each{|sname|
      file.print "  init_lapack_#{sname}(mLapack);\n"
    }
    file.print "}\n"
  }
end





@@debug = ARGV.delete("--debug")

dname = ARGV.shift || raise("Usage: ruby #$0 path_to_lapack_src [name0, name1, ..]")
unless File.directory?(dname)
  raise "the first argument must be directory"
end

unless ARGV.empty?
  names = ARGV
  @@debug = true
else
  names = nil
end

reg = File.join(dname, "[a-z]*[a-z0-9]")
fnames = Dir[reg]

generate_code(fnames, names)

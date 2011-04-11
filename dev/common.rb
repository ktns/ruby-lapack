def get_vars(dim)
  ary = Array.new
  dim.gsub(/MAX\(/,",").gsub(/MIN\(/,",").gsub(/log\(/,",").gsub(/abs\(/,",").gsub(/sqrt\(/,",").gsub(/pow\(/,",").gsub(/LG\(/,",").gsub(/lsame_\(\&([^,]+),[^)]+\)/,'\1').gsub(/ilatrans_\([^)]+\)/,",").gsub(/ilaenv_\(([^,]+),[^,]+/,'\1').gsub(/[\(\)\+\-\*\/:\?=\&\|]+/,",").split(",").each{|d|
    d.strip!
    next if (d == "") || (/^\d+(\.\d+)?$/ =~ d) || /^\"[^\"]+\"$/ =~ d || d=="int" || d=="double"
    ary.push d
  }
  ary
end

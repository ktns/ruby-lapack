require "narray"
require "numru/lapack.so"



class NMatrix

  # to lapack matrix
  def to_lm
    NArray.ref(self.transpose)
  end

  # to lapack band storage
  def to_lb(kl, ku)
    n = shape[0]
    na = NArray.ref(self)
    lb = NArray.new(typecode, 2*kl+ku+1, n)
    n.times do |j|
      i0 = [n-1,j+kl].min
      i1 = [0,j-ku].max
      l = i0 - i1 + 1
      lb[-i1-1..-i0-1,j] = na[j,i0..i1]
    end
    lb
  end
end

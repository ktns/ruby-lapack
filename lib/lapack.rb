require "narray"
require "numru/lapack.so"



class NMatrix
  def to_lm
    NArray.ref(self.transpose)
  end
end

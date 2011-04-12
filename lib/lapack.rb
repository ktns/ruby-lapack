require "narray"
require "numru/lapack.so"



class NMatrix

  # to lapack matrix
  def to_lm
    NArray.ref(self.transpose)
  end

  # to lapack band matrix
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

  # to lapack symmetrix band matrix
  def to_lsb(uplo, kd)
    n = shape[0]
    lsb = NArray.new(typecode, kd+1, n)
    na = NArray.ref(self)
    case uplo
    when /U/i
      n.times do |j|
        i0 = [0,j-kd].max
        i1 = j
        lsb[i0+kd-j..i1+kd-j, j] = na[j,i0..i1]
      end
    when /L/i
      n.times do |j|
        i0 = j
        i1 = [n-1,j+kd].min
        lsb[i0-j..i1-j, j] = na[j,i0..i1]
      end
    else
      raise "uplo is invalid"
    end
    lsb
   end

end

require "test/unit"
require "numru/lapack"

class GelsdTest < Test::Unit::TestCase

  def setup
    @a = NMatrix[[-0.09, -1.56, -1.48, -1.09,  0.08, -1.59],
                 [ 0.14,  0.20, -0.43,  0.84,  0.55, -0.72],
                 [-0.46,  0.29,  0.89,  0.77, -1.13,  1.06],
                 [ 0.68,  1.09, -0.71,  2.11,  0.14,  1.24],
                 [ 1.29,  0.51, -0.96, -1.27,  1.74,  0.34]].to_lm
    @b = NVector[[7.4, 4.3, -8.1, 1.8, 8.7]]
    @rcond = 0.01

    @bout = NArray[[1.5938, -0.1180, -3.1501, 0.1554, 2.5529, -1.6730]]
    @s = NArray[3.9997, 2.9962, 2.0001, 0.9988, 0.0025]
  end

  %w(s d).each do |sd|
    method = "#{sd}gelsd"
    define_method("test_#{method}") do
      s, rank, work, info, b = NumRu::Lapack.send(method, @a, @b, @rcond, -1)
      assert_equal 0, info
      lwork = work[0].to_i
      assert_equal 991, lwork
      s, rank, work, info, b = NumRu::Lapack.send(method, @a, @b, @rcond, lwork)
      assert_equal 0, info
      assert_equal lwork, work[0].to_i
      assert_equal 4, rank
      assert (@bout-b).abs.max < 10e-4
      assert (@s-s).abs.max  < 10e-4
    end
  end

end

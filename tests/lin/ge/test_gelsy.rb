require "test/unit"
require "numru/lapack"

class GelsyTest < Test::Unit::TestCase

  def setup
    @a = NMatrix[[-0.09,  0.14, -0.46,  0.68,  1.29],
                 [-1.56,  0.20,  0.29,  1.09,  0.51],
                 [-1.48, -0.43,  0.89, -0.71, -0.96],
                 [-1.09,  0.84,  0.77,  2.11, -1.27],
                 [ 0.08,  0.55, -1.13,  0.14,  1.74],
                 [-1.59, -0.72,  1.06,  1.24,  0.34]].to_lm
    @b = NVector[[7.4, 4.2, -8.3, 1.8, 8.6, 2.1]]
    @rcond = 0.01
    @n = @a.shape[1]
    @jpvt = NArray.int(@n)

    @lss = NArray[[0.6344, 0.9699, -1.4403, 3.3678, 3.3992]]
  end

  %w(s d).each do |sd|
    method = "#{sd}gelsy"
    define_method("test_#{method}") do
      rank, work, info, a, b, jpvt = NumRu::Lapack.send(method, @a, @b, @jpvt, @rcond, :lwork => -1)
      assert_equal(0, info)
      lwork = work[0].to_i
      assert_equal(27, lwork)
      [lwork, nil].each do |lw|
      rank, work, info, a, b, jpvt = NumRu::Lapack.send(method, @a, @b, @jpvt, @rcond, :lwork => lw)
        assert_equal(0, info)
        assert_equal(lwork, work[0].to_i)
        assert (@lss - b).abs.max < 1e-4
        assert 4, rank
      end
    end
  end

end

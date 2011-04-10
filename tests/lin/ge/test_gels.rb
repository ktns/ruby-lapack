require "test/unit"
require "numru/lapack"

class DelsTest < Test::Unit::TestCase

  def setup
    @a = NMatrix[[-0.57, -1.28, -0.39,  0.25],
                 [-1.93,  1.08, -0.31, -2.14],
                 [ 2.30,  0.24,  0.40, -0.35],
                 [-1.93,  0.64, -0.66,  0.08],
                 [ 0.15,  0.30,  0.15, -2.13],
                 [-0.02,  1.03, -1.43,  0.50]]

    @b = NVector[[-2.67, -0.55, 3.34, -0.77, 0.48, 4.10]]
  end

  def test_dgels
    work, info, a, b = NumRu::Lapack.dgels("N", @a.to_lm, @b)
    assert_equal(0, info)
    lwork = work[0].to_i
    work, info, a, b = NumRu::Lapack.dgels("N", @a.to_lm, @b, :lwork => -1)
    assert_equal(0, info)
    assert_equal(lwork, work[0].to_i)
    work, info, a, b = NumRu::Lapack.dgels("N", @a.to_lm, @b, :lwork => lwork)
    assert_equal(0, info)
    assert_equal(12, work[0])
    n = @a.shape[0]
    assert( (NArray[[1.5339, 1.8707, -1.5241, 0.0392]]-b[0...n,true]).abs.max < 1.0e-4 )
  end
end

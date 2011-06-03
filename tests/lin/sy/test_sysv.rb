$:.push File.dirname(__FILE__) + "/../.."
require "lapack_test"

class SysvTest < Test::Unit::TestCase
  include LapackTest

  def setup
    @a = Hash.new
    @b = Hash.new
    @b_exp = Hash.new
    @ipiv_exp = Hash.new

    @a[:r] = NMatrix[[-1.81,  2.06,  0.63, -1.15],
                     [ 2.06,  1.15,  1.87,  4.20],
                     [ 0.63,  1.87, -0.21,  3.87],
                     [-1.15,  4.20,  3.87,  2.07]].to_lm
    @b[:r] = NVector[[0.96, 6.07, 8.38, 9.50]]
    @b_exp[:r] = NArray[[-5.0, -2.0, 1.0, 4.0]]
    @ipiv_exp[:r] = NArray[1, 2, -2, -2]

    @a[:c] = NMatrix[[-0.56+0.12*I, -1.54-2.86*I,  5.32-1.59*I,  3.80+0.92*I],
                     [-1.54-2.86*I, -2.83-0.03*I, -3.52+0.58*I, -7.86-2.96*I],
                     [ 5.32-1.59*I, -3.52+0.58*I,  8.86+1.81*I,  5.14-0.64*I],
                     [ 3.80+0.92*I, -7.86-2.96*I,  5.14-0.64*I, -0.39-0.71*I]].to_lm
    @b[:c] = NVector[[-6.43+19.24*I, -0.49-1.47*I, -48.18+66.0*I, -55.64+41.22*I]]
    @b_exp[:c] = NArray[[-4.0+3.0*I, 3.0-2.0*I, -2.0+5.0*I, 1.0-1.0*I]]
    @ipiv_exp[:c] = NArray[1, 2, -2, -2]
  end

  %w(s d c z).each do |x|
    method = "#{x}sysv"
    rc = LapackTest.get_rc(x)

    define_method("test_#{method}") do
      ipiv, work, info, a, b = NumRu::Lapack.send(method, "U", @a[rc], @b[rc])
      assert_equal 0, info
      assert_narray @b_exp[rc], b
      assert_narray @ipiv_exp[rc], ipiv
    end

    define_method("test_#{method}_inquiring_lwork") do
      rank, work, info, = NumRu::Lapack.send(method, "U", @a[rc], @b[rc], :lwork => -1)
      assert_equal 0, info
      lwork = get_int(work[0])
      ipiv, work, info, a, b = NumRu::Lapack.send(method, "U", @a[rc], @b[rc], :lwork => lwork)
      assert_equal 0, info
      assert_equal lwork, get_int(work[0])
      assert_narray @b_exp[rc], b
      assert_narray @ipiv_exp[rc], ipiv
    end

    define_method("test_#{method}_inquiring_lwork_oldargstyle") do
      rank, work, info, = NumRu::Lapack.send(method, "U", @a[rc], @b[rc], :lwork => -1)
      assert_equal 0, info
      lwork = get_int(work[0])
      rank, work, info, = NumRu::Lapack.send(method, "U", @a[rc], @b[rc], -1)
      assert_equal 0, info
      assert_equal lwork, get_int(work[0])
    end

  end

end

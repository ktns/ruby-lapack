$:.push File.dirname(__FILE__) + "/../.."
require "lapack_test"

class GesvTest < Test::Unit::TestCase
  include LapackTest

  def setup
    @a = Hash.new
    @b = Hash.new
    @b_exp = Hash.new
    @ipiv_exp = Hash.new

    @a[:r] = NMatrix[[ 1.80,  2.88,  2.05, -0.89],
                     [ 5.25, -2.95, -0.95, -3.80],
                     [ 1.58, -2.69, -2.90, -1.04],
                     [-1.11, -0.66, -0.59,  0.80]].to_lm
    @b[:r] = NVector[[9.52, 24.35, 0.77, -6.22]]
    @b_exp[:r] = NArray[[1.0, -1.0, 3.0, -5.0]]
    @ipiv_exp[:r] = NArray[2, 2, 3, 4]

    @a[:c] = NMatrix[[-1.34+2.55*I,  0.28+3.17*I, -6.39-2.20*I,  0.72-0.92*I],
                     [-0.17-1.41*I,  3.31-0.15*I, -0.15+1.34*I,  1.29+1.38*I],
                     [-3.29-2.39*I, -1.91+4.42*I, -0.14-1.35*I,  1.72+1.35*I],
                     [ 2.41+0.39*I, -0.56+1.47*I, -0.83-0.69*I, -1.96+0.67*I]].to_lm
    @b[:c] = NVector[[26.26+51.78*I, 6.43-8.68*I, -5.75+25.31*I, 1.16+2.57*I]]
    @b_exp[:c] = NArray[[1.0+1.0*I, 2.0-3.0*I, -4.0-5.0*I, 0.0+6.0*I]]
    @ipiv_exp[:c] = NArray[3, 2, 3, 4]
  end

  %w(s d c z).each do |x|
    method = "#{x}gesv"
    rc = LapackTest.get_rc(x)

    define_method("test_#{method}") do 
      ipiv, info, a, b = NumRu::Lapack.send(method, @a[rc], @b[rc])
      assert_equal 0, info
      assert_narray @b_exp[rc], b
      assert_narray @ipiv_exp[rc], ipiv
    end

  end

end

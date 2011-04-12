$:.push File.dirname(__FILE__) + "/../.."
require "lapack_test"

class GgsvdTest < Test::Unit::TestCase
  include LapackTest

  def setup
    @a = Hash.new
    @b = Hash.new
    @k_exp = Hash.new
    @l_exp = Hash.new
    @gsv_exp = Hash.new
    @u_exp = Hash.new
    @v_exp = Hash.new
    @q_exp = Hash.new

    @a[:r] =  NMatrix[[1.0, 2.0, 3.0],
                      [3.0, 2.0, 1.0],
                      [4.0, 5.0, 6.0],
                      [7.0, 8.0, 8.0]].to_lm
    @b[:r] = NMatrix[[-2.0, -3.0, 3.0],
                     [ 4.0,  6.0, 5.0]].to_lm
    @k_exp[:r] = 1
    @l_exp[:r] = 2
    @gsv_exp[:r] = NArray[1.3151, 8.0185e-2]
    @u_exp[:r] = NMatrix[[-1.3484e-1,  5.2524e-1, -2.0924e-1,  8.1373e-1 ],
                         [ 6.7420e-1, -5.2213e-1, -3.8886e-1,  3.4874e-1 ],
                         [ 2.6968e-1,  5.2757e-1, -6.5782e-1, -4.6499e-1 ],
                         [ 6.7420e-1,  4.1615e-1,  6.1014e-1,  1.5127e-15]].to_lm
    @v_exp[:r] = NMatrix[[3.5539e-1, -9.3472e-1],
                         [9.3472e-1,  3.5539e-1]].to_lm
    @q_exp[:r] = NMatrix[[-8.3205e-1, -9.4633e-2, -5.4657e-1],
                         [ 5.5470e-1, -1.4195e-1, -8.1985e-1],
                         [ 0.0000e+0, -9.8534e-1,  1.7060e-1]].to_lm

    @a[:c] = NMatrix[[ 0.96-0.81*I, -0.03+0.96*I, -0.91+2.06*I, -0.05+0.41*I],
                     [-0.98+1.98*I, -1.20+0.19*I, -0.66+0.42*I, -0.81+0.56*I],
                     [ 0.62-0.46*I,  1.01+0.02*I,  0.63-0.17*I, -1.11+0.60*I],
                     [ 0.37+0.38*I,  0.19-0.54*I, -0.98-0.36*I,  0.22-0.20*I],
                     [ 0.83+0.51*I,  0.20+0.01*I, -0.17-0.46*I,  1.47+1.59*I],
                     [ 1.08-0.28*I,  0.20-0.12*I, -0.07+1.23*I,  0.26+0.26*I]].to_lm
    @b[:c] = NMatrix[[ 1.00+0.00*I,  0.00+0.00*I, -1.00+0.00*I,  0.00+0.00*I],
                     [ 0.00+0.00*I,  1.00+0.00*I,  0.00+0.00*I, -1.00+0.00*I]].to_lm
    @k_exp[:c] = 2
    @l_exp[:c] = 2
    @gsv_exp[:c] = NArray[2.0720e+0, 1.1058e+0]
    @u_exp[:c] = NMatrix[[-1.3038e-02-3.2595e-01*I, -1.4039e-01-2.6167e-01*I,  2.5177e-01-7.9789e-01*I, -5.0956e-02-2.1750e-01*I, -4.5947e-02+1.4052e-04*I, -5.2773e-02-2.2492e-01*I],
                         [ 4.2764e-01-6.2582e-01*I,  8.6298e-02-3.8174e-02*I, -3.2188e-01+1.6112e-01*I,  1.1979e-01+1.6319e-01*I, -8.0311e-02-4.3605e-01*I, -3.8117e-02-2.1907e-01*I],
                         [-3.2595e-01+1.6428e-01*I,  3.8163e-01-1.8219e-01*I,  1.3231e-01-1.4565e-02*I, -5.0671e-01+1.8615e-01*I,  5.9714e-02-5.8974e-01*I, -1.3850e-01-9.0941e-02*I],
                         [ 1.5906e-01-5.2151e-03*I, -2.8207e-01+1.9732e-01*I,  2.1598e-01+1.8813e-01*I, -4.0163e-01+2.6787e-01*I, -4.6443e-02+3.0864e-01*I, -3.7354e-01-5.5148e-01*I],
                         [-1.7210e-01-1.3038e-02*I, -5.0942e-01-5.0319e-01*I,  3.6488e-02+2.0316e-01*I,  1.9271e-01+1.5574e-01*I,  5.7843e-01-1.2439e-01*I, -1.8815e-02-5.5686e-02*I],
                         [-2.6336e-01-2.4772e-01*I, -1.0861e-01+2.8474e-01*I,  1.0906e-01-1.2712e-01*I, -8.8159e-02+5.6169e-01*I,  1.5763e-02+4.7130e-02*I,  6.5007e-01+4.9173e-03*I]].to_lm
    @v_exp[:c] = NMatrix[[ 9.8930e-01+1.9041e-19*I, -1.1461e-01+9.0250e-02*I],
                         [-1.1461e-01-9.0250e-02*I, -9.8930e-01+1.9041e-19*I]].to_lm
    @q_exp[:c] = NMatrix[[7.0711e-01+0.0000e+00*I,  0.0000e+00+0.0000e+00*I,  6.9954e-01+4.7274e-19*I,  8.1044e-02-6.3817e-02*I],
                         [0.0000e+00+0.0000e+00*I,  7.0711e-01+0.0000e+00*I, -8.1044e-02-6.3817e-02*I,  6.9954e-01-4.7274e-19*I],
                         [7.0711e-01+0.0000e+00*I,  0.0000e+00+0.0000e+00*I, -6.9954e-01-4.7274e-19*I, -8.1044e-02+6.3817e-02*I],
                         [0.0000e+00+0.0000e+00*I,  7.0711e-01+0.0000e+00*I,  8.1044e-02+6.3817e-02*I, -6.9954e-01+4.7274e-19*I]].to_lm
  end

  %w(s d c z).each do |x|
    method = "#{x}ggsvd"
    rc = LapackTest.get_rc(x)

    define_method("test_#{method}") do
      k, l, alpha, beta, u, v, q, iwork, info, a, b = NumRu::Lapack.send(method, "U", "V", "Q", @a[rc], @b[rc])
        assert_equal 0, info
        assert_narray @gsv_exp[rc], alpha[k...k+l]/beta[k...k+l], 1.0e-4
        assert_narray @u_exp[rc], u, 1.0e-4
        assert_narray @v_exp[rc], v, 1.0e-4
        assert_narray @q_exp[rc], q, 1.0e-4
      end

  end

end

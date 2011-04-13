$:.push File.dirname(__FILE__) + "/../.."
require "lapack_test"

class GesvdTest < Test::Unit::TestCase
  include LapackTest

  def setup
    @a = Hash.new
    @s_exp = Hash.new
    @vt_exp = Hash.new
    @a_exp = Hash.new

    @a[:r] = NMatrix[[ 2.27, -1.54,  1.15, -1.94],
                     [ 0.28, -1.67,  0.94, -0.78],
                     [-0.48, -3.09,  0.99, -0.21],
                     [ 1.07,  1.22,  0.79,  0.63],
                     [-2.35,  2.93, -1.45,  2.30],
                     [ 0.62, -7.39,  1.03, -2.57]].to_lm
    @s_exp[:r] = NArray[9.9966, 3.6831, 1.3569, 0.5000]
    @vt_exp[:r] = NMatrix[[-0.1921,  0.8794, -0.2140,  0.3795],
                          [ 0.8030,  0.3926,  0.2980, -0.3351],
                          [ 0.0041, -0.0752,  0.7827,  0.6178],
                          [-0.5642,  0.2587,  0.5027, -0.6017]].to_lm
    @a_exp[:r] = NMatrix[[-0.2774,  0.6003, -0.1277,  0.1323],
                         [-0.2020,  0.0301,  0.2805,  0.7034],
                         [-0.2918, -0.3348,  0.6453,  0.1906],
                         [ 0.0938,  0.3699,  0.6781, -0.5399],
                         [ 0.4213, -0.5266,  0.0413, -0.0575],
                         [-0.7816, -0.3353, -0.1645, -0.3957] ].to_lm

    @a[:c] = NMatrix[[ 0.96+0.81*I, -0.03-0.96*I, -0.91-2.06*I, -0.05-0.41*I],
                     [-0.98-1.98*I, -1.20-0.19*I, -0.66-0.42*I, -0.81-0.56*I],
                     [ 0.62+0.46*I,  1.01-0.02*I,  0.63+0.17*I, -1.11-0.60*I],
                     [-0.37-0.38*I,  0.19+0.54*I, -0.98+0.36*I,  0.22+0.20*I],
                     [ 0.83-0.51*I,  0.20-0.01*I, -0.17+0.46*I,  1.47-1.59*I],
                     [ 1.08+0.28*I,  0.20+0.12*I, -0.07-1.23*I,  0.26-0.26*I]].to_lm
    @s_exp[:c] = NArray[3.9994, 3.0003, 1.9944, 0.9995]
    @vt_exp[:c] = NMatrix[[ 0.6971+0.0*I,  0.0867-0.3548*I, -0.0560-0.5400*I, 0.1878-0.2253*I],
                          [ 0.2403+0.0*I,  0.0725+0.2336*I, -0.2477+0.5291*I, 0.7026-0.2177*I],
                          [-0.5123+0.0*I, -0.3030+0.1735*I,  0.0678-0.5162*I, 0.4418-0.3864*I],
                          [-0.4403+0.0*I,  0.5294-0.6361*I, -0.3027+0.0346*I, 0.1667-0.0258*I]].to_lm
    @a_exp[:c] = NMatrix[[ 0.5634+0.0016*I, -0.2687+0.2749*I, 0.2451-0.4657*I, 0.3787-0.2987*I],
                         [-0.1205-0.6108*I, -0.2909-0.1085*I,  0.4329+0.1758*I, -0.0182+0.0437*I],
                         [ 0.0816+0.1613*I, -0.1660-0.3885*I, -0.4667-0.3821*I, -0.0800+0.2276*I],
                         [-0.1441-0.1532*I,  0.1984+0.1737*I, -0.0034-0.1555*I,  0.2608+0.5382*I],
                         [ 0.2487-0.0926*I,  0.6253-0.3304*I,  0.2643+0.0194*I,  0.1002-0.0140*I],
                         [ 0.3758+0.0793*I, -0.0307+0.0816*I,  0.1266-0.1747*I, -0.4175+0.4058*I]].to_lm

  end

  %w(s d c z).each do |x|
    method = "#{x}gesvd"
    rc = LapackTest.get_rc(x)

    define_method("test_#{method}") do
      s, u, vt, work, info, a = NumRu::Lapack.send(method, "O", "S", @a[rc])
      assert_equal 0, info
      assert_narray @s_exp[rc], s, 1.0e-4

      vt.shape[0].times do |i|
        vt[i,true] *= -1 if comp_sign(vt[i,0], @vt_exp[rc][i,0])
      end
      assert_narray @vt_exp[rc], u, 1.0e-4
      a.shape[1].times do |i|
        a[true,i] *= -1 if comp_sign(a[0,i], @a_exp[rc][0,i])
      end
      assert_narray @a_exp[rc], a, 1.0e-4
    end

    define_method("test_#{method}_inquireing_lwork") do
      s, u, vt, work, info, = NumRu::Lapack.send(method, "O", "S", @a[rc], :lwork => -1)
      assert_equal 0, info
      lwork = get_int(work[0])
      s, u, vt, work, info, a = NumRu::Lapack.send(method, "O", "S", @a[rc], :lwork => lwork)
      assert_equal 0, info
      assert_equal lwork, get_int(work[0])
      assert_narray @s_exp[rc], s, 1.0e-4
      vt.shape[0].times do |i|
        vt[i,true] *= -1 if comp_sign(vt[i,0], @vt_exp[rc][i,0])
      end
      assert_narray @vt_exp[rc], vt, 1.0e-4
      a.shape[1].times do |i|
        a[true,i] *= -1 if comp_sign(a[0,i], @a_exp[rc][0,i])
      end
      assert_narray @a_exp[rc], a, 1.0e-4
    end

    define_method("test_#{method}_inquireing_lwork_oldargstyle") do
      s, u, vt, work, info, a = NumRu::Lapack.send(method, "O", "S", @a[rc], :lwork => -1)
      assert_equal 0, info
      lwork = get_int(work[0])
      s, u, vt, work, info, = NumRu::Lapack.send(method, "O", "S", @a[rc], -1)
      assert_equal 0, info
      assert_equal lwork, get_int(work[0])
    end

  end

end

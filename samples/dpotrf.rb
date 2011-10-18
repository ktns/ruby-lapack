require "numru/lapack"

uplo = "L"
p a = NArray[[ 4.16, -3.12,  0.56, -0.10],
             [-3.12,  5.03, -0.83,  1.18],
             [ 0.56, -0.83,  0.76,  0.34],
             [-0.10,  1.18,  0.34,  1.18]]
a_org = a.dup
info, a = NumRu::Lapack.dpotrf(uplo, a)

p info
p a

for i in 0...a.shape[0]
  for j in 0...i
    a[j,i] = 0.0
  end
end

a = NMatrix.ref(a)

if (NArray.ref(a.transpose * a) - a_org).abs.gt(1.0e-10).count_true == 0
  p "OK"
else
  p "NG"
end

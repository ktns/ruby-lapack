require "numru/lapack"

jobz = "V"
range = "I"
uplo = "U"
a = NArray[[1,1,0], [1,2,1], [0,1,1]]
vl = vu = 0 # not be used in this example
il = 1
iu = 3
abstol = 0.0
lwork = 78
liwork = 30

m, w, z, isuppz, work, iwork, info, a =
  NumRu::Lapack.dsyevr(jobz, range, uplo,
                       a, vl, vu, il, iu, abstol, :lwork => lwork, :liwork => liwork)

p m
p w
p z
p isuppz
p work
p iwork
p info
p a


# result test
eps = 1.0e-14

flag = m==3

flag &&= w[0].abs < eps
flag &&= (w[1]-1.0).abs < eps
flag &&= (w[2]-3.0).abs < eps

sqrt2 = Math.sqrt(1/2.0)
sqrt3 = Math.sqrt(1/3.0)
sqrt6 = Math.sqrt(1/6.0)
flag &&= (z[0,0]-sqrt3).abs < eps
flag &&= (z[1,0]+sqrt3).abs < eps
flag &&= (z[2,0]-sqrt3).abs < eps
flag &&= (z[0,1]+sqrt2).abs < eps
flag &&= z[1,1].abs < eps
flag &&= (z[2,1]-sqrt2).abs < eps
flag &&= (z[0,2]-sqrt6).abs < eps
flag &&= (z[1,2]-sqrt6*2).abs < eps
flag &&= (z[2,2]-sqrt6).abs < eps

flag &&= (isuppz == NArray[1,3,1,3,1,3])
flag &&= work[0].to_i == 78
flag &&= iwork[0] == 30
flag &&= info == 0

if flag
  print "OK\n"
else
  print "NG\n"
end

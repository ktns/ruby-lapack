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
                       a, vl, vu, il, iu, abstol, lwork, liwork)

p m
p w
p z
p isuppz
p work
p iwork
p info
p a




using ITensors


N = 6
s = siteinds(2,N)
chi = 4
psi = randomMPS(s;linkdims=chi)

println(psi)
orthogonalize(psi, 3)

println(psi)

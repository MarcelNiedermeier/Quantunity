
using ITensors


function getTensorElement(psi, element, N, sites)
    V = ITensor(1.)
    for i in 1:N
        global V *= (psi[i]*state(sites[i], el[i]))
    end
    v = scalar(V)
    return v
end


function getWavefunction(psi)

end


N = 2

sites = siteinds("Qubit", N)
#sites = siteinds("S=1/2", N)
psi = productMPS(sites, "0")

#println(sites)
#println(psi)


A = Complex.(zeros(2, 2))
A[1, 1] = 1.
sites = siteinds(2, 2)
psi = MPS(A, sites)
println("here is A: ", A)




"""
el = [1, 2, 2, 1, 1]

V = ITensor(1.)
for i in 1:N
    #println(psi[i]*state(sites[i], el[i]))
    global V *= (psi[i]*state(sites[i], el[i]))
end

v = scalar(V)
println("tensor element 0 ", v)


println(state(sites[1], 2))

#println("tensor element ", getTensorElement(psi, el, N, sites))
"""


s = ITensors.sample(orthogonalize!(psi, 1))
println(s)

# single-site gates

X = Complex.([0. 1.; 1. 0.])
X_gate = ITensor(X, sites[1], ITensors.prime(sites[1]))
X_gate2 = ITensor(X, sites[2], ITensors.prime(sites[2]))

newA = X_gate*psi[1]
noprime!(newA)
psi[1] = newA

newA = X_gate2*psi[2]
noprime!(newA)
psi[2] = newA


s = ITensors.sample(orthogonalize!(psi, 1))
println("psi after X: ", s)
println(maxlinkdim(psi))

# two-site gates
cnot = Complex.([1. 0. 0. 0.; 0. 1. 0. 0.; 0. 0. 0. 1.; 0. 0. 1. 0.])
swap = Complex.([1. 0. 0. 0.; 0. 0. 1. 0.; 0. 1. 0. 0.; 0. 0. 0. 1.])

CNOT_gate = ITensor(cnot, sites[1], ITensors.prime(sites[1]), sites[2], ITensors.prime(sites[2]))
#CNOT_gate = ITensor(cnot, sites[2], ITensors.prime(sites[2]), sites[3], ITensors.prime(sites[3]))
#SWAP_gate = ITensor(swap, sites[2], ITensors.prime(sites[2]), sites[3], ITensors.prime(sites[3]))

#println("check cnot: ", array(CNOT_gate))
@show reshape(array(CNOT_gate), 4, 4)

orthogonalize!(psi, 1)
wf = (psi[1]*psi[2])*CNOT_gate
noprime!(wf)
ind = uniqueinds(psi[1], psi[2])
U, S, V = svd(wf, ind, cutoff=1E-8)
psi[1] = U
psi[2] = S*V

#println(maxlinkdim(psi))
s = ITensors.sample(orthogonalize!(psi, 1))
#s = ITensors.sample(psi)
println(s)

"""
orthogonalize!(psi, 2)
wf = (psi[2]*psi[3])*SWAP_gate
noprime!(wf)
ind = uniqueinds(psi[2], psi[3])
U, S, V = svd(wf, ind, cutoff=1E-8)
psi[2] = U
psi[3] = S*V

s = ITensors.sample(orthogonalize!(psi, 1))
s_new = s - ones(Int32, length(s))
println(s)
println(s_new)

println(maxlinkdim(psi))
println("here ", typeof(psi))
"""

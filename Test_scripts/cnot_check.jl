using ITensors


N = 6

# set up state |00>
#A = Complex.(zeros(2, 2, 2))
#A[1, 1, 1] = 1.
#sites = siteinds(2, 3)
sites = siteinds("Qubit", N)
psi = productMPS(sites, "0")
#psi = MPS(A, sites)

# check state
s = ITensors.sample(orthogonalize!(psi, 1))
println("psi initialised: ", s)

# single-site gates
X = Complex.([0. 1.; 1. 0.])
X_gate = ITensor(X, sites[1], ITensors.prime(sites[1]))
X_gate2 = ITensor(X, sites[2], ITensors.prime(sites[2]))

# apply X gate to site 1
newA = X_gate*psi[1]
noprime!(newA)
psi[1] = newA

# apply X gate to site 2
newA = X_gate2*psi[2]
noprime!(newA)
psi[2] = newA

# check state
s = ITensors.sample(orthogonalize!(psi, 1))
println("psi after X gates: ", s)

# two-site gates
#cnot = Complex.([1. 0. 0. 0.; 0. 1. 0. 0.; 0. 0. 0. 1.; 0. 0. 1. 0.])
#CNOT_gate = ITensor(cnot, sites[1], sites[2], ITensors.prime(sites[1]), ITensors.prime(sites[2]))


# contruct CNOT
ampo = AutoMPO()
ampo += "ProjUp", 2
ampo += "ProjDn", 2, "X", 3

# construct SWAP
ampo2 = AutoMPO()
ampo2 += "ProjUp", 1, "ProjUp", 5
ampo2 += "S+", 1, "S-", 5
ampo2 += "S-", 1, "S+", 5
ampo2 += "ProjDn", 1, "ProjDn", 5

# construct TOFFOLI
ampo3 = AutoMPO()
ampo3 += "ProjUp", 2, "ProjUp", 3
ampo3 += "ProjUp", 2, "ProjDn", 3
ampo3 += "ProjDn", 2, "ProjUp", 3
ampo3 += "ProjDn", 2, "ProjDn", 3, "X", 6


CNOT_gate_MPO = MPO(ampo, sites)
SWAP_gate_MPO = MPO(ampo2, sites)
TOFFOLI_gate_MPO = MPO(ampo3, sites)


CNOT_psi = contract(CNOT_gate_MPO, psi)


"""
# apply to MPS
orthogonalize!(psi, 1)
wf = (psi[1]*psi[2])*CNOT_gate
noprime!(wf)
ind = uniqueinds(psi[1], psi[2])
U, S, V = svd(wf, ind, cutoff=1E-8)
psi[1] = U
psi[2] = S*V
"""


# check state
s = ITensors.sample(orthogonalize!(CNOT_psi, 1))
println("psi after cnot gate: ", s)


SWAP_CNOT_psi = contract(SWAP_gate_MPO, CNOT_psi)

# check state
s2 = ITensors.sample(orthogonalize!(SWAP_CNOT_psi, 1))
println("psi after swap gate: ", s2)

TOF_SWAP_CNOT_psi = contract(TOFFOLI_gate_MPO, SWAP_CNOT_psi)

# check state
s3 = ITensors.sample(orthogonalize!(TOF_SWAP_CNOT_psi, 1))
println("psi after Tof gate: ", s3)

println(TOF_SWAP_CNOT_psi)



# contruct CNOT
bondDims = []


N = 80

for i in 1:N-11
    local ampo = AutoMPO()
    ampo += "ProjUp", 1, "ProjUp", 5+i
    ampo += "ProjUp", 1, "ProjDn", 5+i
    ampo += "ProjDn", 1, "ProjUp", 5+i
    ampo += "ProjDn", 1, "ProjDn", 5+i, "X", 10+i
    local sites = siteinds("Qubit", N)
    local CNOT = MPO(ampo, sites)
    push!(bondDims, maxlinkdim(CNOT))
end

println(bondDims)

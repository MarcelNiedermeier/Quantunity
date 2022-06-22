
################
## Fidelity test
################


using LinearAlgebra
include("../QSim.jl")

# constants
N = 12
maxdim = 40
#backend = "ED_Julia"
backend = "MPS_ITensor"
lintop = false
contmethod = "naive"
random = true
randombond = 2

# declare sites, get two MPS
sites = siteinds("QCircuit", N)
qc1 = initialise_qcircuit(N, lintop, backend, maxdim, contmethod, random, randombond)
qc2 = initialise_qcircuit(N, lintop, backend, maxdim, contmethod, random, randombond)

MPS1 = deepcopy(qc1.StateVector)
wf1 = get_wavefunction(qc1)

# do some stuff with quantum circuit

# first qubit
#hadamard!(qc1, [1])
#CRn!(qc1, 5, [5, 1])
#CRn!(qc1, 4, [4, 1])
#CRn!(qc1, 3, [3, 1])
#CRn!(qc1, 2, [2, 1])
#
## second qubit
#hadamard!(qc1, [2])
#CRn!(qc1, 4, [5, 2])
#CRn!(qc1, 3, [4, 2])
#CRn!(qc1, 2, [3, 2])
#
## third qubit
#hadamard!(qc1, [3])
#CRn!(qc1, 3, [5, 3])
#CRn!(qc1, 2, [4, 3])
#
## fourth qubit
#hadamard!(qc1, [4])
#CRn!(qc1, 2, [5, 4])
#
## fifth qubit
#hadamard!(qc1, [5])
#
## final swaps
#fullSwap!(qc1, [1, 5])
#fullSwap!(qc1, [2, 4])


# do QFT and inverse QFT
println("QFT")
QFT!(qc1, 1, N)
println("invQFT")
invQFT!(qc1, 1, N)


# calculate fidelity via in-built function
fid_IT = abs(ITensors.dot(qc1.StateVector, MPS1))^2


# calculate fidelity via direct expensive contraction
wf1_after = get_wavefunction(qc1)
wf2 = get_wavefunction(qc2)
fid_ex = abs(dot(wf1, wf1_after))^2

println("in-built fidelity: $fid_IT")
println("exact fidelity: $fid_ex")

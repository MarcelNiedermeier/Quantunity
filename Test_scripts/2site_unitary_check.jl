
######################
## Test 2-site unitary
######################




###########################################
## Test for simple quantum phase estimation
###########################################

using DelimitedFiles
using Statistics
using Plots
include("../QSim.jl")


# set constants
n_prec = 10 # bits of precision
n_probs = [1, 2, 4, 6] # qubits needed to achieve certain success probability
#n_prob = 1
#N = 1 + n_prec + n_prob # total number of qubits in algorithm
#maxdims = [12, 16, 20, 24, 28, 32, 36, 40, 44, 48, 52, 56, 60, 64]# 64, 128]
#maxdims = [2, 3, 4, 6, 8, 10, 12, 14, 16, 20, 24, 28, 32, 36, 40]
maxdim = 2
#N_sample = 4
N_meas = 2000
#backend = "ED_Julia"
backend = "MPS_ITensor"
#initialisation = "eigenstate"
initialisation = "superposition"
contmethod = "naive"
random = false
lintop = false
randombond = 2






# eigenvalues to be estimated, tolerance of measurement function
#θ1 = 1/√3
θ1 = 1/√101
#θ1 = 1*1/2 + 1*1/4 + 1*1/8 + 0*1/16 + 1*1/32 + 1*1/64
#θ2 = 1/√5
θ2 = 1/√10
#θ2 = 0*1/2 + 1*1/4 + 0*1/8 + 1*1/16 + 1*1/32 + 0*1/64
θ3 = 1/√7
θ4 = 1/√5
eps = 0.0001



###############
# "direct" gate
###############


# intialise quantum circuit
qc = initialise_qcircuit(2, lintop, backend, maxdim, contmethod,
random, randombond)

# construct first Bell state
#hadamard!(qc, [1])
#cnot!(qc, [1, 2])

# construct second Bell state
PauliX!(qc, [1])
hadamard!(qc, [1])
cnot!(qc, [1, 2])

# construct third Bell state
#hadamard!(qc, [1])
#PauliX!(qc, [2])
#cnot!(qc, [1, 2])

# construct fourth Bell state
#hadamard!(qc, [1])
#PauliX!(qc, [2])
#PauliZ!(qc, [1, 2])
#cnot!(qc, [1, 2])

# evaluate exact wave function, multiply with eigenvalue
wf_before = get_wavefunction(qc)
println(wf_before*exp(1.0im*2π*θ2))

# check Bell state statistics
sample_measurement(qc, [1, 2], N_meas, eps, true, "ITensor", true)

# apply 2-site operator U
n = 1
U_2site!(qc, U_n_2site(θ1, θ2, θ3, θ4, n), [1, 2])

# compare with wavefunction after application of unitary operator
wf_after = get_wavefunction(qc)
println(wf_after)
println("check norm: ", norm(wf_before*exp(1.0im*2π*θ2) - wf_after))


sample_measurement(qc, [1, 2], N_meas, eps, true, "ITensor", true)

draw(qc)


#################
# controlled gate
#################

# intialise quantum circuit
qc2 = initialise_qcircuit(3, lintop, backend, maxdim, contmethod,
random, randombond)

# switch control qubit on (or not)
PauliX!(qc2, [1])

# construct first Bell state
#hadamard!(qc2, [2])
#cnot!(qc2, [2, 3])

# construct second Bell state
#PauliX!(qc2, [2])
#hadamard!(qc2, [2])
#cnot!(qc2, [2, 3])

# construct third Bell state
BellState!(qc2, 3, [2, 3])
#hadamard!(qc2, [2])
#PauliX!(qc2, [3])
#cnot!(qc2, [2, 3])

# construct fourth Bell state
#hadamard!(qc2, [2])
#PauliX!(qc2, [3])
#PauliZ!(qc2, [2, 3])
#cnot!(qc2, [2, 3])

# evaluate exact wave function, multiply with eigenvalue
wf_before = get_wavefunction(qc2)
println(wf_before*exp(1.0im*2π*θ3))

# check Bell state statistics
sample_measurement(qc2, [2, 3], N_meas, eps, true, "ITensor", true)

# apply 2-site operator U
n = 1
CU_2site!(qc2, U_n_2site(θ1, θ2, θ3, θ4, n), [1, 2, 3])

# compare with wavefunction after application of unitary operator
wf_after = get_wavefunction(qc2)
println(wf_after)
println("check norm: ", norm(wf_before*exp(1.0im*2π*θ3) - wf_after))


sample_measurement(qc2, [2, 3], N_meas, eps, true, "ITensor", true)

draw(qc2)

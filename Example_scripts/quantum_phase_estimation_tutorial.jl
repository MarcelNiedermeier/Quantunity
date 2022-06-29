

####################################
## Quantum Phase Estimation Tutorial
####################################


include("../QSim.jl")

# set constants
n_prec = 8 # bits of precision
n_prob = 4 # number determines success probability
maxdim = 12
N_meas = 5000
backend = "MPS_ITensor"
initialisation = "superposition" # or "-" or "+"
contmethod = "naive"
random = false
lintop = false
randombond = 2
eps = 0.0001


""" Calculate U^n in computational basis with dummy unitary map that has
eigenvalues θ1 and θ2 for the eigenvectors |+⟩, |-⟩. """
function U_n(θ1, θ2, n)
    H = 1/√2 * Complex.([1. 1.; 1. -1.])
    U_tmp = Complex.([exp(1.0im*2π*n*θ1) 0.; 0. exp(1.0im*2π*n*θ2)])
    return H*U_tmp*H
end


""" Calculate estimated phase in QPE algorithm, given a measurement result
of qubits. """
function recover_phase_estimate(phase_array::Array)
    N = length(phase_array)
    phase = 0.
    for j in 1:N
        phase += (1/2)^(j) * phase_array[j]
    end
    return phase
end


""" Function to find the measurement results with the k highest
probabilities. """
function get_highest_prob_measurement(qc, k::Int)

    # get measurement results
    probs = collect(values(qc.ClassicalBitsProportion))
    states = collect(keys(qc.ClassicalBitsProportion))

    # sort probabilities, find k highest
    probs_tmp = sort(probs, rev=true)
    max_probas = probs_tmp[1:k]

    # get positions of these probabilities in original list
    indices = []
    for i in 1:k
        push!(indices, findall(x->x==max_probas[i], probs))
    end

    # find measurement results corresponding to these probabilities
    measurements_max = []
    for i in 1:k
        push!(measurements_max, states[indices[i]][1])
    end

    return measurements_max, max_probas
end


# total number of qubits in algorithm
N = 1 + n_prec + n_prob

# define phases to be determined
θ1 = 1/√2
θ2 = 1/√5

# intialise quantum circuit
qc_MPS = initialise_qcircuit(N, lintop, "MPS_ITensor", maxdim, contmethod,
random, randombond)

# initialise phase estimation register
hadamard!(qc_MPS, [i for i in 1:(N-1)])

# set eigenstate |+⟩ as input state on last qubit
# alternatively: apply X gate before Hadamard to set to |-⟩
if initialisation == "+"
    hadamard!(qc_MPS, [N])
elseif initialisation == "-"
    PauliX!(qc, [N])
    hadamard!(qc_MPS, [N])
elseif initialisation == "superposition"
    println("Initialising to superposition of |+⟩ and |-⟩!")
else
    println("Chose wrong keywork for initialisation, please check.")
end

# do controlled rotations
n = 1
for i in (N-1):-1:1
    CU!(qc_MPS, U_n(θ1, θ2, n), [i, N])
    global n = 2*n
end

# do inverse QFT
invQFT!(qc_MPS, 1, N-1)

# measurement, "DirectSampling" or "SVDbased" or "ITensor"
sample_measurement(qc_MPS, [i for i in 1:(n_prec)], N_meas, eps, true, "ITensor", true)

# check what it looks like
draw(qc_MPS)

# get two measurement results with highest probabilities
# convert back to phases
meas_max, probs_max = get_highest_prob_measurement(qc_MPS, 2)
phases = []
for state in meas_max
    push!(phases, recover_phase_estimate(state))
end

# compare with input phases
for i in 1:length(phases)
    println("Phase(s) found: $(phases[i])")
end
println("Compare to original phases: θ1 = $(θ1), θ2 = $(θ2)")

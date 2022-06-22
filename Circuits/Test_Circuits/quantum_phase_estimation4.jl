




###########################################
## Test for simple quantum phase estimation
###########################################

using DelimitedFiles
using Statistics
using Plots
include("../QSim.jl")


# set constants
n_prec = 10 # bits of precision
#n_probs = [1, 2, 3, 4] # qubits needed to achieve certain success probability
n_prob = 3
#N = 1 + n_prec + n_prob # total number of qubits in algorithm
#maxdims = [12, 16, 20, 24, 28, 32, 36, 40, 44, 48, 52, 56, 60, 64]# 64, 128]
maxdims = [2, 4, 8, 16, 32, 64, 128]
#maxdim = 200
#N_sample = 4
N_meas = 1000
#backend = "ED_Julia"
backend = "MPS_ITensor"
#initialisation = "eigenstate"
initialisation = "superposition"
contmethod = "naive"
random = false
lintop = false
randombond = 2


""" Calculate U^n in computational basis with example unitary map that has
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
        #phase += phase_array[end-j+1] * 2.0^(-N-1+j)
        phase += (1/2)^(j) * phase_array[j]
        #println("$(phase_array[j]), $((1/2)^(j))")
    end
    return phase
end


""" Function to find the measurement result with the highest proportion
in the whole sample. Returns the bitstring and its corresponding
probability. """
function get_highest_prob_measurement(qc)
    p_max = maximum(keys(qc.ClassicalBitsProportion))
    return qc.ClassicalBitsProportion[p_max], p_max
end


function get_highest_prob_measurement(qc, k::Int)

    println(keys(qc.ClassicalBitsProportion)[1])
    keys_tmp = Array(keys(qc.ClassicalBitsProportion))
    sort!(keys_tmp, rev=true)

    max_probas = keys_temp[1:k]
    measurements_max = []

    for i in 1:k
        append!(measurements_max, qc.ClassicalBitsProportion[max_probas[i]])
    end

    return measurements_max, max_probas

end


""" Function to calculate the coefficients of a Fourier-transformed
state directly. """
function invQFT_ex(t, initial_state)

    coeffs = Complex.(zeros(2^t))
    ω = exp(-1.0im*2*π/(2.0^t))

    for j in 0:2^t-1
        b_j = Complex(0.0)
        for k in 0:(2^t-1)
            b_j += initial_state[k+1] * ω^(j*k)
        end
        coeffs[j+1] = b_j/√(2^t)
    end
    return coeffs
end

""" Function to evaluate the coefficients of the quantum phase estimation
exactly for the two qubit case. """
function phase_estimation_register_exact(θ, t)

    coeffs = [1, exp(2π*1.0im * θ * 2^(t-1))]
    for i in t-2:-1:0
        coeffs = kron(coeffs, [1, exp(2π*1.0im * θ * 2^i)])
    end

    return invQFT_ex(t, 1/2^(t/2) * coeffs)
end


""" Function to combine the exactly simulated register of a QPE with
the eigenstate that is input in the algorithm. """
function phase_estimation_state_exact(register, eigenstate, t)
    return kron(register, eigenstate)
end


function get_bitstrings(t)

    bitstrings = []
    for i in 0:2^t-1
        # need reverse?
        push!(bitstrings, reverse(digits(i, base=2, pad=t)))
    end

    return bitstrings
end


function get_measurement_histogram(qc, n_prec)

    probs_tmp = collect(keys(qc.ClassicalBitsProportion))
    states_tmp = collect(values(qc.ClassicalBitsProportion))

    #println(states_tmp)

    states = zeros(2^n_prec)
    probs = zeros(2^n_prec)
    bitstrings = get_bitstrings(n_prec)
    #println(bitstrings)

    for i in 1:2^n_prec
        states[i] = i
        if bitstrings[i] ∈ states_tmp
            ind = findfirst(item -> item == bitstrings[i], states_tmp)
            probs[i] = probs_tmp[ind]
        end
    end

    return Int.(states), probs


end



# eigenvalues to be estimated, tolerance of measurement function
θ1 = 1/√3
#θ1 = 1*1/2 + 1*1/4 + 1*1/8 + 0*1/16 + 1*1/32 + 1*1/64
θ2 = 1/√5
#θ2 = 0*1/2 + 1*1/4 + 0*1/8 + 1*1/16 + 1*1/32 + 0*1/64
eps = 0.0001


for maxdim in maxdims

    # total number of qubits in algorithm
    local N = 1 + n_prec + n_prob

    println("maxdim = $(maxdim)")

    # intialise quantum circuit
    qc = initialise_qcircuit(N, lintop, backend, maxdim, contmethod,
    random, randombond)

    # initialise phase estimation register
    hadamard!(qc, [i for i in 1:(N-1)])


    # leave last qubit as |0⟩ !
    # set eigenstate |+⟩ as input state on last qubit
    #PauliX!(qc, [N])
    if initialisation == "eigenstate"
        hadamard!(qc, [N])
    end

    # do controlled rotations
    local n = 1
    for i in (N-1):-1:1
        CU!(qc, U_n(θ1, θ2, n), [i, N])
        #global n = 2*n
        n = 2*n
    end

    # do inverse QFT
    # might need explicit renormalisation in QFT
    invQFT!(qc, 1, N-1)

    qc.StateVector = 1/norm(qc.StateVector) * qc.StateVector

    println(norm(qc.StateVector))

    # measurement
    if backend == "MPS_ITensor"
        sample_measurement(qc, [i for i in 1:(n_prec)], N_meas, eps, true, "DirectSampling", true)
    else
        sample_measurement(qc, [i for i in 1:(n_prec)], N_meas, eps, true, true)
    end


    states, probs = get_measurement_histogram(qc, n_prec)


    #p = bar(states, probs, yaxis=(:log10))
    p = bar(states, probs)
    xlabel!("states")
    ylabel!("frequency")
    if backend == "MPS_ITensor"
        if initialisation == "eigenstate"
            println("MPS eigenstate")
            title!("MPS, eigenstate, n_prec = $(n_prec), n_prob = $(n_prob), maxdim = $(maxdim)")
            savefig(p,"../Plots/Quantum_Phase_Estimation/QPE_freq_MPS_eigenstate_maxdim_$(maxdim)_nprec_$(n_prec)_nprob_$(n_prob).png")
        else # superposition
            println("MPS superpos")
            title!("MPS, superposition, n_prec = $(n_prec), n_prob = $(n_prob), maxdim = $(maxdim)")
            savefig(p,"../Plots/Quantum_Phase_Estimation/QPE_freq_MPS_superpos_maxdim_$(maxdim)_nprec_$(n_prec)_nprob_$(n_prob).png")
        end
    else # ED
        if initialisation == "eigenstate"
            println("ED eigenstate")
            title!("ED, eigenstate prepared, n_prec = $(n_prec), n_prob = $(n_prob)")
            savefig(p,"../Plots/Quantum_Phase_Estimation/QPE_freq_ED_eigenstate_nprec_$(n_prec)_nprob_$(n_prob).png")
        else # superposition
            println("ED superpos")
            title!("ED, superposition prepared, n_prec = $(n_prec), n_prob = $(n_prob)")
            savefig(p,"../Plots/Quantum_Phase_Estimation/QPE_freq_ED_superpos_nprec_$(n_prec)_nprob_$(n_prob).png")
        end
    end
end

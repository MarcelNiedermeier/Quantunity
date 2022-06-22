


###########################################
## Test for simple quantum phase estimation
###########################################

using DelimitedFiles
using Statistics
include("../QSim.jl")


# set constants
#n_prec = 7 # bits of precision
#n_probs = [2, 3, 4, 5, 6] # qubits needed to achieve certain success probability
#n_prob = 8
N = 2
#N = 1 + n_prec + n_prob # total number of qubits in algorithm
#maxdims = [12, 16, 20, 24, 28, 32, 36, 40, 44, 48, 52, 56, 60, 64]# 64, 128]
#maxdims = [2, 3, 4, 6, 8, 10, 12, 14, 16, 20, 24, 28, 32, 36, 40]
maxdim = 100
N_sample = 4
N_meas = 10000
backend = "ED_Julia"
#backend = "MPS_ITensor"
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


function Kitaev_post_processing(qc_c, qc_s)

    println(collect(values(qc_c.ClassicalBitsProportion)))
    println(collect(values(qc_s.ClassicalBitsProportion)))
    P0_c = [k for (k,v) in qc_c.ClassicalBitsProportion if v==[0]]
    P0_s = [k for (k,v) in qc_s.ClassicalBitsProportion if v==[0]]
    P1_c = [k for (k,v) in qc_c.ClassicalBitsProportion if v==[1]]
    P1_s = [k for (k,v) in qc_s.ClassicalBitsProportion if v==[1]]

    return atan((1 - 2*P0_s[1])/(2*P0_c[1] - 1)), atan((2*P1_s[1] - 1)/(1 - 2*P1_c[1]))

end




# eigenvalues to be estimated, tolerance of measurement function
θ1 = 1/√3
#θ1 = 1*1/2 + 1*1/4 + 1*1/8 + 0*1/16 + 1*1/32 + 1*1/64
#θ2 = 1/√5
θ2 = 0*1/2 + 1*1/4 + 0*1/8 + 1*1/16 + 1*1/32 + 0*1/64
eps = 0.0005



# total number of qubits in algorithm
#N = 1 + n_prec + n_prob


# intialise quantum circuits
qc_c = initialise_qcircuit(N, lintop, backend, maxdim, contmethod,
random, randombond)
qc_s = initialise_qcircuit(N, lintop, backend, maxdim, contmethod,
random, randombond)

# initialise phase estimation register
hadamard!(qc_c, [1])
hadamard!(qc_s, [1])

# prepare eigenstate |+⟩
hadamard!(qc_c, [2])
hadamard!(qc_s, [2])

# sine circuit requires S gate
SGate!(qc_s, [1])

# controlled unitary
CU!(qc_c, U_n(θ1, θ2, 1), [1, N])
CU!(qc_s, U_n(θ1, θ2, 1), [1, N])

# wrap up with hadamard
hadamard!(qc_c, [1])
hadamard!(qc_s, [1])

# sample measurements for both circuits
#sample_measurement(qc_c, [1], N_meas, eps, true, "SVDbased", true)
#sample_measurement(qc_s, [1], N_meas, eps, true, "SVDbased", true)
sample_measurement(qc_c, [1], N_meas, eps, true, true)
sample_measurement(qc_s, [1], N_meas, eps, true, true)

#println(qc_c.ClassicalBitsProportion)
#println(qc_s.ClassicalBitsProportion)

phase_estimate1, phase_estimate2  = Kitaev_post_processing(qc_c, qc_s)

println("estimated deviation from θ1 = $(1/2π * abs(θ1 - phase_estimate1)), $(1/2π * abs(2π*θ1 - phase_estimate2))")









# leave last qubit as |0⟩ !
# set eigenstate |+⟩ as input state on last qubit
#PauliX!(qc, [N])
#hadamard!(qc, [N])

# do controlled rotations
#n = 1
#for i in (N-1):-1:1
#    CU!(qc, U_n(θ1, θ2, n), [i, N])
#    global n = 2*n
#end

# do inverse QFT
#invQFT!(qc, 1, N-1)

# measurement
#sample_measurement(qc, [i for i in 1:(n_prec-1)], N_meas, eps, true, "SVDbased", true)
#sample_measurement(qc, [i for i in 1:(n_prec-1)], N_meas, eps, true, true)

#println(qc.ClassicalBitsProportion)
#highest_pob_state, p_max = get_highest_prob_measurement(qc)
#highest_pob_states, p_max = get_highest_prob_measurement(qc, 2)

# estimated phase
#phase_est1 = recover_phase_estimate(highest_pob_state[1])
#phase_est2 = recover_phase_estimate(highest_pob_state[2])
#println("comparison to θ1: $(abs(θ1 - phase_est1)), $(abs(θ1 - phase_est2))")
#println("comparison to θ2: $(abs(θ2 - phase_est1)), $(abs(θ1 - phase_est2))")




######################
## QPE success probs 2
######################

using DelimitedFiles
using Statistics
include("../QSim.jl")


# set constants
n_prec = 10 # bits of precision
n_probs = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10] # qubits needed to achieve certain success probability
#N = 1 + n_prec + n_prob # total number of qubits in algorithm
#maxdims = [12, 16, 20, 24, 28, 32, 36, 40, 44, 48, 52, 56, 60, 64]# 64, 128]
maxdims = [2, 8, 16, 24, 32, 40]
N_sample = 4
N_meas = 1000
#backend = "ED_Julia"
backend = "MPS_ITensor"
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


# eigenvalues to be estimated, tolerance of measurement function
θ1 = 1/√3
#θ1 = 1/2 + 1/4 + 1*1/8 + 1/16 + 1/32
θ2 = 1/√5
eps = 0.05


# collect data
av_success_probs = zeros(length(maxdims), length(n_probs))


# loop through different maximum bond dimensions
for k in 1:length(maxdims)

    println("doing dim = $(maxdims[k])")

    # total number of qubits in algorithm
    #N = 1 + n_prec + n_probs[k]

    # check for different lengths of qubit register
    for j in 1:length(n_probs)

        println("calculating n_prob = $(n_probs[j])")

        # total number of qubits in algorithm
        N = 1 + n_prec + n_probs[j]

        sucess_prob_temp = zeros(N_sample)

        # average over different samples
        for l in 1:N_sample

            println("sample $l")

            # intialise quantum circuit
            local qc = initialise_qcircuit(N, lintop, backend, maxdims[k], contmethod,
            random, randombond)

            # initialise phase estimation register
            hadamard!(qc, [i for i in 1:(N-1)])

            # set eigenstate |+⟩ as input state on last qubit
            hadamard!(qc, [N])

            # do controlled rotations
            local n = 1
            for i in (N-1):-1:1
                CU!(qc, U_n(θ1, θ2, n), [i, N])
                n = 2*n
            end

            # do inverse QFT
            invQFT!(qc, 1, N-1)

            # measurement
            sample_measurement(qc, [i for i in 1:(n_prec-1)], N_meas, eps, false, "SVDbased", true)

            #println(qc.ClassicalBitsProportion)
            local highest_pob_state, p_max = get_highest_prob_measurement(qc)
            #push!(sucess_prob_temp, p_max)
            sucess_prob_temp[l] = p_max
            #println(highest_pob_state)
            #println("success probability: ", p_max)

            #local estimated_phase = recover_phase_estimate(highest_pob_state)
            #println("true phase: $(round(θ1, digits=N-1)), approx phase: $(round(estimated_phase, digits=N-1))")

        end

        # average over samples
        av_success_probs[k, j] = Statistics.mean(sucess_prob_temp)

    end
end


println(av_success_probs)

# save output in datafile
datafile = zeros(length(maxdims)+3, length(n_probs))

datafile[1, 1] = n_prec
datafile[1, 2] = N_sample
datafile[1, 3] = length(maxdims)
datafile[2, 1:1+length(maxdims)-1] = maxdims[:]
datafile[3, :] = n_probs[:]

for i in 1:length(maxdims)
    datafile[3+i, :] = av_success_probs[i, :]
end

open("../Data/Quantum_Phase_Estimation/QPE_MPS_success_prob_2_n_prec_$(n_prec)_maxbond_$(maxdims[end])_samples_$(N_sample).csv", "w") do io
    writedlm(io, datafile, ", ")
end

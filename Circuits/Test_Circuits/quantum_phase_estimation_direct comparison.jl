
######################################
## Direct comparison QPE exact and MPS
######################################


using DelimitedFiles
using Statistics
using Plots
include("../QSim.jl")


# set constants
n_prec = 7 # bits of precision
n_probs = [1, 2, 3]#, 4] # qubits needed to achieve certain success probability
#n_prob = 1
#N = 1 + n_prec + n_prob # total number of qubits in algorithm
#maxdims = [12, 16, 20, 24, 28, 32, 36, 40, 44, 48, 52, 56, 60, 64]# 64, 128]
#maxdims = [2, 3, 4, 6, 8, 10, 12, 14, 16, 20, 24, 28, 32, 36, 40]
maxdim = 10
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





function get_bitstrings(t)

    bitstrings = []
    for i in 0:2^t-1
        # need reverse?
        push!(bitstrings, reverse(digits(i, base=2, pad=t)))
    end

    return bitstrings
end


""" Function to extract a histogram from the dictionary ClassicalBitsProportion,
where the different measurement results of the quantum circuit are saved. Returns
a list of the states (identified by their decimal number) and the corresponding
measurement frequencies. """
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
#θ1 = 1/√101
#θ1 = 1*1/2 + 1*1/4 + 1*1/8 + 0*1/16 + 1*1/32 + 1*1/64
#θ2 = 1/√5
#θ2 = 1/√10
θ2 = 0*1/2 + 1*1/4 + 0*1/8 + 1*1/16 + 1*1/32 + 0*1/64
eps = 0.0001


for n_prob in n_probs

    # total number of qubits in algorithm
    local N = 1 + n_prec + n_prob

    println("n_prob = $(n_prob)")

    # intialise quantum circuit
    qc_ED = initialise_qcircuit(N, lintop, "ED_Julia", maxdim, contmethod,
    random, randombond)
    qc_MPS = initialise_qcircuit(N, lintop, "MPS_ITensor", maxdim, contmethod,
    random, randombond)

    # initialise phase estimation register
    hadamard!(qc_ED, [i for i in 1:(N-1)])
    hadamard!(qc_MPS, [i for i in 1:(N-1)])


    # leave last qubit as |0⟩ !
    # set eigenstate |+⟩ as input state on last qubit
    #PauliX!(qc, [N])
    if initialisation == "eigenstate"
        hadamard!(qc_ED, [N])
        hadamard!(qc_MPS, [N])
    end

    # do controlled rotations
    local n = 1
    for i in (N-1):-1:1
        CU!(qc_ED, U_n(θ1, θ2, n), [i, N])
        CU!(qc_MPS, U_n(θ1, θ2, n), [i, N])
        #global n = 2*n
        n = 2*n
    end

    # do inverse QFT
    invQFT!(qc_ED, 1, N-1)
    invQFT!(qc_MPS, 1, N-1)

    # measurement, "DirectSampling" or "SVDbased" or "ITensor"
    sample_measurement(qc_ED, [i for i in 1:(n_prec)], N_meas, eps, true, true)
    sample_measurement(qc_MPS, [i for i in 1:(n_prec)], N_meas, eps, true, "ITensor", true)
    states_ED, probs_ED = get_measurement_histogram(qc_ED, n_prec)
    states_MPS, probs_MPS = get_measurement_histogram(qc_MPS, n_prec)

    # compare state vectors
    println("checking coefficients of ED and MPS statevectors")
    wf_ED = qc_ED.StateVector
    wf_MPS = get_wavefunction(qc_MPS)
    for p in 1:length(wf_ED)
        diff = abs(wf_ED[p])^2 - abs(wf_MPS[p])^2
        if abs(diff) > 0.0000000001
            println("discrepancy detected!")
        end
    end
    println("norm diff of ED and MPS: ", norm(wf_ED - wf_MPS))

    #println(collect(keys(qc.ClassicalBitsProportion)))
    #println(collect(values(qc.ClassicalBitsProportion)))

    #println(qc.ClassicalBitsProportion)
    #highest_pob_state, p_max = get_highest_prob_measurement(qc)
    #highest_pob_states, p_max = get_highest_prob_measurement(qc, 2)

    # estimated phase
    #phase_est1 = recover_phase_estimate(highest_pob_state)
    #phase_est2 = recover_phase_estimate(highest_pob_state[2])
    #println("comparison to θ1: $(abs(θ1 - phase_est1))")
    #println("comparison to θ1: $(abs(θ1 - phase_est1)), $(abs(θ1 - phase_est2))")
    #println("comparison to θ2: $(abs(θ2 - phase_est1)), $(abs(θ1 - phase_est2))")

    #t = n_prec + n_prob
    #eigenstate1 = 1/√2 * [1., 1.]
    #eigenstate2 = 1/√2 * [1., -1.]
    #eigenstate = [0., 1.]
    #register1 = phase_estimation_register_exact(θ1, t)
    #register2 = phase_estimation_register_exact(θ2, t)
    #phase_est_ex = recover_phase_estimate(register)
    #println("comparison ex to θ1: $(abs(θ1 - phase_est_ex))")

    #if initialisation == "eigenstate"
    #    exact_state = phase_estimation_state_exact(register1, eigenstate1, t)
    #else
    #    exact_state1 = phase_estimation_state_exact(register1, eigenstate1, t)
    #    exact_state2 = phase_estimation_state_exact(register2, eigenstate2, t)
    #    exact_state = 1/√2 * (exact_state1 + exact_state2)
    #end

    #println("exact norm ", norm(exact_state))
    #println("qc norm ", norm(qc.StateVector))
    #if backend == "MPS_ITensor"
    #    println("norm diff: ", norm(get_wavefunction(qc) - exact_state))
    #else
    #    println("norm diff: ", norm(qc.StateVector - exact_state))
    #end

    #println("checking differences in measurement results")
    #if backend == "MPS_ITensor"
    #    for p in 1:length(exact_state)
    #        diff = abs(abs.(exact_state[p]).^2 - abs.(get_wavefunction(qc)[p]).^2)
    #        if diff > 0.001
    #            println("here")
    #        end
    #    end
    #else
    #    for p in 1:length(exact_state)
    #        diff = abs(abs.(exact_state[p]).^2 - abs.(qc.StateVector[p]).^2)
    #        if diff > 0.001
    #            println("here")
    #        end
    #    end
    #end
    #println("exact wf squared", abs.(exact_state).^2)

    #p = bar(states, probs, yaxis=(:log10))
    p_ED = bar(states_ED, probs_ED)
    xlabel!("states")
    ylabel!("frequency")
    if initialisation == "eigenstate"
        title!("ED, eigenstate prepared, n_prec = $(n_prec), n_prob = $(n_prob)")
        savefig(p_ED,"../Plots/Quantum_Phase_Estimation/test_QPE_freq_ED_eigenstate_nprec_$(n_prec)_nprob_$(n_prob).png")
    else # superposition
        title!("ED, superposition prepared, n_prec = $(n_prec), n_prob = $(n_prob)")
        savefig(p_ED,"../Plots/Quantum_Phase_Estimation/test_QPE_freq_ED_superpos_nprec_$(n_prec)_nprob_$(n_prob).png")
    end

    p_MPS = bar(states_MPS, probs_MPS)
    xlabel!("states")
    ylabel!("frequency")
    if initialisation == "eigenstate"
        title!("MPS, eigenstate prepared, n_prec = $(n_prec), n_prob = $(n_prob)")
        savefig(p_MPS,"../Plots/Quantum_Phase_Estimation/test_QPE_freq_MPS_eigenstate_nprec_$(n_prec)_nprob_$(n_prob).png")
    else # superposition
        title!("MPS, superposition prepared, n_prec = $(n_prec), n_prob = $(n_prob)")
        savefig(p_MPS,"../Plots/Quantum_Phase_Estimation/test_QPE_freq_MPS_superpos_nprec_$(n_prec)_nprob_$(n_prob).png")
    end

    #if backend == "MPS_ITensor"
    #    if initialisation == "eigenstate"
    #        title!("MPS, eigenstate prepared, n_prec = $(n_prec), n_prob = $(n_prob)")
    #        savefig(p,"../Plots/Quantum_Phase_Estimation/QPE_freq_MPS_eigenstate_nprec_$(n_prec)_nprob_$(n_prob).png")
    #    else # superposition
    #        title!("MPS, superposition prepared, n_prec = $(n_prec), n_prob = $(n_prob)")
    #        savefig(p,"../Plots/Quantum_Phase_Estimation/test_QPE_freq_MPS_superpos_nprec_$(n_prec)_nprob_$(n_prob).png")
    #    end
    #else # ED
    #    if initialisation == "eigenstate"
    #        title!("ED, eigenstate prepared, n_prec = $(n_prec), n_prob = $(n_prob)")
    #        savefig(p,"../Plots/Quantum_Phase_Estimation/QPE_freq_ED_eigenstate_nprec_$(n_prec)_nprob_$(n_prob).png")
    #    else # superposition
    #        title!("ED, superposition prepared, n_prec = $(n_prec), n_prob = $(n_prob)")
    #        savefig(p,"../Plots/Quantum_Phase_Estimation/test_QPE_freq_ED_superpos_nprec_$(n_prec)_nprob_$(n_prob).png")
    #    end
    #end
end

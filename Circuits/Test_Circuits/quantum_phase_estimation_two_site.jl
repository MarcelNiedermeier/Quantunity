
###########################################################
## Test for simple quantum phase estimation for 2-site gate
###########################################################

using DelimitedFiles
using Statistics
using Plots
include("../QSim.jl")


# set constants
n_prec = 7 # bits of precision
n_probs = [1, 2, 3, 4]#, 4, 6] # qubits needed to achieve certain success probability
#n_prob = 1
#N = 1 + n_prec + n_prob # total number of qubits in algorithm
#maxdims = [12, 16, 20, 24, 28, 32, 36, 40, 44, 48, 52, 56, 60, 64]# 64, 128]
#maxdims = [2, 3, 4, 6, 8, 10, 12, 14, 16, 20, 24, 28, 32, 36, 40]
maxdim = 32
#N_sample = 4
N_meas = 2000
#backend = "ED_Julia"
backend = "MPS_ITensor"
initialisation = "eigenstate"
#initialisation = "superposition"
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


for n_prob in n_probs

    # total number of qubits in algorithm
    local N = 2 + n_prec + n_prob

    println("n_prob = $(n_prob)")

    # intialise quantum circuit
    qc = initialise_qcircuit(N, lintop, backend, maxdim, contmethod,
    random, randombond)

    # initialise phase estimation register
    hadamard!(qc, [i for i in 1:(N-2)])


    # leave last qubit as |0⟩ !
    # set eigenstate Bell state as input state on last qubit
    #PauliX!(qc, [N])
    if initialisation == "eigenstate"
        BellState!(qc, 1, [N-1, N])
        #hadamard!(qc, [N-1])
        #cnot!(qc, [N-1, N])
    else # superposition (all Bell states equal)
        hadamard!(qc, [N])
        #hadamard!(qc, [N-1, N])
    end

    U = diagm([θ1, θ2, θ3, θ4])
    println(U)

    # do controlled rotations
    local n = 1
    for i in (N-2):-1:1
        #CU_2site!(qc, U_n_2site(θ1, θ2, θ3, θ4, n), [i, N-1, N])
        CU_general!(qc, U_n_2site(θ1, θ2, θ3, θ4, n), [i], [N-1, N])
        #CU_general!(qc, U^n, [i], [N-1, N])
        #global n = 2*n
        n = 2*n
    end


    # CU_general!(qc::QC_IT_MPS, U, control_qubits, action_qubits, update_rep=true)


    # do inverse QFT
    invQFT!(qc, 1, N-2)

    # measurement, "DirectSampling" or "SVDbased"
    if backend == "MPS_ITensor"
        println("here")
        sample_measurement(qc, [i for i in 1:(n_prec)], N_meas, eps, true, "ITensor", true)
    else
        sample_measurement(qc, [i for i in 1:(n_prec)], N_meas, eps, true, true)
    end

    draw(qc)


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

    if backend == "MPS_ITensor"
        if initialisation == "eigenstate"
            title = "MPS, eigenstate prepared, n_prec = $(n_prec), n_prob = $(n_prob)"
            path = "../Plots/Quantum_Phase_Estimation/two_site/QPE_freq_MPS_eigenstate_nprec_$(n_prec)_nprob_$(n_prob)_maxdim_$(maxdim).png"
        else # superposition
            title = "MPS, superposition prepared, n_prec = $(n_prec), n_prob = $(n_prob)"
            path = "../Plots/Quantum_Phase_Estimation/two_site/QPE_freq_MPS_superpos_nprec_$(n_prec)_nprob_$(n_prob)_maxdim_$(maxdim).png"
        end
        make_histogram(qc, n_prec, title, path, maxdim)
    else # ED
        if initialisation == "eigenstate"
            title = "ED, eigenstate prepared, n_prec = $(n_prec), n_prob = $(n_prob)"
            path = "../Plots/Quantum_Phase_Estimation/two_site/QPE_freq_ED_eigenstate_nprec_$(n_prec)_nprob_$(n_prob).png"
        else # superposition
            title = "ED, superposition prepared, n_prec = $(n_prec), n_prob = $(n_prob)"
            path = "../Plots/Quantum_Phase_Estimation/two_site/QPE_freq_ED_superpos_nprec_$(n_prec)_nprob_$(n_prob).png"
        end
        make_histogram(qc, n_prec, title, path)
    end
end

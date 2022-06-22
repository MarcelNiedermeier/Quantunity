
#######################
## QPE Phase Transition
#######################

using DelimitedFiles
using Statistics
using Plots
include("../../QSim.jl")


# set constants
n_prec = 10 # bits of precision
#n_probs = [1, 2, 3, 4]#, 4, 6] # qubits needed to achieve certain success probability
n_probs = [1, 2]#, 4, 6] # qubits needed to achieve certain success probability
#n_prob = 4 # 2, 3, 4
#N_operator = [1, 2, 3, 4, 5, 6, 7, 8]
#N_operator = [4]
N_op = 4
#N_operator = [5]
#maxdims = [12, 16, 20, 24, 28, 32, 36, 40, 44, 48, 52, 56, 60, 64]# 64, 128]
#global maxdim = 10
maxdims = [i for i in (2^N_op - 1):(2^N_op + 1)]
#maxdims = [i for i in 1:(2^N_op + 10)]
N_meas = 5000
#backend = "ED_Julia"
backend = "MPS_ITensor"
#initialisation = "eigenstate"
initialisation = "superposition"
contmethod = "naive"
random = false
lintop = false
truncateQFT = false
randombond = 2
eps = 0.0001
ϵ = 0.1


###########
# Functions
###########


""" Function to build a diagonal (hermitian) matrix as
diag(1, 2, 3, ..., 2^N). """
function get_diagonal_matrix(N)

    n = 2^N
    entries = Complex.(ones(n))
    for i in 1:n
        #ϕ = (i-1)/n
        #entries[i] = exp(1.0im*2π*ϕ)
        #entries[i] = (i-1)/n * 2π
        entries[i] = i
    end
    return diagm(entries)
end


""" Function to rescale an operator A, such that its spectrum is mapped
from the interval [E0, Emax] to the interval [0, 2π]. Returns the
rescaled operator and bandwidth (Emax, E0). """
function rescale_operator(A, ϵ=0.1)

    # find bandwidth of spectrum
    E0 = eigmin(A)
    Emax = eigmax(A)
    bandwidth = (Emax, E0)

    id_mat = diagm(ones(size(A)[1]))

    # do rescaling of matrix
    #return (A - E0*diagm(ones(size(A)[1])))/(Emax - E0) * (2π - ϵ) .+ ϵ, bandwidth
    return (A - E0*id_mat)/(Emax - E0) * ((2π-ϵ) - ϵ) + ϵ*id_mat, bandwidth
end


""" Function to map the rescaled spectrum back to the original
interval. """
function inverse_rescaling_spectrum(y, bandwidth, ϵ=0.1)

    # get bandwidth
    E0 = bandwidth[2]
    Emax = bandwidth[1]

    #return (y - ϵ*ones(length(y)))/(2π - ϵ) * (Emax-E0) .+ E0
    return (y .- ϵ)/((2π-ϵ) - ϵ) * (Emax - E0) .+ E0
end



""" Function to get phase of complex number mapped to [0, 2π]. """
function get_phase(c)
    ϕ = angle(c)
    if ϕ > 0 # phase in [0, π]
        return ϕ
    else # phase in [π, 2π]
        return 2π + ϕ
    end
end



""" Function to extract a histogram from the dictionary ClassicalBitsProportion,
where the different measurement results of the quantum circuit are saved. Returns
a list of the states (identified by their decimal number) and the corresponding
measurement frequencies. """
function get_measurement_histogram_phases(qc, t)

    # get measurement results, set up containers
    states_tmp = collect(keys(qc.ClassicalBitsProportion))
    probs_tmp = collect(values(qc.ClassicalBitsProportion))
    phases = zeros(2^t)
    probs = zeros(2^t)
    bitstrings = get_bitstrings(t)

    for i in 1:2^t
        phases[i] = i * 1/2^t

        # if a bitstring has been measured: look up corresponding probability
        if bitstrings[i] ∈ states_tmp
            ind = findfirst(item -> item == bitstrings[i], states_tmp)
            probs[i] = probs_tmp[ind]
        end
    end
    return phases, probs
end




""" Function to create a simple histogram with measurement results from
a given quantum circuit object. """
function make_histogram_phases(qc::QC_IT_MPS, t, title::String, path::String, maxdim, bandwidth, ϵ=0.1)

    # get measurement results and corresponding frequencies
    phases_res, freqs = get_measurement_histogram_phases(qc, t)

    # rescale phases to original bandwidth, need phases in interval [0, 2π]
    phases = inverse_rescaling_spectrum(2π*phases_res, bandwidth, ϵ)

    # make (rudimentary) plot; save
    p = bar(phases, freqs, label="measurement")
    xlabel!("spectrum")
    ylabel!("frequency")
    title!(title)
    annotate!((100, 0.2, "χ = $(maxdim)"))
    savefig(p, path)

end


##################
# quantum circuits
##################

# check for diffent success probabilities
for k in 1:length(n_probs)

    println("n_prob = $(n_probs[k])")

    # loop bond dimensions "over phase transition"
    for l in 1:length(maxdims)

        # total number of qubits in algorithm
        local N = N_op + n_prec + n_probs[k]

        println("maxdim = $(maxdims[l])")


        # intialise quantum circuit
        qc = initialise_qcircuit(N, lintop, backend, maxdims[l], contmethod,
        random, randombond)


        # leave last qubit as |0⟩ !
        # set eigenstate Bell state as input state on last qubit
        #PauliX!(qc, [N])
        if initialisation == "eigenstate"
            #BellState!(qc, 1, [N-1, N])
            #hadamard!(qc, [N-1])
            #cnot!(qc, [N-1, N])
        else # superposition (all Bell states equal)
            hadamard!(qc, [j for j in (n_prec+n_probs[k]+1):N])
        end


        # construct diagonal Hamiltonian
        operator = get_diagonal_matrix(N_op)

        # rescale operator
        operator_res, bandwidth = rescale_operator(operator)

        # construct unitary operator
        U_operator_res = exp(1.0im*operator_res)

        # do quantum phase estimation
        QPE!(qc, U_operator_res, 1, n_prec+n_probs[k])
        invQPE!(qc, U_operator_res, 1, n_prec+n_probs[k])


        # measurement, "DirectSampling" or "SVDbased"
        if backend == "MPS_ITensor"
            sample_measurement(qc, [i for i in 1:(n_prec)], N_meas, eps, true, "ITensor", true)
        else
            sample_measurement(qc, [i for i in 1:(n_prec)], N_meas, eps, true, true)
        end

        draw(qc)

        phases_MPS_res, probs_MPS = get_measurement_histogram_phases(qc, n_prec)
        phases = inverse_rescaling_spectrum(2π*phases_MPS_res, bandwidth, ϵ)

        #states_tmp = collect(keys(qc.ClassicalBitsProportion))
        #probs_tmp = collect(values(qc.ClassicalBitsProportion))
        println("probs = ", probs_MPS)
        println("phases = ", phases)
        println("length probs = ", length(probs_MPS))



        # get histogram

        # calculate deviation from even distribution

        # save for given bond dimension and size of operator



        #if backend == "MPS_ITensor"
        #    if truncateQFT
        #        title = "MPS, $(n_H)-qubit U, n_prec = $(n_prec), n_prob = $(n_prob)"
        #        path = "../Plots/Quantum_Phase_Estimation/diagonal_matrix/poster/test_QPE_freq_MPS_n_H_$(n_H)_nprec_$(n_prec)_nprob_$(n_prob)_maxdim_$(maxdim)_QFTdim_$(maxdim2).png"
        #    else # QFT at same bond dim as rest
        #        title = "MPS, $(n_H)-qubit U, n_prec = $(n_prec), n_prob = $(n_prob)"
        #        path = "../Plots/Quantum_Phase_Estimation/diagonal_matrix/poster/test_QPE_freq_MPS_n_H_$(n_H)_nprec_$(n_prec)_nprob_$(n_prob)_maxdim_$(maxdim).png"
        #    end
        #    #make_histogram(qc, n_prec, title, path, maxdim)
        #    #println("bandwidth: ", bandwidth)
        #    make_histogram_phases(qc, n_prec, title, path, maxdim, bandwidth)
        #else # ED
        #    if initialisation == "eigenstate"
        #        title = "ED, eigenstate prepared, n_prec = $(n_prec), n_prob = $(n_prob)"
        #        path = "../Plots/Quantum_Phase_Estimation/diagonal_matrix/QPE_freq_ED_eigenstate_nprec_$(n_prec)_nprob_$(n_prob).png"
        #    else # superposition
        #        title = "ED, superposition prepared, n_prec = $(n_prec), n_prob = $(n_prob)"
        #        path = "../Plots/Quantum_Phase_Estimation/two_sdiagonal_matrixite/QPE_freq_ED_superpos_nprec_$(n_prec)_nprob_$(n_prob).png"
        #    end
        #    make_histogram(qc, n_prec, title, path)
        #end

    end
end

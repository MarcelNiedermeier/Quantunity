
######################
## QPE subroutine test
######################




#####################################
## QPE of multisite diagonal operator
#####################################

using DelimitedFiles
using Statistics
using Plots
include("../QSim.jl")


# set constants
n_prec = 8 # bits of precision
#n_probs = [1, 2, 3, 4]#, 4, 6] # qubits needed to achieve certain success probability
n_prob = 4 # 2, 3, 4
#N_operator = [1, 2, 3, 4, 5, 6, 7, 8]
N_operator = [5]
#N_operator = [5]
#maxdims = [12, 16, 20, 24, 28, 32, 36, 40, 44, 48, 52, 56, 60, 64]# 64, 128]
global maxdim = 150
#global maxdim2 = 2 # [2, 4, 8, 16, 32]
N_meas = 5000
#backend = "ED_Julia"
backend = "MPS_ITensor"
initialisation = "eigenstate"
#initialisation = "superposition"
contmethod = "naive"
random = true
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


for n_H in N_operator

    # total number of qubits in algorithm
    local N = n_H + n_prec + n_prob

    println("N_operator = $(n_H)")

    #if n_H == 6
    #    global maxdim = 64
    #elseif n_H == 7
    #    global maxdim = 128
    #elseif n_H == 8
    #    global maxdim = 160
    #end


    # intialise quantum circuit
    qc = initialise_qcircuit(N, lintop, backend, maxdim, contmethod,
    random, randombond)

    # initialise phase estimation register
    #hadamard!(qc, [i for i in 1:(N-n_H)])


    # leave last qubit as |0⟩ !
    # set eigenstate Bell state as input state on last qubit
    #PauliX!(qc, [N])
    if initialisation == "eigenstate"
        #BellState!(qc, 1, [N-1, N])
        #hadamard!(qc, [N-1])
        #cnot!(qc, [N-1, N])
    else # superposition (all Bell states equal)
        hadamard!(qc, [j for j in (n_prec+n_prob+1):N])
    end


    # construct diagonal Hamiltonian
    operator = get_diagonal_matrix(n_H)

    # rescale operator
    operator_res, bandwidth = rescale_operator(operator)

    #println("operator res: ", operator_res)

    # construct unitary operator
    U_operator_res = exp(1.0im*operator_res)
    U_operator_res_conj = exp(-1.0im*operator_res)

    initial_state = deepcopy(qc.StateVector)

    #hadamard!(qc, [i for i in 1:N])
    QPE!(qc, U_operator_res, 1, n_prec+n_prob, false)
    invQPE!(qc, U_operator_res_conj, 1, n_prec+n_prob, false)

    final_state = qc.StateVector
    norm_diff = norm(initial_state - final_state)
    println("norm diff = $(norm_diff)")

    #PauliX!(qc, [i for i in 1:N])
    #QFT!(qc, 1, N)
    #invQFT!(qc, 1, N)

    #CU_general!(qc, U_operator_res^n, [2], [j for j in (n_prec+n_prob+1):
    #    (N)])
    #CU_general!(qc, U_operator_res_conj^n, [2], [j for j in (n_prec+n_prob+1):
    #    (N)])




    # find exact spectrum of rescaled propagator
    #spectrum_res = eigvals(U_operator_res)
    #phases_res = get_phase.(spectrum_res)

    #println("phases res: ", phases_res)

    # invert rescaling of spectrum and get spectrum of initial operator
    #phases = inverse_rescaling_spectrum(phases_res, bandwidth)

    #println("phases  ", phases)

    #println(get_phase.(spectrum))
    #h = histogram(phases, bins=2^n_prec)
    #savefig(h, "../Plots/Quantum_Phase_Estimation/diagonal_matrix/exact_histogram/QPE_freq_exact_n_H_$(n_H)_nprec_$(n_prec)_nprob_$(n_prob).png")

    # do controlled rotations
    #local n = 1
    #for i in (N-n_H):-1:1

        # get exponentiated Hamiltonian
        #U = exp(1.0im*H*n)
        #U = prop

        #CU_2site!(qc, U_n_2site(θ1, θ2, θ3, θ4, n), [i, N-1, N])
        #CU_general!(qc, U_n_2site(θ1, θ2, θ3, θ4, n), [i], [N-1, N])

        # construct controlled operator
        #CU_general!(qc, U_operator_res^n, [i], [j for j in (n_prec+n_prob+1):N])
        #global n = 2*n
        #n = 2*n
    #end

    # if desired, truncate bond dim before Fourier transform
    #if truncateQFT
    #    qc.MaxBondDim = maxdim2
    #end

    #println("bond dim now: ", qc.MaxBondDim)

    # do inverse QFT
    #invQFT!(qc, 1, N-n_H)

    # measurement, "DirectSampling" or "SVDbased"
    if backend == "MPS_ITensor"
        sample_measurement(qc, [i for i in 1:(N)], N_meas, eps, true, "ITensor", true)
    else
        sample_measurement(qc, [i for i in 1:(n_prec)], N_meas, eps, true, true)
    end

    draw(qc)

    #phases_MPS_res, probs_MPS = get_measurement_histogram_phases(qc, n_prec)

    #states_tmp = collect(keys(qc.ClassicalBitsProportion))
    #probs_tmp = collect(values(qc.ClassicalBitsProportion))
    #println("probs = ", probs_MPS)
    #println("phases_res = ", phases_MPS_res)



    #if backend == "MPS_ITensor"
    #    if truncateQFT
    #        title = "MPS, $(n_H)-qubit U, n_prec = $(n_prec), n_prob = $(n_prob)"
    #        path = "../Plots/Quantum_Phase_Estimation/diagonal_matrix/poster/QPE_freq_MPS_n_H_$(n_H)_nprec_$(n_prec)_nprob_$(n_prob)_maxdim_$(maxdim)_QFTdim_$(maxdim2).png"
    #    else # QFT at same bond dim as rest
    #        title = "MPS, $(n_H)-qubit U, n_prec = $(n_prec), n_prob = $(n_prob)"
    #        path = "../Plots/Quantum_Phase_Estimation/diagonal_matrix/poster/QPE_freq_MPS_n_H_$(n_H)_nprec_$(n_prec)_nprob_$(n_prob)_maxdim_$(maxdim).png"
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

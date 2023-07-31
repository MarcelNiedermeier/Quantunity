
###########################
## Random Circuit MPS vs DM
###########################


#include("/scratch/work/niederm1/apps/QuantumSimulator/QSim.jl")
include("../../Functions/density_matrix_simulator.jl")
include("../../QSim.jl")

#using DelimitedFiles



function random_circuit(ind)


    # hard parameters
    bitstring = 200
    #N_sample = 100
    p1 = 0
    γ = 0

    # get parameter space
    N_qubits = [14]#, 12, 14, 16, 18, 20]
    Ds = [24]#, 20, 40, 80]
    bond_dims = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130]
    #bond_dims = [10]
    p_depols = [10 .^range(-4, stop=0, length=13)]
    N_samples = [i for i in 1:150]
    parameter_space = collect(Iterators.product(N_qubits, Ds, bond_dims, p_depols, N_samples))

    ps = parameter_space[ind]
    N_qubit = ps[1]
    D = ps[2]
    bond_dim = ps[3]
    p_depol = ps[4]
    N_sample = ps[5]

    #bitstring_probs_MPS = zeros(N_sample)
    #bitstring_probs_DM = zeros(N_sample)

    #for sample in 1:N_sample

    #println("Doing sample $sample")

    # simulate random circuit
    qc_MPS = initialise_qcircuit(N_qubit, "MPS_ITensor", maxdim=bond_dim)
    qc_DM = initialise_dmcircuit(N_qubit)
    randomCircuit!(qc_MPS, 1, N_qubit, D)
    randomCircuit!(qc_DM, 1, N_qubit, D, p1, γ, p_depol)

    # calculate probability of bitstring
    #bitstring_probs_MPS[sample] = abs(get_bitstring_coefficient(qc_MPS, bitstring))^2
    #bitstring_probs_DM[sample] = get_bitstring_probability(qc_DM, bitstring)
    bitstring_prob_MPS = abs(get_bitstring_coefficient(qc_MPS, bitstring))^2
    bitstring_prob_DM = get_bitstring_probability(qc_DM, bitstring)

    #end

    #println("bitstring prob MPS: ", bitstring_probs_MPS)
    #println("bitstring prob DM: ", bitstring_probs_DM)

    # get cumulative distribution
    #bitstring_probs_cum_MPS = cum_bitstring_dist(bitstring_probs_MPS, N_qubit)
    #bitstring_probs_cum_DM = cum_bitstring_dist(bitstring_probs_DM, N_qubit)
    #println("check probs sum MPS: ", bitstring_probs_cum_MPS[end])
    #println("check probs sum DM: ", bitstring_probs_cum_DM[end])

    #cum_distributions[i, :] = bitstring_probs_cum

    # datafile for output
    #datafile = zeros(2, N_sample)
    #datafile[1, :] = bitstring_probs_cum_MPS
    #datafile[2, :] = bitstring_probs_cum_DM

    # save data
    #open("random_circuit_MPS_DM_test.csv", "w") do io
    #    writedlm(io, datafile, ", ")
    #end

    return (bitstring, p1, γ, N_qubit, D, bond_dim, p_depol, N_sample, bitstring_prob_MPS, bitstring_prob_DM)
end


#@time random_circuit(1)

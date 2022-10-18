
###################################
## Random Quantum Circuits Tutorial
###################################


# check convergence of cumulative distribution of bitstring
# probabilities to Porter-Thomas distribution

#using DelimitedFiles
include("../QSim.jl")


""" Function to apply a layer of CNOT gates on neighbouring qubit lines
to a quantum circuit. Can indicate the starting position. """
function cnot_layer!(qc, start_pos)
    for i in start_pos:2:(N-1)
        cnot!(qc, [i, i+1])
    end
end




""" Cumulative Porter-Thomas distribution, as given in PRX "What limits
the simulation of quantum computers?" paper. Describes distribution of
bitstring probabilities sampled over random unitary circuits (assuming
"perfect" sampling of unitaries according to the Haar measure). """
function porter_thomas(x, N)
    out = zeros(length(x))
    for i in 1:length(x)
        out[i] = 1 - (1 - x[i])^(2^N - 1)
    end
    return out
end


""" Porter-Thomas distribution as given in Sean Mullane's article. """
function porter_thomas2(x, N)
    out = zeros(length(x))
    for i in 1:length(x)
        out[i] = N*exp(-N*x[i])
    end
    return out
end


""" Function to calculate the discrete cumulative distribution of
a set of bitstring probabilities. The bitstring probabilities are
statistically distributed around 1/2^N. Here, they are automatically
rescaled by 2^N, such that they are distributed around 1. The cumulative
probabilities are calculated in the interval [0, max_prob]. """
function cum_bitstring_dist(bitstring_probs, N, max_prob=7)

    # get number of qubits, rescale bitstring probabilities
    bitstring_probs = 2^N * bitstring_probs

    # get x-range
    N_samples = length(bitstring_probs)
    x_prob = LinRange(0, max_prob, N_samples)

    # get cumulative distribution
    bitstring_probs_cum = zeros(N_samples)
    for i in 1:N_samples
        bitstring_probs_cum[i] = length(bitstring_probs[bitstring_probs .<= x_prob[i]])/N_samples
    end

    return bitstring_probs_cum
end


backend = "MPS_ITensor"
#backend = "ED_Julia"

# constants
N = 8
pos = 1
num = 8
lintop = false
maxdim = 32
bitstring = 200
N_samples = 2

# different circuit depths
#depths = [100]#, 24, 48, 100]
D = 6

# save cumulative distributions
#cum_distributions = zeros(length(depths), N_samples)

# x for cumulative distribution
max_prob = 10
x_prob = LinRange(0, max_prob, N_samples)


# initialise circuit object
qc = initialise_qcircuit(N, lintop, backend, maxdim)

# build random circuit
#for j in 1:D

    # apply random unitaries to each site
#    random_single_site_gates!(qc, [i for i in 1:N])

    # apply vertically stacked CNOT gates
#    if isodd(j)
#        cnot_layer!(qc, 1)
#    else
#        cnot_layer!(qc, 2)
#    end
#end

randomCircuit!(qc, pos, num, D, true)

draw(qc)


"""
# investigate different circuit depths
for i in 1:length(depths)

    # save results for given bitstring
    bitstring_probs = zeros(N_samples)

    println("Now depth $(depths[i])")


    # sample N_samples different random quantum circuits
    for sample in 1:N_samples

        if sample%10 == 0
            println("Calculating sample $sample")
        end

        # initialise circuit object
        local qc = initialise_qcircuit(N, lintop, backend, maxdim)

        # build random circuit
        for j in 1:(depths[i]÷2)

            # apply random unitaries to each site
            random_single_site_gates!(qc, [i for i in 1:N])

            # apply vertically stacked CNOT gates
            if isodd(j)
                sequential_cnot!(qc, [1])
            else
                sequential_cnot!(qc, [2])
            end

            # renormalise every four layers
            if j%4 == 0
                qc.StateVector = (1/real(norm(qc.StateVector)))*qc.StateVector
            end
        end

        # calculate probability of bitstring
        coeff = get_bitstring_coefficient(qc, bitstring)
        bitstring_probs[sample] = abs(coeff)^2

        #println("Norm: ", norm(qc.StateVector))
        println("2-Fid: ", qc.TwoQubitFidelity)
        #println("N-Fid: ", qc.NQubitFidelity[end])
        #println("Av-2-Fid: ", qc.AverageTwoQubitFidelity[end])

    end


    # get cumulative distribution
    bitstring_probs_cum = cum_bitstring_dist(bitstring_probs, N)
    println("check probs sum: ", bitstring_probs_cum[end])
    cum_distributions[i, :] = bitstring_probs_cum

end

# get Porter Thomas
PT_prob = porter_thomas(2.0^(-N)*x_prob, N)

# datafile for output
datafile = zeros(length(depths)+2, N_samples)

# general information ("header")
datafile[1, 1] = N
datafile[1, 2] = maxdim
datafile[1, 3] = bitstring
datafile[1, 4] = length(depths)
datafile[1, 5] = N_samples
for i in 1:length(depths)
    datafile[1, 5+i] = depths[i]
end

# data
datafile[2, :] = PT_prob
for i in 1:length(depths)
    datafile[2+i, :] = cum_distributions[i, :]
end

# save data
open("Data/Random_circuits/random_circuit_N_$(N)_bond_$(maxdim)_samples_$(N_samples).csv", "w") do io
    writedlm(io, datafile, ", ")
end




# plotting
#plot(x_prob, PT_prob, label="Porter-Thomas")
#
#for i in 1:length(depths)
#    plot!(x_prob, cum_distributions[i, :], label="N = $(N), χ = $(maxdim), D = $(depths[i])")
#end
#xlabel!("2^N × ρ")
#ylabel!("P(p_U(x) <= ρ)")
#savefig("Plots/random_circuit_N_$(N)_bond_$(maxdim)_samples_$(N_samples).png")

"""

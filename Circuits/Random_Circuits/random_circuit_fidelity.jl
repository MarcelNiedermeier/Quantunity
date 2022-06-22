
#######################################
## Investigate a random quantum circuit
#######################################

# check fidelities for fixed bond dimension and depth, but variable
# qubit numbers

using DelimitedFiles
using Statistics
include("../QSim.jl")


backend = "MPS_ITensor"
#backend = "ED_Julia"

# constants
Ns = [20, 40]#, 40, 60]
lintop = false
maxdim = 32 # 16, 32, 64
bitstring = 200
N_samples = 20
depth = 100
calculate_ex = false # compare with exactly calculated fidelities
random = false
randombond = 32
contmethod = "naive"


# save data for different qubit numbers
twofidelities = zeros(length(Ns), depth÷2+1)
Nfidelities = zeros(length(Ns), depth÷2+1)
Avtwofidelities = zeros(length(Ns), depth÷2+1)

if calculate_ex
    twofidelities_ex = zeros(length(Ns), depth÷2+1)
    Nfidelities_ex = zeros(length(Ns), depth÷2+1)
    Avtwofidelities_ex = zeros(length(Ns), depth÷2+1)
end

twofidelities_std = zeros(length(Ns), depth÷2+1)
Nfidelities_std = zeros(length(Ns), depth÷2+1)
Avtwofidelities_std = zeros(length(Ns), depth÷2+1)

# investigate different numbers of qubits
for i in 1:length(Ns)

    # save results for fidelities per sample
    twofidelities_tmp = zeros(N_samples, depth÷2+1)
    Nfidelities_tmp = zeros(N_samples, depth÷2+1)
    Avtwofidelities_tmp = zeros(N_samples, depth÷2+1)

    if calculate_ex
        twofidelities_ex_tmp = zeros(N_samples, depth÷2+1)
        Nfidelities_ex_tmp = zeros(N_samples, depth÷2+1)
        Avtwofidelities_ex_tmp = zeros(N_samples, depth÷2+1)
    end

    println("Now N =  $(Ns[i])")

    # sample N_samples different random quantum circuits
    for samp in 1:N_samples

        if samp%4 == 0
            println("Calculating sample $samp")
        end

        # initialise circuit object
        local qc = initialise_qcircuit(Ns[i], lintop, backend, maxdim, contmethod, random, randombond)

        # build random circuit until desired depth
        for j in 1:(depth÷2)

            #if j%20 == 0
            #    println("applying layers $j and $(j+1)")
            #end

            # apply random unitaries to each site
            random_single_site_gates!(qc, [i for i in 1:Ns[i]])

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

        # record fidelities
        twofidelities_tmp[samp, :] = qc.TwoQubitFidelitySelected
        Nfidelities_tmp[samp, :] = qc.NQubitFidelitySelected
        Avtwofidelities_tmp[samp, :] = qc.AverageTwoQubitFidelitySelected

        if calculate_ex
            twofidelities_ex_tmp[samp, :] = qc.TwoQubitFidelityExactSelected
            Nfidelities_ex_tmp[samp, :] = qc.NQubitFidelityExactSelected
            Avtwofidelities_ex_tmp[samp, :] = qc.AverageTwoQubitFidelityExactSelected
        end

    end

    # average recorded fidelities over samples and find std
    twofidelities[i, :] = Statistics.mean(twofidelities_tmp, dims=1)
    Nfidelities[i, :] = Statistics.mean(Nfidelities_tmp, dims=1)
    Avtwofidelities[i, :] = Statistics.mean(Avtwofidelities_tmp, dims=1)

    if calculate_ex
        twofidelities_ex[i, :] = Statistics.mean(twofidelities_ex_tmp, dims=1)
        Nfidelities_ex[i, :] = Statistics.mean(Nfidelities_ex_tmp, dims=1)
        Avtwofidelities_ex[i, :] = Statistics.mean(Avtwofidelities_ex_tmp, dims=1)
    end

    twofidelities_std[i, :] = std(twofidelities_tmp, dims=1)
    Nfidelities_std[i, :] = std(Nfidelities_tmp, dims=1)
    Avtwofidelities_std[i, :] = std(Avtwofidelities_tmp, dims=1)

end


# datafile for output
datafile = zeros(6*length(Ns)+1, depth÷2+1)
if calculate_ex
    datafile_ex = zeros(3*length(Ns), depth÷2+1)
end

# general information ("header")
datafile[1, 1] = N_samples
datafile[1, 2] = maxdim
datafile[1, 3] = bitstring
datafile[1, 4] = length(Ns)
datafile[1, 5] = depth
for i in 1:length(Ns)
    datafile[1, 5+i] = Ns[i]
end

# data
for i in 1:length(Ns)
    datafile[1+i, :] = twofidelities[i, :]
    datafile[1+length(Ns)+i, :] = Nfidelities[i, :]
    datafile[1+2*length(Ns)+i, :] = Avtwofidelities[i, :]
    datafile[1+3*length(Ns)+i, :] = twofidelities_std[i, :]
    datafile[1+4*length(Ns)+i, :] = Nfidelities_std[i, :]
    datafile[1+5*length(Ns)+i, :] = Avtwofidelities_std[i, :]
end

# save data for exact fidelities
if calculate_ex
    for i in 1:length(Ns)
        datafile_ex[i, :] = twofidelities_ex[i, :]
        datafile_ex[length(Ns)+i, :] = Nfidelities_ex[i, :]
        datafile_ex[2*length(Ns)+i, :] = Avtwofidelities_ex[i, :]
    end
end

# save data
if random
    open("Data/Random_circuits/random_circuit_rand_MPS_fid_D_$(depth)_bond_$(maxdim)_samples_$(N_samples).csv", "w") do io
        writedlm(io, datafile, ", ")
    end

    if calculate_ex
        open("Data/Random_circuits/random_circuit_rand_MPS_fid_ex_D_$(depth)_bond_$(maxdim)_samples_$(N_samples).csv", "w") do io
            writedlm(io, datafile_ex, ", ")
        end
    end
else
    open("Data/Random_circuits/random_circuit_fid_D_$(depth)_bond_$(maxdim)_samples_$(N_samples).csv", "w") do io
        writedlm(io, datafile, ", ")
    end

    if calculate_ex
        open("Data/Random_circuits/random_circuit_fid_ex_D_$(depth)_bond_$(maxdim)_samples_$(N_samples).csv", "w") do io
            writedlm(io, datafile_ex, ", ")
        end
    end
end

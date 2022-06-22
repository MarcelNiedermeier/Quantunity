

#######################################
## Investigate a random quantum circuit
#######################################

# check fidelities for fixed qubit numbers and depth, but variable
# bond dimensions


using DelimitedFiles
using Statistics
include("../QSim.jl")


backend = "MPS_ITensor"
#backend = "ED_Julia"

# constants
N = 60
lintop = false
#maxdim = 16 # 32, 64
bond_dims = [8, 16, 32, 64, 128, 256]#, 512]
bitstring = 200
N_samples = 3
depth = 100


# save data for different qubit numbers
twofidelities = zeros(length(bond_dims), depth÷2+1)
Nfidelities = zeros(length(bond_dims), depth÷2+1)
Avtwofidelities = zeros(length(bond_dims), depth÷2+1)

twofidelities_std = zeros(length(bond_dims), depth÷2+1)
Nfidelities_std = zeros(length(bond_dims), depth÷2+1)
Avtwofidelities_std = zeros(length(bond_dims), depth÷2+1)

# investigate different bond dimensions
for i in 1:length(bond_dims)

    # save results for fidelities per sample
    twofidelities_tmp = zeros(N_samples, depth÷2+1)
    Nfidelities_tmp = zeros(N_samples, depth÷2+1)
    Avtwofidelities_tmp = zeros(N_samples, depth÷2+1)

    println("Now bond_dim =  $(bond_dims[i])")

    # sample N_samples different random quantum circuits
    for samp in 1:N_samples

        #if samp%4 == 0
        #    println("Calculating sample $samp")
        #end
        println("Calculating sample $samp")

        # initialise circuit object
        local qc = initialise_qcircuit(N, lintop, backend, bond_dims[i])

        # build random circuit until desired depth
        for j in 1:(depth÷2)

            #if j%20 == 0
            #    println("applying layers $j and $(j+1)")
            #end

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

        # record fidelities
        twofidelities_tmp[samp, :] = qc.TwoQubitFidelity
        Nfidelities_tmp[samp, :] = qc.NQubitFidelity
        Avtwofidelities_tmp[samp, :] = qc.AverageTwoQubitFidelity

    end

    # average recorded fidelities over samples and find std
    twofidelities[i, :] = Statistics.mean(twofidelities_tmp, dims=1)
    Nfidelities[i, :] = Statistics.mean(Nfidelities_tmp, dims=1)
    Avtwofidelities[i, :] = Statistics.mean(Avtwofidelities_tmp, dims=1)

    twofidelities_std[i, :] = std(twofidelities_tmp, dims=1)
    Nfidelities_std[i, :] = std(Nfidelities_tmp, dims=1)
    Avtwofidelities_std[i, :] = std(Avtwofidelities_tmp, dims=1)


end


# datafile for output
datafile = zeros(6*length(bond_dims)+1, depth÷2+1)

# general information ("header")
datafile[1, 1] = N_samples
datafile[1, 2] = N
datafile[1, 3] = bitstring
datafile[1, 4] = length(bond_dims)
datafile[1, 5] = depth
for i in 1:length(bond_dims)
    datafile[1, 5+i] = bond_dims[i]
end

# data
for i in 1:length(bond_dims)
    datafile[1+i, :] = twofidelities[i, :]
    datafile[1+length(bond_dims)+i, :] = Nfidelities[i, :]
    datafile[1+2*length(bond_dims)+i, :] = Avtwofidelities[i, :]
    datafile[1+3*length(bond_dims)+i, :] = twofidelities_std[i, :]
    datafile[1+4*length(bond_dims)+i, :] = Nfidelities_std[i, :]
    datafile[1+5*length(bond_dims)+i, :] = Avtwofidelities_std[i, :]
end

# save data
open("Data/Random_circuits/random_circuit_fid_D_$(depth)_N_$(N)_samples_$(N_samples).csv", "w") do io
    writedlm(io, datafile, ", ")
end

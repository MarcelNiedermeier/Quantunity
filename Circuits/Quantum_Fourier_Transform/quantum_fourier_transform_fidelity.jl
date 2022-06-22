

#######################
# Study Fidelity of QFT
#######################


using DelimitedFiles
using Statistics
include("../QSim.jl")


# set constants
N = 40 # 16, 24, 32, 40
#maxdim = 10
#maxdims =  [1, 2, 3]
maxdims = [4, 8, 12, 16, 20, 24, 28, 32, 36, 40, 44, 48, 52, 56, 60, 64]# 64, 128]
N_sample = 4
#N_meas = 10
#backend = "ED_Julia"
backend = "MPS_ITensor"
contmethod = "naive"
random = true
lintop = false
#randombond = 2
randombonds = [1, 2, 4]#, 8, 16]


# collect data
av_fidelities = zeros(length(randombonds), length(maxdims))


for k in 1:length(randombonds)

    println("doing randombond = $(randombonds[k])")

    for j in 1:length(maxdims)

        println("maxdim = $(maxdims[j])")

        # save different fidelities before averaging
        fid_temp = zeros(N_sample)

        for i in 1:N_sample

            # initialise circuit, save intial state
            local qc = initialise_qcircuit(N, lintop, backend, maxdims[j], contmethod, random, randombond)
            #qc = initialise_qcircuit(N, lintop, backend)
            #println(qc.MaxBondDim)
            local initial_state = deepcopy(qc.StateVector)

            # do QFT and inverse QFT
            println("QFT")
            QFT!(qc, 1, N)
            println("invQFT")
            invQFT!(qc, 1, N)

            #draw(qc)

            # compare truncated state with initial state
            local fid = abs(ITensors.dot(initial_state, qc.StateVector))^2
            #println("fidelity: $fid")

            fid_temp[i] = fid

        end

        av_fidelities[k, j] = Statistics.mean(fid_temp)
    end
end

println(av_fidelities)

# save output in datafile
datafile = zeros(length(randombonds)+3, length(maxdims))

datafile[1, 1] = N
datafile[1, 2] = N_sample
datafile[1, 3] = length(randombonds)
datafile[2, 1:1+length(randombonds)-1] = randombonds[:]
datafile[3, :] = maxdims[:]

for i in 1:length(randombonds)
    datafile[3+i, :] = av_fidelities[i, :]
end

open("../Data/Quantum_Fourier_Transform/QFT_MPS_fid_N_$(N)_maxbond_$(maxdims[end])_samples_$(N_sample).csv", "w") do io
    writedlm(io, datafile, ", ")
end

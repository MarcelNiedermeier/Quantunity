
###############################
## Compare MPS QFT to exact QFT
###############################


using DelimitedFiles
using Statistics
include("../QSim.jl")


# set constants
#N = 10 # 16, 24, 32, 40
Ns = [10, 12, 14, 16]
D = 20 # depth of random circuit
maxdims =  [10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80]
#maxdims = [4, 8, 12, 16, 20]#, 24, 28, 32, 36, 40, 44, 48, 52, 56, 60, 64]# 64, 128]
N_sample = 4
#backend = "ED_Julia"
backend = "MPS_ITensor"
contmethod = "naive"
random = false
lintop = false
#randombonds = [1, 2, 4]#, 8, 16]


""" Function to calculate the coefficients of a Fourier-transformed
state directly. """
function QFT_coeff(j, N, initial_state)
    ω = exp(1.0im*2*π/(2.0^N))
    b_j = Complex(0.0)
    for k in 0:(2^N-1)
        b_j += initial_state[k+1] * ω^(j*k)
    end
    return b_j/√(2^N)
end


# collect data
av_fidelities = zeros(length(Ns), length(maxdims))

# evaluate different system sizes
for p in 1:length(Ns)

    println("system length $(Ns[p])")

    # loop through different bond dimensions
    for j in 1:length(maxdims)

        println("maxdim = $(maxdims[j])")

        # save different fidelities before averaging
        fid_temp = zeros(N_sample)

        # average over multiple samples
        for i in 1:N_sample

            println("doing sample $i")

            # initialise circuit, save intial state
            local qc = initialise_qcircuit(Ns[p], lintop, backend, maxdims[j], contmethod, random, randombond)
            #local initial_state = deepcopy(qc.StateVector)
            #local initial_state = deepcopy(get_wavefunction(qc))

            # set up random state
            println("setting up random state")
            randomCircuit!(qc, 1, Ns[p], D, true)

            # save initial state
            local initial_state = deepcopy(get_wavefunction(qc))

            # do QFT and inverse QFT
            println("QFT")
            QFT!(qc, 1, Ns[p])

            # calculate Fourier coefficients by hand
            coeffs = []
            k = 0
            for k in 0:(2^Ns[p]-1)
                #println("coeff for j=$j: $(round(QFT_coeff_had(j, N), digits=5))")
                #println("coeff for j=$j: $(round(QFT_coeff(j, N, initial_state), digits=5))")
                push!(coeffs, round(QFT_coeff(k, Ns[p], initial_state), digits=16))
            end

            # calculate fidelity
            fid_temp[i] = abs(dot(get_wavefunction(qc), coeffs))^2


            # compare truncated state with initial state
            #local fid = abs(ITensors.dot(initial_state, qc.StateVector))^2
        end

        # average over samples
        av_fidelities[p, j] = Statistics.mean(fid_temp)

    end
end


println(av_fidelities)

# save output in datafile
datafile = zeros(length(Ns)+3, length(maxdims))

datafile[1, 1] = D
datafile[1, 2] = N_sample
datafile[1, 3] = length(Ns)
datafile[2, 1:1+length(Ns)-1] = Ns[:]
datafile[3, :] = maxdims[:]

for i in 1:length(Ns)
    datafile[3+i, :] = av_fidelities[i, :]
end

open("../Data/Quantum_Fourier_Transform/QFT_MPS_exact_fid_maxN_$(Ns[end])_maxbond_$(maxdims[end])_samples_$(N_sample).csv", "w") do io
    writedlm(io, datafile, ", ")
end


############################
## Quantum Fourier Transform
############################

using DelimitedFiles
using Statistics
using LinearAlgebra
include("../QSim.jl")


# set constants
N = 5
N2 = 12
maxdim = 30
N_meas = 100
#backend = "ED_Julia"
backend = "MPS_ITensor"
lintop = false
contmethod = "naive"
random = true
#random = false
randombond = 3

# set up circuit
qc = initialise_qcircuit(N, lintop, backend)
qc2 = initialise_qcircuit(N2, lintop, backend, maxdim, contmethod, random, randombond)
#qc3 = initialise_qcircuit(2, lintop, "ED_Julia")


###############################
# build easy QFT on five qubits
###############################

#PauliX!(qc, [1])

#sites = qc.IndexSet
#n = 3
#gate = op("SWAP", sites[1], sites[2])
#println(gate)

#println(qc.StateVector)
#orthogonalize(qc.StateVector, 3)
#wf = (qc.StateVector[1]*qc.StateVector[2])*gate
#noprime!(wf)
#ind = uniqueinds(qc.StateVector[1], qc.StateVector[2])
#U, S, V = svd(wf, ind, maxdim=qc.MaxBondDim)

# first qubit
hadamard!(qc, [1])
CRn!(qc, 5, [5, 1])
CRn!(qc, 4, [4, 1])
CRn!(qc, 3, [3, 1])
CRn!(qc, 2, [2, 1])

# second qubit
hadamard!(qc, [2])
CRn!(qc, 4, [5, 2])
CRn!(qc, 3, [4, 2])
CRn!(qc, 2, [3, 2])

# third qubit
hadamard!(qc, [3])
CRn!(qc, 3, [5, 3])
CRn!(qc, 2, [4, 3])

# fourth qubit
hadamard!(qc, [4])
CRn!(qc, 2, [5, 4])

# fifth qubit
hadamard!(qc, [5])

# final swaps
fullSwap!(qc, [1, 5])
fullSwap!(qc, [2, 4])

draw(qc, false)
#sample_measurement(qc, [i for i in 1:N], N_meas)


""" Apply num-qubit QFT to quantum circuit qc, starting position pos. """
function QFT_test!(qc, pos, num)

    # Hadamard gates and controlled rotations
    for i in 1:num
        hadamard!(qc, [i])
        k = 0
        for j in pos+num-i:-1:pos+1
            #println("j+1: $(j+1)")
            local state_tmp = qc.StateVector
            #CRn!(qc, j, [pos+num-1-k, i])
            cnot!(qc, [pos+num-1-k, i])

            # other CRn scheme
            #CRn!(qc, j, [j+1, i])
            k = k+1
            fid = abs(ITensors.dot(state_tmp, qc.StateVector))^2
            #norm = abs(ITensors.dot(qc.StateVector, qc.StateVector))^2
            #println("fidelity: $fid")
            #println("norm: $norm")
        end
    end

    # Swaps in the end
    for i in 0:(num÷2-1)
        fullSwap!(qc, [pos+i, pos+num-1-i])
    end
end


""" Inverse Quantum Fourier transform. """
function invQFT_test!(qc, pos, num)

    # Swaps in the beginning
    for i in (num÷2-1):-1:0
        fullSwap!(qc, [pos+i, pos+num-1-i])
    end

    # Hadamard gates and hermitian conjugate of controlled rotations
    for i in num:-1:1
        #println("i: $i")
        #hadamard!(qc, [i])
        #k = 0
        #for j in pos+num-i:-1:pos+1
        #for j in pos+1:pos+num-i
        for j in i+1:pos+num-1
            #println("j: $j")
            #println("j-1+1: $(j-i+1)")
            #println("pos+num-1-k: $(pos+num-1-k)")
            local state_tmp = qc.StateVector
            #CRn!(qc, -(j-i+1), [j, i])
            cnot!(qc, [j, i])

            # other CRn scheme
            #CRn!(qc, -j, [pos+num-1-k, i])
            #k = k+1
            fid = abs(ITensors.dot(state_tmp, qc.StateVector))^2
            #norm = abs(ITensors.dot(qc.StateVector, qc.StateVector))^2
            #println("fidelity: $fid")
            #println("norm: $norm")
        end
        hadamard!(qc, [i])
    end
end




# different input state
PauliX!(qc2, [1, 3])#, 5, 7])#, 6, 8])
#sample_measurement(qc2, [i for i in 1:N2], N_meas)
#println("initial coefficient of [1, 0, 1, 0, 1, 0]: $(get_coefficient(qc2, [1, 0, 1, 0, 1, 0, 1, 0]))")
wf_in = get_wavefunction(qc2)
#println("initial wave function: $(get_wavefunction(qc2))")

initial_state = qc2.StateVector
println("norm initial state: $(norm(initial_state))")
#println("norm initial state: $(norm(initial_state))")
fid = abs(ITensors.dot(initial_state, initial_state))^2
println("fidelity initial state: $fid")

# do QFT on quantum circuit
println("QFT")
QFT_test!(qc2, 1, N2)
println("invQFT")
invQFT_test!(qc2, 1, N2)


#println("final coefficient of [1, 0, 1, 0, 1, 0]: $(get_coefficient(qc2, [1, 0, 1, 0, 1, 0, 1, 0]))")
wf_fin = get_wavefunction(qc2)
#println("final wave function: $(get_wavefunction(qc2))")
println("compare wave functions")
#for i in 1:length(wf_fin)
#    #println("$(round(wf_in[i], digits=10)),   $(round(wf_fin[i], digits=10))")
#    println("$(round(wf_in[i], digits=10)-round(wf_fin[i], digits=10))")
#end
println("check norm: $(norm(wf_fin-wf_in))")
println("check overlap: $(dot(wf_in, wf_fin)^2)")
println("final state: ", qc2.StateVector)


#println("norm final state: $(norm(qc2.StateVector))")
println("norm initial state: $(ITensors.dot(initial_state, initial_state))")
println("norm final state: $(ITensors.dot(qc2.StateVector, qc2.StateVector))")
println(ITensors.dot(initial_state, qc2.StateVector))
println(norm(initial_state + (-1)*qc2.StateVector)^2)
println(abs(ITensors.dot(initial_state + (-1)*qc2.StateVector, initial_state + (-1)*qc2.StateVector)))
fidelity = abs(ITensors.dot(initial_state, qc2.StateVector))^2
#fidelity = abs(ITensors.dot(initial_state, initial_state))^2
println("fidelity: $fidelity")


#println("initial state: ", initial_state)
#println("final state: ", qc2.StateVector)


sample_measurement(qc2, [i for i in 1:N2], N_meas)
draw(qc2, false)


#PauliX!(qc3, [1, 2])
#hadamard!(qc3, [1])
#CRn!(qc3, 3, [2, 1])
#println(qc3.StateVector)

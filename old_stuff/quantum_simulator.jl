
###############################
# Quantum Simulator header file
###############################

#using ITensors
using IterTools
using StatsBase
using Random
using DataStructures
using LinearAlgebra
using DelimitedFiles
using CSV
using DataFrames


####################
# struct for circuit
####################

mutable struct QC

    # principle properties
    StateVector::Vector{ComplexF64}
    NumQubits::Int32
    CircuitDepth::Int32

    # if (single) measurements are performed, sve results as classical bits
    ClassicalBits::SortedDict

    # keep track of gates in circuit
    Representation::SortedDict
    RepresentationFull::SortedDict # including swaps to realise 2-site gates
    ConversionTable::Dict # holds representation of gates

    # topology
    LinearTopology::Bool


    # coherence times of qubits/error rates
    # other?

end


# struct for gates as well?




#################
# basic functions
#################

""" Fast way to generate the Pauli matrices X, Y and Z as well as
the identity matrix E in two dimensions. """
function get_Pauli_matrices()
    E = Complex.([1. 0.; 0. 1.])
    X = Complex.([0. 1.; 1. 0.])
    Y = [0. -1.0im; 1.0im 0.]
    Z = Complex.([1. 0.; 0. -1.])
    return E, X, Y, Z
end

""" Function to build the initial state |00...0> of a quantum
circuit as an exact state vector. Needs the desired number of qubits
as input. """
function initialise_state(N)

    # initialise complex d^N dim vector of zeros, set first element to 1.0
    #psi = zeros(ntuple(i->d, N)) + zeros(ntuple(i->d, N))*im
    #psi = zeros(d^N) + zeros(d^N)*im
    psi = Complex.(zeros(2^N))
    psi[1] = 1.0

    # create ITensor object from this vector
    #sites = siteinds(d, N)
    #initial_state = ITensor(psi, sites)
    #return initial_state, sites
    return psi
end

""" Function to initialise a quantum circuit object with N qubits. """
function initialise_qcircuit(N, lintop=false)

    # initialise state as exact vector and construct circuit object
    psi = initialise_state(N)
    qcircuit = QC(psi, N, 0, SortedDict(), SortedDict(), SortedDict(), Dict(), lintop)

    # set up representation of circuit
    for i in 1:N
        qcircuit.Representation[i] = []
        qcircuit.RepresentationFull[i] = []
    end

    # initialise conversion table to print out gates
    qcircuit.ConversionTable[0] = "---"
    qcircuit.ConversionTable[1] = "-H-"
    qcircuit.ConversionTable[2] = "-X-"
    qcircuit.ConversionTable[3] = "-○-"
    qcircuit.ConversionTable[-3] = "-+-"
    qcircuit.ConversionTable[4] = "-⊗-"
    qcircuit.ConversionTable[5] = "-Z-"
    qcircuit.ConversionTable[6] = "-Y-"
    qcircuit.ConversionTable[7] = "-S-"
    qcircuit.ConversionTable[8] = "-T-"
    qcircuit.ConversionTable[9] = "-Rˣ-"
    qcircuit.ConversionTable[10] = "-Rʸ-"
    qcircuit.ConversionTable[11] = "-Rᶻ-"
    qcircuit.ConversionTable[12] = "-P-"
    qcircuit.ConversionTable[13] = "-○-"
    qcircuit.ConversionTable[-13] = "-U-"
    qcircuit.ConversionTable[14] = "-○-"
    qcircuit.ConversionTable[-14] = "-+-"
    qcircuit.ConversionTable[15] = "-|-"
    qcircuit.ConversionTable[16] = "-M-"
    qcircuit.ConversionTable[17] = "--M"

    return qcircuit
end

""" Function that takes in the meta-information about the quantum
circuit (after it has been constructed) and prints out a very schematic
representation in the console.
TO DO: correct alignment of gate, maybe extend to be (pretty-) printed
into a text file. """
function draw_circuit(qc::QC, full=false)

    # prompt for file name is desire to save in txt file!

    # check desired printing format
    if full
        qc_rep = qc.RepresentationFull
        println("Schematic representation of full quantum circuit, depth = $(qc.CircuitDepth): ")
    else
        qc_rep = qc.Representation
        println("Schematic representation of quantum circuit, depth = $(qc.CircuitDepth): ")
    end

    # loop through gate matrix and look up corresponding representation
    println(" ")
    for (qubit, gates) in qc_rep
        print(qubit, " |0⟩ ")
        for i in 1:length(gates)
            if i == length(gates)
                println(qc.ConversionTable[gates[i]])
            else
                print(qc.ConversionTable[gates[i]])
            end
        end
    end
    println("\n")
end





## in development

function save_circuit(qc, name)

    # check if path is string
    if typeof(name) != String
        error("Wrong data type for file name (need string).")
    end

    # get data from circuit
    N = qc.NumQubits
    rep = qc.Representation
    repFull = qc.RepresentationFull
    len = length(collect(values(rep))[1])
    lenFull = length(collect(values(repFull))[1])

    matQC = zeros(N, len)
    matQCFull = zeros(N, lenFull)

    # write to matrix
    for i in 1:N
        for j in 1:len
            matQC[i, j] = Int(rep[i][j])
        end
        for k in 1:lenFull
            matQCFull[i, k] = Int(repFull[i][k])
        end
    end

    # save in CSV file
    # header line with information
    # matrices with representations



    open(name, "w") do io
        writedlm(io, matQCFull)
    end
    #CSV.write(name, DataFrame(matQC, :auto))
end


function load_circuit(name)

    matQC = readdlm(name)
    println(matQC)

    N = size(matQC)[1]

end

## - end development




#####################
# auxiliary functions
#####################

function filter_binary_numbers(list_binaries, list_values, list_indices)
    out = []
    for num in list_binaries
        keep = true
        for i in 1:length(list_indices)
            if num[list_indices[i]] != list_values[i]
                keep = false
            end
        end
        if keep
            push!(out, num)
        end
    end
    return out
end


function bit_array_to_int(arr)
    s = 0
    N = length(arr)
    for i in 1:N
        s += arr[i] * 2^(N - i)
    end
    return s
end


""" Function to extract the 2-permutations from a general array of poisitions,
pos = [pos[1], pos[2]], such that the permutation represented by pos is
equivalent to the sequence of 2-permutations. """
function get_2_perms_old(pos)
    perms = []
    for i in (pos[1]+1):pos[2]
        push!(perms, [i-1, i])
    end
    return perms
end

""" Function to extract the 2-permutations from a general array of poisitions,
pos = [pos[1], pos[2]], such that the permutation represented by pos is
equivalent to the sequence of 2-permutations. """
function get_2_perms(pos)
    perms_forward = []
    perms_backward = []
    for i in (pos[1]+1):pos[2]
        push!(perms, [i-1, i])
    end
    return perms
end


function get_2_perms_from_cycle(cycle)

    perms = []

    for i in 1:(length(cycle)-1)

        if cycle[i+1] - cycle[i] == 1
            push!(perms, [cycle[i], cycle[i+1]])

        elseif cycle[i+1] - cycle[i] == -1
            push!(perms, [cycle[i+1], cycle[i]])
        else
            perms_tmp = get_2_perms_old([cycle[i], cycle[i+1]])
            perms_tmp_rev = reverse(perms_tmp[1:end-1])
            for p in perms_tmp
                push!(perms, p)
            end
            for p in perms_tmp_rev
                push!(perms, p)
            end
        end
    end

    return perms
end



#############
# measurement
#############

""" Function to sample the measurements of the wave function psi N_meas times.
Register specifies the sequence of qubits to be measured; can be an arbitrary
subregister of the full state. Prints out the different outcomes with their
corresponding frequency of occurrence ("empirical probability"). To be used at
the end of a quantum circuit to evaluate the result; doesn't perform collapse
of the wave function (as multiple samples are usually desired). """
function sample_measurement(qc::QC, register::Array{Int64, 1}, N_meas=100)

    # check correct ordering of values in register
    for i in 1:length(register)-1
        if register[i] >= register[i+1]
            error("Wrong ordering of qubit indices in register.")
        end
    end

    # obtain full probability distribution from wavefunction
    N_qubits = qc.NumQubits
    states = [reverse(digits(i, base=2, pad=N_qubits)) for i in 0:(2^N_qubits-1)]
    probabilities = abs.(qc.StateVector).^2

    # create dictionary with different states and corresponding probabilities
    probs = Dict()
    for i in 1:2^N_qubits
        probs[states[i]] = probabilities[i]
    end

    # if only measuring subregister: caculate marginalised probability distribution
    N_register = 2^length(register)
    states_marg = [reverse(digits(i, base=2, pad=length(register))) for i in 0:N_register-1]
    prob_list = []
    register_marg = Dict()

    # create a lists of equivalent states with same qubits values in selected registers
    for st in states_marg
        equiv_states = filter_binary_numbers(states, st, register)

        # calculate probability of corresponding subspaces
        prob_marg = 0
        for e_st in equiv_states
            prob_marg += probs[e_st]
        end
        push!(prob_list, prob_marg)
        register_marg[st] = (prob_marg, equiv_states)
    end

    # calculate weights corresponding to marginalised probs, sample number
    # of measurements
    weights = Weights(Vector{Float64}(collect(values(prob_list))))
    samp = StatsBase.sample(Random.GLOBAL_RNG, states_marg, weights, N_meas)
    samp = proportionmap(samp)

    # summarise result of measurements
    println("\n")
    println("Full qubit register: $([i for i in 1:N_qubits])")
    println("Doing $N_meas measurements of the register $register: ")
    println("Obtain the following states with their corresponding frequencies: ")
    for key in sort!(collect(keys(samp)))
        item = samp[key]
        println("State: ", key, ", p = ", item)
    end
    println("\n")

    # update representing matrix of quantum circuit
    for i in 1:qc.NumQubits
        if i in register
            push!(qc.Representation[i], 17)
            push!(qc.RepresentationFull[i], 17)
        else
            push!(qc.Representation[i], 0)
            push!(qc.RepresentationFull[i], 0)
        end
    end
end

""" Function to perform a measurement of the wave function psi.
Register specifies the sequence of qubits to be measured; can be an arbitrary
subregister of the full state. In contrast to the function sample_measurement(),
this function only evaluates a single sample and collapses the wavefunction
accordingly. To be used for a measurement within the circuit. The results of the
measurement are saved as classical bits in the ClassicalBits dictionary of the
quantum circuit object. """
function measure!(qc::QC, register::Array{Int64, 1})

    # check correct ordering of values in register
    for i in 1:length(register)-1
        if register[i] >= register[i+1]
            error("Wrong ordering of qubit indices in register.")
        end
    end

    # obtain full probability distribution from wavefunction
    N_qubits = qc.NumQubits
    states = [reverse(digits(i, base=2, pad=N_qubits)) for i in 0:(2^N_qubits-1)]
    probabilities = abs.(qc.StateVector).^2

    # create dictionary with different states and corresponding probabilities
    probs = Dict()
    for i in 1:2^N_qubits
        probs[states[i]] = probabilities[i]
    end

    # if only measuring subregister: caculate marginalised probability distribution
    N_register = 2^length(register)
    states_marg = [reverse(digits(i, base=2, pad=length(register))) for i in 0:N_register-1]
    prob_list = []
    register_marg = Dict()

    # create a lists of equivalent states with same qubits values in selected registers
    for st in states_marg
        equiv_states = filter_binary_numbers(states, st, register)

        # calculate probability of corresponding subspaces
        prob_marg = 0
        for e_st in equiv_states
            prob_marg += probs[e_st]
        end
        push!(prob_list, prob_marg)
        register_marg[st] = (prob_marg, equiv_states)
    end

    # calculate weights corresponding to marginalised probs, sample state
    weights = Weights(Vector{Float64}(prob_list))
    samp = StatsBase.sample(Random.GLOBAL_RNG, states_marg, weights, 1)

    # contruct projector onto measured subspace
    proj = Complex.(zeros(2^N_qubits, 2^N_qubits))
    for eq_st in register_marg[samp[1]][2]
        ind = bit_array_to_int(eq_st) + 1
        proj[ind, ind] = 1.0
    end

    # apply projector to state, renormalise
    qc.StateVector .= proj*qc.StateVector
    qc.StateVector .= qc.StateVector/sqrt(sum(abs.(qc.StateVector).^2))

    # save classical bits obtained from measurement
    for i in 1:length(register)
        qc.ClassicalBits[register[i]] = samp[1][i]
    end

    #sort!(collect(keys(qc.ClassicalBits)))

    # summarise result of measurements
    println("\n")
    println("Full qubit register: $([i for i in 1:N_qubits])")
    println("Measured register $register: ")
    println("The qubit register has collapsed to $(samp[1]). ")
    println("\n")

    # update representing matrix of quantum circuit
    for i in 1:qc.NumQubits
        if i in register
            push!(qc.Representation[i], 16)
            push!(qc.RepresentationFull[i], 16)
        else
            push!(qc.Representation[i], 0)
            push!(qc.RepresentationFull[i], 0)
        end
    end
end











####################
# single-qubit gates
####################

""" Function to apple the (single-qubit) Hadamard gate to a state psi
(with N qubits) in every position specified in the array pos. """
function hadamard!(qc::QC, pos::Array{Int64, 1})

   # get matrices
   h = 1/sqrt(2) * Complex.([1. 1.; 1. -1.])
   E, _, _, _ = get_Pauli_matrices()

   # initialise first matrix
   if pos[1] == 1
       gate = h
   else
       gate = E
   end

   # construct whole operator
   for i in 2:qc.NumQubits
       if i in pos
           gate = kron(gate, h)
       else
           gate = kron(gate, E)
       end
   end

   # create a new set of outgoing indices and define ITensor from Hadamard array
   #sites_out = siteinds(d, N)
   #gate = ITensor(gate, sites, sites_out)

   #return gate*psi
   qc.StateVector .= gate*qc.StateVector
   #psi .= gate*psi

   # update representing matrix of quantum circuit
   for i in 1:qc.NumQubits
       if i in pos
           push!(qc.Representation[i], 1)
           push!(qc.RepresentationFull[i], 1)
       else
           push!(qc.Representation[i], 0)
           push!(qc.RepresentationFull[i], 0)
       end
   end

   # update circuit depth
   qc.CircuitDepth += 1
end

""" Function to apple the (single-qubit) Pauli-X gate to a state psi
(with N qubits) in every position specified in the array pos. """
function PauliX!(qc::QC, pos::Array{Int64, 1})

   # get matrices
   E, X, _, _ = get_Pauli_matrices()

   # initialise first matrix
   if pos[1] == 1
       gate = X
   else
       gate = E
   end

   # construct whole operator
   for i in 2:qc.NumQubits
       if i in pos
           gate = kron(gate, X)
       else
           gate = kron(gate, E)
       end
   end

   # create a new set of outgoing indices and define ITensor from Hadamard array
   #sites_out = siteinds(d, N)
   #gate = ITensor(gate, sites, sites_out)
   #return gate*psi
   #psi = gate*psi
   qc.StateVector .= gate*qc.StateVector

   # update representing matrix of quantum circuit
   for i in 1:qc.NumQubits
       if i in pos
           push!(qc.Representation[i], 2)
           push!(qc.RepresentationFull[i], 2)
       else
           push!(qc.Representation[i], 0)
           push!(qc.RepresentationFull[i], 0)
       end
   end

   # update circuit depth
   qc.CircuitDepth += 1
end

""" Function to apple the (single-qubit) Pauli-Y gate to a state psi
(with N qubits) in every position specified in the array pos. """
function PauliY!(qc::QC, pos::Array{Int64, 1})

   # get matrices
   E, _, Y, _ = get_Pauli_matrices()

   # initialise first matrix
   if pos[1] == 1
       gate = Y
   else
       gate = E
   end

   # construct whole operator
   for i in 2:qc.NumQubits
       if i in pos
           gate = kron(gate, Y)
       else
           gate = kron(gate, E)
       end
   end

   # create a new set of outgoing indices and define ITensor from Hadamard array
   #sites_out = siteinds(d, N)
   #gate = ITensor(gate, sites, sites_out)
   #return gate*psi
   qc.StateVector .= gate*qc.StateVector

   # update representing matrix of quantum circuit
   for i in 1:qc.NumQubits
       if i in pos
           push!(qc.Representation[i], 6)
           push!(qc.RepresentationFull[i], 6)
       else
           push!(qc.Representation[i], 0)
           push!(qc.RepresentationFull[i], 0)
       end
   end

   # update circuit depth
   qc.CircuitDepth += 1
end

""" Function to apple the (single-qubit) Pauli-Z gate to a state psi
(with N qubits) in every position specified in the array pos. """
function PauliZ!(qc::QC, pos::Array{Int64, 1})

   # get matrices
   E, _, _, Z = get_Pauli_matrices()

   # initialise first matrix
   if pos[1] == 1
       gate = Z
   else
       gate = E
   end

   # construct whole operator
   for i in 2:qc.NumQubits
       if i in pos
           gate = kron(gate, Z)
       else
           gate = kron(gate, E)
       end
   end

   # create a new set of outgoing indices and define ITensor from Hadamard array
   #sites_out = siteinds(d, N)
   #gate = ITensor(gate, sites, sites_out)
   #return gate*psi
   qc.StateVector .= gate*qc.StateVector

   # update representing matrix of quantum circuit
   for i in 1:qc.NumQubits
       if i in pos
           push!(qc.Representation[i], 5)
           push!(qc.RepresentationFull[i], 5)
       else
           push!(qc.Representation[i], 0)
           push!(qc.RepresentationFull[i], 0)
       end
   end

   # update circuit depth
   qc.CircuitDepth += 1
end

""" Function to apple the (single-qubit) S-gate to a state psi
(with N qubits) in every position specified in the array pos. """
function SGate!(qc::QC, pos::Array{Int64, 1})

   # get matrices
   E, _, _, _ = get_Pauli_matrices()
   S = Complex.([1. 0.; 0. 1.0im])

   # initialise first matrix
   if pos[1] == 1
       gate = S
   else
       gate = E
   end

   # construct whole operator
   for i in 2:qc.NumQubits
       if i in pos
           gate = kron(gate, S)
       else
           gate = kron(gate, E)
       end
   end

   # create a new set of outgoing indices and define ITensor from Hadamard array
   #sites_out = siteinds(d, N)
   #gate = ITensor(gate, sites, sites_out)
   #return gate*psi
   qc.StateVector .= gate*qc.StateVector

   # update representing matrix of quantum circuit
   for i in 1:qc.NumQubits
       if i in pos
           push!(qc.Representation[i], 7)
           push!(qc.RepresentationFull[i], 7)
       else
           push!(qc.Representation[i], 0)
           push!(qc.RepresentationFull[i], 0)
       end
   end

   # update circuit depth
   qc.CircuitDepth += 1
end

""" Function to apple the (single-qubit) T-gate to a state psi
(with N qubits) in every position specified in the array pos. """
function TGate!(qc::QC, pos::Array{Int64, 1})

   # get matrices
   E, _, _, _ = get_Pauli_matrices()
   T = Complex.([1. 0.; 0. exp(1.0im*pi/4)])

   # initialise first matrix
   if pos[1] == 1
       gate = T
   else
       gate = E
   end

   # construct whole operator
   for i in 2:qc.NumQubits
       if i in pos
           gate = kron(gate, T)
       else
           gate = kron(gate, E)
       end
   end

   # create a new set of outgoing indices and define ITensor from Hadamard array
   #sites_out = siteinds(d, N)
   #gate = ITensor(gate, sites, sites_out)
   #return gate*psi
   qc.StateVector .= gate*qc.StateVector

   # update representing matrix of quantum circuit
   for i in 1:qc.NumQubits
       if i in pos
           push!(qc.Representation[i], 8)
           push!(qc.RepresentationFull[i], 8)
       else
           push!(qc.Representation[i], 0)
           push!(qc.RepresentationFull[i], 0)
       end
   end

   # update circuit depth
   qc.CircuitDepth += 1
end

""" Function to apple the (single-qubit) Rx-gate to a state psi
(with N qubits) in every position specified in the array pos. Performs
a rotation around the x-axis specified by theta."""
function RXGate!(qc::QC, pos::Array{Int64, 1}, theta)

   # get matrices
   E, X, _, _ = get_Pauli_matrices()
   RX = exp(-1.0im * X * theta/2)

   # initialise first matrix
   if pos[1] == 1
       gate = RX
   else
       gate = E
   end

   # construct whole operator
   for i in 2:qc.NumQubits
       if i in pos
           gate = kron(gate, RX)
       else
           gate = kron(gate, E)
       end
   end

   # create a new set of outgoing indices and define ITensor from Hadamard array
   #sites_out = siteinds(d, N)
   #gate = ITensor(gate, sites, sites_out)
   #return gate*psi
   qc.StateVector .= gate*qc.StateVector

   # update representing matrix of quantum circuit
   for i in 1:qc.NumQubits
       if i in pos
           push!(qc.Representation[i], 9)
           push!(qc.RepresentationFull[i], 9)
       else
           push!(qc.Representation[i], 0)
           push!(qc.RepresentationFull[i], 0)
       end
   end

   # update circuit depth
   qc.CircuitDepth += 1
end

""" Function to apple the (single-qubit) Ry-gate to a state psi
(with N qubits) in every position specified in the array pos. Performs
a rotation around the y-axis specified by theta."""
function RYGate!(qc::QC, pos::Array{Int64, 1}, theta)

   # get matrices
   E, _, Y, _ = get_Pauli_matrices()
   RY = exp(-1.0im * Y * theta/2)

   # initialise first matrix
   if pos[1] == 1
       gate = RY
   else
       gate = E
   end

   # construct whole operator
   for i in 2:qc.NumQubits
       if i in pos
           gate = kron(gate, RY)
       else
           gate = kron(gate, E)
       end
   end

   # create a new set of outgoing indices and define ITensor from Hadamard array
   #sites_out = siteinds(d, N)
   #gate = ITensor(gate, sites, sites_out)
   #return gate*psi
   qc.StateVector .= gate*qc.StateVector

   # update representing matrix of quantum circuit
   for i in 1:qc.NumQubits
       if i in pos
           push!(qc.Representation[i], 10)
           push!(qc.RepresentationFull[i], 10)
       else
           push!(qc.Representation[i], 0)
           push!(qc.RepresentationFull[i], 0)
       end
   end

   # update circuit depth
   qc.CircuitDepth += 1
end

""" Function to apple the (single-qubit) Rz-gate to a state psi
(with N qubits) in every position specified in the array pos. Performs
a rotation around the z-axis specified by theta."""
function RZGate!(qc::QC, pos::Array{Int64, 1}, theta)

   # get matrices
   E, _, _, Z = get_Pauli_matrices()
   RZ = exp(-1.0im * Z * theta/2)

   # initialise first matrix
   if pos[1] == 1
       gate = RZ
   else
       gate = E
   end

   # construct whole operator
   for i in 2:qc.NumQubits
       if i in pos
           gate = kron(gate, RZ)
       else
           gate = kron(gate, E)
       end
   end

   # create a new set of outgoing indices and define ITensor from Hadamard array
   #sites_out = siteinds(d, N)
   #gate = ITensor(gate, sites, sites_out)
   #return gate*psi
   qc.StateVector .= gate*qc.StateVector

   # update representing matrix of quantum circuit
   for i in 1:qc.NumQubits
       if i in pos
           push!(qc.Representation[i], 11)
           push!(qc.RepresentationFull[i], 11)
       else
           push!(qc.Representation[i], 0)
           push!(qc.RepresentationFull[i], 0)
       end
   end

   # update circuit depth
   qc.CircuitDepth += 1
end

""" Function to apply a phase shift by the angle theta to a state psi
(with N qubits) in every position specified in the array pos. """
function PhaseShift!(qc::QC, pos::Array{Int64, 1}, theta)

   # get matrices
   E, _, _, _ = get_Pauli_matrices()

   # initialise first matrix
   if pos[1] == 1
       gate = E * exp(1.0im * theta)
   else
       gate = E
   end

   # construct whole operator
   for i in 2:qc.NumQubits
       if i in pos
           gate = kron(gate, E * exp(1.0im * theta))
       else
           gate = kron(gate, E)
       end
   end

   # create a new set of outgoing indices and define ITensor from Hadamard array
   #sites_out = siteinds(d, N)
   #gate = ITensor(gate, sites, sites_out)
   #return gate*psi
   qc.StateVector .= gate*qc.StateVector

   # update representing matrix of quantum circuit
   for i in 1:qc.NumQubits
       if i in pos
           push!(qc.Representation[i], 12)
           push!(qc.RepresentationFull[i], 12)
       else
           push!(qc.Representation[i], 0)
           push!(qc.RepresentationFull[i], 0)
       end
   end

   # update circuit depth
   qc.CircuitDepth += 1
end


############################
# two-qubit gates (adjacent)
############################

""" Function to apply the CNOT-gate to two adjacent qubits in positions
pos = [pos[1], pos[2]] = [pos[1], pos[1]+1]. The "upper" qubit (pos[1])
is the control qubit, the "lower" qubit (pos[2]) the one that is acted
on. This function serves mostly as an auxiliary function in the more
general cnot! function, see below. """
function cnot2!(qc::QC, pos::Array{Int64, 1})

   if length(pos) != 2
       error("Incorrect specification of indices (need 2 positions).")
   elseif abs(pos[1]-pos[2]) != 1
       error("Incorrect specification of indices (need adjacent positions).")
   end

   # get matrices
   E, _, _, _ = get_Pauli_matrices()
   cnot = Complex.([1. 0. 0. 0.; 0. 1. 0. 0.; 0. 0. 0. 1.; 0. 0. 1. 0.])

   if pos[1] == 1
       gate = cnot

       for i in 3:qc.NumQubits
           gate = kron(gate, E)
       end
   else
       gate = E
       for i in 2:pos[1]-1
           gate = kron(gate, E)
       end
       gate = kron(gate, cnot)
       for i in (pos[2]+1):qc.NumQubits
           gate = kron(gate, E)
       end
   end

   # create a new set of outgoing indices and define ITensor from Hadamard array
   #sites_out = siteinds(d, N)
   #gate = ITensor(gate, sites, sites_out)
   #return gate*psi
   qc.StateVector .= gate*qc.StateVector

   # update representing matrix of quantum circuit
   for i in 1:qc.NumQubits
       if i == pos[1]
           push!(qc.RepresentationFull[i], 3) # control qubit
       elseif i == pos[2]
           push!(qc.RepresentationFull[i], -3) # action qubit
       else
           push!(qc.RepresentationFull[i], 0)
       end
   end

   # update circuit depth
   qc.CircuitDepth += 1
end

""" Function to swap two adjacent qubits in positions pos = [pos[1], pos[2]]
= [pos[1], pos[1]+1]. """
function swap2!(qc::QC, pos::Array{Int64, 1})

   # check structure of gate indices
   if length(pos) != 2
       error("Incorrect specification of indices (need 2 positions).")
   elseif abs(pos[1]-pos[2]) != 1
       error("Incorrect specification of indices (need adjacent positions).")
   end

   # get matrices
   E, _, _, _ = get_Pauli_matrices()
   swap = Complex.([1. 0. 0. 0.; 0. 0. 1. 0.; 0. 1. 0. 0.; 0. 0. 0. 1.])

   if pos[1] == 1
       gate = swap

       for i in 3:qc.NumQubits
           gate = kron(gate, E)
       end
   else
       gate = E
       for i in 2:pos[1]-1
           gate = kron(gate, E)
       end
       gate = kron(gate, swap)
       for i in (pos[2]+1):qc.NumQubits
           gate = kron(gate, E)
       end
   end

   # create a new set of outgoing indices and define ITensor from Hadamard array
   #sites_out = siteinds(d, N)
   #gate = ITensor(gate, sites, sites_out)
   #return gate*psi
   #psi = gate*psi
   qc.StateVector .= gate*qc.StateVector

   # update representing matrix of quantum circuit
   for i in 1:qc.NumQubits
       if i in pos
           push!(qc.RepresentationFull[i], 4)
       else
           push!(qc.RepresentationFull[i], 0)
       end
   end

   # update circuit depth
   qc.CircuitDepth += 1
end

""" Function to apply a controlled-U gate to two adjacent qubits in positions
pos = [pos[1], pos[2]] = [pos[1], pos[1]+1], where U is a unitary single-qubit
gate. The "upper" qubit (pos[1]) is the control qubit, the "lower" qubit
(pos[2]) the one that is acted on. This function serves mostly as an auxiliary
function in the more general CU! function, see below. """
function CU2!(qc::QC, pos::Array{Int64, 1})

   # check if always gets two indices!
   if length(pos) != 2
       error("Incorrect specification of indices (need 2 positions).")
   elseif abs(pos[1]-pos[2]) != 1
       error("Incorrect specification of indices (need adjacent positions).")
   end

   # get matrices
   E, _, _, _ = get_Pauli_matrices()
   CU = Complex.(zeros(4,4))
   CU[1:2, 1:2] = E
   CU[3:4, 3:4] = U

   if pos[1] == 1
       gate = CU

       for i in 3:qc.NumQubits
           gate = kron(gate, E)
       end
   else
       gate = E
       for i in 2:pos[1]-1
           gate = kron(gate, E)
       end
       gate = kron(gate, CU)
       for i in (pos[2]+1):qc.NumQubits
           gate = kron(gate, E)
       end
   end

   # create a new set of outgoing indices and define ITensor from Hadamard array
   #sites_out = siteinds(d, N)
   #gate = ITensor(gate, sites, sites_out)
   #return gate*psi
   qc.StateVector .= gate*qc.StateVector

   # update representing matrix of quantum circuit
   for i in 1:qc.NumQubits
       if i == pos[1]
           push!(qc.RepresentationFull[i], 13) # control qubit
       elseif i == pos[2]
           push!(qc.RepresentationFull[i], -13) # action qubit
       else
           push!(qc.RepresentationFull[i], 0)
       end
   end

   # update circuit depth
   qc.CircuitDepth += 1
end


################################
# two-qubit gates (not adjacent)
################################

""" Function to perfom a swap of two qubits over an arbitrary separation,
which is given by the array pos. pos[1] is understood to be smaller than
pos[2]. The swap is obtained through an equivalent sequence of 2-swaps. """
function swap!(qc::QC, pos::Array{Int64, 1})

    # check correct format of indices
    if length(pos) != 2
        error("Incorrect specification of indices (need 2 positions).")
    elseif pos[1] >= pos[2]
        error("Incorrect specification of indices (pos[1] should be smaller than pos[2]).")
    end

    if pos[1] == (pos[2]-1)
        swap2!(qc::QC, pos)
    else
        positions = get_2_perms_old(pos)
        for p in positions
            swap2!(qc::QC, p)
        end
    end
end

""" Function to reverse the arbitrary swap, performed by the swap! function."""
function unswap!(qc::QC, pos::Array{Int64, 1})

    # check correct format of indices
    if length(pos) != 2
        error("Incorrect specification of indices (need 2 positions).")
    elseif pos[1] >= pos[2]
        error("Incorrect specification of indices (pos[1] should be smaller than pos[2]).")
    end

    if pos[1] == (pos[2]-1)
        swap2!(qc::QC, pos)
    else
        positions = get_2_perms_old(pos)
        for p in reverse(positions)
            swap2!(qc::QC, p)
        end
    end
end


# the two functions above perform "half-swaps", still need function
# for full arbitrary swap!




""" Function to apply the CNOT gate to two arbitrary qubits, defined by
positions in pos = [pos[1], pos[2]]. The first qubits in is the control
qubit, the second the qubit that is acted on. The control qubit can be
"above" or "below" the action qubit, the function implements the desired
arrangement. """
function cnot!(qc::QC, pos::Array{Int64, 1})

    # check correct format of indices
    if length(pos) != 2
        error("Incorrect specification of indices (need 2 positions).")
    elseif pos[1] == pos[2]
        error("Incorrect specification of indices.")
    end

    # apply general CNOT via SWAP gates for linear topology
    if qc.LinearTopology == true


        # control qubit "above" qubit that is acted on
        if pos[2] > pos[1]

            # qubits on two adjacent lines
            if pos[1] == (pos[2]-1)
                cnot2!(qc::QC, pos)

            # qubits not on adjacent lines
            else
                swap!(qc::QC, [pos[1], pos[2]-1])
                cnot2!(qc::QC, [pos[2]-1, pos[2]])
                unswap!(qc::QC, [pos[1], pos[2]-1])
            end

        # control qubit "below" qubit that is acted on
        else

            # qubits on two adjacent lines
            if pos[1] == (pos[2]+1)
                swap2!(qc::QC, pos)
                cnot2!(qc::QC, pos)
                swap2!(qc::QC, pos)

            # qubits on two adjacent lines
            else
                swap!(qc::QC, [pos[2], pos[1]])
                cnot2!(qc::QC, [pos[1]-1, pos[1]])
                unswap!(qc::QC, [pos[2], pos[1]])
            end
        end

        # update representing matrix of quantum circuit
        for i in 1:qc.NumQubits
            if i == pos[1]
                push!(qc.Representation[i], 3) # control qubit
            elseif i == pos[2]
                push!(qc.Representation[i], -3) # action qubit
            elseif minimum(pos) < i && i < maximum(pos)
                push!(qc.Representation[i], 15) # vertical line
            else
                push!(qc.Representation[i], 0)
            end
        end

    # apply general CNOT via direct construction of matrix for master topology
    else

        println("am here")
        println("position: ", pos)

        N = qc.NumQubits
        E, X, _, - = get_Pauli_matrices()

        # projectors |0><0| and |1><1|
        proj0 = Complex.([1. 0.; 0. 0.])
        proj1 = Complex.([0. 0.; 0. 1.])
        println(proj0)
        println(proj1)

        # control qubit "above" qubit that is acted on
        if pos[2] > pos[1]

            # qubit 1 is control bit
            if pos[1] == 1
                gate_tmp0 = proj0
                gate_tmp1 = proj1

                # fill with identities
                println("now here")
                println("position: ", pos)
                println(typeof(pos))
                for i in 2:pos[2]-1
                    println(i)
                end
                println("now here")
                println("position: ", pos)
                for i in 2:pos[2]-1
                    println("now here")
                    gate_tmp0 = kron(gate_tmp0, E)
                    gate_tmp1 = kron(gate_tmp1, E)
                end

                # apply X gate in action space
                gate_tmp0 = kron(gate_tmp0, E)
                gate_tmp1 = kron(gate_tmp1, X)

                # fill with identities
                for i in pos[2]+1:N
                    gate_tmp0 = kron(gate_tmp0, E)
                    gate_tmp1 = kron(gate_tmp1, E)
                end

                # complete CNOT gate
                gate = gate_tmp0 + gate_tmp1

            # qubit 1 is NOT control bit
            else

                # initialise with identities
                gate_tmp0 = E
                gate_tmp1 = E

                # fill with identities
                for i in 2:pos[1]-1
                    gate_tmp0 = kron(gate_tmp0, E)
                    gate_tmp1 = kron(gate_tmp1, E)
                end

                gate_tmp0 = kron(gate_tmp0, proj0)
                gate_tmp1 = kron(gate_tmp1, proj1)

                # fill with identities
                for i in pos[1]:pos[2]-1
                    gate_tmp0 = kron(gate_tmp0, E)
                    gate_tmp1 = kron(gate_tmp1, E)
                end

                # apply X gate in action space
                gate_tmp0 = kron(gate_tmp0, E)
                gate_tmp1 = kron(gate_tmp1, X)

                # fill with identities
                for i in pos[2]+1:N
                    gate_tmp0 = kron(gate_tmp0, E)
                    gate_tmp1 = kron(gate_tmp1, E)
                end

                # complete CNOT gate
                gate = gate_tmp0 + gate_tmp1
            end

        # control qubit "below" qubit that is acted on
        else

            # qubit 1 is action bit
            if pos[1] == 1
                gate_tmp0 = E
                gate_tmp1 = X

                # fill with identities
                for i in 2:pos[2]-1
                    gate_tmp0 = kron(gate_tmp0, E)
                    gate_tmp1 = kron(gate_tmp1, E)
                end

                # apply projectors in control space
                gate_tmp0 = kron(gate_tmp0, proj0)
                gate_tmp1 = kron(gate_tmp1, proj1)

                # fill with identities
                for i in pos[2]+1:N
                    gate_tmp0 = kron(gate_tmp0, E)
                    gate_tmp1 = kron(gate_tmp1, E)
                end

                # complete CNOT gate
                gate = gate_tmp0 + gate_tmp1

            # qubit 1 is NOT action bit
            else

                # initialise with identities
                gate_tmp0 = E
                gate_tmp1 = E

                # fill with identities
                for i in 2:pos[1]-1
                    gate_tmp0 = kron(gate_tmp0, E)
                    gate_tmp1 = kron(gate_tmp1, E)
                end

                gate_tmp0 = kron(gate_tmp0, E)
                gate_tmp1 = kron(gate_tmp1, X)

                # fill with identities
                for i in pos[1]:pos[2]-1
                    gate_tmp0 = kron(gate_tmp0, E)
                    gate_tmp1 = kron(gate_tmp1, E)
                end

                # apply projectors in control space
                gate_tmp0 = kron(gate_tmp0, proj0)
                gate_tmp1 = kron(gate_tmp1, proj1)

                # fill with identities
                for i in pos[2]+1:N
                    gate_tmp0 = kron(gate_tmp0, E)
                    gate_tmp1 = kron(gate_tmp1, E)
                end

                # complete CNOT gate
                gate = gate_tmp0 + gate_tmp1
            end
        end

        # update statevector
        qc.StateVector .= gate*qc.StateVector

        # update representing matrix of quantum circuit
        for i in 1:qc.NumQubits
            if i == pos[1]
                push!(qc.Representation[i], 3) # control qubit
                push!(qc.RepresentationFull[i], 3) # control qubit
            elseif i == pos[2]
                push!(qc.Representation[i], -3) # action qubit
                push!(qc.RepresentationFull[i], -3) # action qubit
            elseif minimum(pos) < i && i < maximum(pos)
                push!(qc.Representation[i], 15) # vertical line
                push!(qc.RepresentationFull[i], 15) # vertical line
            else
                push!(qc.Representation[i], 0)
                push!(qc.RepresentationFull[i], 0)
            end
        end

        # update cicruit depth
        qc.CircuitDepth += 1

    end
end

""" Function to apply the CU gate to two arbitrary qubits, defined by
positions in pos = [pos[1], pos[2]]. U is a unitary single-qubit gate.
The first qubits in is the control qubit, the second the qubit that is
acted on. The control qubit can be "above" or "below" the action qubit,
the function implements the desired arrangement. """
function CU!(qc::QC, pos::Array{Int64, 1}, U)

    # check correct format of indices
    if length(pos) != 2
        error("Incorrect specification of indices (need 2 positions).")
    elseif pos[1] == pos[2]
        error("Incorrect specification of indices.")
    end

    # apply general CNOT via SWAP gates for linear topology
    if qc.LinearTopology == true

        # control qubit "above" qubit that is acted on
        if pos[2] > pos[1]

            # qubits on two adjacent lines
            if pos[1] == (pos[2]-1)
                CU2!(qc::QC, pos, U)

            # qubits not on adjacent lines
            else
                swap!(qc::QC, [pos[1], pos[2]-1])
                CU2!(qc::QC, [pos[2]-1, pos[2]], U)
                unswap!(qc::QC, [pos[1], pos[2]-1])
            end

        # control qubit "below" qubit that is acted on
        else

            # qubits on two adjacent lines
            if pos[1] == (pos[2]+1)
                swap2!(qc::QC, pos)
                CU2!(qc::QC, pos, U)
                swap2!(qc::QC, pos)

            # qubits on two adjacent lines
            else
                swap!(qc::QC, [pos[2], pos[1]])
                CU2!(qc::QC, [pos[1]-1, pos[1]], U)
                unswap!(qc::QC, [pos[2], pos[1]])
            end
        end

        # update representing matrix of quantum circuit
        for i in 1:qc.NumQubits
            if i == pos[1]
                push!(qc.Representation[i], 13) # control qubit
            elseif i == pos[2]
                push!(qc.Representation[i], -13) # action qubit
            elseif minimum(pos) < i && i < maximum(pos)
                push!(qc.Representation[i], 15) # vertical line
            else
                push!(qc.Representation[i], 0)
            end
        end

    # apply general CU via direct construction of matrix for master topology
    else

        N = qc.NumQubits
        E, _, _, _  = get_Pauli_matrices()

        # projectors |0><0| and |1><1|
        proj0 = Complex.([1. 0.; 0. 0.])
        proj1 = Complex.([0. 0.; 0. 1.])

        # control qubit "above" qubit that is acted on
        if pos[2] > pos[1]

            # qubit 1 is control bit
            if pos[1] == 1
                gate_tmp0 = proj0
                gate_tmp1 = proj1

                # fill with identities
                for i in 2:pos[2]-1
                    gate_tmp0 = kron(gate_tmp0, E)
                    gate_tmp1 = kron(gate_tmp1, E)
                end

                # apply X gate in action space
                gate_tmp0 = kron(gate_tmp0, E)
                gate_tmp1 = kron(gate_tmp1, U)

                # fill with identities
                for i in pos[2]+1:N
                    gate_tmp0 = kron(gate_tmp0, E)
                    gate_tmp1 = kron(gate_tmp1, E)
                end

                # complete CU gate
                gate = gate_tmp0 + gate_tmp1

            # qubit 1 is NOT control bit
            else

                # initialise with identities
                gate_tmp0 = E
                gate_tmp1 = E

                # fill with identities
                for i in 2:pos[1]-1
                    gate_tmp0 = kron(gate_tmp0, E)
                    gate_tmp1 = kron(gate_tmp1, E)
                end

                gate_tmp0 = kron(gate_tmp0, proj0)
                gate_tmp1 = kron(gate_tmp1, proj1)

                # fill with identities
                for i in pos[1]:pos[2]-1
                    gate_tmp0 = kron(gate_tmp0, E)
                    gate_tmp1 = kron(gate_tmp1, E)
                end

                # apply X gate in action space
                gate_tmp0 = kron(gate_tmp0, E)
                gate_tmp1 = kron(gate_tmp1, U)

                # fill with identities
                for i in pos[2]+1:N
                    gate_tmp0 = kron(gate_tmp0, E)
                    gate_tmp1 = kron(gate_tmp1, E)
                end

                # complete CU gate
                gate = gate_tmp0 + gate_tmp1
            end

        # control qubit "below" qubit that is acted on
        else

            # qubit 1 is action bit
            if pos[1] == 1
                gate_tmp0 = E
                gate_tmp1 = U

                # fill with identities
                for i in 2:pos[2]-1
                    gate_tmp0 = kron(gate_tmp0, E)
                    gate_tmp1 = kron(gate_tmp1, E)
                end

                # apply projectors in control space
                gate_tmp0 = kron(gate_tmp0, proj0)
                gate_tmp1 = kron(gate_tmp1, proj1)

                # fill with identities
                for i in pos[2]+1:N
                    gate_tmp0 = kron(gate_tmp0, E)
                    gate_tmp1 = kron(gate_tmp1, E)
                end

                # complete CU gate
                gate = gate_tmp0 + gate_tmp1

            # qubit 1 is NOT action bit
            else

                # initialise with identities
                gate_tmp0 = E
                gate_tmp1 = E

                # fill with identities
                for i in 2:pos[1]-1
                    gate_tmp0 = kron(gate_tmp0, E)
                    gate_tmp1 = kron(gate_tmp1, E)
                end

                gate_tmp0 = kron(gate_tmp0, E)
                gate_tmp1 = kron(gate_tmp1, U)

                # fill with identities
                for i in pos[1]:pos[2]-1
                    gate_tmp0 = kron(gate_tmp0, E)
                    gate_tmp1 = kron(gate_tmp1, E)
                end

                # apply projectors in control space
                gate_tmp0 = kron(gate_tmp0, proj0)
                gate_tmp1 = kron(gate_tmp1, proj1)

                # fill with identities
                for i in pos[2]+1:N
                    gate_tmp0 = kron(gate_tmp0, E)
                    gate_tmp1 = kron(gate_tmp1, E)
                end

                # complete CU gate
                gate = gate_tmp0 + gate_tmp1
            end
        end

        # update statevector
        qc.StateVector .= gate*qc.StateVector

        # update representing matrix of quantum circuit
        for i in 1:qc.NumQubits
            if i == pos[1]
                push!(qc.Representation[i], 13) # control qubit
                push!(qc.RepresentationFull[i], 13) # control qubit
            elseif i == pos[2]
                push!(qc.Representation[i], -13) # action qubit
                push!(qc.RepresentationFull[i], -13) # action qubit
            elseif minimum(pos) < i && i < maximum(pos)
                push!(qc.Representation[i], 15) # vertical line
                push!(qc.RepresentationFull[i], 15) # vertical line
            else
                push!(qc.Representation[i], 0)
                push!(qc.RepresentationFull[i], 0)
            end
        end

        # update circuit depth
        qc.CircuitDepth += 1

    end
end


##############################
# three-qubit gates (adjacent)
##############################

""" Function to apply the Toffoli gate to three adjacent qubits in positions
pos = [pos[1], pos[2], pos[3]] = [pos[1], pos[1]+1, pos[1]+2]. The "upper"
qubits (pos[1], pos[2]) are the control qubits, the "lower" qubit (pos[3])
the one that is acted on. This function serves mostly as an auxiliary function
in the more general toffoli! function, see below. """
function toffoli3!(qc::QC, pos::Array{Int64, 1})

   # check correct format of indices and size of circuit
   if qc.NumQubits < 3
       error("Trying to apply a 3-qubit gate to less than 3 qubits!")
   elseif length(pos) != 3
       error("Incorrect specification of indices (need 3 positions).")
   elseif pos[3]-pos[2] != 1 && pos[2]-pos[1] != 1
       error("Incorrect specification of indices (need adjecent positions).")
   end

   # get matrices
   E, _, _, _ = get_Pauli_matrices()
   toffoli = Complex.(Matrix{Float64}(I, 8, 8))
   toffoli[7:8, 7:8] = Complex.([0. 1.; 1. 0.])

   if pos[1] == 1
       gate = toffoli
       for i in 4:qc.NumQubits
           gate = kron(gate, E)
       end
   else
       gate = E
       for i in 2:pos[1]-1
           gate = kron(gate, E)
       end
       gate = kron(gate, toffoli)
       for i in (pos[3]+1):qc.NumQubits
           gate = kron(gate, E)
       end
   end

   # create a new set of outgoing indices and define ITensor from Hadamard array
   #sites_out = siteinds(d, N)
   #gate = ITensor(gate, sites, sites_out)
   #return gate*psi
   qc.StateVector .= gate*qc.StateVector

   # update representing matrix of quantum circuit
   for i in 1:qc.NumQubits
       if i == pos[1]
           #push!(qc.Representation[i], 12)
           push!(qc.RepresentationFull[i], 14) # control qubit
       elseif i == pos[2]
           push!(qc.RepresentationFull[i], -14) # action qubit
       else
           #push!(qc.Representation[i], 0)
           push!(qc.RepresentationFull[i], 0)
       end
   end

   # update circuit depth
   qc.CircuitDepth += 1
end


##################################
# three-qubit gates (not adjacent)
##################################

""" Function to apply the Toffoli gate in an arbitrary configuration
of the input positions. """
function toffoli!(qc::QC, pos::Array{Int64, 1})

    # check correct format of indices and size of circuit
    if qc.NumQubits < 3
        error("Trying to apply a 3-qubit gate to less than 3 qubits!")
    elseif length(pos) != 3
        error("Incorrect specification of indices (need 3 positions).")
    elseif pos[1] >= pos[2]
        error("Incorrect specification of indices (check control-qubits in ascending order).")
    end

    # apply general CNOT via SWAP gates for linear topology
    if qc.LinearTopology == true

        # reference array of qubits
        qubits = [i for i in 1:qc.NumQubits]

        # get permutation of qubits from position indication
        perm = [i for i in 1:(minimum(pos)-1)]
        append!(perm, pos)
        append!(perm, setdiff(qubits, perm))

        # get non-trivial cycles that the permutation above defines
        cycs = collect(cycles(Perm(perm)))
        for cyc in cycs
            if length(cyc) == 1
                deleteat!(cycs, findall(x->x==cyc, cycs))
            end
        end

        # convert into sequence of 2-swaps on neighbouring qubits
        swaps = []
        for cyc in cycs
            swap_tmp = get_2_perms_from_cycle(cyc)
            for sw in swap_tmp
                push!(swaps, sw)
            end
        end

        # swap qubits, apply Toffoli, swap back
        for sw in swaps
            swap2!(qc::QC, sw)
        end
        toffoli3!(qc::QC, [minimum(pos), minimum(pos)+1, minimum(pos)+2])
        for sw in reverse(swaps)
            swap2!(qc::QC, sw)
        end

        # update representing matrix of quantum circuit
        for i in 1:qc.NumQubits
            if i == pos[1]
                push!(qc.Representation[i], 14) # control qubit
            elseif i == pos[2]
                    push!(qc.Representation[i], 14) # control qubit
            elseif i == pos[3]
                push!(qc.Representation[i], -14) # action qubit
            elseif minimum(pos) < i && i < maximum(pos)
                push!(qc.Representation[i], 15) # vertical line
            else
                push!(qc.Representation[i], 0)
            end
        end

    # apply general Toffoli gate via direct construction of matrix for master topology
    else

        N = qc.NumQubits
        E, X, _, _  = get_Pauli_matrices()

        # projectors |0><0| and |1><1|
        proj0 = Complex.([1. 0.; 0. 0.])
        proj1 = Complex.([0. 0.; 0. 1.])

        # action qubit below both control bits
        if pos[3] > pos[1] && if pos[3] > pos[2]

            # first control bit in first position 1
            if pos[1] == 1

                # initialise as projectors
                gate_tmp00 = proj0
                gate_tmp01 = proj0
                gate_tmp10 = proj1
                gate_tmp11 = proj1

                # fill with identities
                for i in 2:pos[2]-1
                    gate_tmp00 = kron(gate_tmp00, E)
                    gate_tmp01 = kron(gate_tmp01, E)
                    gate_tmp10 = kron(gate_tmp10, E)
                    gate_tmp11 = kron(gate_tmp11, E)
                end

                # projection space
                gate_tmp00 = kron(gate_tmp00, proj0)
                gate_tmp01 = kron(gate_tmp01, proj1)
                gate_tmp10 = kron(gate_tmp10, proj0)
                gate_tmp11 = kron(gate_tmp11, proj1)

                # fill with identities
                for i in pos[2]+1:pos[3]-1
                    gate_tmp00 = kron(gate_tmp00, E)
                    gate_tmp01 = kron(gate_tmp01, E)
                    gate_tmp10 = kron(gate_tmp10, E)
                    gate_tmp11 = kron(gate_tmp11, E)
                end

                # apply X gate in action space
                gate_tmp00 = kron(gate_tmp00, E)
                gate_tmp01 = kron(gate_tmp01, E)
                gate_tmp10 = kron(gate_tmp10, E)
                gate_tmp11 = kron(gate_tmp11, X)

                # fill with identities
                for i in pos[3]+1:N
                    gate_tmp00 = kron(gate_tmp00, E)
                    gate_tmp01 = kron(gate_tmp01, E)
                    gate_tmp10 = kron(gate_tmp10, E)
                    gate_tmp11 = kron(gate_tmp11, E)
                end

                # complete Toffoli gate
                gate = gate_tmp00 + gate_tmp01 + gate_tmp10 + gate_tmp11

            # first control bit NOT in first position 1
            else

                # initialise as identities
                gate_tmp00 = E
                gate_tmp01 = E
                gate_tmp10 = E
                gate_tmp11 = E

                # fill with identities
                for i in 2:pos[1]-1
                    gate_tmp00 = kron(gate_tmp00, E)
                    gate_tmp01 = kron(gate_tmp01, E)
                    gate_tmp10 = kron(gate_tmp10, E)
                    gate_tmp11 = kron(gate_tmp11, E)
                end

                # projection space
                gate_tmp00 = kron(gate_tmp00, proj0)
                gate_tmp01 = kron(gate_tmp01, proj0)
                gate_tmp10 = kron(gate_tmp10, proj1)
                gate_tmp11 = kron(gate_tmp11, proj1)

                # fill with identities
                for i in pos[1]+1:pos[2]-1
                    gate_tmp00 = kron(gate_tmp00, E)
                    gate_tmp01 = kron(gate_tmp01, E)
                    gate_tmp10 = kron(gate_tmp10, E)
                    gate_tmp11 = kron(gate_tmp11, E)
                end

                # projection space
                gate_tmp00 = kron(gate_tmp00, proj0)
                gate_tmp01 = kron(gate_tmp01, proj1)
                gate_tmp10 = kron(gate_tmp10, proj0)
                gate_tmp11 = kron(gate_tmp11, proj1)

                # fill with identities
                for i in pos[2]+1:pos[3]-1
                    gate_tmp00 = kron(gate_tmp00, E)
                    gate_tmp01 = kron(gate_tmp01, E)
                    gate_tmp10 = kron(gate_tmp10, E)
                    gate_tmp11 = kron(gate_tmp11, E)
                end

                # apply X gate in action space
                gate_tmp00 = kron(gate_tmp00, E)
                gate_tmp01 = kron(gate_tmp01, E)
                gate_tmp10 = kron(gate_tmp10, E)
                gate_tmp11 = kron(gate_tmp11, X)

                # fill with identities
                for i in pos[3]+1:N
                    gate_tmp00 = kron(gate_tmp00, E)
                    gate_tmp01 = kron(gate_tmp01, E)
                    gate_tmp10 = kron(gate_tmp10, E)
                    gate_tmp11 = kron(gate_tmp11, E)
                end

                # complete Toffoli gate
                gate = gate_tmp00 + gate_tmp01 + gate_tmp10 + gate_tmp11
            end

        # action qubit between control bit
        elseif pos[3] > pos[1] && if pos[3] < pos[2]

            # first control bit in first position 1
            if pos[1] == 1

                # initialise as projectors
                gate_tmp00 = proj0
                gate_tmp01 = proj0
                gate_tmp10 = proj1
                gate_tmp11 = proj1

                # fill with identities
                for i in 2:pos[2]-1
                    gate_tmp00 = kron(gate_tmp00, E)
                    gate_tmp01 = kron(gate_tmp01, E)
                    gate_tmp10 = kron(gate_tmp10, E)
                    gate_tmp11 = kron(gate_tmp11, E)
                end

                # action space
                gate_tmp00 = kron(gate_tmp00, E)
                gate_tmp01 = kron(gate_tmp01, E)
                gate_tmp10 = kron(gate_tmp10, E)
                gate_tmp11 = kron(gate_tmp11, X)


                # fill with identities
                for i in pos[2]+1:pos[3]-1
                    gate_tmp00 = kron(gate_tmp00, E)
                    gate_tmp01 = kron(gate_tmp01, E)
                    gate_tmp10 = kron(gate_tmp10, E)
                    gate_tmp11 = kron(gate_tmp11, E)
                end

                # projection space
                gate_tmp00 = kron(gate_tmp00, proj0)
                gate_tmp01 = kron(gate_tmp01, proj1)
                gate_tmp10 = kron(gate_tmp10, proj0)
                gate_tmp11 = kron(gate_tmp11, proj1)

                # fill with identities
                for i in pos[3]+1:N
                    gate_tmp00 = kron(gate_tmp00, E)
                    gate_tmp01 = kron(gate_tmp01, E)
                    gate_tmp10 = kron(gate_tmp10, E)
                    gate_tmp11 = kron(gate_tmp11, E)
                end

                # complete Toffoli gate
                gate = gate_tmp00 + gate_tmp01 + gate_tmp10 + gate_tmp11

            # first control qubit NOT in position 1
            else

                # initialise as identities
                gate_tmp00 = E
                gate_tmp01 = E
                gate_tmp10 = E
                gate_tmp11 = E

                # fill with identities
                for i in 2:pos[1]-1
                    gate_tmp00 = kron(gate_tmp00, E)
                    gate_tmp01 = kron(gate_tmp01, E)
                    gate_tmp10 = kron(gate_tmp10, E)
                    gate_tmp11 = kron(gate_tmp11, E)
                end

                # projection space
                gate_tmp00 = kron(gate_tmp00, proj0)
                gate_tmp01 = kron(gate_tmp01, proj0)
                gate_tmp10 = kron(gate_tmp10, proj1)
                gate_tmp11 = kron(gate_tmp11, proj1)

                # fill with identities
                for i in pos[1]+1:pos[2]-1
                    gate_tmp00 = kron(gate_tmp00, E)
                    gate_tmp01 = kron(gate_tmp01, E)
                    gate_tmp10 = kron(gate_tmp10, E)
                    gate_tmp11 = kron(gate_tmp11, E)
                end

                # action space
                gate_tmp00 = kron(gate_tmp00, E)
                gate_tmp01 = kron(gate_tmp01, E)
                gate_tmp10 = kron(gate_tmp10, E)
                gate_tmp11 = kron(gate_tmp11, X)

                # fill with identities
                for i in pos[2]+1:pos[3]-1
                    gate_tmp00 = kron(gate_tmp00, E)
                    gate_tmp01 = kron(gate_tmp01, E)
                    gate_tmp10 = kron(gate_tmp10, E)
                    gate_tmp11 = kron(gate_tmp11, E)
                end

                # projection space
                gate_tmp00 = kron(gate_tmp00, proj0)
                gate_tmp01 = kron(gate_tmp01, proj1)
                gate_tmp10 = kron(gate_tmp10, proj0)
                gate_tmp11 = kron(gate_tmp11, proj1)

                # fill with identities
                for i in pos[3]+1:N
                    gate_tmp00 = kron(gate_tmp00, E)
                    gate_tmp01 = kron(gate_tmp01, E)
                    gate_tmp10 = kron(gate_tmp10, E)
                    gate_tmp11 = kron(gate_tmp11, E)
                end

                # complete Toffoli gate
                gate = gate_tmp00 + gate_tmp01 + gate_tmp10 + gate_tmp11
            end

        # action qubit above both control bits
        else

            # action bit in first position 1
            if pos[1] == 1

                # initialise as projectors
                gate_tmp00 = proj0
                gate_tmp01 = proj0
                gate_tmp10 = proj1
                gate_tmp11 = proj1

                # fill with identities
                for i in 2:pos[2]-1
                    gate_tmp00 = kron(gate_tmp00, E)
                    gate_tmp01 = kron(gate_tmp01, E)
                    gate_tmp10 = kron(gate_tmp10, E)
                    gate_tmp11 = kron(gate_tmp11, E)
                end

                # projection space
                gate_tmp00 = kron(gate_tmp00, proj0)
                gate_tmp01 = kron(gate_tmp01, proj1)
                gate_tmp10 = kron(gate_tmp10, proj0)
                gate_tmp11 = kron(gate_tmp11, proj1)

                # fill with identities
                for i in pos[2]+1:pos[3]-1
                    gate_tmp00 = kron(gate_tmp00, E)
                    gate_tmp01 = kron(gate_tmp01, E)
                    gate_tmp10 = kron(gate_tmp10, E)
                    gate_tmp11 = kron(gate_tmp11, E)
                end

                # action space
                gate_tmp00 = kron(gate_tmp00, E)
                gate_tmp01 = kron(gate_tmp01, E)
                gate_tmp10 = kron(gate_tmp10, E)
                gate_tmp11 = kron(gate_tmp11, X)

                # fill with identities
                for i in pos[3]+1:N
                    gate_tmp00 = kron(gate_tmp00, E)
                    gate_tmp01 = kron(gate_tmp01, E)
                    gate_tmp10 = kron(gate_tmp10, E)
                    gate_tmp11 = kron(gate_tmp11, E)
                end

                # complete Toffoli gate
                gate = gate_tmp00 + gate_tmp01 + gate_tmp10 + gate_tmp11

            # action bit NOT in position 1
            else

                # initialise as identities
                gate_tmp00 = E
                gate_tmp01 = E
                gate_tmp10 = E
                gate_tmp11 = E

                # fill with identities
                for i in 2:pos[1]-1
                    gate_tmp00 = kron(gate_tmp00, E)
                    gate_tmp01 = kron(gate_tmp01, E)
                    gate_tmp10 = kron(gate_tmp10, E)
                    gate_tmp11 = kron(gate_tmp11, E)
                end

                # action space
                gate_tmp00 = kron(gate_tmp00, E)
                gate_tmp01 = kron(gate_tmp01, E)
                gate_tmp10 = kron(gate_tmp10, E)
                gate_tmp11 = kron(gate_tmp11, X)

                # fill with identities
                for i in pos[1]+1:pos[2]-1
                    gate_tmp00 = kron(gate_tmp00, E)
                    gate_tmp01 = kron(gate_tmp01, E)
                    gate_tmp10 = kron(gate_tmp10, E)
                    gate_tmp11 = kron(gate_tmp11, E)
                end

                # projection space
                gate_tmp00 = kron(gate_tmp00, proj0)
                gate_tmp01 = kron(gate_tmp01, proj0)
                gate_tmp10 = kron(gate_tmp10, proj1)
                gate_tmp11 = kron(gate_tmp11, proj1)

                # fill with identities
                for i in pos[2]+1:pos[3]-1
                    gate_tmp00 = kron(gate_tmp00, E)
                    gate_tmp01 = kron(gate_tmp01, E)
                    gate_tmp10 = kron(gate_tmp10, E)
                    gate_tmp11 = kron(gate_tmp11, E)
                end

                # projection space
                gate_tmp00 = kron(gate_tmp00, proj0)
                gate_tmp01 = kron(gate_tmp01, proj1)
                gate_tmp10 = kron(gate_tmp10, proj0)
                gate_tmp11 = kron(gate_tmp11, proj1)

                # fill with identities
                for i in pos[3]+1:N
                    gate_tmp00 = kron(gate_tmp00, E)
                    gate_tmp01 = kron(gate_tmp01, E)
                    gate_tmp10 = kron(gate_tmp10, E)
                    gate_tmp11 = kron(gate_tmp11, E)
                end

                # complete Toffoli gate
                gate = gate_tmp00 + gate_tmp01 + gate_tmp10 + gate_tmp11
            end
        end

        # update state vector
        qc.StateVector .= gate*qc.StateVector

        # update representing matrix of quantum circuit
        for i in 1:qc.NumQubits
            if i == pos[1]
                push!(qc.Representation[i], 14) # control qubit
                push!(qc.RepresentationFull[i], 14) # control qubit
            elseif i == pos[2]
                push!(qc.Representation[i], -14) # action qubit
                push!(qc.RepresentationFull[i], -14) # action qubit
            else
                push!(qc.Representation[i], 0)
                push!(qc.RepresentationFull[i], 0)
            end
        end

        # update circuit depth
        qc.CircuitDepth += 1
    end
end
end
end

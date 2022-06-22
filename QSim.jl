
##############################
## Quantum Simulator Main File
##############################


######################################
# dependencies on other Julia packages
######################################

using ITensors # most important backbone
using IterTools
using StatsBase
using Random
using DataStructures
using LinearAlgebra
using DelimitedFiles
#using CSV
using DataFrames
using AbstractAlgebra






####################################################
# Structs defining quantum circuits and custom gates
####################################################


""" Struct defining the basic quantum circuit object that is manipulated.
Holds information about the exact state vector (as a Julia array), the number
of qubits and the circuit depth. In addition, the sequences of gates applied
to the circuit are kept track of in two representation (sorted) dictionaries
(one just accounting for high-level gates, the other saving information about
the low-level gates potentially needed to compose a high-level gate). If a
measurement is performed within the circuit, classical bits of information
are saved in the corresponding dictionary. Finally, one can choose between
working with a linear topology (requiring additional swaps for non-adjacent
2-qubit gates) and a master topology (every qubit connected to every other
qubit, i.e. 2-site gates can be applied directly through the corresponding
matrix representation). """
mutable struct QC

    # principle properties
    StateVector::Vector{ComplexF64}
    NumQubits::Int32
    CircuitDepth::Int32

    # if (single) measurements are performed, save results as classical bits
    ClassicalBits::SortedDict
    ClassicalBitsProportion::SortedDict

    # keep track of gates in circuit
    Representation::SortedDict
    RepresentationFull::SortedDict # including swaps to realise 2-site gates
    ConversionTable::Dict # holds representation of gates
    HasCustomGate::Bool

    # topology
    LinearTopology::Bool

    # coherence times of qubits/error rates
    # other?

end


""" Struct defining the basic quantum circuit object that is manipulated.
Holds information about the state vector (as an ITensor MPS). Further parameters
which are recorded are the number of qubits, circuit depth, the ITensor index
set of the wave function and the bond dimension after each gate layer (as an
array with a number of entries equal to the circuit depth). One may specify a
maximum bond dimension up to which contractions with MPOs are performed exactly,
as well as the contraction method. To assess the quality of approximation to the
true state, different fidelities are monitored and updated throughout the quantum
circuit. In addition, the sequences of gates applied to the circuit are kept
track of in two representation (sorted) dictionaries (one just accounting for
high-level gates, the other saving information about the low-level gates potentially
needed to compose a high-level gate). If a measurement is performed within the
circuit, classical bits of information are saved in the corresponding dictionary.
Finally, one can choose between working with a linear topology (requiring
additional swaps for non-adjacent 2-qubit gates) and a master topology (every
qubit connected to every other qubit, i.e. 2-site gates can be applied directly. """
mutable struct QC_IT_MPS

    # principle properties
    StateVector::MPS
    NumQubits::Int32
    CircuitDepth::Int32
    IndexSet::Vector{Index{Int64}}
    MaxBondDim::Int32
    BondDim::Array # obtain list of bond dims (for each layer)
    ContractionMethod::String # method for MPS-MPO contractions

    # fidelities
    TwoQubitFidelity::Array
    TwoQubitFidelitySelected::Array
    TwoQubitFidelityExactSelected::Array
    NQubitFidelity::Array
    NQubitFidelitySelected::Array
    NQubitFidelityExactSelected::Array
    AverageTwoQubitFidelity::Array
    AverageTwoQubitFidelitySelected::Array
    AverageTwoQubitFidelityExactSelected::Array

    # if (single) measurements are performed, save results as classical bits
    ClassicalBits::SortedDict
    ClassicalBitsProportion::SortedDict

    # keep track of gates in circuit
    Representation::SortedDict
    RepresentationFull::SortedDict # including swaps to realise 2-site gates
    ConversionTable::Dict # holds representation of gates
    HasCustomGate::Bool

    # topology
    LinearTopology::Bool

    # coherence times of qubits/error rates
    # other?

end


""" Struct allowing to group several gates into a custom gate, if
the repetition of a more complex object (such as an oracle, time-evolution
gate etc.) is required. """
mutable struct custom_gate

    # principle properties
    NumQubitsGate::Int64
    GateName::String

    # keep track of gates in circuit
    Representation::SortedDict

end


# struct for gates as well?
# might be useful for saving extra information with gates, such as error rates

#mutable struct gate
#    Matrix::Matrix{Complex}
#    ErrorRate::Float32
#end



#########################################
# basic functions: initilisation, drawing
#########################################


""" Function to build the initial state |00...0> of a quantum
circuit as an exact state vector. Needs the desired number of qubits
as input. Alternatively, create (normalised) random state. """
function initialise_state(N, random=false)

    if random
        psi = rand(Complex{Float64}, 2^N)
        psi = psi/norm(psi)
    else
        psi = Complex.(zeros(2^N))
        psi[1] = 1.0
    end

    return psi
end


""" Function to initialise a custom gate object.
Only gets meta-information about gate object, doesn't do construction! """
function initialise_custom_gate(qc, N_reg::Int64, gate_name::String)#, backend="MPS_ITensor", maxdim=100, contmethod="naive")

    if length(gate_name) > 5
        error("Chosen gate name is too long (max. 5 characters)!")
    end

    # initialise custom gate object
    CGate = custom_gate(N_reg, gate_name, SortedDict())

    # set up representation of custom gate
    for i in 1:N_reg
        CGate.Representation[i] = []
    end
    return CGate
end


""" Function to initialise a quantum circuit object with N qubits. By setting
the lintop parameter, one may specify a linear topology or a master topology
of the qubits. The backend parameter determines whether the circuit is built
from Julia arrays or ITensor MPS objects. """
function initialise_qcircuit(N, lintop=false, backend="MPS_ITensor", maxdim=100,
    contmethod="naive", random=false, randombond=16)

    # check inputs
    if typeof(N) != Int
        error("Wrong data type for number of qubits (must be Int). ")
    elseif typeof(lintop) != Bool
        error("Wrong data type for lintop parameter (need Bool). ")
    elseif backend ∉ ["ED_Julia", "MPS_ITensor"]
        error("Incorrect backend (needs ED_Julia or MPS_ITensor). ")
    elseif typeof(maxdim) != Int
        error("Wrong data type for maximal bond dimension (must be Int). ")
    elseif typeof(contmethod) != String
        error("Wrong data type for MPS-MPO contraction method (must be String). ")
    end

    # initialise state vector and construct circuit object, depending on backend
    if backend == "ED_Julia"
        psi = initialise_state(N, random)
        qcircuit = QC(psi, N, 0, SortedDict(), SortedDict(), SortedDict(), SortedDict(),
        Dict(), false, lintop)
    else
        sites = siteinds("QCircuit", N)
        if random
            psi = randomMPS(sites, linkdims=randombond)
        else
            psi = productMPS(sites, "0")
        end
        qcircuit = QC_IT_MPS(psi, N, 0, sites, maxdim, [maxlinkdim(psi)], contmethod,
        [1.0], [1.0], [1.0], [1.0], [1.0], [1.0], [1.0], [1.0], [1.0], SortedDict(), SortedDict(),
        SortedDict(), SortedDict(), Dict(), false, lintop)
    end

    # set up representation of circuit
    for i in 1:N
        qcircuit.Representation[i] = []
        qcircuit.RepresentationFull[i] = []
    end

    # initialise conversion table to print out gates
    qcircuit.ConversionTable[0] = "----"
    qcircuit.ConversionTable[1] = "-H--"
    qcircuit.ConversionTable[2] = "-X--"
    qcircuit.ConversionTable[3] = "-○--"
    qcircuit.ConversionTable[-3] = "-+--"
    qcircuit.ConversionTable[4] = "-⊗--"
    qcircuit.ConversionTable[-4] = "-⊗--"
    qcircuit.ConversionTable[5] = "-Z--"
    qcircuit.ConversionTable[6] = "-Y--"
    qcircuit.ConversionTable[7] = "-S--"
    qcircuit.ConversionTable[8] = "-T--"
    qcircuit.ConversionTable[9] = "-Rˣ-"
    qcircuit.ConversionTable[10] = "-Rʸ-"
    qcircuit.ConversionTable[11] = "-Rᶻ-"
    qcircuit.ConversionTable[12] = "-P--"
    qcircuit.ConversionTable[13] = "-○--"
    qcircuit.ConversionTable[-13] = "-U--"
    qcircuit.ConversionTable[14] = "-○--"
    qcircuit.ConversionTable[-14] = "-+--"
    qcircuit.ConversionTable[15] = "-|--"
    qcircuit.ConversionTable[16] = "-M--"
    qcircuit.ConversionTable[17] = "---M"
    qcircuit.ConversionTable[18] = "-√X-"
    qcircuit.ConversionTable[19] = "-U--"
    qcircuit.ConversionTable[20] = "-|  "
    qcircuit.ConversionTable[21] = " |--"
    qcircuit.ConversionTable[22] = "-+++"
    qcircuit.ConversionTable[23] = "++--"
    qcircuit.ConversionTable[24] = "-○--"
    qcircuit.ConversionTable[-24] = "-Rn-"
    qcircuit.ConversionTable[25] = "-U--"
    qcircuit.ConversionTable[-25] = "-U--"

    return qcircuit
end


""" Function that takes in the meta-information about the quantum
circuit (after it has been constructed) and prints out a (very) schematic
representation in the console. Works for the different backends, as only the
metainformation about the circuit is needed is draw it on the terminal.
TO DO: maybe extend to be (pretty-) printed into a file, see matplotlib. """
function draw(qc::QC, full=false)

    # prompt for file name if desired to save in txt file!

    # set maximum length of terminal up to which the circuit can be printed
    maxprintlength = 34

    # check desired printing format
    if full
        qc_rep = qc.RepresentationFull
        println("Schematic representation of full quantum circuit, depth = $(qc.CircuitDepth): ")
        if length(collect(values(qc_rep))[1]) > maxprintlength
            println("Circuit too long to be printed in terminal!")
            println("Cutting after $maxprintlength gates.")
        end
    else
        qc_rep = qc.Representation
        println("Schematic representation of quantum circuit, depth = $(qc.CircuitDepth): ")
        if length(collect(values(qc_rep))[1]) > maxprintlength
            println("Circuit too long to be printed in terminal!")
            println("Cutting after $maxprintlength gates.")
        end
    end

    # loop through gate matrix and look up corresponding representation
    println(" ")
    for (qubit, gates) in qc_rep
        if qubit < 10
            print(qubit, "   |0⟩ ")
        elseif qubit >= 10 && qubit < 100
            print(qubit, "  |0⟩ ")
        else
            print(qubit, " |0⟩ ")
        end
        if length(gates) < maxprintlength
            for i in 1:length(gates)
                if i == length(gates)
                    println(qc.ConversionTable[gates[i]])
                else
                    print(qc.ConversionTable[gates[i]])
                end
            end
        else
            for i in 1:maxprintlength
                if i == maxprintlength
                    println(qc.ConversionTable[gates[i]])
                else
                    print(qc.ConversionTable[gates[i]])
                end
            end
        end
    end
    println("\n")
end


""" Function that takes in the meta-information about the quantum
circuit (after it has been constructed) and prints out a (very) schematic
representation in the console. Works for the different backends, as only the
metainformation about the circuit is needed is draw it on the terminal.
TO DO: maybe extend to be (pretty-) printed into a file, see matplotlib. """
function draw(qc::QC_IT_MPS, full=false)

    # prompt for file name if desired to save in txt file!

    # set maximum length of terminal up to which the circuit can be printed
    maxprintlength = 34

    # check desired printing format
    if full
        qc_rep = qc.RepresentationFull
        println("\n")
        println("Schematic representation of full quantum circuit, depth = $(qc.CircuitDepth): ")
        println("Maximum possible bond dimension: $(2^(qc.NumQubits÷2))")
        println("Maximum allowed bond dimension before truncation: $(qc.MaxBondDim)")
        if length(collect(values(qc_rep))[1]) > maxprintlength
            println("Circuit too long to be printed in terminal!")
            println("Cutting after $maxprintlength gates.")
        end
    else
        qc_rep = qc.Representation
        println("\n")
        println("Schematic representation of quantum circuit, depth = $(qc.CircuitDepth): ")
        println("Maximum possible bond dimension: $(2^(qc.NumQubits÷2))")
        println("Maximum allowed bond dimension before truncation: $(qc.MaxBondDim)")
        if length(collect(values(qc_rep))[1]) > maxprintlength
            println("Circuit too long to be printed in terminal!")
            println("Cutting after $maxprintlength gates.")
        end
    end

    # loop through gate matrix and look up corresponding representation
    println(" ")
    for (qubit, gates) in qc_rep
        if qubit < 10
            print(qubit, "   |0⟩ ")
        elseif qubit >= 10 && qubit < 100
            print(qubit, "  |0⟩ ")
        else
            print(qubit, " |0⟩ ")
        end
        if length(gates) < maxprintlength
            for i in 1:length(gates)
                if i == length(gates)
                    if qubit < qc.NumQubits
                        print(qc.ConversionTable[gates[i]])
                        println("  $(ITensors.dim(linkind(qc.StateVector, qubit)))")
                    else
                        println(qc.ConversionTable[gates[i]])
                    end
                else
                    print(qc.ConversionTable[gates[i]])
                end
            end
        else
            for i in 1:maxprintlength
                if i == maxprintlength
                    if qubit < qc.NumQubits
                        print(qc.ConversionTable[gates[i]])
                        println("  $(ITensors.dim(linkind(qc.StateVector, qubit)))")
                    else
                        println(qc.ConversionTable[gates[i]])
                    end
                else
                    print(qc.ConversionTable[gates[i]])
                end
            end
        end
    end
    if qc.HasCustomGate == false
        print("Chi:  ")
        if length(qc.BondDim) < maxprintlength
            for i in 1:length(qc.BondDim)
                if qc.BondDim[i] < 10
                    print(" $(qc.BondDim[i])  ")
                elseif qc.BondDim[i] >= 10 && qc.BondDim[i] < 100
                    print(" $(qc.BondDim[i]) ")
                else
                    print(" $(qc.BondDim[i])")
                end
            end
        else
            for i in 1:maxprintlength
                if qc.BondDim[i] < 10
                    print(" $(qc.BondDim[i])  ")
                elseif qc.BondDim[i] >= 10 && qc.BondDim[i] < 100
                    print(" $(qc.BondDim[i]) ")
                else
                    print(" $(qc.BondDim[i])")
                end
            end
        end
    end
    println("\n")
    println("Vertical rightmost column: bond dimension throughout final MPS (bond between current line and next line).")
    println("Lowest horizontal line: max bond dimension after each layer of gates.\n")
end



#####################
# auxiliary functions
#####################


""" Function to filter all the binary number from a list of binary numbers,
which have a certain configuration in given positions. E.g. for the list
[[000], [001], [010], [011], [100], [101], [110], [111]], the filtered
numbers with configuration [11] in positions [1, 3] are [101] and [111].
Note that the binary numbers should be lists themselves (i.e. [1, 1, 1]
instead of [111]) and have only been shortened for reasons of brevity. """
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


""" Function to convert an array of bits into the corresponding base-10
integer. The array represents the base-2 number from left to right, i.e.
[1, 0, 1] corresponds to 1*2² + 0*2¹ + 1*2⁰ = 5 in base 10."""
function bit_array_to_int(arr)
    s = 0
    N = length(arr)
    for i in 1:N
        s += arr[i] * 2^(N - i)
    end
    return s
end


""" Function to extract the 2-permutations from a general array of positions,
pos = [pos[1], pos[2]], such that the permutation represented by pos is
equivalent to the composition of those 2-permutations. Example:
(2, 5) = (2, 3)(3, 4)(4, 5). Note: here given in permutation (i.e. tuple)
notation. """
function get_2_perms_old(pos)
    perms = []
    for i in (pos[1]+1):pos[2]
        push!(perms, [i-1, i])
    end
    return perms
end


""" """
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


""" Function to create any given state in the computational basis as
MPS. To generate the state |01101>, specify the configuration
[0, 1, 1, 0, 1]. Need also specify the site indices the MPS is
constructed with. """
function MPS_computationalBasis(sites, configuration::Array)

    # initialise state
    basis_state = productMPS(sites, "0")

    # create desired state in comp. basis by applying X gates in resp. positions
    for i in 1:length(configuration)
        if configuration[i] == 1
            gate = op("X", sites[i])
            updatedTensor = gate*basis_state[i]
            noprime!(updatedTensor)
            basis_state[i] = updatedTensor
        end
    end
    return basis_state
end


#######################################
# Basic plotting of measurement results
#######################################


""" Auxiliary function to get a list of bitstrings representing a t-qubit
register. """
function get_bitstrings(t)
    bitstrings = []
    for i in 0:2^t-1
        push!(bitstrings, reverse(digits(i, base=2, pad=t)))
    end
    return bitstrings
end


""" Function to extract a histogram from the dictionary ClassicalBitsProportion,
where the different measurement results of the quantum circuit are saved. Returns
a list of the states (identified by their decimal number) and the corresponding
measurement frequencies. """
function get_measurement_histogram(qc, t)

    # get measurement results, set up containers
    #probs_tmp = collect(keys(qc.ClassicalBitsProportion))
    #states_tmp = collect(values(qc.ClassicalBitsProportion))
    states_tmp = collect(keys(qc.ClassicalBitsProportion))
    probs_tmp = collect(values(qc.ClassicalBitsProportion))
    states = zeros(2^t)
    probs = zeros(2^t)
    bitstrings = get_bitstrings(t)

    for i in 1:2^t
        states[i] = i
        if bitstrings[i] ∈ states_tmp
            ind = findfirst(item -> item == bitstrings[i], states_tmp)
            probs[i] = probs_tmp[ind]
        end
    end
    return Int.(states), probs
end


""" Function to create a simple histogram with measurement results from
a given quantum circuit object. """
function make_histogram(qc::QC, t, title::String, path::String)

    # get measurement results and corresponding frequencies
    states, freqs = get_measurement_histogram(qc, t)

    # make (rudimentary) plot; save
    p = bar(states, freqs, label="measurement")
    xlabel!("states (decimal expansion)")
    ylabel!("frequency")
    title!(title)
    savefig(p, path)

end


""" Function to create a simple histogram with measurement results from
a given quantum circuit object. MPS version also prints max bond dimension
on plot. """
function make_histogram(qc::QC_IT_MPS, t, title::String, path::String, maxdim::Int)

    # get measurement results and corresponding frequencies
    states, freqs = get_measurement_histogram(qc, t)

    # make (rudimentary) plot; save
    p = bar(states, freqs, label="measurement")
    xlabel!("states (decimal expansion)")
    ylabel!("frequency")
    title!(title)
    annotate!((100, 0.2, "χ = $(maxdim)"))
    savefig(p, path)

end


###############################################
# Imports of quantum simulation functionalities
###############################################

###############
# Hilbert space
###############

include("Hilbert_space_QC.jl")


##############
# Measurements
##############

include("measurement.jl")


####################
# Single-qubit gates
####################

include("Functions/single_qubit.jl")
include("Functions/single_qubit_MPS.jl")


#################
# Two-qubit gates
#################

include("Functions/two_qubit.jl")
include("Functions/two_qubit_MPS.jl")


###################
# Three-qubit gates
###################

include("Functions/three_qubit.jl")
include("Functions/three_qubit_MPS.jl")


############################
# Arbitrary controlled gates
############################

include("Functions/arbitrary_CU_MPS.jl")


##############
# Custom gates
##############

include("Functions/custom_gate.jl")


#############
# Subroutines
#############

include("Functions/subroutines.jl")


#################################
# Exact evaluations of algorithms
#################################

include("exact_algorithms.jl")

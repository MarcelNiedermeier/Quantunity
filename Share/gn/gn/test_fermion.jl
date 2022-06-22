
######################################
## Test script for fermionic operators
######################################

include("../../QSim.jl")



N = 8+1 # total number of qubits
maxdim = 20 # maximum allowed bond dimension for MPS
contmethod = "naive" # contraction method for MPS (not so important right now)
random = false # no random intial state for circuit
lintop = false # possibility to use linear qubit topology, also don't touch for now :)
randombond = 10 # bond dimension of random initial MPS, if needed'
dt = 0.5 #timestep
last_qubit = N

#some constants
r = 1
k = 0.5
m = 1
g = 1.01

backend_MPS = "MPS_ITensor" # calculate with MPS (any higher number of qubits)




""" Function implementing a fermionic creation operator, encoded by
a Jordan Wigner transformation. Specify position by giving a single
integer pos. """
function creation_operator!(qc, pos::Int, update_rep=true)

    mpo_array1 = []
    mpo_array2 = []

    push!(mpo_array1, 1/2)
    push!(mpo_array2, -1.0im/2)

    for i in 1:pos-1
        push!(mpo_array1, "Z")
        push!(mpo_array2, "Z")
        push!(mpo_array1, i)
        push!(mpo_array2, i)
    end
    push!(mpo_array1, "X")
    push!(mpo_array2, "Y")
    push!(mpo_array1, pos)
    push!(mpo_array2, pos)

    mpo_tuple1 = tuple(mpo_array1...)
    mpo_tuple2 = tuple(mpo_array2...)

    # initialise creation operator
    creation_op = OpSum()
    creation_op += mpo_tuple1
    creation_op += mpo_tuple2

    # get sites
    sites = qc.IndexSet

    # turn into MPO
    gate = MPO(creation_op, sites)

    # apply to state vector, unprime
    qc.StateVector = contract(gate, qc.StateVector, maxdim=qc.MaxBondDim)#, method="naive")
    noprime!(qc.StateVector)

    # update representing matrix of quantum circuit
    #if update_rep
    #   update_representation_single_site!(qc, pos, 1)
    #end

    # update cicruit depth and bond dimension
    qc.CircuitDepth += 1
    push!(qc.BondDim, maxlinkdim(qc.StateVector))
end


""" Function implementing a fermionic annhilation operator, encoded by
a Jordan Wigner transformation. Specify position by giving a single
integer pos. """
function annihilation_operator!(qc, pos::Int, update_rep=true)

    mpo_array1 = []
    mpo_array2 = []

    push!(mpo_array1, 1/2)
    push!(mpo_array2, 1.0im/2)

    for i in 1:pos-1
        push!(mpo_array1, "Z")
        push!(mpo_array2, "Z")
        push!(mpo_array1, i)
        push!(mpo_array2, i)
    end
    push!(mpo_array1, "X")
    push!(mpo_array2, "Y")
    push!(mpo_array1, pos)
    push!(mpo_array2, pos)

    mpo_tuple1 = tuple(mpo_array1...)
    mpo_tuple2 = tuple(mpo_array2...)

    # initialise annihilation operator
    annihilation_op = OpSum()
    annihilation_op += mpo_tuple1
    annihilation_op += mpo_tuple2

    # get sites
    sites = qc.IndexSet

    # turn into MPO
    gate = MPO(annihilation_op, sites)

    # apply to state vector, unprime
    qc.StateVector = contract(gate, qc.StateVector, maxdim=qc.MaxBondDim)#, method="naive")
    noprime!(qc.StateVector)

    # update representing matrix of quantum circuit
    #if update_rep
    #   update_representation_single_site!(qc, pos, 1)
    #end

    # update cicruit depth and bond dimension
    qc.CircuitDepth += 1
    push!(qc.BondDim, maxlinkdim(qc.StateVector))
end


# test functions with some random settings

qc_MPS = initialise_qcircuit(N, lintop, backend_MPS, maxdim,
contmethod, random, randombond)

pos = 4
creation_operator!(qc_MPS, pos)
creation_operator!(qc_MPS, pos+2)
annihilation_operator!(qc_MPS, pos)


# check that we still have properly normalised states
println(ITensors.norm(qc_MPS.StateVector))

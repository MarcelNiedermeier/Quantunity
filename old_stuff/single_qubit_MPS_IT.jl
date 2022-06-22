
#################################################################
## Single qubit gates for Quantum Simulator - MPS ITensor version
#################################################################


#####################
# Auxiliary functions
#####################


""" Auxiliary function to apply a given single-site gate "gateName" (as
specified in the Hilbert space) to a quantum circuit in the positions
given in the array pos. """
function apply_single_site_gates!(qc::QC_IT_MPS, pos::Array{Int64, 1}, gateName::String)

    # loop through given positions and apply gate
    for i in pos
        gate = op(gateName, qc.IndexSet[i])
        updatedTensor = gate*qc.StateVector[i]
        noprime!(updatedTensor)
        qc.StateVector[i] = updatedTensor
    end
end


""" Auxiliary function to apply a given single-site gate "gateName" (as
specified in the Hilbert space, with additional parameter theta) to a
quantum circuit in the positions given in the array pos. """
function apply_single_site_gates!(qc::QC_IT_MPS, pos::Array{Int64, 1}, gateName::String, theta::Float64)

    # loop through given positions and apply gate
    for i in pos
        gate = op(gateName, qc.IndexSet[i]; theta)
        updatedTensor = gate*qc.StateVector[i]
        noprime!(updatedTensor)
        qc.StateVector[i] = updatedTensor
    end
end


""" Function to generate an MPO representing a single-site gate applied
on arbitrary positions given in the array pos. E.g. specifying "X" on
pos[1, 3] for a 4-qubit system constructs the MPO X ⊗ E ⊗ X ⊗ E. """
function single_site_MPO(qc::QC_IT_MPS, gate::String, pos::Array{Int64, 1})

    # get indices in QC Hilbert space, build MPO for single-site operator
    MPO_list = []
    for i in 1:qc.NumQubits
        if i in pos
            push!(MPO_list, gate)
        else
            push!(MPO_list, "E")
        end
    end

    return MPO(qc.IndexSet, String.(MPO_list))
end


""" Function to generate an MPO representing a single-site gate applied
on arbitrary positions given in the array pos. E.g. specifying "X" on
pos[1, 3] for a 4-qubit system constructs the MPO X ⊗ E ⊗ X ⊗ E. """
#function single_site_MPO(qc::QC_IT_MPS, gate::String, pos::Array{Int64, 1}, theta::Float64)
#
#    # get indices in QC Hilbert space, build MPO for single-site operator
#    MPO_list = Vector{ITensor}([])
#    for i in 1:qc.NumQubits
#        if i in pos
#            push!(MPO_list, op(gate, qc.IndexSet[i]; theta))
#        else
#            push!(MPO_list, op("E", qc.IndexSet[i]))
#        end
#    end
#    println("MPO list: ", typeof(MPO_list))
#
#    return MPO(MPO_list)
#end


""" Function to update the representation of a quantum circuit qc for a
gate (specified by num) applied positions given through the array pos. """
function update_representation_single_site!(qc::QC_IT_MPS, pos::Array{Int64, 1}, num::Int64)
    for i in 1:qc.NumQubits
        if i in pos
            push!(qc.Representation[i], num)
            push!(qc.RepresentationFull[i], num)
        else
            push!(qc.Representation[i], 0)
            push!(qc.RepresentationFull[i], 0)
        end
    end
end


#######
# Gates
#######


""" Function to apply the (single-qubit) Hadamard gate to a quantum circuit qc
in every position specified in the array pos. """
function hadamard!(qc::QC_IT_MPS, pos::Array{Int64, 1})

   # if quantum circuit has linear topology: apply single-site operators
   #if qc.LinearTopology == true

   # apply Hadamard gate(s) to quantum circuit
   apply_single_site_gates!(qc, pos, "H")

   # if quantum circuit has master topology: apply in MPO form
   #else

   #   # get MPO of single site gate, apply to state
   #   gate = single_site_MPO(qc, "H", pos)
   #   qc.StateVector = contract(gate, qc.StateVector, method="naive")

   #end

   # update representing matrix of quantum circuit
   update_representation_single_site!(qc, pos, 1)

   # update circuit depth and bond dimension
   qc.CircuitDepth += 1
   push!(qc.BondDim, maxlinkdim(qc.StateVector))
end


""" Function to apply the (single-qubit) Pauli X gate to a quantum circuit qc
in every position specified in the array pos. """
function PauliX!(qc::QC_IT_MPS, pos::Array{Int64, 1})

   # if quantum circuit has linear topology: apply single-site operators
   #if qc.LinearTopology == true

   # apply Pauli X gate(s) to quantum circuit
   apply_single_site_gates!(qc, pos, "X")

   # if quantum circuit has master topology: apply in MPO form
   #else

   #   # get MPO of single site gate, apply to state
   #   gate = single_site_MPO(qc, "X", pos)
   #   qc.StateVector = contract(gate, qc.StateVector, method="naive")

   #end

   # update representing matrix of quantum circuit
   update_representation_single_site!(qc, pos, 2)

   # update circuit depth and bond dimension
   qc.CircuitDepth += 1
   push!(qc.BondDim, maxlinkdim(qc.StateVector))
end


""" Function to apply the (single-qubit) Pauli Y gate to a quantum circuit qc
in every position specified in the array pos. """
function PauliY!(qc::QC_IT_MPS, pos::Array{Int64, 1})

   # if quantum circuit has linear topology: apply single-site operators
   #if qc.LinearTopology == true

   # apply Pauli Y gate(s) to quantum circuit
   apply_single_site_gates!(qc, pos, "Y")

   # if quantum circuit has master topology: apply in MPO form
   #else

   #   # get MPO of single site gate, apply to state
   #   gate = single_site_MPO(qc, "Y", pos)
   #   qc.StateVector = contract(gate, qc.StateVector, method="naive")

   #end

   # update representing matrix of quantum circuit
   update_representation_single_site!(qc, pos, 6)

   # update circuit depth and bond dimension
   qc.CircuitDepth += 1
   push!(qc.BondDim, maxlinkdim(qc.StateVector))
end


""" Function to apply the (single-qubit) Pauli Z gate to a quantum circuit qc
in every position specified in the array pos. """
function PauliZ!(qc::QC_IT_MPS, pos::Array{Int64, 1})

   # if quantum circuit has linear topology: apply single-site operators
   #if qc.LinearTopology == true

   # apply Pauli Z gate(s) to quantum circuit
   apply_single_site_gates!(qc, pos, "Z")

   # if quantum circuit has master topology: apply in MPO form
   #else

   #   # get MPO of single site gate, apply to state
   #   gate = single_site_MPO(qc, "Z", pos)
   #   qc.StateVector = contract(gate, qc.StateVector, method="naive")

   #end

   # update representing matrix of quantum circuit
   update_representation_single_site!(qc, pos, 5)

   # update circuit depth and bond dimension
   qc.CircuitDepth += 1
   push!(qc.BondDim, maxlinkdim(qc.StateVector))
end


""" Function to apply the (single-qubit) √X gate to a quantum circuit qc
in every position specified in the array pos. """
function SqrtX!(qc::QC_IT_MPS, pos::Array{Int64, 1})

   # if quantum circuit has linear topology: apply single-site operators
   #if qc.LinearTopology == true

   # apply √X gate(s) to quantum circuit
   apply_single_site_gates!(qc, pos, "√X")

   # if quantum circuit has master topology: apply in MPO form
   #else

   #   # get MPO of single site gate, apply to state
   #   gate = single_site_MPO(qc, "√X", pos)
   #   qc.StateVector = contract(gate, qc.StateVector, method="naive")

   #end

   # update representing matrix of quantum circuit
   update_representation_single_site!(qc, pos, 18)

   # update circuit depth and bond dimension
   qc.CircuitDepth += 1
   push!(qc.BondDim, maxlinkdim(qc.StateVector))
end


""" Function to apply the (single-qubit) S gate to a quantum circuit qc
in every position specified in the array pos. """
function SGate!(qc::QC_IT_MPS, pos::Array{Int64, 1})

   # if quantum circuit has linear topology: apply single-site operators
   #if qc.LinearTopology == true

   # apply S gate(s) to quantum circuit
   apply_single_site_gates!(qc, pos, "S")

   # if quantum circuit has master topology: apply in MPO form
   #else

   #   # get MPO of single site gate, apply to state
   #   gate = single_site_MPO(qc, "S", pos)
   #   qc.StateVector = contract(gate, qc.StateVector, method="naive")

   #end

   # update representing matrix of quantum circuit
   update_representation_single_site!(qc, pos, 7)

   # update circuit depth and bond dimension
   qc.CircuitDepth += 1
   push!(qc.BondDim, maxlinkdim(qc.StateVector))
end


""" Function to apply the (single-qubit) T gate to a quantum circuit qc
in every position specified in the array pos. """
function TGate!(qc::QC_IT_MPS, pos::Array{Int64, 1})

   # if quantum circuit has linear topology: apply single-site operators
   #if qc.LinearTopology == true

   # apply T gate(s) to quantum circuit
   apply_single_site_gates!(qc, pos, "T")

   # if quantum circuit has master topology: apply in MPO form
   #else

   #   # get MPO of single site gate, apply to state
   #   gate = single_site_MPO(qc, "T", pos)
   #   qc.StateVector = contract(gate, qc.StateVector, method="naive")

   #end

   # update representing matrix of quantum circuit
   update_representation_single_site!(qc, pos, 8)

   # update circuit depth and bond dimension
   qc.CircuitDepth += 1
   push!(qc.BondDim, maxlinkdim(qc.StateVector))
end


""" Function to apply the (single-qubit) Rx gate to a quantum circuit qc
in every position specified in the array pos. Performs a rotation around
the x-axis specified by theta."""
function RXGate!(qc::QC_IT_MPS, pos::Array{Int64, 1}, theta::Float64)

   # if quantum circuit has linear topology: apply single-site operators
   #if qc.LinearTopology == true

   # apply Rx gate(s) to quantum circuit
   apply_single_site_gates!(qc, pos, "Rx", theta)

   # if quantum circuit has master topology: apply in MPO form
   #else

   #   # apply Hadamard gate(s) to quantum circuit
   #   apply_single_site_gates!(qc, pos, "Rx", theta)

   #   # get MPO of single site gate, apply to state
   #   #gate = single_site_MPO(qc, "Rx", pos, theta)
   #   #println("check MPO: ", gate)
   #   #qc.StateVector = contract(gate, qc.StateVector)#, method="naive")

   #end

   # update representing matrix of quantum circuit
   update_representation_single_site!(qc, pos, 9)

   # update circuit depth and bond dimension
   qc.CircuitDepth += 1
   push!(qc.BondDim, maxlinkdim(qc.StateVector))
end


""" Function to apply the (single-qubit) Ry gate to a quantum circuit qc
in every position specified in the array pos. Performs a rotation around
the x-axis specified by theta."""
function RYGate!(qc::QC_IT_MPS, pos::Array{Int64, 1}, theta::Float64)

   # if quantum circuit has linear topology: apply single-site operators
   #if qc.LinearTopology == true

   # apply Ry gate(s) to quantum circuit
   apply_single_site_gates!(qc, pos, "Ry", theta)

   # if quantum circuit has master topology: apply in MPO form
   #else

   #   # apply Hadamard gate(s) to quantum circuit
   #   apply_single_site_gates!(qc, pos, "Ry", theta)

   #   # get MPO of single site gate, apply to state
   #   #gate = single_site_MPO(qc, "Ry", pos)
   #   #qc.StateVector = contract(gate, qc.StateVector, method="naive")

   #end

   # update representing matrix of quantum circuit
   update_representation_single_site!(qc, pos, 10)

   # update circuit depth and bond dimension
   qc.CircuitDepth += 1
   push!(qc.BondDim, maxlinkdim(qc.StateVector))
end


""" Function to apply the (single-qubit) Rz gate to a quantum circuit qc
in every position specified in the array pos. Performs a rotation around
the x-axis specified by theta."""
function RZGate!(qc::QC_IT_MPS, pos::Array{Int64, 1}, theta::Float64)

   # if quantum circuit has linear topology: apply single-site operators
   #if qc.LinearTopology == true

   # apply Rz gate(s) to quantum circuit
   apply_single_site_gates!(qc, pos, "Rz", theta)

   # if quantum circuit has master topology: apply in MPO form
   #else

   #   # apply Hadamard gate(s) to quantum circuit
   #   apply_single_site_gates!(qc, pos, "Rz", theta)

   #   # get MPO of single site gate, apply to state
   #   #gate = single_site_MPO(qc, "Rz", pos)
   #   #qc.StateVector = contract(gate, qc.StateVector, method="naive")

   #end

   # update representing matrix of quantum circuit
   update_representation_single_site!(qc, pos, 11)

   # update circuit depth and bond dimension
   qc.CircuitDepth += 1
   push!(qc.BondDim, maxlinkdim(qc.StateVector))
end


""" Function to apply a phase shift by the angle theta to a quantum circuit
qc in every position specified in the array pos. """
function PhaseShift!(qc::QC_IT_MPS, pos::Array{Int64, 1}, theta::Float64)

   # if quantum circuit has linear topology: apply single-site operators
   #if qc.LinearTopology == true

   # apply P gate(s) to quantum circuit
   apply_single_site_gates!(qc, pos, "P", theta)

   # if quantum circuit has master topology: apply in MPO form
   #else

   #   # apply Hadamard gate(s) to quantum circuit
   #   apply_single_site_gates!(qc, pos, "P", theta)

   #   # get MPO of single site gate, apply to state
   #   #gate = single_site_MPO(qc, "P", pos)
   #   #qc.StateVector = contract(gate, qc.StateVector, method="naive")

   #end

   # update representing matrix of quantum circuit
   update_representation_single_site!(qc, pos, 12)

   # update circuit depth and bond dimension
   qc.CircuitDepth += 1
   push!(qc.BondDim, maxlinkdim(qc.StateVector))
end

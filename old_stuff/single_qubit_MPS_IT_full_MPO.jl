
################################################################################
## Single qubit gates for Quantum Simulator - MPS ITensor version with full MPOs
################################################################################


""" Function to generate an MPO representing a single-site gate applied
on arbitrary positions given in the array pos. E.g. specifying "X" on
pos[1, 3] for a 4-qubit system constructs the MPO X ⊗ E ⊗ X ⊗ E. """
function single_site_MPO(qc::QC_IT_MPS, gate::String, pos::Array{Int64, 1})

    # get indices in QC Hilbert space, build MPO for single-site operator
    #sites = siteinds("QCircuit", N)
    #sites = qc.IndexSet
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



#######
# Gates
#######


""" Function to apply the (single-qubit) Hadamard gate to a quantum circuit qc
in every position specified in the array pos. """
function hadamard!(qc::QC_IT_MPS, pos::Array{Int64, 1})

   # get MPO of single site gate, apply to state
   gate = single_site_MPO(qc, "H", pos)
   qc.StateVector = contract(gate, qc.StateVector, method="naive")

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

   # update circuit depth and bond dimension
   qc.CircuitDepth += 1
   push!(qc.BondDim, maxlinkdim(qc.StateVector))
end


""" Function to apply the (single-qubit) Pauli X gate to a quantum circuit qc
in every position specified in the array pos. """
function PauliX!(qc::QC_IT_MPS, pos::Array{Int64, 1})

   # get MPO of single site gate, apply to state
   gate = single_site_MPO(qc, "X", pos)
   qc.StateVector = contract(qc.StateVector, gate, method="naive")

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

   # update circuit depth and bond dimension
   qc.CircuitDepth += 1
   push!(qc.BondDim, maxlinkdim(qc.StateVector))
end


""" Function to apply the (single-qubit) Pauli Y gate to a quantum circuit qc
in every position specified in the array pos. """
function PauliY!(qc::QC_IT_MPS, pos::Array{Int64, 1})

   gate = single_site_MPO(qc, "Y", pos)
   qc.StateVector = contract(qc.StateVector, gate, method="naive")

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

   # update circuit depth and bond dimension
   qc.CircuitDepth += 1
   push!(qc.BondDim, maxlinkdim(qc.StateVector))
end


""" Function to apply the (single-qubit) Pauli Z gate to a quantum circuit qc
in every position specified in the array pos. """
function PauliZ!(qc::QC_IT_MPS, pos::Array{Int64, 1})

   gate = single_site_MPO(qc, "Z", pos)
   qc.StateVector = contract(qc.StateVector, gate, method="naive")

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

   # update circuit depth and bond dimension
   qc.CircuitDepth += 1
   push!(qc.BondDim, maxlinkdim(qc.StateVector))
end


""" Function to apply the (single-qubit) √X gate to a quantum circuit qc
in every position specified in the array pos. """
function SqrtX!(qc::QC_IT_MPS, pos::Array{Int64, 1})

   gate = single_site_MPO(qc, "√X", pos)
   qc.StateVector = contract(qc.StateVector, gate, method="naive")

   # update representing matrix of quantum circuit
   for i in 1:qc.NumQubits
       if i in pos
           push!(qc.Representation[i], 18)
           push!(qc.RepresentationFull[i], 18)
       else
           push!(qc.Representation[i], 0)
           push!(qc.RepresentationFull[i], 0)
       end
   end

   # update circuit depth and bond dimension
   qc.CircuitDepth += 1
   push!(qc.BondDim, maxlinkdim(qc.StateVector))
end


""" Function to apply the (single-qubit) S gate to a quantum circuit qc
in every position specified in the array pos. """
function SGate!(qc::QC_IT_MPS, pos::Array{Int64, 1})

   gate = single_site_MPO(qc, "S", pos)
   qc.StateVector = contract(qc.StateVector, gate, method="naive")

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

   # update circuit depth and bond dimension
   qc.CircuitDepth += 1
   push!(qc.BondDim, maxlinkdim(qc.StateVector))
end


""" Function to apply the (single-qubit) T gate to a quantum circuit qc
in every position specified in the array pos. """
function TGate!(qc::QC_IT_MPS, pos::Array{Int64, 1})

   gate = single_site_MPO(qc, "T", pos)
   qc.StateVector = contract(qc.StateVector, gate, method="naive")

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

   # update circuit depth and bond dimension
   qc.CircuitDepth += 1
   push!(qc.BondDim, maxlinkdim(qc.StateVector))
end


""" Function to apply the (single-qubit) Rx gate to a quantum circuit qc
in every position specified in the array pos. Performs a rotation around
the x-axis specified by theta."""
function RXGate!(qc::QC_IT_MPS, pos::Array{Int64, 1}, theta)

   gate = single_site_MPO(qc, "Rx", pos)
   qc.StateVector = contract(qc.StateVector, gate, method="naive")

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

   # update circuit depth and bond dimension
   qc.CircuitDepth += 1
   push!(qc.BondDim, maxlinkdim(qc.StateVector))
end


""" Function to apply the (single-qubit) Ry gate to a quantum circuit qc
in every position specified in the array pos. Performs a rotation around
the x-axis specified by theta."""
function RYGate!(qc::QC_IT_MPS, pos::Array{Int64, 1}, theta)

   gate = single_site_MPO(qc, "Ry", pos)
   qc.StateVector = contract(qc.StateVector, gate, method="naive")

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

   # update circuit depth and bond dimension
   qc.CircuitDepth += 1
   push!(qc.BondDim, maxlinkdim(qc.StateVector))
end


""" Function to apply the (single-qubit) Rz gate to a quantum circuit qc
in every position specified in the array pos. Performs a rotation around
the x-axis specified by theta."""
function RZGate!(qc::QC_IT_MPS, pos::Array{Int64, 1}, theta)

   gate = single_site_MPO(qc, "Rz", pos)
   qc.StateVector = contract(qc.StateVector, gate, method="naive")

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

   # update circuit depth and bond dimension
   qc.CircuitDepth += 1
   push!(qc.BondDim, maxlinkdim(qc.StateVector))
end


""" Function to apply a phase shift by the angle theta to a quantum circuit
qc in every position specified in the array pos. """
function PhaseShift!(qc::QC_IT_MPS, pos::Array{Int64, 1}, theta)

   gate = single_site_MPO(qc, "P", pos)
   qc.StateVector = contract(qc.StateVector, gate, method="naive")

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

   # update circuit depth and bond dimension
   qc.CircuitDepth += 1
   push!(qc.BondDim, maxlinkdim(qc.StateVector))
end

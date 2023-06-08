
#################################################################
## Single qubit gates for Quantum Simulator - MPS ITensor version
#################################################################


#####################
# Auxiliary functions
#####################


""" Auxiliary function to apply a given single-site gate "gateName" (as
specified in the Hilbert space) to a quantum circuit in the positions
given in the array pos. """
function apply_single_site_gates!(qc::QC_IT_MPS, pos, gateName::String)

    # loop through given positions and apply gate
    for i in pos
        gate = op(gateName, qc.IndexSet[i])
        updatedTensor = gate*qc.StateVector[i]
        noprime!(updatedTensor)
        qc.StateVector[i] = updatedTensor
    end
end


""" Auxiliary function to apply a given single-site gate "gateName" (as
specified in the Hilbert space, with additional parameter θ) to a
quantum circuit in the positions given in the array pos. """
function apply_single_site_gates!(qc::QC_IT_MPS, pos,
   gateName::String, θ::Number)

    # loop through given positions and apply gate
    for i in pos
        gate = op(gateName, qc.IndexSet[i]; θ=θ)
        updatedTensor = gate*qc.StateVector[i]
        noprime!(updatedTensor)
        qc.StateVector[i] = updatedTensor
    end
end


""" Auxiliary function to apply a given single-site gate "gateName" (as
specified in the Hilbert space, with additional parameters α, β, γ, δ) to a
quantum circuit in the positions given in the array pos. """
function apply_single_site_gates!(qc::QC_IT_MPS, pos, gateName::String,
   α::Number, β::Number, γ::Number, δ::Number)

    # check unitarity
    if gateName == "U"
       s = Index(2, "QCircuit")
       E = array(op("E", s))
       U = array(op("U", s; α=α, β=β, γ=γ, δ=δ))
       if norm(U*adjoint(U) - E) > 1e-14
           error("The matrix you have entered is not unitary, check parameters.")
       end
    end

    # loop through given positions and apply gate
    for i in pos
        gate = op(gateName, qc.IndexSet[i]; α=α, β=β, γ=γ, δ=δ)
        updatedTensor = gate*qc.StateVector[i]
        noprime!(updatedTensor)
        qc.StateVector[i] = updatedTensor
    end
end


""" Function to update the representation of a quantum circuit qc for a
gate (specified by num) applied positions given through the array pos. """
function update_representation_single_site!(qc::QC_IT_MPS, pos,
   num::Int64)
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
function Hadamard!(qc::QC_IT_MPS, pos; update_rep=true)

   # apply Hadamard gate(s) to quantum circuit
   apply_single_site_gates!(qc, pos, "H")

   # update representing matrix of quantum circuit
   if update_rep
      update_representation_single_site!(qc, pos, 1)
   end

   # update circuit depth and bond dimension
   qc.CircuitDepth += 1
   push!(qc.BondDim, maxlinkdim(qc.StateVector))
end


""" Function to apply the (single-qubit) Pauli X gate to a quantum circuit qc
in every position specified in the array pos. """
function PauliX!(qc::QC_IT_MPS, pos; update_rep=true)

   # apply Pauli X gate(s) to quantum circuit
   apply_single_site_gates!(qc, pos, "X")

   # update representing matrix of quantum circuit
   if update_rep
      update_representation_single_site!(qc, pos, 2)
   end

   # update circuit depth and bond dimension
   qc.CircuitDepth += 1
   push!(qc.BondDim, maxlinkdim(qc.StateVector))
end


""" Function to apply the (single-qubit) Pauli Y gate to a quantum circuit qc
in every position specified in the array pos. """
function PauliY!(qc::QC_IT_MPS, pos; update_rep=true)

   # apply Pauli Y gate(s) to quantum circuit
   apply_single_site_gates!(qc, pos, "Y")

   # update representing matrix of quantum circuit
   if update_rep
      update_representation_single_site!(qc, pos, 6)
   end

   # update circuit depth and bond dimension
   qc.CircuitDepth += 1
   push!(qc.BondDim, maxlinkdim(qc.StateVector))
end


""" Function to apply the (single-qubit) Pauli Z gate to a quantum circuit qc
in every position specified in the array pos. """
function PauliZ!(qc::QC_IT_MPS, pos; update_rep=true)

   # apply Pauli Z gate(s) to quantum circuit
   apply_single_site_gates!(qc, pos, "Z")

   # update representing matrix of quantum circuit
   if update_rep
      update_representation_single_site!(qc, pos, 5)
   end

   # update circuit depth and bond dimension
   qc.CircuitDepth += 1
   push!(qc.BondDim, maxlinkdim(qc.StateVector))
end


""" Function to apply the (single-qubit) √X gate to a quantum circuit qc
in every position specified in the array pos. """
function SqrtX!(qc::QC_IT_MPS, pos; update_rep=true)

   # apply √X gate(s) to quantum circuit
   apply_single_site_gates!(qc, pos, "√X")

   # update representing matrix of quantum circuit
   if update_rep
      update_representation_single_site!(qc, pos, 18)
   end

   # update circuit depth and bond dimension
   qc.CircuitDepth += 1
   push!(qc.BondDim, maxlinkdim(qc.StateVector))
end


""" Function to apply the (single-qubit) S gate to a quantum circuit qc
in every position specified in the array pos. """
function SGate!(qc::QC_IT_MPS, pos; update_rep=true)

   # apply S gate(s) to quantum circuit
   apply_single_site_gates!(qc, pos, "S")

   # update representing matrix of quantum circuit
   if update_rep
      update_representation_single_site!(qc, pos, 7)
   end

   # update circuit depth and bond dimension
   qc.CircuitDepth += 1
   push!(qc.BondDim, maxlinkdim(qc.StateVector))
end


""" Function to apply the (single-qubit) S† gate to a quantum circuit qc
in every position specified in the array pos. """
function S_dagGate!(qc::QC_IT_MPS, pos; update_rep=true)

   # apply S gate(s) to quantum circuit
   apply_single_site_gates!(qc, pos, "S†")

   # update representing matrix of quantum circuit
   if update_rep
      update_representation_single_site!(qc, pos, 7)
   end

   # update circuit depth and bond dimension
   qc.CircuitDepth += 1
   push!(qc.BondDim, maxlinkdim(qc.StateVector))
end


""" Function to apply the (single-qubit) T gate to a quantum circuit qc
in every position specified in the array pos. """
function TGate!(qc::QC_IT_MPS, pos; update_rep=true)

   # apply T gate(s) to quantum circuit
   apply_single_site_gates!(qc, pos, "T")

   # update representing matrix of quantum circuit
   if update_rep
      update_representation_single_site!(qc, pos, 8)
   end

   # update circuit depth and bond dimension
   qc.CircuitDepth += 1
   push!(qc.BondDim, maxlinkdim(qc.StateVector))
end


""" Function to apply the (single-qubit) Ry gate to a quantum circuit qc
in every position specified in the array pos. Performs a rotation around
the x-axis specified by θ."""
function RxGate!(qc::QC_IT_MPS, pos, θ::Number; update_rep=true)

   # apply Rx gate(s) to quantum circuit
   apply_single_site_gates!(qc, pos, "Rx", θ)

   # update representing matrix of quantum circuit
   if update_rep
      update_representation_single_site!(qc, pos, 9)
   end

   # update circuit depth and bond dimension
   qc.CircuitDepth += 1
   push!(qc.BondDim, maxlinkdim(qc.StateVector))
end


""" Function to apply the (single-qubit) Ry gate to a quantum circuit qc
in every position specified in the array pos. Performs a rotation around
the y-axis specified by θ."""
function RyGate!(qc::QC_IT_MPS, pos, θ::Number; update_rep=true)

   # apply Ry gate(s) to quantum circuit
   apply_single_site_gates!(qc, pos, "Ry", θ)

   # update representing matrix of quantum circuit
   if update_rep
      update_representation_single_site!(qc, pos, 10)
   end

   # update circuit depth and bond dimension
   qc.CircuitDepth += 1
   push!(qc.BondDim, maxlinkdim(qc.StateVector))
end


""" Function to apply the (single-qubit) Rz gate to a quantum circuit qc
in every position specified in the array pos. Performs a rotation around
the z-axis specified by θ."""
function RzGate!(qc::QC_IT_MPS, pos, θ::Number; update_rep=true)

   # apply Rz gate(s) to quantum circuit
   apply_single_site_gates!(qc, pos, "Rz", θ)

   # update representing matrix of quantum circuit
   if update_rep
      update_representation_single_site!(qc, pos, 11)
   end

   # update circuit depth and bond dimension
   qc.CircuitDepth += 1
   push!(qc.BondDim, maxlinkdim(qc.StateVector))
end


""" Function to apply a phase shift by the angle θ to a quantum circuit
qc in every position specified in the array pos. """
function PhaseShift!(qc::QC_IT_MPS, pos, θ::Number; update_rep=true)

   # apply P gate(s) to quantum circuit
   apply_single_site_gates!(qc, pos, "P", θ)

   # update representing matrix of quantum circuit
   if update_rep
      update_representation_single_site!(qc, pos, 12)
   end

   # update circuit depth and bond dimension
   qc.CircuitDepth += 1
   push!(qc.BondDim, maxlinkdim(qc.StateVector))
end


""" Function to apply custom-defined unitary operator U to a quantum circuit
qc in every position specified in the array pos. """
function UGate!(qc::QC_IT_MPS, U, pos; update_rep=true)

   # apply U gate(s) to quantum circuit
   apply_single_site_gates!(qc, pos, "U", U[1, 1], U[1, 2], U[2, 1], U[2, 2])

   # update representing matrix of quantum circuit
   if update_rep
      update_representation_single_site!(qc, pos, 19)
   end

   # update circuit depth and bond dimension
   qc.CircuitDepth += 1
   push!(qc.BondDim, maxlinkdim(qc.StateVector))
end

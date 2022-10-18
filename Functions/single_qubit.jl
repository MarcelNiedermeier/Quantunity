
##############################################################
## Single qubit gates for Quantum Simulator - ED Julia version
##############################################################


#####################
# Auxiliary functions
#####################


""" Auxiliary function to generate an operator O representing a single-site gate
applied on arbitrary positions given in the array pos as a Kronecker product.
E.g. specifying X on pos[1, 3] for a 4-qubit system constructs the operator
X ⊗ E ⊗ X ⊗ E. """
#function get_gate_single_site(qc::QC, O::Matrix{ComplexF64}, pos::Array{Int64, 1})
function get_gate_single_site(qc::QC, O::Matrix, pos::Array{Int64, 1})

    # get identity
    s = Index(2, "QCircuit")
    E = array(op("E", s))

    # initialise first matrix
    if pos[1] == 1
        gate = O
    else
        gate = E
    end

    # construct whole operator
    for i in 2:qc.NumQubits
        if i in pos
            gate = kron(gate, O)
        else
            gate = kron(gate, E)
        end
    end

    return gate
end


""" Function to update the representation of a quantum circuit qc for a
gate (specified by num) applied positions given through the array pos. """
function update_representation_single_site!(qc::QC, pos::Array{Int64, 1}, num::Int64)
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
function hadamard!(qc::QC, pos::Array{Int64, 1}, update_rep=true)

   # get matrices
   s = Index(2, "QCircuit")
   H = array(op("H", s))

   # obtain tensor product operator to act on state vector
   gate = get_gate_single_site(qc, H, pos)

   # update state vector
   qc.StateVector .= gate*qc.StateVector

   # update representing matrix of quantum circuit
   if update_rep
      update_representation_single_site!(qc, pos, 1)
   end

   # update circuit depth
   qc.CircuitDepth += 1
end


""" Function to apply the (single-qubit) Pauli X gate to a quantum circuit qc
in every position specified in the array pos. """
function PauliX!(qc::QC, pos::Array{Int64, 1}, update_rep=true)

   # get matrices
   s = Index(2, "QCircuit")
   X = array(op("X", s))

   # obtain tensor product operator to act on state vector
   gate = get_gate_single_site(qc, X, pos)

   # update state vector
   qc.StateVector .= gate*qc.StateVector

   # update representing matrix of quantum circuit
   if update_rep
      update_representation_single_site!(qc, pos, 2)
   end

   # update circuit depth
   qc.CircuitDepth += 1
end


""" Function to apply the (single-qubit) Pauli Y gate to a quantum circuit qc
in every position specified in the array pos. """
function PauliY!(qc::QC, pos::Array{Int64, 1}, update_rep=true)

   # get matrices
   s = Index(2, "QCircuit")
   Y = array(op("Y", s))

   # obtain tensor product operator to act on state vector
   gate = get_gate_single_site(qc, Y, pos)

   # update state vector
   qc.StateVector .= gate*qc.StateVector

   # update representing matrix of quantum circuit
   if update_rep
      update_representation_single_site!(qc, pos, 6)
   end

   # update circuit depth
   qc.CircuitDepth += 1
end


""" Function to apply the (single-qubit) Pauli Z gate to a quantum circuit qc
in every position specified in the array pos. """
function PauliZ!(qc::QC, pos::Array{Int64, 1}, update_rep=true)

   # get matrices
   s = Index(2, "QCircuit")
   Z = array(op("Z", s))

   # obtain tensor product operator to act on state vector
   gate = get_gate_single_site(qc, Z, pos)

   # update state vector
   qc.StateVector .= gate*qc.StateVector

   # update representing matrix of quantum circuit
   if update_rep
      update_representation_single_site!(qc, pos, 5)
   end

   # update circuit depth
   qc.CircuitDepth += 1
end


""" Function to apply the (single-qubit) √X gate to a quantum circuit qc
in every position specified in the array pos. """
function SqrtX!(qc::QC, pos::Array{Int64, 1}, update_rep=true)

   # get matrices
   s = Index(2, "QCircuit")
   sqrtX = array(op("√X", s))

   # obtain tensor product operator to act on state vector
   gate = get_gate_single_site(qc, sqrtX, pos)

   # update state vector
   qc.StateVector .= gate*qc.StateVector

   # update representing matrix of quantum circuit
   if update_rep
      update_representation_single_site!(qc, pos, 18)
   end

   # update circuit depth
   qc.CircuitDepth += 1
end


""" Function to apply the (single-qubit) S gate to a quantum circuit qc
in every position specified in the array pos. """
function SGate!(qc::QC, pos::Array{Int64, 1}, update_rep=true)

   # get matrices
   s = Index(2, "QCircuit")
   S = array(op("S", s))

   # obtain tensor product operator to act on state vector
   gate = get_gate_single_site(qc, S, pos)

   # update state vector
   qc.StateVector .= gate*qc.StateVector

   # update representing matrix of quantum circuit
   if update_rep
      update_representation_single_site!(qc, pos, 7)
   end

   # update circuit depth
   qc.CircuitDepth += 1
end


""" Function to apply the (single-qubit) T gate to a quantum circuit qc
in every position specified in the array pos. """
function TGate!(qc::QC, pos::Array{Int64, 1}, update_rep=true)

   # get matrices
   s = Index(2, "QCircuit")
   T = array(op("T", s))

   # obtain tensor product operator to act on state vector
   gate = get_gate_single_site(qc, T, pos)

   # update state vector
   qc.StateVector .= gate*qc.StateVector

   # update representing matrix of quantum circuit
   if update_rep
      update_representation_single_site!(qc, pos, 8)
   end

   # update circuit depth
   qc.CircuitDepth += 1
end


""" Function to apply the (single-qubit) Rx gate to a quantum circuit qc
in every position specified in the array pos. Performs a rotation around
the x-axis specified by θ."""
function RXGate!(qc::QC, pos::Array{Int64, 1}, θ::Number, update_rep=true)

   # get matrices
   s = Index(2, "QCircuit")
   Rx = array(op("Rx", s; θ=θ))

   # obtain tensor product operator to act on state vector
   gate = get_gate_single_site(qc, Rx, pos)

   # update state vector
   qc.StateVector .= gate*qc.StateVector

   # update representing matrix of quantum circuit
   if update_rep
      update_representation_single_site!(qc, pos, 9)
   end

   # update circuit depth
   qc.CircuitDepth += 1
end


""" Function to apply the (single-qubit) Ry gate to a quantum circuit qc
in every position specified in the array pos. Performs a rotation around
the x-axis specified by θ."""
function RYGate!(qc::QC, pos::Array{Int64, 1}, θ::Number, update_rep=true)

   # get matrices
   s = Index(2, "QCircuit")
   Ry = array(op("Ry", s; θ=θ))

   # obtain tensor product operator to act on state vector
   gate = get_gate_single_site(qc, Ry, pos)

   # update state vector
   qc.StateVector .= gate*qc.StateVector

   # update representing matrix of quantum circuit
   if update_rep
      update_representation_single_site!(qc, pos, 10)
   end

   # update circuit depth
   qc.CircuitDepth += 1
end


""" Function to apply the (single-qubit) Rz gate to a quantum circuit qc
in every position specified in the array pos. Performs a rotation around
the x-axis specified by θ."""
function RZGate!(qc::QC, pos::Array{Int64, 1}, θ::Number, update_rep=true)

   # get matrices
   s = Index(2, "QCircuit")
   Rz = array(op("Rz", s; θ=θ))

   # obtain tensor product operator to act on state vector
   gate = get_gate_single_site(qc, Rz, pos)

   # update state vector
   qc.StateVector .= gate*qc.StateVector

   # update representing matrix of quantum circuit
   if update_rep
      update_representation_single_site!(qc, pos, 11)
   end

   # update circuit depth
   qc.CircuitDepth += 1
end


""" Function to apply a phase shift by the angle θ to a quantum circuit
qc in every position specified in the array pos. """
function PhaseShift!(qc::QC, pos::Array{Int64, 1}, θ::Number, update_rep=true)

   # get matrices
   s = Index(2, "QCircuit")
   P = array(op("P", s; θ=θ))

   # obtain tensor product operator to act on state vector
   gate = get_gate_single_site(qc, P, pos)

   # update state vector
   qc.StateVector .= gate*qc.StateVector

   # update representing matrix of quantum circuit
   if update_rep
      update_representation_single_site!(qc, pos, 12)
   end

   # update circuit depth
   qc.CircuitDepth += 1
end


""" Function to apply custom-defined unitary operator U to a quantum circuit
qc in every position specified in the array pos. The matrix is defined as
U = [α β; γ δ] with the corresponding parameters. """
function UGate!(qc::QC, pos::Array{Int64, 1}, α::Number, β::Number,
   γ::Number, δ::Number, update_rep=true)

   # get matrices
   s = Index(2, "QCircuit")
   E = array(op("E", s))
   U = array(op("U", s; α=α, β=β, γ=γ, δ=δ))

   # check unitarity
   if norm(U*adjoint(U) - E) > 1e-14
      error("The matrix you have entered is not unitary, check parameters.")
   end

   # obtain tensor product operator to act on state vector
   gate = get_gate_single_site(qc, U, pos)

   # update state vector
   qc.StateVector .= gate*qc.StateVector

   # update representing matrix of quantum circuit
   if update_rep
      update_representation_single_site!(qc, pos, 19)
   end

   # update circuit depth
   qc.CircuitDepth += 1
end


#################
# Random Circuits
#################


""" Function to apply a (different) random unitary single-site gate to a
quantum circuit in the positions given in the array pos. """
function random_single_site_gates!(qc::QC, pos::Array{Int64, 1}; update_rep=true)

    # get identity
    s = Index(2, "QCircuit")
    E = array(op("E", s))

    # initialise first matrix
    if pos[1] == 1

       # get random angles and random unitary 2x2 matrix
       α = 2π * rand()
       β = 2π * rand()
       γ = 2π * rand()
       δ = 2π * rand()
       U = array(op("RandU", s; α=α, β=β, γ=γ, δ=δ))
       gate = U

    else
       gate = E
    end

    # construct whole operator
    for i in 2:qc.NumQubits
       if i in pos

          # get random angles and random unitary 2x2 matrix
          α = 2π * rand()
          β = 2π * rand()
          γ = 2π * rand()
          δ = 2π * rand()
          U = array(op("RandU", s; α=α, β=β, γ=γ, δ=δ))
          gate = kron(gate, U)

       else
          gate = kron(gate, E)
       end
    end

    # update state vector
    qc.StateVector .= gate*qc.StateVector

    # update representing matrix of quantum circuit
    if update_rep
      update_representation_single_site!(qc, pos, 19)
    end

    # update circuit depth
    qc.CircuitDepth += 1

end

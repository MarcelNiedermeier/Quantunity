
####################################################################
## Single qubit gates for Quantum Simulator - Density matrix version
####################################################################


#####################
# Auxiliary functions
#####################

""" Auxiliary function to generate an operator O representing a single-site gate
applied on arbitrary positions given in the array pos as a Kronecker product.
E.g. specifying X on pos[1, 3] for a 4-qubit system constructs the operator
X ⊗ E ⊗ X ⊗ E. """
function get_single_site_gate(qc::QC_DM, matrix::SparseMatrixCSC{ComplexF64, Int64},
   pos::Array{Int64, 1})

    # get identity and number of qubits
    s = Index(2, "QCircuit")
    E = sparse(array(op("E", s)))
    N = qc.NumQubits

    if 1 ∈ pos
        gate = matrix
    else
        gate = E
    end
    for i in 2:N
        if i ∈ pos
            gate = kron(gate, matrix)
        else
            gate = kron(gate, E)
        end
    end
    return gate
end


""" Function to apply a single-site gate to the density matrix representing
the quantum circuit. """
function apply_single_site_gate!(qc::QC_DM, matrix::SparseMatrixCSC{ComplexF64, Int64},
   pos::Array{Int64, 1})
    gate = get_single_site_gate(qc, matrix, pos)
    qc.StateVector = gate * qc.StateVector * gate'
end


""" Function to update the representation of a quantum circuit qc for a
gate (specified by num) applied positions given through the array pos. """
function update_representation_single_site!(qc::QC_DM, pos::Array{Int64, 1},
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


#####################
# Single-qubit errors
#####################


# bit flip, phase flip?

""" Function to apply the one-qubit amplitude damping noise channel in
position pos to a quantum circuit, with damping constant γ. """
function one_qubit_amplitude_damping!(qc::QC_DM, pos)

    # get error rate
    γ = qc.p_amp_damp

    # construct Kraus operators and corresponding gates
    K0 = sparse(Complex.([1. 0.; 0. sqrt(1-γ)]))
    K1 = sparse(Complex.([0. sqrt(γ); 0. 0.]))
    K0_gate = get_single_site_gate(qc, K0, pos)
    K1_gate = get_single_site_gate(qc, K1, pos)

    qc.StateVector =  K0_gate*qc.StateVector*K0_gate' + K1_gate*qc.StateVector*K1_gate'
end


""" Function to apply the one-qubit depolarisation noise channel in
position pos to a quantum circuit, with probability p. """
function one_qubit_depolarisation!(qc::QC_DM, pos)

    # get identity and number of qubits
    s = Index(2, "QCircuit")
    E = sparse(array(op("E", s)))
    X = sparse(array(op("X", s)))
    Y = sparse(array(op("Y", s)))
    Z = sparse(array(op("Z", s)))

    # get error rate
    p = qc.p_depol1

    # case where no noise happens
    qc_no_noise = (1-p)*qc.StateVector

    # 1-qubit Pauli group without E
    Pauli_ops = collect(Iterators.product([E, X, Y, Z]))[2:end]

    # get gates for depolarisation error
    error_gates = []
    for i in 1:length(Pauli_ops)
        #push!(error_gates, get_single_site_gate(qc, Pauli_ops[i]..., pos[1]))
        push!(error_gates, get_single_site_gate(qc, Pauli_ops[i]..., pos))
    end

    # apply depolarisation error
    qc_noise = p/3 * error_gates[1]*qc.StateVector*error_gates[1]
    for i in 2:length(error_gates)
        qc_noise += p/3 * error_gates[i]*qc.StateVector*error_gates[i]
    end

    qc.StateVector = qc_no_noise + qc_noise
end


#######
# Gates
#######


""" Function to apply the (single-qubit) Hadamard gate to a quantum circuit qc
in every position specified in the array pos. """
function Hadamard!(qc::QC_DM, pos::Array{Int64, 1}; update_rep=true)

    # get matrices
    s = Index(2, "QCircuit")
    H = sparse(array(op("H", s)))

    # apply gate
    apply_single_site_gate!(qc, H, pos)

    # do error simulation
    one_qubit_amplitude_damping!(qc, pos)
    one_qubit_depolarisation!(qc, pos)

    # update representing matrix of quantum circuit
    if update_rep
       update_representation_single_site!(qc, pos, 1)
    end

    # update circuit depth
    qc.CircuitDepth += 1
end


""" Function to apply the (single-qubit) Pauli X gate to a quantum circuit qc
in every position specified in the array pos. """
function PauliX!(qc::QC_DM, pos::Array{Int64, 1}; update_rep=true)

   # get matrices
   s = Index(2, "QCircuit")
   X = sparse(array(op("X", s)))

   # apply gate
   apply_single_site_gate!(qc, X, pos)

   # do error simulation
   one_qubit_amplitude_damping!(qc, pos)
   one_qubit_depolarisation!(qc, pos)

   # update representing matrix of quantum circuit
   if update_rep
      update_representation_single_site!(qc, pos, 2)
   end

   # update circuit depth
   qc.CircuitDepth += 1
end


""" Function to apply the (single-qubit) Pauli Y gate to a quantum circuit qc
in every position specified in the array pos. """
function PauliY!(qc::QC_DM, pos::Array{Int64, 1}; update_rep=true)

   # get matrices
   s = Index(2, "QCircuit")
   Y = sparse(array(op("Y", s)))

   # apply gate
   apply_single_site_gate!(qc, Y, pos)

   # do error simulation
   one_qubit_amplitude_damping!(qc, pos)
   one_qubit_depolarisation!(qc, pos)

   # update representing matrix of quantum circuit
   if update_rep
      update_representation_single_site!(qc, pos, 6)
   end

   # update circuit depth
   qc.CircuitDepth += 1
end


""" Function to apply the (single-qubit) Pauli Z gate to a quantum circuit qc
in every position specified in the array pos. """
function PauliZ!(qc::QC_DM, pos::Array{Int64, 1}; update_rep=true)

   # get matrices
   s = Index(2, "QCircuit")
   Z = sparse(array(op("Z", s)))

   # apply gate
   apply_single_site_gate!(qc, Z, pos)

   # do error simulation
   one_qubit_amplitude_damping!(qc, pos)
   one_qubit_depolarisation!(qc, pos)

   # update representing matrix of quantum circuit
   if update_rep
      update_representation_single_site!(qc, pos, 5)
   end

   # update circuit depth
   qc.CircuitDepth += 1
end


""" Function to apply the (single-qubit) √X gate to a quantum circuit qc
in every position specified in the array pos. """
function SqrtX!(qc::QC_DM, pos::Array{Int64, 1}; update_rep=true)

   # get matrices
   s = Index(2, "QCircuit")
   sqrt_X = sparse(array(op("√X", s)))

   # apply gate
   apply_single_site_gate!(qc, sqrt_X, pos)

   # do error simulation
   one_qubit_amplitude_damping!(qc, pos)
   one_qubit_depolarisation!(qc, pos)

   # update representing matrix of quantum circuit
   if update_rep
      update_representation_single_site!(qc, pos, 18)
   end

   # update circuit depth
   qc.CircuitDepth += 1
end


""" Function to apply the (single-qubit) S gate to a quantum circuit qc
in every position specified in the array pos. """
function SGate!(qc::QC_DM, pos::Array{Int64, 1}; update_rep=true)

   # get matrices
   s = Index(2, "QCircuit")
   S = sparse(array(op("S", s)))

   # apply gate
   apply_single_site_gate!(qc, S, pos)

   # do error simulation
   one_qubit_amplitude_damping!(qc, pos)
   one_qubit_depolarisation!(qc, pos)

   # update representing matrix of quantum circuit
   if update_rep
      update_representation_single_site!(qc, pos, 7)
   end

   # update circuit depth
   qc.CircuitDepth += 1
end


""" Function to apply the (single-qubit) T gate to a quantum circuit qc
in every position specified in the array pos. """
function TGate!(qc::QC_DM, pos::Array{Int64, 1}; update_rep=true)

   # get matrices
   s = Index(2, "QCircuit")
   T = sparse(array(op("T", s)))

   # apply gate
   apply_single_site_gate!(qc, T, pos)

   # do error simulation
   one_qubit_amplitude_damping!(qc, pos)
   one_qubit_depolarisation!(qc, pos)

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
function RXGate!(qc::QC_DM, pos::Array{Int64, 1}, θ::Number; update_rep=true)

   # get matrices
   s = Index(2, "QCircuit")
   Rx = sparse(array(op("Rx", s; θ=θ)))

   # apply gate
   apply_single_site_gate!(qc, Rx, pos)

   # do error simulation
   one_qubit_amplitude_damping!(qc, pos)
   one_qubit_depolarisation!(qc, pos)

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
function RYGate!(qc::QC_DM, pos::Array{Int64, 1}, θ::Number; update_rep=true)

   # get matrices
   s = Index(2, "QCircuit")
   Ry = sparse(array(op("Ry", s; θ=θ)))

   # apply gate
   apply_single_site_gate!(qc, Ry, pos)

   # do error simulation
   one_qubit_amplitude_damping!(qc, pos)
   one_qubit_depolarisation!(qc, pos)

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
function RZGate!(qc::QC_DM, pos::Array{Int64, 1}, θ::Number; update_rep=true)

   # get matrices
   s = Index(2, "QCircuit")
   Rz = sparse(array(op("Rz", s; θ=θ)))

   # apply gate
   apply_single_site_gate!(qc, Rz, pos)

   # do error simulation
   one_qubit_amplitude_damping!(qc, pos)
   one_qubit_depolarisation!(qc, pos)

   # update representing matrix of quantum circuit
   if update_rep
      update_representation_single_site!(qc, pos, 11)
   end

   # update circuit depth
   qc.CircuitDepth += 1
end


""" Function to apply a phase shift by the angle θ to a quantum circuit
qc in every position specified in the array pos. """
function PhaseShift!(qc::QC_DM, pos::Array{Int64, 1}, θ::Number; update_rep=true)

   # get matrices
   s = Index(2, "QCircuit")
   P = sparse(array(op("P", s; θ=θ)))

   # apply gate
   apply_single_site_gate!(qc, P, pos)

   # do error simulation
   one_qubit_amplitude_damping!(qc, pos)
   one_qubit_depolarisation!(qc, pos)

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
function UGate!(qc::QC_DM, U, pos::Array{Int64, 1}; update_rep=true)

   α = U[1, 1]
   β = U[1, 2]
   γ = U[2, 1]
   δ = U[2, 2]

   # get matrices
   s = Index(2, "QCircuit")
   E = sparse(array(op("E", s)))
   U = sparse(array(op("U", s; α=α, β=β, γ=γ, δ=δ)))

   # check unitarity
   if norm(U*adjoint(U) - E) > 1e-14
      error("The matrix you have entered is not unitary, check parameters.")
   end

   # apply gate
   apply_single_site_gate!(qc, U, pos)

   # do error simulation
   one_qubit_amplitude_damping!(qc, pos)
   one_qubit_depolarisation!(qc, pos)

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


""" Function to apply a (the same) random unitary single-site gate to a
quantum circuit in the positions given in the array pos. """
function RandomU!(qc::QC_DM, pos::Array{Int64, 1}; update_rep=true)

    # get number of qubits and random unitary 2x2 matrix
    N = qc.NumQubits
    s = Index(2, "QCircuit")
    α_rand = 2π * rand()
    β_rand = 2π * rand()
    γ_rand = 2π * rand()
    δ_rand = 2π * rand()
    U = sparse(array(op("RandU", s; α=α_rand, β=β_rand, γ=γ_rand, δ=δ_rand)))

    # apply gate
    apply_single_site_gate!(qc, U, pos)

    # do error simulation
    one_qubit_amplitude_damping!(qc, pos)
    one_qubit_depolarisation!(qc, pos)

    # update representing matrix of quantum circuit
    if update_rep
      update_representation_single_site!(qc, pos, 19)
    end

    # update circuit depth
    qc.CircuitDepth += 1
end


""" Function to apply a (different) random unitary single-site gate to a
quantum circuit in the positions given in the array pos. """
function random_single_site_gates!(qc::QC, pos::Array{Int64, 1}; update_rep=true)

   # initialise
   N = qc.NumQubits
   s = Index(2, "QCircuit")

   for i in [i for i in pos:pos+num-1]

        # get different random unitary for each position
        α_rand = 2π * rand()
        β_rand = 2π * rand()
        γ_rand = 2π * rand()
        δ_rand = 2π * rand()
        U = sparse(array(op("RandU", s; α=α, β=β, γ=γ, δ=δ)))

        # apply gate, do error simulation
        apply_single_site_gate!(qc, U, pos)
        one_qubit_amplitude_damping!(qc, pos)
        one_qubit_depolarisation!(qc, pos)
   end

   # update representing matrix of quantum circuit
   if update_rep
     update_representation_single_site!(qc, pos, 19)
   end

   # update circuit depth
   qc.CircuitDepth += 1
end

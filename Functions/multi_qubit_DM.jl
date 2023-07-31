
#################################
## Multi-Qubit Gates - DM Version
#################################


#####################
# Auxiliary Functions
#####################

# all in arbitrary_CU_DM.jl !


##################
# Two-qubit errors
##################


""" Function to apply the two-qubit depolarisation noise channel in
positions pos[1], pos[2] to a quantum circuit, with probability p. """
function two_qubit_depolarisation!(qc::QC_DM, pos)

    # get identity and number of qubits
    s = Index(2, "QCircuit")
    E = sparse(array(op("E", s)))
    X = sparse(array(op("X", s)))
    Y = sparse(array(op("Y", s)))
    Z = sparse(array(op("Z", s)))

    # get error rate
    p = qc.p_depol2

    # case where no noise happens
    qc_no_noise = (1-p)*qc.StateVector

    # 2-qubit Pauli group without E ⊗ E
    Pauli_ops = collect(Iterators.product([E, X, Y, Z], [E, X, Y, Z]))[2:end]

    # get gates for depolarisation error
    error_gates = []
    for i in 1:length(Pauli_ops)
        ops = Pauli_ops[i]
        error_gate1 = get_single_site_gate(qc, ops[1], [pos[1]])
        error_gate2 = get_single_site_gate(qc, ops[2], [pos[2]])
        push!(error_gates, error_gate1*error_gate2)
    end

    # apply depolarisation error
    qc_noise = p/15 * error_gates[1]*qc.StateVector*error_gates[1]
    for i in 2:length(error_gates)
        qc_noise += p/15 * error_gates[i]*qc.StateVector*error_gates[i]
    end

    qc.StateVector = qc_no_noise + qc_noise
end


###########################################
# Shorthands for Standard Multi-Qubit Gates
###########################################


""" Function to apply the CNOT gate to a quantum circuit qc
in position specified in the array pos; pos[1] = control qubit,
pos[2] = target qubit. """
function Cnot!(qc::QC_DM, pos; update_rep=true)

   # check correct specification of inputs
   if length(pos) ≠ 2
      error("Incorrect specification of position (need array of length 2).")
   end

   # get matrices
   s = Index(2, "QCircuit")
   X = sparse(array(op("X", s)))

   # apply gate
   multiply_controlled_single_site!(qc, X, [pos[1]], [pos[2]],
      update_rep=update_rep, num=4)

   # do error simulation
   two_qubit_depolarisation!(qc, pos)

end


""" Function to apply the TOFFOLI gate to a quantum circuit qc
in position specified in the array pos; pos[1] = control qubit 1,
pos[2] = control qubit 2, pos[3] = target qubit. """
function Toffoli!(qc::QC, pos; update_rep=true)

   # check correct specification of inputs
   if length(pos) ≠ 3
      error("Incorrect specification of position (need array of length 3).")
   end

   # get matrices
   s = Index(2, "QCircuit")
   X = sparse(array(op("X", s)))

   # apply Pauli X gate to quantum circuit
   multiply_controlled_single_site!(qc, X, [pos[1], pos[2]], [pos[3]],
      update_rep=update_rep, num=4)
end


""" Function to apply the Rn gate to a quantum circuit qc
in position specified in the array pos; pos[1] = control qubit,
pos[2] = target qubit. n specifies the rotation angle 2π/2^n. """
function CRn!(qc::QC_DM, n, pos; update_rep=true)

   # check correct specification of inputs
   if length(pos) ≠ 2
      error("Incorrect specification of position (need array of length 2).")
   end

   # get matrices
   s = Index(2, "QCircuit")
   Rn = sparse(array(op("Rn", s; n=n)))

   # apply X gate to quantum circuit
   multiply_controlled_single_site!(qc, Rn, [pos[1]], [pos[2]],
      update_rep=update_rep, num=24)

   # do error simulation
   two_qubit_depolarisation!(qc, pos)
end


""" Function to apply the SWAP gate to a quantum circuit qc;
switches the two qubits given in the array pos. """
function Swap!(qc::QC_DM, pos; update_rep=true)

   # check correct specification of inputs
   if length(pos) ≠ 2
      error("Incorrect specification of position (need array of length 2).")
   end

   Cnot!(qc, pos; update_rep=false)
   Cnot!(qc, reverse(pos); update_rep=false)
   Cnot!(qc, pos; update_rep=false)

   # application via CNOTs seems to have similar consequences (doesn't preserve norm)
   #Cnot!(qc, pos, update_rep=update_rep, recordEE=recordEE, α=α, cutoff=cutoff)
   #Cnot!(qc, reverse(pos), update_rep=update_rep, recordEE=recordEE, α=α, cutoff=cutoff)
   #Cnot!(qc, pos, update_rep=update_rep, recordEE=recordEE, α=α, cutoff=cutoff)

   # apply standard, uncontrolled SWAP gate
   #multiply_controlled_general_SWAP!(qc, [], pos, update_rep=update_rep, num=4,
   #   recordEE=recordEE, α=α, cutoff=cutoff)

   # explicit renormalisation
   #qc.StateVector = 1/norm(qc.StateVector) * qc.StateVector

   # update representing matrix of quantum circuit
   if update_rep
      update_representation_arbitrary_CU!(qc, control_qubits, action_qubits, num=4)
   end

   # do error simulation
   two_qubit_depolarisation!(qc, pos)
end


###########################
# Currently not implemented
###########################


""" Function to apply the FREDKIN gate to a quantum circuit qc;
switches the two qubits pos[2] and pos[3] if the control qubit
in pos[1] is set.
CURRENTLY NOT IMPLEMENTED! """
function Fredkin!(qc::QC_DM, pos; update_rep=true)

   # check correct specification of inputs
   if length(pos) ≠ 3
      error("Incorrect specification of position (need array of length 3).")
   end

   # apply standard, controlled SWAP gate
   #multiply_controlled_general_SWAP!(qc, [pos[1]], [pos[2], pos[3]],
   #   update_rep=update_rep, num=4, recordEE=recordEE, α=α, cutoff=cutoff)

end


""" Function to apply the √SWAP gate to a quantum circuit qc;
switches the two qubits given in the array action_qubits; admits
arbitrary number of control qubits.
CURRENTLY NOT IMPLEMENTED!"""
function SqrtSwap!(qc::QC_DM, control_qubits, action_qubits; update_rep=true, β=1/2)

   # check correct specification of inputs
   if length(action_qubits) ≠ 2
      error("Incorrect specification of position (need array of length 2).")
   end

   # apply controlled √SWAP gate
   #multiply_controlled_general_SWAP!(qc, control_qubits, action_qubits,
   #   update_rep=update_rep, num=26, recordEE=recordEE, α=α, cutoff=cutoff, β=β)
end


""" Function to apply the SWAP^β gate to a quantum circuit qc;
switches the two qubits given in the array action_qubits; admits
arbitrary number of control qubits.
CURRENTLY NOT IMPLEMENTED!"""
function GeneralSwap!(qc::QC_DM, control_qubits, action_qubits; update_rep=true, β=1/2)

   if β == 1/2
      println("Check if √SWAP gate could be more convenient.")
   end

   # check correct specification of inputs
   if length(action_qubits) ≠ 2
      error("Incorrect specification of position (need array of length 2).")
   end

   # apply controlled SWAP gate to arbitrary power β
   #multiply_controlled_general_SWAP!(qc, control_qubits, action_qubits,
   #   update_rep=update_rep, num=27, recordEE=recordEE, α=α, cutoff=cutoff, β=β)
end



##########################################
# Arbitrarily controlled single-site gates
##########################################


""" Function to apply the (single-qubit) Hadamard gate to a quantum circuit qc
in every position specified in the array pos; an arbitrary number of control
qubits can be given. """
function C_Hadamard!(qc::QC_DM, control_qubits, pos; update_rep=true)

   # get matrices
   s = Index(2, "QCircuit")
   H = sparse(array(op("H", s)))

   # apply Hadamard gate to quantum circuit
   multiply_controlled_single_site!(qc, H, control_qubits, pos,
      update_rep=update_rep, num=1)

   # do error simulation (for 2 qubits)
   if length(control_qubits) + length(pos) == 2
      two_qubit_depolarisation!(qc, pos)
   end
end


""" Function to apply the (single-qubit) Pauli X gate gate to a quantum circuit qc
in every position specified in the array pos; an arbitrary number of control
qubits can be given. """
function C_PauliX!(qc::QC_DM, control_qubits, pos; update_rep=true)

   # get matrices
   s = Index(2, "QCircuit")
   X = sparse(array(op("X", s)))

   # apply Pauli X gate to quantum circuit
   multiply_controlled_single_site!(qc, "X", control_qubits, pos,
      update_rep=update_rep, num=2)

   # do error simulation (for 2 qubits)
   if length(control_qubits) + length(pos) == 2
      two_qubit_depolarisation!(qc, pos)
   end
end


""" Function to apply the (single-qubit) Pauli Y gate gate to a quantum circuit qc
in every position specified in the array pos; an arbitrary number of control
qubits can be given. """
function C_PauliY!(qc::QC_DM, control_qubits, pos; update_rep=true)

   # get matrices
   s = Index(2, "QCircuit")
   Y = sparse(array(op("Y", s)))

   # apply Pauli Y gate to quantum circuit
   multiply_controlled_single_site!(qc, Y, control_qubits, pos,
      update_rep=update_rep, num=6)

   # do error simulation (for 2 qubits)
   if length(control_qubits) + length(pos) == 2
      two_qubit_depolarisation!(qc, pos)
   end
end


""" Function to apply the (single-qubit) Pauli Z gate gate to a quantum circuit qc
in every position specified in the array pos; an arbitrary number of control
qubits can be given. """
function C_PauliZ!(qc::QC_DM, control_qubits, pos; update_rep=true)

   # get matrices
   s = Index(2, "QCircuit")
   Z = sparse(array(op("Z", s)))

   # apply Pauli Z gate to quantum circuit
   multiply_controlled_single_site!(qc, Z, control_qubits, pos,
      update_rep=update_rep, num=5)

   # do error simulation (for 2 qubits)
   if length(control_qubits) + length(pos) == 2
      two_qubit_depolarisation!(qc, pos)
   end
end


""" Function to apply the (single-qubit) √X gate gate to a quantum circuit qc
in every position specified in the array pos; an arbitrary number of control
qubits can be given. """
function C_SqrtX!(qc::QC_DM, control_qubits, pos; update_rep=true)

   # get matrices
   s = Index(2, "QCircuit")
   sqrt_X = sparse(array(op("√X", s)))

   # apply √X gate to quantum circuit
   multiply_controlled_single_site!(qc, sqrt_X, control_qubits, pos,
      update_rep=update_rep, num=18)

   # do error simulation (for 2 qubits)
   if length(control_qubits) + length(pos) == 2
      two_qubit_depolarisation!(qc, pos)
   end
end


""" Function to apply the (single-qubit) S gate gate to a quantum circuit qc
in every position specified in the array pos; an arbitrary number of control
qubits can be given. """
function C_SGate!(qc::QC_DM, control_qubits, pos; update_rep=true)

   # get matrices
   s = Index(2, "QCircuit")
   S = sparse(array(op("S", s)))

   # apply S gate to quantum circuit
   multiply_controlled_single_site!(qc, S, control_qubits, pos,
      update_rep=update_rep, num=7)

   # do error simulation (for 2 qubits)
   if length(control_qubits) + length(pos) == 2
      two_qubit_depolarisation!(qc, pos)
   end
end


""" Function to apply the (single-qubit) T gate gate to a quantum circuit qc
in every position specified in the array pos; an arbitrary number of control
qubits can be given. """
function C_TGate!(qc::QC_DM, control_qubits, pos; update_rep=true)

   # get matrices
   s = Index(2, "QCircuit")
   T = sparse(array(op("T", s)))

   # apply T gate to quantum circuit
   multiply_controlled_single_site!(qc, T, control_qubits, pos,
      update_rep=update_rep, num=8)

   # do error simulation (for 2 qubits)
   if length(control_qubits) + length(pos) == 2
      two_qubit_depolarisation!(qc, pos)
   end
end


##################
# Parametric Gates
##################


""" Function to apply the (single-qubit) Rx gate gate to a quantum circuit qc
in a position specified in the array pos; an arbitrary number of control
qubits can be given.
Parameters:
- θ: rotation angle.
- update_rep: save graphical representation of gate. """
function C_RxGate!(qc::QC_DM, control_qubits, pos, θ::Number; update_rep=true)

   # get matrices
   s = Index(2, "QCircuit")
   Rx = sparse(array(op("Rx", s; θ=θ)))

   # apply Rx gate to quantum circuit
   multiply_controlled_single_site!(qc, Rx, control_qubits, pos,
      update_rep=update_rep, num=9)

   # do error simulation (for 2 qubits)
   if length(control_qubits) + length(pos) == 2
      two_qubit_depolarisation!(qc, pos)
   end
end


""" Function to apply the (single-qubit) Ry gate gate to a quantum circuit qc
in a position specified in the array pos; an arbitrary number of control
qubits can be given.
Parameters:
- θ: rotation angle.
- update_rep: save graphical representation of gate. """
function C_RyGate!(qc::QC_DM, control_qubits, pos, θ::Number; update_rep=true)

   # get matrices
   s = Index(2, "QCircuit")
   Ry = sparse(array(op("Ry", s; θ=θ)))

   # apply Ry gate to quantum circuit
   multiply_controlled_single_site!(qc, Ry, control_qubits, pos,
      update_rep=update_rep, num=10)

   # do error simulation (for 2 qubits)
   if length(control_qubits) + length(pos) == 2
      two_qubit_depolarisation!(qc, pos)
   end
end


""" Function to apply the (single-qubit) Rz gate gate to a quantum circuit qc
in a position specified in the array pos; an arbitrary number of control
qubits can be given.
Parameters:
- θ: rotation angle.
- update_rep: save graphical representation of gate. """
function C_RzGate!(qc::QC_DM, control_qubits, pos, θ::Number; update_rep=true)

   # get matrices
   s = Index(2, "QCircuit")
   Rz = sparse(array(op("Rz", s; θ=θ)))

   # apply Rz gate to quantum circuit
   multiply_controlled_single_site!(qc, Rz, control_qubits, pos,
      update_rep=update_rep, num=11)

   # do error simulation (for 2 qubits)
   if length(control_qubits) + length(pos) == 2
      two_qubit_depolarisation!(qc, pos)
   end
end


""" Function to apply the (single-qubit) phase shift gate to a quantum circuit qc
in a position specified in the array pos; an arbitrary number of control
qubits can be given.
Parameters:
- θ: phase angle.
- update_rep: save graphical representation of gate. """
function C_PhaseShift!(qc::QC_DM, control_qubits, pos, θ::Number; update_rep=true)

   # get matrices
   s = Index(2, "QCircuit")
   P = sparse(array(op("P", s; θ=θ)))

   # apply P gate to quantum circuit
   multiply_controlled_single_site!(qc, P, control_qubits, pos,
      update_rep=update_rep, num=12)

   # do error simulation (for 2 qubits)
   if length(control_qubits) + length(pos) == 2
      two_qubit_depolarisation!(qc, pos)
   end
end


""" Function to apply the (single-qubit) Rn gate to a quantum circuit qc
in a position specified in the array pos; an arbitrary number of control
qubits can be given.
Parameters:
- θ: phase angle.
- update_rep: save graphical representation of gate. """
function C_Rn!(qc::QC_DM, control_qubits, pos, n::Number; update_rep=true)

   # get matrices
   s = Index(2, "QCircuit")
   Rn = sparse(array(op("Rn", s; n=n)))

   # apply Rn gate to quantum circuit
   multiply_controlled_single_site!(qc, Rn, control_qubits, pos,
      update_rep=update_rep, num=24)

   # do error simulation (for 2 qubits)
   if length(control_qubits) + length(pos) == 2
      two_qubit_depolarisation!(qc, pos)
   end
end


""" Function to apply the (single-qubit) custom-defined unitary operator U to
a quantum circuit qc in a position specified in the array pos. An arbitrary
number of control qubits can be given.
Parameters:
- update_rep: save graphical representation of gate. """
function C_UGate!(qc::QC_DM, U, control_qubits, pos; update_rep=true)

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

   # apply U gate to quantum circuit
   multiply_controlled_single_site!(qc, U, control_qubits, pos,
      update_rep=update_rep, num=19, recordEE=recordEE)

   # do error simulation (for 2 qubits)
   if length(control_qubits) + length(pos) == 2
      two_qubit_depolarisation!(qc, pos)
   end
end

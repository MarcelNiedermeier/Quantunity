
##########################################
## Multi-Qubit Gates - MPS ITensor Version
##########################################


#####################
# Auxiliary Functions
#####################

# all in arbitrary_CU_MPS.jl !


###########################################
# Shorthands for Standard Multi-Qubit Gates
###########################################


""" Function to apply the CNOT gate to a quantum circuit qc
in position specified in the array pos; pos[1] = control qubit,
pos[2] = target qubit.
Parameters:
- update_rep: save graphical representation of gate.
- recordEE: calculate entanglement entropy after gate.
- α: parameter in Rényi entropy (α=1 for von Neumann).
- cutoff: smallest singular to be kept for entropy calulation. """
function Cnot!(qc::QC_IT_MPS, pos; update_rep=true,
   recordEE=true, α=1, cutoff=1E-6)

   # check correct specification of inputs
   if length(pos) ≠ 2
      error("Incorrect specification of position (need array of length 2).")
   end

   # apply Pauli X gate to quantum circuit
   multiply_controlled_single_site!(qc, "X", [pos[1]], [pos[2]],
      update_rep=update_rep, num=4, recordEE=recordEE, α=α, cutoff=cutoff)
end


""" Function to apply the TOFFOLI gate to a quantum circuit qc
in position specified in the array pos; pos[1] = control qubit 1,
pos[2] = control qubit 2, pos[3] = target qubit.
Parameters:
- update_rep: save graphical representation of gate.
- recordEE: calculate entanglement entropy after gate.
- α: parameter in Rényi entropy (α=1 for von Neumann).
- cutoff: smallest singular to be kept for entropy calulation. """
function Toffoli!(qc::QC_IT_MPS, pos; update_rep=true,
   recordEE=true, α=1, cutoff=1E-6)

   # check correct specification of inputs
   if length(pos) ≠ 3
      error("Incorrect specification of position (need array of length 3).")
   end

   # apply Pauli X gate to quantum circuit
   multiply_controlled_single_site!(qc, "X", [pos[1], pos[2]], [pos[3]],
      update_rep=update_rep, num=4, recordEE=recordEE, α=α, cutoff=cutoff)
end


""" Function to apply the Rn gate to a quantum circuit qc
in position specified in the array pos; pos[1] = control qubit,
pos[2] = target qubit. n specifies the rotation angle 2π/2^n.
Parameters:
- update_rep: save graphical representation of gate.
- recordEE: calculate entanglement entropy after gate.
- α: parameter in Rényi entropy (α=1 for von Neumann).
- cutoff: smallest singular to be kept for entropy calulation. """
function CRn!(qc::QC_IT_MPS, n, pos; update_rep=true, recordEE=true,
   α=1, cutoff=1E-6)

   # check correct specification of inputs
   if length(pos) ≠ 2
      error("Incorrect specification of position (need array of length 2).")
   end

   # apply Rn gate to quantum circuit
   C_Rn!(qc, [pos[1]], [pos[2]], n, update_rep=update_rep,
      recordEE=recordEE, α=α, cutoff=cutoff)

end


""" Function to apply the SWAP gate to a quantum circuit qc;
switches the two qubits given in the array pos.
Parameters:
- update_rep: save graphical representation of gate.
- recordEE: calculate entanglement entropy after gate.
- α: parameter in Rényi entropy (α=1 for von Neumann).
- cutoff: smallest singular to be kept for entropy calulation. """
function Swap!(qc::QC_IT_MPS, pos; update_rep=true,
   recordEE=true, α=1, cutoff=1E-6)

   # check correct specification of inputs
   if length(pos) ≠ 2
      error("Incorrect specification of position (need array of length 2).")
   end

   # application via CNOTs seems to have similar consequences (doesn't preserve norm)
   #Cnot!(qc, pos, update_rep=update_rep, recordEE=recordEE, α=α, cutoff=cutoff)
   #Cnot!(qc, reverse(pos), update_rep=update_rep, recordEE=recordEE, α=α, cutoff=cutoff)
   #Cnot!(qc, pos, update_rep=update_rep, recordEE=recordEE, α=α, cutoff=cutoff)

   # apply standard, uncontrolled SWAP gate
   multiply_controlled_general_SWAP!(qc, [], pos, update_rep=update_rep, num=4,
      recordEE=recordEE, α=α, cutoff=cutoff)

   # explicit renormalisation
   qc.StateVector = 1/norm(qc.StateVector) * qc.StateVector

end


""" Function to apply the FREDKIN gate to a quantum circuit qc;
switches the two qubits pos[2] and pos[3] if the control qubit
in pos[1] is set.
Parameters:
- update_rep: save graphical representation of gate.
- recordEE: calculate entanglement entropy after gate.
- α: parameter in Rényi entropy (α=1 for von Neumann).
- cutoff: smallest singular to be kept for entropy calulation. """
function Fredkin!(qc::QC_IT_MPS, pos; update_rep=true,
   recordEE=true, α=1, cutoff=1E-6)

   # check correct specification of inputs
   if length(pos) ≠ 3
      error("Incorrect specification of position (need array of length 3).")
   end

   # apply standard, controlled SWAP gate
   multiply_controlled_general_SWAP!(qc, [pos[1]], [pos[2], pos[3]],
      update_rep=update_rep, num=4, recordEE=recordEE, α=α, cutoff=cutoff)
end


""" Function to apply the √SWAP gate to a quantum circuit qc;
switches the two qubits given in the array action_qubits; admits
arbitrary number of control qubits.
Parameters:
- update_rep: save graphical representation of gate.
- recordEE: calculate entanglement entropy after gate.
- α: parameter in Rényi entropy (α=1 for von Neumann).
- cutoff: smallest singular to be kept for entropy calulation. """
function SqrtSwap!(qc::QC_IT_MPS, control_qubits, action_qubits; update_rep=true,
   recordEE=true, α=1, cutoff=1E-6, β=1/2)

   # check correct specification of inputs
   if length(action_qubits) ≠ 2
      error("Incorrect specification of position (need array of length 2).")
   end

   # apply controlled √SWAP gate
   multiply_controlled_general_SWAP!(qc, control_qubits, action_qubits,
      update_rep=update_rep, num=26, recordEE=recordEE, α=α, cutoff=cutoff, β=β)
end


""" Function to apply the SWAP^β gate to a quantum circuit qc;
switches the two qubits given in the array action_qubits; admits
arbitrary number of control qubits.
Parameters:
- update_rep: save graphical representation of gate.
- recordEE: calculate entanglement entropy after gate.
- α: parameter in Rényi entropy (α=1 for von Neumann).
- cutoff: smallest singular to be kept for entropy calulation. """
function GeneralSwap!(qc::QC_IT_MPS, control_qubits, action_qubits; update_rep=true,
   recordEE=true, α=1, cutoff=1E-6, β=1/2)

   if β == 1/2
      println("Check if √SWAP gate could be more convenient.")
   end

   # check correct specification of inputs
   if length(action_qubits) ≠ 2
      error("Incorrect specification of position (need array of length 2).")
   end

   # apply controlled SWAP gate to arbitrary power β
   multiply_controlled_general_SWAP!(qc, control_qubits, action_qubits,
      update_rep=update_rep, num=27, recordEE=recordEE, α=α, cutoff=cutoff, β=β)
end


#############
# Fixed Gates
#############


""" Function to apply the (single-qubit) Hadamard gate to a quantum circuit qc
in every position specified in the array pos; an arbitrary number of control
qubits can be given.
Parameters:
- update_rep: save graphical representation of gate.
- recordEE: calculate entanglement entropy after gate.
- α: parameter in Rényi entropy (α=1 for von Neumann).
- cutoff: smallest singular to be kept for entropy calulation. """
function C_Hadamard!(qc::QC_IT_MPS, control_qubits, pos;
   update_rep=true, recordEE=true, α=1, cutoff=1E-6)

   # apply Hadamard gate to quantum circuit
   multiply_controlled_single_site!(qc, "H", control_qubits, pos,
      update_rep=update_rep, num=1, recordEE=recordEE, α=α, cutoff=cutoff)
end


""" Function to apply the (single-qubit) Pauli X gate gate to a quantum circuit qc
in every position specified in the array pos; an arbitrary number of control
qubits can be given.
Parameters:
- update_rep: save graphical representation of gate.
- recordEE: calculate entanglement entropy after gate.
- α: parameter in Rényi entropy (α=1 for von Neumann).
- cutoff: smallest singular to be kept for entropy calulation. """
function C_PauliX!(qc::QC_IT_MPS, control_qubits, pos;
   update_rep=true, recordEE=true, α=1, cutoff=1E-6)

   # apply Pauli X gate to quantum circuit
   multiply_controlled_single_site!(qc, "X", control_qubits, pos,
      update_rep=update_rep, num=2, recordEE=recordEE, α=α, cutoff=cutoff)
end


""" Function to apply the (single-qubit) Pauli Y gate gate to a quantum circuit qc
in every position specified in the array pos; an arbitrary number of control
qubits can be given.
Parameters:
- update_rep: save graphical representation of gate.
- recordEE: calculate entanglement entropy after gate.
- α: parameter in Rényi entropy (α=1 for von Neumann).
- cutoff: smallest singular to be kept for entropy calulation. """
function C_PauliY!(qc::QC_IT_MPS, control_qubits, pos;
   update_rep=true, recordEE=true, α=1, cutoff=1E-6)

   # apply Pauli Y gate to quantum circuit
   multiply_controlled_single_site!(qc, "Y", control_qubits, pos,
      update_rep=update_rep, num=6, recordEE=recordEE, α=α, cutoff=cutoff)
end


""" Function to apply the (single-qubit) Pauli Z gate gate to a quantum circuit qc
in every position specified in the array pos; an arbitrary number of control
qubits can be given.
Parameters:
- update_rep: save graphical representation of gate.
- recordEE: calculate entanglement entropy after gate.
- α: parameter in Rényi entropy (α=1 for von Neumann).
- cutoff: smallest singular to be kept for entropy calulation. """
function C_PauliZ!(qc::QC_IT_MPS, control_qubits, pos;
   update_rep=true, recordEE=true, α=1, cutoff=1E-6)

   # apply Pauli Z gate to quantum circuit
   multiply_controlled_single_site!(qc, "Z", control_qubits, pos,
      update_rep=update_rep, num=5, recordEE=recordEE, α=α, cutoff=cutoff)
end


""" Function to apply the (single-qubit) √X gate gate to a quantum circuit qc
in every position specified in the array pos; an arbitrary number of control
qubits can be given.
Parameters:
- update_rep: save graphical representation of gate.
- recordEE: calculate entanglement entropy after gate.
- α: parameter in Rényi entropy (α=1 for von Neumann).
- cutoff: smallest singular to be kept for entropy calulation. """
function C_SqrtX!(qc::QC_IT_MPS, control_qubits, pos;
   update_rep=true, recordEE=true, α=1, cutoff=1E-6)

   # apply √X gate to quantum circuit
   multiply_controlled_single_site!(qc, "√X", control_qubits, pos,
      update_rep=update_rep, num=18, recordEE=recordEE, α=α, cutoff=cutoff)
end


""" Function to apply the (single-qubit) S gate gate to a quantum circuit qc
in every position specified in the array pos; an arbitrary number of control
qubits can be given.
Parameters:
- update_rep: save graphical representation of gate.
- recordEE: calculate entanglement entropy after gate.
- α: parameter in Rényi entropy (α=1 for von Neumann).
- cutoff: smallest singular to be kept for entropy calulation. """
function C_SGate!(qc::QC_IT_MPS, control_qubits, pos;
   update_rep=true, recordEE=true, α=1, cutoff=1E-6)

   # apply S gate to quantum circuit
   multiply_controlled_single_site!(qc, "S", control_qubits, pos,
      update_rep=update_rep, num=7, recordEE=recordEE, α=α, cutoff=cutoff)
end


""" Function to apply the (single-qubit) T gate gate to a quantum circuit qc
in every position specified in the array pos; an arbitrary number of control
qubits can be given.
Parameters:
- update_rep: save graphical representation of gate.
- recordEE: calculate entanglement entropy after gate.
- α: parameter in Rényi entropy (α=1 for von Neumann).
- cutoff: smallest singular to be kept for entropy calulation. """
function C_TGate!(qc::QC_IT_MPS, control_qubits, pos;
   update_rep=true, recordEE=true, α=1, cutoff=1E-6)

   # apply T gate to quantum circuit
   multiply_controlled_single_site!(qc, "T", control_qubits, pos,
      update_rep=update_rep, num=8, recordEE=recordEE, α=α, cutoff=cutoff)
end


##################
# Parametric Gates
##################


""" Function to apply the (single-qubit) Rx gate gate to a quantum circuit qc
in a position specified in the array pos; an arbitrary number of control
qubits can be given.
Parameters:
- θ: rotation angle.
- update_rep: save graphical representation of gate.
- recordEE: calculate entanglement entropy after gate.
- α: parameter in Rényi entropy (α=1 for von Neumann).
- cutoff: smallest singular to be kept for entropy calulation. """
function C_RxGate!(qc::QC_IT_MPS, control_qubits, pos, θ::Number;
   update_rep=true, recordEE=true, α=1, cutoff=1E-6)

   # apply Rx gate to quantum circuit
   multiply_controlled_single_site!(qc, "Rx", control_qubits, pos,
      update_rep=update_rep, num=9, recordEE=recordEE, params=[θ],
      α=α, cutoff=cutoff)
end


""" Function to apply the (single-qubit) Ry gate gate to a quantum circuit qc
in a position specified in the array pos; an arbitrary number of control
qubits can be given.
Parameters:
- θ: rotation angle.
- update_rep: save graphical representation of gate.
- recordEE: calculate entanglement entropy after gate.
- α: parameter in Rényi entropy (α=1 for von Neumann).
- cutoff: smallest singular to be kept for entropy calulation. """
function C_RyGate!(qc::QC_IT_MPS, control_qubits, pos, θ::Number;
   update_rep=true, recordEE=true, α=1, cutoff=1E-6)

   # apply Ry gate to quantum circuit
   multiply_controlled_single_site!(qc, "Ry", control_qubits, pos,
      update_rep=update_rep, num=10, recordEE=recordEE, params=[θ],
      α=α, cutoff=cutoff)
end


""" Function to apply the (single-qubit) Rz gate gate to a quantum circuit qc
in a position specified in the array pos; an arbitrary number of control
qubits can be given.
Parameters:
- θ: rotation angle.
- update_rep: save graphical representation of gate.
- recordEE: calculate entanglement entropy after gate.
- α: parameter in Rényi entropy (α=1 for von Neumann).
- cutoff: smallest singular to be kept for entropy calulation. """
function C_RzGate!(qc::QC_IT_MPS, control_qubits, pos, θ::Number;
   update_rep=true, recordEE=true, α=1, cutoff=1E-6)

   # apply Rz gate to quantum circuit
   multiply_controlled_single_site!(qc, "Rz", control_qubits, pos,
      update_rep=update_rep, num=11, recordEE=recordEE, params=[θ], α=α,
      cutoff=cutoff)
end


""" Function to apply the (single-qubit) phase shift gate to a quantum circuit qc
in a position specified in the array pos; an arbitrary number of control
qubits can be given.
Parameters:
- θ: phase angle.
- update_rep: save graphical representation of gate.
- recordEE: calculate entanglement entropy after gate.
- α: parameter in Rényi entropy (α=1 for von Neumann).
- cutoff: smallest singular to be kept for entropy calulation. """
function C_PhaseShift!(qc::QC_IT_MPS, control_qubits, pos, θ::Number;
   update_rep=true, recordEE=true, α=1, cutoff=1E-6)

   # apply P gate to quantum circuit
   multiply_controlled_single_site!(qc, "P", control_qubits, pos,
      update_rep=update_rep, num=12, recordEE=recordEE, params=[θ],
      α=α, cutoff=cutoff)
end


""" Function to apply the (single-qubit) Rn gate to a quantum circuit qc
in a position specified in the array pos; an arbitrary number of control
qubits can be given.
Parameters:
- θ: phase angle.
- update_rep: save graphical representation of gate.
- recordEE: calculate entanglement entropy after gate.
- α: parameter in Rényi entropy (α=1 for von Neumann).
- cutoff: smallest singular to be kept for entropy calulation. """
function C_Rn!(qc::QC_IT_MPS, control_qubits, pos, n::Number;
   update_rep=true, recordEE=true, α=1, cutoff=1E-6)

   # apply Rn gate to quantum circuit
   multiply_controlled_single_site!(qc, "Rn", control_qubits, pos,
      update_rep=update_rep, num=24, recordEE=recordEE, params=[n],
      α=α, cutoff=cutoff)
end


""" Function to apply the (single-qubit) custom-defined unitary operator U to
a quantum circuit qc in a position specified in the array pos. An arbitrary
number of control qubits can be given.
Parameters:
- update_rep: save graphical representation of gate.
- recordEE: calculate entanglement entropy after gate.
- α: parameter in Rényi entropy (α=1 for von Neumann).
- cutoff: smallest singular to be kept for entropy calulation. """
function C_UGate!(qc::QC_IT_MPS, U, control_qubits, pos; update_rep=true,
   recordEE=true, α=1, cutoff=1E-6)

   # apply U gate to quantum circuit
   multiply_controlled_single_site!(qc, "U", control_qubits, pos,
      update_rep=update_rep, num=19, recordEE=recordEE,
      params=[U[1, 1], U[1, 2], U[2, 1], U[2, 2]], α=α, cutoff=cutoff)
end

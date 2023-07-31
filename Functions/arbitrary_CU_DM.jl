

######################################################
## Arbitrary controlled unitary operators - DM version
######################################################


#####################
# Auxiliary functions
#####################


""" Function to update the representation of a quantum circuit qc for an
arbitrary controlled gate. The gate is generically represented by "U". """
function update_representation_arbitrary_CU!(qc::QC_DM, control_qubits,
    action_qubits, num)

    for i in 1:qc.NumQubits
        if i ∈ control_qubits
            push!(qc.Representation[i], 3) # control qubit
        elseif i ∈ action_qubits
            push!(qc.Representation[i], num) # action qubit
        elseif (i ∉ control_qubits) && (i ∉ action_qubits) && (i > minimum(control_qubits)) && (i < minimum(action_qubits))
            push!(qc.Representation[i], 15) # vertical line
        elseif (i ∉ control_qubits) && (i ∉ action_qubits) && (i < maximum(control_qubits)) && (i > maximum(action_qubits))
            push!(qc.Representation[i], 15) # vertical line
        else
            push!(qc.Representation[i], 0)
        end
    end
end


""" Function to construct the N-qubit matrix representing an arbitrarily
controlled quantum gate, potentially on multiple sites. """
function multiply_controlled_single_site!(qc::QC_DM, matrix, control_pos, action_pos;
    update_rep=true, num=19)

    # get matrices and number of qubits
    N = qc.NumQubits
    s = Index(2, "QCircuit")
    E = sparse(array(op("E", s)))
    proj11 = sparse(array(op("Proj11", s))) # |1><1|

    # build trivial part
    gate_triv = E
    for i in 2:N
        gate_triv = kron(gate_triv, E)
    end

    # build non-trivial part
    if 1 ∈ control_pos
        gate_non_triv = proj11
    elseif 1 ∈ action_pos
        gate_non_triv = matrix - E
    else
        gate_non_triv = E
    end
    for i in 2:N
        if i ∈ control_pos
            gate_non_triv = kron(gate_non_triv, proj11)
        elseif i ∈ action_pos
            gate_non_triv = kron(gate_non_triv, matrix - E)
        else
            gate_non_triv = kron(gate_non_triv, E)
        end
    end

    # build controlled gate
    gate =  gate_triv+gate_non_triv

    # update state vector
    qc.StateVector = gate * qc.StateVector * gate'

    # update representing matrix of quantum circuit
    if update_rep
        update_representation_arbitrary_CU!(qc, control_pos, action_pos, num)
    end

    # update circuit depth
    qc.CircuitDepth += 1
end


""" Function to apply multiply-controlled (general) SWAP gates.
CURRENTLY NOT IMPLEMENTED! """
function multiply_controlled_general_SWAP!(qc::QC_DM, control_pos, action_pos;
    update_rep=true, num=4, β=1)

    #
    #
    #
    #

end

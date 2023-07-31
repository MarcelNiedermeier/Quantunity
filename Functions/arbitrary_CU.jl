
#######################################################
## Arbitrary controlled unitary operators - ED version
#######################################################


#####################
# Auxiliary functions
#####################


""" Function to update the representation of a quantum circuit qc for an
arbitrary controlled gate. The gate is generically represented by "U". """
function update_representation_arbitrary_CU!(qc::QC, control_qubits,
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
function multiply_controlled_single_site!(qc::QC, matrix, control_pos, action_pos;
    update_rep=true, num=19)

    # get matrices and number of qubits
    N = qc.NumQubits
    s = Index(2, "QCircuit")
    E = array(op("E", s))
    proj11 = array(op("Proj11", s)) # |1><1|

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
    qc.StateVector .= gate*qc.StateVector

    # update representing matrix of quantum circuit
    if update_rep
        update_representation_arbitrary_CU!(qc, control_pos, action_pos, num)
    end

    # update circuit depth
    qc.CircuitDepth += 1
end


""" Function to apply multiply-controlled (general) SWAP gates.
CURRENTLY NOT IMPLEMENTED! """
function multiply_controlled_general_SWAP!(qc::QC, control_qubits, action_qubits;
    update_rep=true, num=4, β=1)

    #
    #
    #
    #

end






###############
# Old functions
###############



""" Function to get a list of all possible projectors for the control
space for a given number of num_control control qubits.
E.g. for a doubly-controlled gate: [[00, 00], [00, 11], [11, 00]]. """
function get_control_space_projectors(num_control)

    projectors = ["Proj00", "Proj11"]
    projectors_tmp = [projectors for i in 1:num_control]
    proj_product = Iterators.product(projectors_tmp...)
    proj_product = reshape(collect(proj_product), (2^num_control, 1))

    # only keep up to second last element (last one projector for action space!)
    proj = []
    for i in 1:2^num_control-1
        push!(proj, proj_product[i])
    end

    return proj
end


""" Function to create a single term in the sum needed to implemented
general controlled gates. Need to specify a list of projectors, their
positions as well as the "action operator" and its position. """
function get_partial_gate(qc, U, projector_list, control_qubits, action_qubits)

    # check inputs
    if maximum(control_qubits) > qc.NumQubits
        error("Invalid choice of positions for given circuit!")
    elseif maximum(action_qubits) > qc.NumQubits
        error("Invalid choice of positions for given circuit!")
    elseif length(intersect(control_qubits, action_qubits)) != 0
        error("Control qubits must be different from action qubits")
    elseif diff(action_qubits) != ones(length(action_qubits)-1)
        error("Incorrect specification of qubits U acts on (must be adjacent).")
    elseif length(projector_list) != length(control_qubits)
        error("List of projectors doesn't coincide with given positions.")
    end

    num_control = length(control_qubits)
    num_action = length(action_qubits)
    N = qc.NumQubits

    # get matrices
    s = Index(2, "QCircuit")
    E = array(op("E", s))
    proj0 = array(op("Proj00", s)) # |0><0|
    proj1 = array(op("Proj11", s)) # |1><1|



    # gate initialisation if U is applied on qubit 1!



    # initialise gate
    if control_qubits[1] == 1
        gate = array(op(projector_list[1], s))
        proj_counter = 2
    else
        gate = E
        proj_counter = 1
    end

    #println("projector list ", projector_list)

    # qubit lines before U: check if projector needed, insert if yes
    #proj_counter = 2
    for i in 2:(action_qubits[1]-1)
        if i ∈ control_qubits
            gate = kron(gate, array(op(projector_list[proj_counter], s)))
            proj_counter += 1
        else
            gate = kron(gate, E)
        end
    end

    # tensor product with U in "action qubits"
    gate = kron(gate, U)

    # qubit lines after U: check if projector needed, insert if yes
    for i in (action_qubits[end]+1):N
        if i ∈ control_qubits
            gate = kron(gate, array(op(projector_list[proj_counter], s)))
            proj_counter += 1
        else
            gate = kron(gate, E)
        end
    end
    return gate
end


#######
# Gates
#######

""" Function to implement a controlled unitary operator (on adjacent
sites). """
function CU_general!(qc::QC, U, control_qubits, action_qubits, update_rep=true)

    num_control = length(control_qubits)
    num_action = length(action_qubits)
    N = qc.NumQubits

    # get matrices
    s = Index(2, "QCircuit")
    E = array(op("E", s))

    # get projectors for control and action space
    c_proj = get_control_space_projectors(num_control)
    a_proj =  ["Proj11" for i in 1:num_control]

    #println("c_proj ", c_proj)
    #println("a_proj ", a_proj)

    # initialise gate by constructing action space
    gate = get_partial_gate(qc, U, a_proj, control_qubits, action_qubits)

    # get identities in different dimensions
    if num_action == 1
        identity = E
    elseif num_action == 2
        identity = kron(E, E)
    elseif num_action == 3
        identity = kron(kron(E, E), E)
    elseif num_action == 4
        identity = kron(kron(kron(E, E), E), E)
    elseif num_action == 5
        identity = kron(kron(kron(kron(E, E), E), E), E)
    elseif num_action == 6
        identity = kron(kron(kron(kron(kron(E, E), E), E), E), E)
    elseif num_action == 7
        identity = kron(kron(kron(kron(kron(kron(E, E), E), E), E), E), E)
    elseif num_action == 8
        identity = kron(kron(kron(kron(kron(kron(kron(E, E), E), E), E), E), E), E)
    elseif num_action > 8
        error("More than 8 action qubits currently not supported.")
    end

    # finish gate by adding the corresponding control space gates
    for proj in c_proj
        gate += get_partial_gate(qc, identity, proj, control_qubits, action_qubits)
    end

    # update state vector
    qc.StateVector .= gate*qc.StateVector

    # update representing matrix of quantum circuit
    if update_rep
        update_representation_arbitrary_CU!(qc, control_qubits, action_qubits)
    end

    # update circuit depth
    qc.CircuitDepth += 1

end




######
# Test
######

#s = Index(2, "QCircuit")
#E = array(op("E", s))

#X = [0. 1.;
#     1. 0.]

#H = 1/√2 .* [1. 1.;
#            1. -1.]


#U = kron(kron(X, X), X)
#U2 = kron(kron(H, H), H)


# set constants
#N = 8
#maxdim = 32
#N_meas = 5000
#backend = "ED_Julia"
#eps = 0.01
#backend = "MPS_ITensor"
#contmethod = "naive"
#random = false
#lintop = false
#randombond = 2

#qc = initialise_qcircuit(N, lintop, backend)

#num_qubits = 8
#control_qubits = [1, 3]
#action_qubits = [4, 5, 6]
#c_proj = get_control_space_projectors(length(control_qubits))
#println(c_proj)

#projector_list = c_proj[1]

#gate = get_partial_gate(qc, U, projector_list, control_qubits, action_qubits)

#PauliX!(qc, [1, 3])
#sample_measurement(qc, [1, 2, 3, 4, 5, 6, 7, 8], N_meas, eps)
#CU_general!(qc, U2, control_qubits, action_qubits)
#sample_measurement(qc, [1, 2, 3, 4, 5, 6, 7, 8], N_meas, eps)

#draw(qc)

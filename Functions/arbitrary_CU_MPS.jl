
#########################################
## Arbitrary controlled unitary operators
#########################################


#####################
# Auxiliary functions
#####################


""" Function to update the representation of a quantum circuit qc for an
arbitrary controlled gate. The parameter "num" specifies the representation
of the gate. """
function update_representation_arbitrary_CU!(qc::QC_IT_MPS, control_qubits,
    action_qubits, num)

    for i in 1:qc.NumQubits
        if i ∈ control_qubits
            push!(qc.Representation[i], 3) # control qubit
        elseif i ∈ action_qubits
            push!(qc.Representation[i], num) # action qubit
        elseif (i ∉ control_qubits) && (i ∉ action_qubits) && (i > minimum(union(control_qubits, action_qubits))) && (i < maximum(union(control_qubits, action_qubits)))
            push!(qc.Representation[i], 15) # vertical line
        #elseif (i ∉ control_qubits) && (i ∉ action_qubits) && (i < maximum(control_qubits)) && (i > maximum(action_qubits))
        #elseif (i ∉ control_qubits) && (i ∉ action_qubits) && (i < maximum(union(control_qubits, action_qubits)))
        #    push!(qc.Representation[i], 15) # vertical line
        else
            push!(qc.Representation[i], 0)
        end
    end
end


""" Function to update the representation of a quantum circuit qc for an
arbitrary controlled gate. The gate is generically represented by "U". """
function update_representation_arbitrary_CU!(qc::QC_IT_MPS, control_qubits, action_qubits)

    for i in 1:qc.NumQubits
        if i ∈ control_qubits
            push!(qc.Representation[i], 3) # control qubit
        elseif i ∈ action_qubits
            push!(qc.Representation[i], 19) # action qubit
        elseif (i ∉ control_qubits) && (i ∉ action_qubits) && (i > minimum(control_qubits)) && (i < minimum(action_qubits))
            push!(qc.Representation[i], 15) # vertical line
        elseif (i ∉ control_qubits) && (i ∉ action_qubits) && (i < maximum(control_qubits)) && (i > maximum(action_qubits))
            push!(qc.Representation[i], 15) # vertical line
        else
            push!(qc.Representation[i], 0)
        end
    end
end


""" Function to obtain a list of (single-qubit) projectors needed to
describe the matrix element at position index of a N-qubit operator. """
function get_projector_list(N, index)

    # check that index structure matches size of operator!
    if index[1] > 2^N || index[1] > 2^N
        error("Incorrect indices [$(index[1]), $(index[1])] for an $(N)-qubit operator.")
    end

    # get position of element in matrix (decimal)
    position = 2^N * (index[1]-1) + (index[2]-1)

    # convert to binary
    position_bin = reverse(digits(position, base=2, pad=2*N))

    proj_numbers = []
    for i in 1:N
        push!(proj_numbers, [position_bin[i], position_bin[i+N]])
    end

    proj_list = []
    for i in 1:N
        if proj_numbers[i] == [0, 0]
            push!(proj_list, "Proj00")
        elseif proj_numbers[i] == [0, 1]
            push!(proj_list, "Proj01")
        elseif proj_numbers[i] == [1, 0]
            push!(proj_list, "Proj10")
        else # proj_numbers[i] == [1, 1]
            push!(proj_list, "Proj11")
        end
    end

    return proj_list
end


""" Function to get the full list of projectors needed to construct a controlled
unitary operator acting on length(action_qubits) action qubits, with
length(control_qubits) control qubits, at the position specified by index."""
function get_projector_list_full(num_qubits, control_qubits, action_qubits, index)

    N = length(action_qubits)

    # check that index structure matches action qubits!
    if index[1] > 2^N || index[1] > 2^N
        error("Incorrect indices [$(index[1]), $(index[1])] for an $(N)-qubit operator.")
    end

    projector_list_tmp = get_projector_list(N, index)
    projector_list = []

    counter = 1
    for i in 1:num_qubits
        if i ∈ control_qubits
            push!(projector_list, "Proj11")
        elseif i ∈ action_qubits
            push!(projector_list, projector_list_tmp[counter])
            counter += 1
        end
    end
    return projector_list
end


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


#######
# Gates
#######

""" Function to apply multiply-controlled quantum gates on a single
site.  """
function multiply_controlled_single_site!(qc::QC_IT_MPS, U,
    control_qubits, action_qubits; update_rep=true, num=19,
    params=[], recordEE=true, α=1, cutoff=1E-6)

    # check inputs
    if maximum(control_qubits) > qc.NumQubits || maximum(action_qubits) > qc.NumQubits
        error("Invalid choice of positions for given circuit!")
    elseif length(action_qubits)+length(control_qubits) > qc.NumQubits
        error("Number of desired qubits greater than circuit size!")
    elseif length(intersect(control_qubits, action_qubits)) != 0
        error("Control qubits must be different from action qubits")
    elseif U ∈ ["Rx", "Ry", "Rz", "P", "Rn"] && length(action_qubits) > 1
        error("Application of controlled multiple parametric single-site
            gates (Rx, Ry, Rz, P, UGate) not implemented")
    elseif U ∈ ["U"] && length(action_qubits) > 1
        error("Application of controlled multiple parametric single-site
            gates (Rx, Ry, Rz, P, UGate) not implemented")
    end

    # initialise CU gate
    cu = OpSum()

    # layer of identities
    unities = []
    for i in 1:2*qc.NumQubits
        if isodd(i)
            push!(unities, "E")
        else # iseven
            push!(unities, i÷2)
        end
    end
    unities = Tuple(unities) # OpSum() object needs Tuple
    cu += unities

    # subtract final element
    activated_control_space = []
    push!(activated_control_space, -1) # -1 is overall coefficient
    for i in 1:2*length(control_qubits)
        if isodd(i)
            push!(activated_control_space, "Proj11")
        else # iseven
            push!(activated_control_space, control_qubits[i÷2])
        end
    end
    activated_control_space = Tuple(activated_control_space)
    cu += activated_control_space

    # non-parametric gates
    if U ∈ ["X", "Y", "Z", "H", "S", "√X", "T"]

        # set up target space
        target_space = []

        # add (multi-) single-site operators
        for i in 1:2*length(action_qubits)
            if isodd(i)
                push!(target_space, U)
            else # iseven
                push!(target_space, action_qubits[i÷2])
            end
        end

        # activated control space
        for i in 1:2*length(control_qubits)
            if isodd(i)
                push!(target_space, "Proj11")
            else # iseven
                push!(target_space, control_qubits[i÷2])
            end
        end

        target_space = Tuple(target_space)
        cu += target_space

    # gates with single parameter
elseif U ∈ ["Rx", "Ry", "Rz", "P", "Rn"]

        # get gate matrix
        if U ∈ ["Rx", "Ry", "Rz", "P"]
            s = Index(2, "QCircuit")
            gate_matrix = array(op(U, s; θ=params[1]))
        else # Rn gate
            s = Index(2, "QCircuit")
            gate_matrix = array(op(U, s; n=params[1]))
        end

        # set up target space
        target_space00 = [gate_matrix[1, 1], "Proj00", action_qubits[1]]
        target_space01 = [gate_matrix[1, 2], "Proj01", action_qubits[1]]
        target_space10 = [gate_matrix[2, 1], "Proj10", action_qubits[1]]
        target_space11 = [gate_matrix[2, 2], "Proj11", action_qubits[1]]

        # add (multi-) single-site operators
        for i in 1:2*length(action_qubits)
            if isodd(i)
                push!(target_space00, "Proj11")
                push!(target_space01, "Proj11")
                push!(target_space10, "Proj11")
                push!(target_space11, "Proj11")
            else # iseven
                push!(target_space00, control_qubits[i÷2])
                push!(target_space01, control_qubits[i÷2])
                push!(target_space10, control_qubits[i÷2])
                push!(target_space11, control_qubits[i÷2])
            end
        end

        # make Tuple and add to MPO under construction
        target_space00 = Tuple(target_space00)
        target_space01 = Tuple(target_space01)
        target_space10 = Tuple(target_space10)
        target_space11 = Tuple(target_space11)
        cu += target_space00
        cu += target_space01
        cu += target_space10
        cu += target_space11

    # gates with multiple parameters
    elseif U ∈ ["U"]

        # get gate matrix
        s = Index(2, "QCircuit")
        gate_matrix = array(op(U, s; α=params[1], β=params[2], γ=params[3],
            δ=params[4]))

        # set up target space
        target_space00 = [gate_matrix[1, 1], "Proj00", action_qubits[1]]
        target_space01 = [gate_matrix[1, 2], "Proj01", action_qubits[1]]
        target_space10 = [gate_matrix[2, 1], "Proj10", action_qubits[1]]
        target_space11 = [gate_matrix[2, 2], "Proj11", action_qubits[1]]

        # add (multi-) single-site operators
        for i in 1:2*length(action_qubits)
            if isodd(i)
                push!(target_space00, "Proj11")
                push!(target_space01, "Proj11")
                push!(target_space10, "Proj11")
                push!(target_space11, "Proj11")
            else # iseven
                push!(target_space00, control_qubits[i÷2])
                push!(target_space01, control_qubits[i÷2])
                push!(target_space10, control_qubits[i÷2])
                push!(target_space11, control_qubits[i÷2])
            end
        end

        # make Tuple and add to MPO under construction
        target_space00 = Tuple(target_space00)
        target_space01 = Tuple(target_space01)
        target_space10 = Tuple(target_space10)
        target_space11 = Tuple(target_space11)
        cu += target_space00
        cu += target_space01
        cu += target_space10
        cu += target_space11
    end

    # get sites
    sites = qc.IndexSet

    # turn into MPO
    gate = MPO(cu, sites)
    #println("maxlinkdim = ", maxlinkdim(gate))

    # apply to state vector, unprime
    qc.StateVector = contract(gate, qc.StateVector, maxdim=qc.MaxBondDim)#, method=qc.ContractionMethod)
    noprime!(qc.StateVector)

    # update representing matrix of quantum circuit
    if update_rep
        update_representation_arbitrary_CU!(qc, control_qubits, action_qubits, num)
    end

    # update circuit depth and bond dimension
    qc.CircuitDepth += 1
    push!(qc.BondDim, maxlinkdim(qc.StateVector))

    # measure entanglement entropy after application of gate
    if recordEE
        push!(qc.EntanglementEntropy, entanglement_entropy(qc, α=α, cutoff=cutoff))
    end
end




""" Function to apply multiply-controlled (general) SWAP gates.  """
function multiply_controlled_general_SWAP!(qc::QC_IT_MPS,
    control_qubits, action_qubits; update_rep=true, num=4,
    β=1, recordEE=true, α=1, cutoff=1E-6)

    # check inputs
    if length(control_qubits) > 0
        if maximum(control_qubits) > qc.NumQubits || maximum(action_qubits) > qc.NumQubits
            error("Invalid choice of positions for given circuit!")
        elseif length(action_qubits)+length(control_qubits) > qc.NumQubits
            error("Number of desired qubits greater than circuit size!")
        end
    end
    if length(intersect(control_qubits, action_qubits)) != 0
        error("Control qubits must be different from action qubits")
    elseif length(action_qubits) ≠ 2
        error("SWAP gate switches only two qubits (wrong number given).")
    end


    # regular SWAP
    if β == 1

        # SWAP gate not controlled
        if length(control_qubits) == 0

            # initialise SWAP gate
            cswap = OpSum()

            # layer of identities
            unities = []
            for i in 1:2*qc.NumQubits
                if isodd(i)
                    push!(unities, "E")
                else # iseven
                    push!(unities, i÷2)
                end
            end
            prepend!(unities, 1/2)
            cswap += Tuple(unities) # OpSum() object needs Tuple

            # add (multi-) single-site operators
            for U in ["X", "Y", "Z"]
                action_space = []
                for i in 1:2*length(action_qubits)
                    if isodd(i)
                        push!(action_space, U)
                    else # iseven
                        push!(action_space, action_qubits[i÷2])
                    end
                end
                prepend!(action_space, 1/2)
                cswap += Tuple(action_space)
            end

        # controlled SWAP gate
        else

            # initialise SWAP gate
            cswap = OpSum()

            # layer of identities
            unities = []
            for i in 1:2*qc.NumQubits
                if isodd(i)
                    push!(unities, "E")
                else # iseven
                    push!(unities, i÷2)
                end
            end
            prepend!(unities, 1/2)
            cswap += Tuple(unities) # OpSum() object needs Tuple

            # build activated target space
            for U in ["X", "Y", "Z"]
                action_space = []
                for i in 1:2*length(action_qubits)
                    if isodd(i)
                        push!(action_space, U)
                    else # iseven
                        push!(action_space, action_qubits[i÷2])
                    end
                end

                # add activated control space projectors
                for i in 1:2*length(control_qubits)
                    if isodd(i)
                        push!(action_space, "Proj11")
                    else # iseven
                        push!(action_space, control_qubits[i÷2])
                    end
                end
                prepend!(action_space, 1/2)
                cswap += Tuple(action_space)
            end
        end

    else # β ≠ 1, general (SWAP)^β

        # SWAP gate not controlled
        if length(control_qubits) == 0

            # initialise SWAP gate
            cswap = OpSum()

            # add terms
            cswap += "Proj00", action_qubits[1], "Proj00", action_qubits[2]
            cswap += "Proj11", action_qubits[1], "Proj11", action_qubits[2]
            cswap += (exp(1.0im*π*β)+1)/2, "Proj00", action_qubits[1], "Proj11", action_qubits[2]
            cswap += (exp(1.0im*π*β)+1)/2, "Proj11", action_qubits[1], "Proj00", action_qubits[2]
            cswap += (exp(1.0im*π*β)-1)/2, "Proj01", action_qubits[1], "Proj10", action_qubits[2]
            cswap += (exp(1.0im*π*β)-1)/2, "Proj10", action_qubits[1], "Proj01", action_qubits[2]

        # controlled SWAP gate
        else

            # initialise SWAP gate
            cswap = OpSum()

            # layer of identities
            unities = []
            for i in 1:2*qc.NumQubits
                if isodd(i)
                    push!(unities, "E")
                else # iseven
                    push!(unities, i÷2)
                end
            end
            cswap += Tuple(unities) # OpSum() object needs Tuple

            # subtract final element
            activated_control_space = []
            push!(activated_control_space, -1) # -1 is overall coefficient
            for i in 1:2*length(control_qubits)
                if isodd(i)
                    push!(activated_control_space, "Proj11")
                else # iseven
                    push!(activated_control_space, control_qubits[i÷2])
                end
            end
            cswap += Tuple(activated_control_space)

            # add non-trivial coefficients
            action_space1 = ["Proj00", action_qubits[1], "Proj00", action_qubits[2]]
            action_space2 = ["Proj11", action_qubits[1], "Proj11", action_qubits[2]]
            action_space3 = [(exp(1.0im*π*β)+1)/2, "Proj00", action_qubits[1], "Proj11", action_qubits[2]]
            action_space4 = [(exp(1.0im*π*β)+1)/2, "Proj11", action_qubits[1], "Proj00", action_qubits[2]]
            action_space5 = [(exp(1.0im*π*β)-1)/2, "Proj01", action_qubits[1], "Proj10", action_qubits[2]]
            action_space6 = [(exp(1.0im*π*β)-1)/2, "Proj10", action_qubits[1], "Proj01", action_qubits[2]]
            action_spaces = [action_space1, action_space2, action_space3, action_space4, action_space5, action_space5]

            for action_space in action_spaces
                for i in 1:2*length(control_qubits)
                    if isodd(i)
                        push!(action_space, "Proj11")
                    else # iseven
                        push!(action_space, control_qubits[i÷2])
                    end
                end
                cswap += Tuple(action_space)
            end
        end
    end

    # get sites
    sites = qc.IndexSet

    # turn into MPO
    gate = MPO(cswap, sites)
    #println("maxlinkdim = ", maxlinkdim(gate))

    # apply to state vector, unprime
    qc.StateVector = contract(gate, qc.StateVector, maxdim=qc.MaxBondDim)#, method=qc.ContractionMethod)
    noprime!(qc.StateVector)

    # update representing matrix of quantum circuit
    if update_rep
        update_representation_arbitrary_CU!(qc, control_qubits, action_qubits, num)
    end

    # update circuit depth and bond dimension
    qc.CircuitDepth += 1
    push!(qc.BondDim, maxlinkdim(qc.StateVector))

    # measure entanglement entropy after application of gate
    if recordEE
        push!(qc.EntanglementEntropy, entanglement_entropy(qc, α=α, cutoff=cutoff))
    end
end


###############
# Old Functions
###############



""" Function to implement a controlled single-site operator """
function CU_single_site_general!(qc::QC_IT_MPS, U, control_qubits, action_qubits,
    update_rep=true, num=19)

    # check inputs
    if maximum(control_qubits) > qc.NumQubits
        error("Invalid choice of positions for given circuit!")
    elseif length(action_qubits) != 1
        error("Too many action qubits specified (need only 1)!")
    elseif action_qubits[1] > qc.NumQubits
        error("Invalid choice of positions for given circuit!")
    elseif length(intersect(control_qubits, action_qubits)) != 0
        error("Control qubits must be different from action qubits")
    end

    num_control = length(control_qubits)
    num_action = length(action_qubits) # = 1

    # apply general CNOT via SWAP gates for linear topology
    if qc.LinearTopology == true

        println("Linear topology currently not implemented!")


    # apply general controlled unitary gate via direct construction of MPO for master topology
    else

        # initialise CU gate
        cu = OpSum()

        # get list of control space projectors
        c_proj = get_control_space_projectors(num_control)

        # build control space depending on number of control qubits
        if num_control == 1
            for i in 1:length(c_proj)
                cu += c_proj[i][1], control_qubits[1]
            end
        elseif num_control == 2
            for i in 1:length(c_proj)
                #println("c-projectors: ", c_proj[i])
                cu += c_proj[i][1], control_qubits[1], c_proj[i][2], control_qubits[2]
            end
        elseif num_control == 3
            for i in 1:length(c_proj)
                cu += c_proj[i][1], control_qubits[1], c_proj[i][2], control_qubits[2],
                c_proj[i][3], control_qubits[3]
            end
        elseif num_control == 4
            for i in 1:length(c_proj)
                cu += c_proj[i][1], control_qubits[1], c_proj[i][2], control_qubits[2],
                c_proj[i][3], control_qubits[3], c_proj[i][4], control_qubits[4]
            end
        elseif num_control == 5
            for i in 1:length(c_proj)
                cu += c_proj[i][1], control_qubits[1], c_proj[i][2], control_qubits[2],
                c_proj[i][3], control_qubits[3], c_proj[i][4], control_qubits[4],
                c_proj[i][5], control_qubits[5]
            end
        elseif num_control == 6
            for i in 1:length(c_proj)
                cu += c_proj[i][1], control_qubits[1], c_proj[i][2], control_qubits[2],
                c_proj[i][3], control_qubits[3], c_proj[i][4], control_qubits[4],
                c_proj[i][5], control_qubits[5], c_proj[i][6], control_qubits[6]
            end
        elseif num_control == 7
            for i in 1:length(c_proj)
                cu += c_proj[i][1], control_qubits[1], c_proj[i][2], control_qubits[2],
                c_proj[i][3], control_qubits[3], c_proj[i][4], control_qubits[4],
                c_proj[i][5], control_qubits[5], c_proj[i][6], control_qubits[6],
                c_proj[i][7], control_qubits[7]
            end
        elseif num_control == 8
            for i in 1:length(c_proj)
                cu += c_proj[i][1], control_qubits[1], c_proj[i][2], control_qubits[2],
                c_proj[i][3], control_qubits[3], c_proj[i][4], control_qubits[4],
                c_proj[i][5], control_qubits[5], c_proj[i][6], control_qubits[6],
                c_proj[i][7], control_qubits[7], c_proj[i][8], control_qubits[8]
            end
        elseif num_control == 9
            for i in 1:length(c_proj)
                cu += c_proj[i][1], control_qubits[1], c_proj[i][2], control_qubits[2],
                c_proj[i][3], control_qubits[3], c_proj[i][4], control_qubits[4],
                c_proj[i][5], control_qubits[5], c_proj[i][6], control_qubits[6],
                c_proj[i][7], control_qubits[7], c_proj[i][8], control_qubits[8],
                c_proj[i][9], control_qubits[9]
            end
        elseif num_control == 10
            for i in 1:length(c_proj)
                cu += c_proj[i][1], control_qubits[1], c_proj[i][2], control_qubits[2],
                c_proj[i][3], control_qubits[3], c_proj[i][4], control_qubits[4],
                c_proj[i][5], control_qubits[5], c_proj[i][6], control_qubits[6],
                c_proj[i][7], control_qubits[7], c_proj[i][8], control_qubits[8],
                c_proj[i][9], control_qubits[9], c_proj[i][10], control_qubits[10]
            end
        elseif num_control == 11
            for i in 1:length(c_proj)
                cu += c_proj[i][1], control_qubits[1], c_proj[i][2], control_qubits[2],
                c_proj[i][3], control_qubits[3], c_proj[i][4], control_qubits[4],
                c_proj[i][5], control_qubits[5], c_proj[i][6], control_qubits[6],
                c_proj[i][7], control_qubits[7], c_proj[i][8], control_qubits[8],
                c_proj[i][9], control_qubits[9], c_proj[i][10], control_qubits[10],
                c_proj[i][11], control_qubits[11]
            end
        elseif num_control == 12
            for i in 1:length(c_proj)
                cu += c_proj[i][1], control_qubits[1], c_proj[i][2], control_qubits[2],
                c_proj[i][3], control_qubits[3], c_proj[i][4], control_qubits[4],
                c_proj[i][5], control_qubits[5], c_proj[i][6], control_qubits[6],
                c_proj[i][7], control_qubits[7], c_proj[i][8], control_qubits[8],
                c_proj[i][9], control_qubits[9], c_proj[i][10], control_qubits[10],
                c_proj[i][11], control_qubits[11], c_proj[i][12], control_qubits[12]
            end
        elseif num_control == 13
            for i in 1:length(c_proj)
                cu += c_proj[i][1], control_qubits[1], c_proj[i][2], control_qubits[2],
                c_proj[i][3], control_qubits[3], c_proj[i][4], control_qubits[4],
                c_proj[i][5], control_qubits[5], c_proj[i][6], control_qubits[6],
                c_proj[i][7], control_qubits[7], c_proj[i][8], control_qubits[8],
                c_proj[i][9], control_qubits[9], c_proj[i][10], control_qubits[10],
                c_proj[i][11], control_qubits[11], c_proj[i][12], control_qubits[12],
                c_proj[i][13], control_qubits[13]
            end
        elseif num_control == 14
            for i in 1:length(c_proj)
                cu += c_proj[i][1], control_qubits[1], c_proj[i][2], control_qubits[2],
                c_proj[i][3], control_qubits[3], c_proj[i][4], control_qubits[4],
                c_proj[i][5], control_qubits[5], c_proj[i][6], control_qubits[6],
                c_proj[i][7], control_qubits[7], c_proj[i][8], control_qubits[8],
                c_proj[i][9], control_qubits[9], c_proj[i][10], control_qubits[10],
                c_proj[i][11], control_qubits[11], c_proj[i][12], control_qubits[12],
                c_proj[i][13], control_qubits[13], c_proj[i][14], control_qubits[14]
            end
        elseif num_control == 15
            for i in 1:length(c_proj)
                cu += c_proj[i][1], control_qubits[1], c_proj[i][2], control_qubits[2],
                c_proj[i][3], control_qubits[3], c_proj[i][4], control_qubits[4],
                c_proj[i][5], control_qubits[5], c_proj[i][6], control_qubits[6],
                c_proj[i][7], control_qubits[7], c_proj[i][8], control_qubits[8],
                c_proj[i][9], control_qubits[9], c_proj[i][10], control_qubits[10],
                c_proj[i][11], control_qubits[11], c_proj[i][12], control_qubits[12],
                c_proj[i][13], control_qubits[13], c_proj[i][14], control_qubits[14],
                c_proj[i][15], control_qubits[15]
            end
        elseif num_control == 16
            for i in 1:length(c_proj)
                cu += c_proj[i][1], control_qubits[1], c_proj[i][2], control_qubits[2],
                c_proj[i][3], control_qubits[3], c_proj[i][4], control_qubits[4],
                c_proj[i][5], control_qubits[5], c_proj[i][6], control_qubits[6],
                c_proj[i][7], control_qubits[7], c_proj[i][8], control_qubits[8],
                c_proj[i][9], control_qubits[9], c_proj[i][10], control_qubits[10],
                c_proj[i][11], control_qubits[11], c_proj[i][12], control_qubits[12],
                c_proj[i][13], control_qubits[13], c_proj[i][14], control_qubits[14],
                c_proj[i][15], control_qubits[15], c_proj[i][16], control_qubits[16]
            end
        elseif num_control == 17
            for i in 1:length(c_proj)
                cu += c_proj[i][1], control_qubits[1], c_proj[i][2], control_qubits[2],
                c_proj[i][3], control_qubits[3], c_proj[i][4], control_qubits[4],
                c_proj[i][5], control_qubits[5], c_proj[i][6], control_qubits[6],
                c_proj[i][7], control_qubits[7], c_proj[i][8], control_qubits[8],
                c_proj[i][9], control_qubits[9], c_proj[i][10], control_qubits[10],
                c_proj[i][11], control_qubits[11], c_proj[i][12], control_qubits[12],
                c_proj[i][13], control_qubits[13], c_proj[i][14], control_qubits[14],
                c_proj[i][15], control_qubits[15], c_proj[i][16], control_qubits[16],
                c_proj[i][17], control_qubits[17]
            end
        elseif num_control == 18
            for i in 1:length(c_proj)
                cu += c_proj[i][1], control_qubits[1], c_proj[i][2], control_qubits[2],
                c_proj[i][3], control_qubits[3], c_proj[i][4], control_qubits[4],
                c_proj[i][5], control_qubits[5], c_proj[i][6], control_qubits[6],
                c_proj[i][7], control_qubits[7], c_proj[i][8], control_qubits[8],
                c_proj[i][9], control_qubits[9], c_proj[i][10], control_qubits[10],
                c_proj[i][11], control_qubits[11], c_proj[i][12], control_qubits[12],
                c_proj[i][13], control_qubits[13], c_proj[i][14], control_qubits[14],
                c_proj[i][15], control_qubits[15], c_proj[i][16], control_qubits[16],
                c_proj[i][17], control_qubits[17], c_proj[i][18], control_qubits[18]
            end
        elseif num_control == 19
            for i in 1:length(c_proj)
                cu += c_proj[i][1], control_qubits[1], c_proj[i][2], control_qubits[2],
                c_proj[i][3], control_qubits[3], c_proj[i][4], control_qubits[4],
                c_proj[i][5], control_qubits[5], c_proj[i][6], control_qubits[6],
                c_proj[i][7], control_qubits[7], c_proj[i][8], control_qubits[8],
                c_proj[i][9], control_qubits[9], c_proj[i][10], control_qubits[10],
                c_proj[i][11], control_qubits[11], c_proj[i][12], control_qubits[12],
                c_proj[i][13], control_qubits[13], c_proj[i][14], control_qubits[14],
                c_proj[i][15], control_qubits[15], c_proj[i][16], control_qubits[16],
                c_proj[i][17], control_qubits[17], c_proj[i][18], control_qubits[18],
                c_proj[i][19], control_qubits[19]
            end
        elseif num_control == 20
            for i in 1:length(c_proj)
                cu += c_proj[i][1], control_qubits[1], c_proj[i][2], control_qubits[2],
                c_proj[i][3], control_qubits[3], c_proj[i][4], control_qubits[4],
                c_proj[i][5], control_qubits[5], c_proj[i][6], control_qubits[6],
                c_proj[i][7], control_qubits[7], c_proj[i][8], control_qubits[8],
                c_proj[i][9], control_qubits[9], c_proj[i][10], control_qubits[10],
                c_proj[i][11], control_qubits[11], c_proj[i][12], control_qubits[12],
                c_proj[i][13], control_qubits[13], c_proj[i][14], control_qubits[14],
                c_proj[i][15], control_qubits[15], c_proj[i][16], control_qubits[16],
                c_proj[i][17], control_qubits[17], c_proj[i][18], control_qubits[18],
                c_proj[i][19], control_qubits[19], c_proj[i][20], control_qubits[20]
            end
        elseif num_control == 21
            for i in 1:length(c_proj)
                cu += c_proj[i][1], control_qubits[1], c_proj[i][2], control_qubits[2],
                c_proj[i][3], control_qubits[3], c_proj[i][4], control_qubits[4],
                c_proj[i][5], control_qubits[5], c_proj[i][6], control_qubits[6],
                c_proj[i][7], control_qubits[7], c_proj[i][8], control_qubits[8],
                c_proj[i][9], control_qubits[9], c_proj[i][10], control_qubits[10],
                c_proj[i][11], control_qubits[11], c_proj[i][12], control_qubits[12],
                c_proj[i][13], control_qubits[13], c_proj[i][14], control_qubits[14],
                c_proj[i][15], control_qubits[15], c_proj[i][16], control_qubits[16],
                c_proj[i][17], control_qubits[17], c_proj[i][18], control_qubits[18],
                c_proj[i][19], control_qubits[19], c_proj[i][20], control_qubits[20],
                c_proj[i][21], control_qubits[21]
            end
        elseif num_control > 21
            error("More than 21 control qubits currently not supported.")
        end


        # loop through matrix entries and construct action space

        if num_control == 0
            for i in 1:size(U)[1]
                for j in 1:size(U)[2]
                    if U[i, j] != Complex(0.)
                        #proj_list = get_projector_list_full(qc.NumQubits, control_qubits, action_qubits, [i,j])
                        proj_list = get_projector_list(num_action, [i,j])
                        cu += U[i,j], proj_list[1], action_qubits[1]
                    end
                end
            end
        elseif num_control == 1
            for i in 1:size(U)[1]
                for j in 1:size(U)[2]
                    if U[i, j] != Complex(0.)
                        proj_list = get_projector_list(num_action, [i,j])
                        cu += U[i,j], "Proj11", control_qubits[1], proj_list[1], action_qubits[1]
                    end
                end
            end
        elseif num_control == 2
            #println("I am here!")
            #println("control ", control_qubits)
            #println("action ", action_qubits)
            for i in 1:size(U)[1]
                for j in 1:size(U)[2]
                    if U[i, j] != Complex(0.)
                        proj_list = get_projector_list(num_action, [i,j])
                        cu += U[i,j], "Proj11", control_qubits[1], "Proj11", control_qubits[2], proj_list[1], action_qubits[1]
                    end
                end
            end
        elseif num_control == 3
            for i in 1:size(U)[1]
                for j in 1:size(U)[2]
                    if U[i, j] != Complex(0.)
                        proj_list = get_projector_list(num_action, [i,j])
                        cu += U[i,j], "Proj11", control_qubits[1], "Proj11", control_qubits[2], "Proj11", control_qubits[3],
                        proj_list[1], action_qubits[1]
                    end
                end
            end
        elseif num_control == 4
            for i in 1:size(U)[1]
                for j in 1:size(U)[2]
                    if U[i, j] != Complex(0.)
                        proj_list = get_projector_list(num_action, [i,j])
                        cu += U[i,j], "Proj11", control_qubits[1], "Proj11", control_qubits[2], "Proj11", control_qubits[3],
                        "Proj11", control_qubits[4], proj_list[1], action_qubits[1]
                    end
                end
            end
        elseif num_control == 5
            for i in 1:size(U)[1]
                for j in 1:size(U)[2]
                    if U[i, j] != Complex(0.)
                        proj_list = get_projector_list(num_action, [i,j])
                        cu += U[i,j], "Proj11", control_qubits[1], "Proj11", control_qubits[2], "Proj11", control_qubits[3],
                        "Proj11", control_qubits[4], "Proj11", control_qubits[5], proj_list[1], action_qubits[1]
                    end
                end
            end
        elseif num_control == 6
            for i in 1:size(U)[1]
                for j in 1:size(U)[2]
                    if U[i, j] != Complex(0.)
                        proj_list = get_projector_list(num_action, [i,j])
                        cu += U[i,j], "Proj11", control_qubits[1], "Proj11", control_qubits[2], "Proj11", control_qubits[3],
                        "Proj11", control_qubits[4], "Proj11", control_qubits[5], "Proj11", control_qubits[6],
                        proj_list[1], action_qubits[1]
                    end
                end
            end
        elseif num_control == 7
            for i in 1:size(U)[1]
                for j in 1:size(U)[2]
                    if U[i, j] != Complex(0.)
                        proj_list = get_projector_list(num_action, [i,j])
                        cu += U[i,j], "Proj11", control_qubits[1], "Proj11", control_qubits[2], "Proj11", control_qubits[3],
                        "Proj11", control_qubits[4], "Proj11", control_qubits[5], "Proj11", control_qubits[6],
                        "Proj11", control_qubits[7], proj_list[1], action_qubits[1]
                    end
                end
            end
        elseif num_control == 8
            for i in 1:size(U)[1]
                for j in 1:size(U)[2]
                    if U[i, j] != Complex(0.)
                        proj_list = get_projector_list(num_action, [i,j])
                        cu += U[i,j], "Proj11", control_qubits[1], "Proj11", control_qubits[2], "Proj11", control_qubits[3],
                        "Proj11", control_qubits[4], "Proj11", control_qubits[5], "Proj11", control_qubits[6],
                        "Proj11", control_qubits[7], "Proj11", control_qubits[8], proj_list[1], action_qubits[1]
                    end
                end
            end
        elseif num_control == 9
            for i in 1:size(U)[1]
                for j in 1:size(U)[2]
                    if U[i, j] != Complex(0.)
                        proj_list = get_projector_list(num_action, [i,j])
                        cu += U[i,j], "Proj11", control_qubits[1], "Proj11", control_qubits[2], "Proj11", control_qubits[3],
                        "Proj11", control_qubits[4], "Proj11", control_qubits[5], "Proj11", control_qubits[6],
                        "Proj11", control_qubits[7], "Proj11", control_qubits[8], "Proj11", control_qubits[9], proj_list[1], action_qubits[1]
                    end
                end
            end
        elseif num_control == 10
            for i in 1:size(U)[1]
                for j in 1:size(U)[2]
                    if U[i, j] != Complex(0.)
                        proj_list = get_projector_list(num_action, [i,j])
                        cu += U[i,j], "Proj11", control_qubits[1], "Proj11", control_qubits[2], "Proj11", control_qubits[3],
                        "Proj11", control_qubits[4], "Proj11", control_qubits[5], "Proj11", control_qubits[6],
                        "Proj11", control_qubits[7], "Proj11", control_qubits[8], "Proj11", control_qubits[9],
                        "Proj11", control_qubits[10], proj_list[1], action_qubits[1]
                    end
                end
            end
        elseif num_control == 11
            for i in 1:size(U)[1]
                for j in 1:size(U)[2]
                    if U[i, j] != Complex(0.)
                        proj_list = get_projector_list(num_action, [i,j])
                        cu += U[i,j], "Proj11", control_qubits[1], "Proj11", control_qubits[2], "Proj11", control_qubits[3],
                        "Proj11", control_qubits[4], "Proj11", control_qubits[5], "Proj11", control_qubits[6],
                        "Proj11", control_qubits[7], "Proj11", control_qubits[8], "Proj11", control_qubits[9],
                        "Proj11", control_qubits[10], "Proj11", control_qubits[11], proj_list[1], action_qubits[1]
                    end
                end
            end
        elseif num_control == 12
            for i in 1:size(U)[1]
                for j in 1:size(U)[2]
                    if U[i, j] != Complex(0.)
                        proj_list = get_projector_list(num_action, [i,j])
                        cu += U[i,j], "Proj11", control_qubits[1], "Proj11", control_qubits[2], "Proj11", control_qubits[3],
                        "Proj11", control_qubits[4], "Proj11", control_qubits[5], "Proj11", control_qubits[6],
                        "Proj11", control_qubits[7], "Proj11", control_qubits[8], "Proj11", control_qubits[9],
                        "Proj11", control_qubits[10], "Proj11", control_qubits[11], "Proj11", control_qubits[12],
                        proj_list[1], action_qubits[1]
                    end
                end
            end
        elseif num_control == 13
            for i in 1:size(U)[1]
                for j in 1:size(U)[2]
                    if U[i, j] != Complex(0.)
                        proj_list = get_projector_list(num_action, [i,j])
                        cu += U[i,j], "Proj11", control_qubits[1], "Proj11", control_qubits[2], "Proj11", control_qubits[3],
                        "Proj11", control_qubits[4], "Proj11", control_qubits[5], "Proj11", control_qubits[6],
                        "Proj11", control_qubits[7], "Proj11", control_qubits[8], "Proj11", control_qubits[9],
                        "Proj11", control_qubits[10], "Proj11", control_qubits[11], "Proj11", control_qubits[12],
                        "Proj11", control_qubits[13], proj_list[1], action_qubits[1]
                    end
                end
            end
        elseif num_control == 14
            for i in 1:size(U)[1]
                for j in 1:size(U)[2]
                    if U[i, j] != Complex(0.)
                        proj_list = get_projector_list(num_action, [i,j])
                        cu += U[i,j], "Proj11", control_qubits[1], "Proj11", control_qubits[2], "Proj11", control_qubits[3],
                        "Proj11", control_qubits[4], "Proj11", control_qubits[5], "Proj11", control_qubits[6],
                        "Proj11", control_qubits[7], "Proj11", control_qubits[8], "Proj11", control_qubits[9],
                        "Proj11", control_qubits[10], "Proj11", control_qubits[11], "Proj11", control_qubits[12],
                        "Proj11", control_qubits[13], "Proj11", control_qubits[14], proj_list[1], action_qubits[1]
                    end
                end
            end
        elseif num_control == 15
            for i in 1:size(U)[1]
                for j in 1:size(U)[2]
                    if U[i, j] != Complex(0.)
                        proj_list = get_projector_list(num_action, [i,j])
                        cu += U[i,j], "Proj11", control_qubits[1], "Proj11", control_qubits[2], "Proj11", control_qubits[3],
                        "Proj11", control_qubits[4], "Proj11", control_qubits[5], "Proj11", control_qubits[6],
                        "Proj11", control_qubits[7], "Proj11", control_qubits[8], "Proj11", control_qubits[9],
                        "Proj11", control_qubits[10], "Proj11", control_qubits[11], "Proj11", control_qubits[12],
                        "Proj11", control_qubits[13], "Proj11", control_qubits[14], "Proj11", control_qubits[15],
                        proj_list[1], action_qubits[1]
                    end
                end
            end
        elseif num_control == 16
            for i in 1:size(U)[1]
                for j in 1:size(U)[2]
                    if U[i, j] != Complex(0.)
                        proj_list = get_projector_list(num_action, [i,j])
                        cu += U[i,j], "Proj11", control_qubits[1], "Proj11", control_qubits[2], "Proj11", control_qubits[3],
                        "Proj11", control_qubits[4], "Proj11", control_qubits[5], "Proj11", control_qubits[6],
                        "Proj11", control_qubits[7], "Proj11", control_qubits[8], "Proj11", control_qubits[9],
                        "Proj11", control_qubits[10], "Proj11", control_qubits[11], "Proj11", control_qubits[12],
                        "Proj11", control_qubits[13], "Proj11", control_qubits[14], "Proj11", control_qubits[15],
                        "Proj11", control_qubits[16], proj_list[1], action_qubits[1]
                    end
                end
            end
        elseif num_control == 17
            for i in 1:size(U)[1]
                for j in 1:size(U)[2]
                    if U[i, j] != Complex(0.)
                        proj_list = get_projector_list(num_action, [i,j])
                        cu += U[i,j], "Proj11", control_qubits[1], "Proj11", control_qubits[2], "Proj11", control_qubits[3],
                        "Proj11", control_qubits[4], "Proj11", control_qubits[5], "Proj11", control_qubits[6],
                        "Proj11", control_qubits[7], "Proj11", control_qubits[8], "Proj11", control_qubits[9],
                        "Proj11", control_qubits[10], "Proj11", control_qubits[11], "Proj11", control_qubits[12],
                        "Proj11", control_qubits[13], "Proj11", control_qubits[14], "Proj11", control_qubits[15],
                        "Proj11", control_qubits[16], "Proj11", control_qubits[17], proj_list[1], action_qubits[1]
                    end
                end
            end
        elseif num_control == 18
            for i in 1:size(U)[1]
                for j in 1:size(U)[2]
                    if U[i, j] != Complex(0.)
                        proj_list = get_projector_list(num_action, [i,j])
                        cu += U[i,j], "Proj11", control_qubits[1], "Proj11", control_qubits[2], "Proj11", control_qubits[3],
                        "Proj11", control_qubits[4], "Proj11", control_qubits[5], "Proj11", control_qubits[6],
                        "Proj11", control_qubits[7], "Proj11", control_qubits[8], "Proj11", control_qubits[9],
                        "Proj11", control_qubits[10], "Proj11", control_qubits[11], "Proj11", control_qubits[12],
                        "Proj11", control_qubits[13], "Proj11", control_qubits[14], "Proj11", control_qubits[15],
                        "Proj11", control_qubits[16], "Proj11", control_qubits[17], "Proj11", control_qubits[18],
                        proj_list[1], action_qubits[1]
                    end
                end
            end
        elseif num_control == 19
            for i in 1:size(U)[1]
                for j in 1:size(U)[2]
                    if U[i, j] != Complex(0.)
                        proj_list = get_projector_list(num_action, [i,j])
                        cu += U[i,j], "Proj11", control_qubits[1], "Proj11", control_qubits[2], "Proj11", control_qubits[3],
                        "Proj11", control_qubits[4], "Proj11", control_qubits[5], "Proj11", control_qubits[6],
                        "Proj11", control_qubits[7], "Proj11", control_qubits[8], "Proj11", control_qubits[9],
                        "Proj11", control_qubits[10], "Proj11", control_qubits[11], "Proj11", control_qubits[12],
                        "Proj11", control_qubits[13], "Proj11", control_qubits[14], "Proj11", control_qubits[15],
                        "Proj11", control_qubits[16], "Proj11", control_qubits[17], "Proj11", control_qubits[18],
                        "Proj11", control_qubits[19], proj_list[1], action_qubits[1]
                    end
                end
            end
        elseif num_control == 20
            for i in 1:size(U)[1]
                for j in 1:size(U)[2]
                    if U[i, j] != Complex(0.)
                        proj_list = get_projector_list(num_action, [i,j])
                        cu += U[i,j], "Proj11", control_qubits[1], "Proj11", control_qubits[2], "Proj11", control_qubits[3],
                        "Proj11", control_qubits[4], "Proj11", control_qubits[5], "Proj11", control_qubits[6],
                        "Proj11", control_qubits[7], "Proj11", control_qubits[8], "Proj11", control_qubits[9],
                        "Proj11", control_qubits[10], "Proj11", control_qubits[11], "Proj11", control_qubits[12],
                        "Proj11", control_qubits[13], "Proj11", control_qubits[14], "Proj11", control_qubits[15],
                        "Proj11", control_qubits[16], "Proj11", control_qubits[17], "Proj11", control_qubits[18],
                        "Proj11", control_qubits[19], "Proj11", control_qubits[20], proj_list[1], action_qubits[1]
                    end
                end
            end
        elseif num_control == 21
            for i in 1:size(U)[1]
                for j in 1:size(U)[2]
                    if U[i, j] != Complex(0.)
                        proj_list = get_projector_list(num_action, [i,j])
                        cu += U[i,j], "Proj11", control_qubits[1], "Proj11", control_qubits[2], "Proj11", control_qubits[3],
                        "Proj11", control_qubits[4], "Proj11", control_qubits[5], "Proj11", control_qubits[6],
                        "Proj11", control_qubits[7], "Proj11", control_qubits[8], "Proj11", control_qubits[9],
                        "Proj11", control_qubits[10], "Proj11", control_qubits[11], "Proj11", control_qubits[12],
                        "Proj11", control_qubits[13], "Proj11", control_qubits[14], "Proj11", control_qubits[15],
                        "Proj11", control_qubits[16], "Proj11", control_qubits[17], "Proj11", control_qubits[18],
                        "Proj11", control_qubits[19], "Proj11", control_qubits[20], "Proj11", control_qubits[21],
                        proj_list[1], action_qubits[1]
                    end
                end
            end
        elseif num_control > 21
            error("More than 21 control qubits currently not supported.")
        end

        # get sites
        sites = qc.IndexSet

        # turn into MPO
        gate = MPO(cu, sites)
        #println("maxlinkdim = ", maxlinkdim(gate))
        #println(gate)

        # apply to state vector, unprime
        qc.StateVector = contract(gate, qc.StateVector, maxdim=qc.MaxBondDim)#, method=qc.ContractionMethod)
        noprime!(qc.StateVector)

        # update representing matrix of quantum circuit
        if update_rep
            update_representation_arbitrary_CU!(qc, control_qubits, action_qubits, num)
        end

        # update circuit depth and bond dimension
        qc.CircuitDepth += 1
        push!(qc.BondDim, maxlinkdim(qc.StateVector))
    end
end


""" Function to implement a controlled unitary operator (on adjacent
sites). """
function CU_general!(qc::QC_IT_MPS, U, control_qubits, action_qubits,
    update_rep=true, num=19)

    # check inputs
    if maximum(control_qubits) > qc.NumQubits
        error("Invalid choice of positions for given circuit!")
    elseif maximum(action_qubits) > qc.NumQubits
        error("Invalid choice of positions for given circuit!")
    elseif length(intersect(control_qubits, action_qubits)) != 0
        error("Control qubits must be different from action qubits")
    elseif diff(action_qubits) != ones(length(action_qubits)-1)
        error("Incorrect specification of qubits U acts on (must be adjacent).")
    end

    num_control = length(control_qubits)
    num_action = length(action_qubits)

    # apply general CNOT via SWAP gates for linear topology
    if qc.LinearTopology == true

        println("Linear topology currently not implemented!")


    # apply general controlled unitary gate via direct construction of MPO for master topology
    else

        # initialise CU gate
        cu = OpSum()

        # get list of control space projectors
        c_proj = get_control_space_projectors(num_control)

        # build control space depending on number of control qubits
        if num_control == 1
            for i in 1:length(c_proj)
                cu += c_proj[i][1], control_qubits[1]
            end
        elseif num_control == 2
            for i in 1:length(c_proj)
                #println("c-projectors: ", c_proj[i])
                cu += c_proj[i][1], control_qubits[1], c_proj[i][2], control_qubits[2]
            end
        elseif num_control == 3
            for i in 1:length(c_proj)
                cu += c_proj[i][1], control_qubits[1], c_proj[i][2], control_qubits[2],
                c_proj[i][3], control_qubits[3]
            end
        elseif num_control == 4
            for i in 1:length(c_proj)
                cu += c_proj[i][1], control_qubits[1], c_proj[i][2], control_qubits[2],
                c_proj[i][3], control_qubits[3], c_proj[i][4], control_qubits[4]
            end
        elseif num_control == 5
            for i in 1:length(c_proj)
                cu += c_proj[i][1], control_qubits[1], c_proj[i][2], control_qubits[2],
                c_proj[i][3], control_qubits[3], c_proj[i][4], control_qubits[4],
                c_proj[i][5], control_qubits[5]
            end
        elseif num_control > 5
            error("More than 5 control qubits currently not supported.")
        end


        # loop through matrix entries and construct action space

        ###################
        if num_control == 0
        ###################

            if num_action == 1
                for i in 1:size(U)[1]
                    for j in 1:size(U)[2]
                        #proj_list = get_projector_list_full(qc.NumQubits, control_qubits, action_qubits, [i,j])
                        proj_list = get_projector_list(num_action, [i,j])
                        cu += U[i,j], proj_list[1], action_qubits[1]
                    end
                end

            elseif num_action == 2
                for i in 1:size(U)[1]
                    for j in 1:size(U)[2]
                        #proj_list = get_projector_list_full(qc.NumQubits, control_qubits, action_qubits, [i,j])
                        proj_list = get_projector_list(num_action, [i,j])
                        cu += U[i,j], proj_list[1], action_qubits[1], proj_list[2], action_qubits[2]
                    end
                end

            elseif num_action == 3
                for i in 1:size(U)[1]
                    for j in 1:size(U)[2]
                        #proj_list = get_projector_list_full(qc.NumQubits, control_qubits, action_qubits, [i,j])
                        proj_list = get_projector_list(num_action, [i,j])
                        cu += U[i,j], proj_list[1], action_qubits[1], proj_list[2], action_qubits[2],
                        proj_list[3], action_qubits[3]
                    end
                end

            elseif num_action == 4
                for i in 1:size(U)[1]
                    for j in 1:size(U)[2]
                        #proj_list = get_projector_list_full(qc.NumQubits, control_qubits, action_qubits, [i,j])
                        proj_list = get_projector_list(num_action, [i,j])
                        cu += U[i,j], proj_list[1], action_qubits[1], proj_list[2], action_qubits[2],
                        proj_list[3], action_qubits[3], proj_list[4], action_qubits[4]
                    end
                end

            elseif num_action == 5
                for i in 1:size(U)[1]
                    for j in 1:size(U)[2]
                        #proj_list = get_projector_list_full(qc.NumQubits, control_qubits, action_qubits, [i,j])
                        proj_list = get_projector_list(num_action, [i,j])
                        cu += U[i,j], proj_list[1], action_qubits[1], proj_list[2], action_qubits[2],
                        proj_list[3], action_qubits[3], proj_list[4], action_qubits[4], proj_list[5], action_qubits[5]
                    end
                end

            elseif num_action == 6
                for i in 1:size(U)[1]
                    for j in 1:size(U)[2]
                        #proj_list = get_projector_list_full(qc.NumQubits, control_qubits, action_qubits, [i,j])
                        proj_list = get_projector_list(num_action, [i,j])
                        cu += U[i,j], proj_list[1], action_qubits[1], proj_list[2], action_qubits[2],
                        proj_list[3], action_qubits[3], proj_list[4], action_qubits[4], proj_list[5], action_qubits[5],
                        proj_list[6], action_qubits[6]
                    end
                end

            elseif num_action == 7
                for i in 1:size(U)[1]
                    for j in 1:size(U)[2]
                        #proj_list = get_projector_list_full(qc.NumQubits, control_qubits, action_qubits, [i,j])
                        proj_list = get_projector_list(num_action, [i,j])
                        cu += U[i,j], proj_list[1], action_qubits[1], proj_list[2], action_qubits[2],
                        proj_list[3], action_qubits[3], proj_list[4], action_qubits[4], proj_list[5], action_qubits[5],
                        proj_list[6], action_qubits[6], proj_list[7], action_qubits[7]
                    end
                end

            elseif num_action == 8
                for i in 1:size(U)[1]
                    for j in 1:size(U)[2]
                        #proj_list = get_projector_list_full(qc.NumQubits, control_qubits, action_qubits, [i,j])
                        proj_list = get_projector_list(num_action, [i,j])
                        cu += U[i,j], proj_list[1], action_qubits[1], proj_list[2], action_qubits[2],
                        proj_list[3], action_qubits[3], proj_list[4], action_qubits[4], proj_list[5], action_qubits[5],
                        proj_list[6], action_qubits[6], proj_list[7], action_qubits[7], proj_list[8], action_qubits[8]
                    end
                end

            elseif num_action > 8
                error("More than 8 action qubits currently not supported.")
            end

        #######################
        elseif num_control == 1
        #######################

            if num_action == 1
                for i in 1:size(U)[1]
                    for j in 1:size(U)[2]
                        proj_list = get_projector_list(num_action, [i,j])
                        cu += U[i,j], "Proj11", control_qubits[1], proj_list[1], action_qubits[1]
                    end
                end

            elseif num_action == 2
                #println("I am here")
                for i in 1:size(U)[1]
                    for j in 1:size(U)[2]
                        proj_list = get_projector_list(num_action, [i,j])
                        #println("projectors: ", proj_list)
                        #println("U[$i, $j] = $(U[i,j])")
                        cu += U[i,j], "Proj11", control_qubits[1], proj_list[1], action_qubits[1],
                        proj_list[2], action_qubits[2]
                    end
                end

            elseif num_action == 3
                for i in 1:size(U)[1]
                    for j in 1:size(U)[2]
                        proj_list = get_projector_list(num_action, [i,j])
                        cu += U[i,j], "Proj11", control_qubits[1], proj_list[1], action_qubits[1],
                        proj_list[2], action_qubits[2], proj_list[3], action_qubits[3]
                    end
                end

            elseif num_action == 4
                for i in 1:size(U)[1]
                    for j in 1:size(U)[2]
                        proj_list = get_projector_list(num_action, [i,j])
                        cu += U[i,j], "Proj11", control_qubits[1], proj_list[1], action_qubits[1],
                        proj_list[2], action_qubits[2], proj_list[3], action_qubits[3], proj_list[4], action_qubits[4]
                    end
                end

            elseif num_action == 5
                for i in 1:size(U)[1]
                    for j in 1:size(U)[2]
                        proj_list = get_projector_list(num_action, [i,j])
                        cu += U[i,j], "Proj11", control_qubits[1], proj_list[1], action_qubits[1],
                        proj_list[2], action_qubits[2], proj_list[3], action_qubits[3], proj_list[4], action_qubits[4],
                        proj_list[5], action_qubits[5]
                    end
                end

            elseif num_action == 6
                for i in 1:size(U)[1]
                    for j in 1:size(U)[2]
                        proj_list = get_projector_list(num_action, [i,j])
                        cu += U[i,j], "Proj11", control_qubits[1], proj_list[1], action_qubits[1],
                        proj_list[2], action_qubits[2], proj_list[3], action_qubits[3], proj_list[4], action_qubits[4],
                        proj_list[5], action_qubits[5], proj_list[6], action_qubits[6]
                    end
                end

            elseif num_action == 7
                for i in 1:size(U)[1]
                    for j in 1:size(U)[2]
                        proj_list = get_projector_list(num_action, [i,j])
                        cu += U[i,j], "Proj11", control_qubits[1], proj_list[1], action_qubits[1],
                        proj_list[2], action_qubits[2], proj_list[3], action_qubits[3], proj_list[4], action_qubits[4],
                        proj_list[5], action_qubits[5], proj_list[6], action_qubits[6], proj_list[7], action_qubits[7]
                    end
                end

            elseif num_action == 8
                for i in 1:size(U)[1]
                    for j in 1:size(U)[2]
                        proj_list = get_projector_list(num_action, [i,j])
                        cu += U[i,j], "Proj11", control_qubits[1], proj_list[1], action_qubits[1],
                        proj_list[2], action_qubits[2], proj_list[3], action_qubits[3], proj_list[4], action_qubits[4],
                        proj_list[5], action_qubits[5], proj_list[6], action_qubits[6], proj_list[7], action_qubits[7],
                        proj_list[8], action_qubits[8]
                    end
                end

            elseif num_action > 8
                error("More than 8 action qubits currently not supported.")
            end

        #######################
        elseif num_control == 2
        #######################

            if num_action == 1
                #println("I am here!")
                #println("control ", control_qubits)
                #println("action ", action_qubits)
                for i in 1:size(U)[1]
                    for j in 1:size(U)[2]
                        proj_list = get_projector_list(num_action, [i,j])
                        cu += U[i,j], "Proj11", control_qubits[1], "Proj11", control_qubits[2], proj_list[1], action_qubits[1]
                    end
                end

            elseif num_action == 2
                for i in 1:size(U)[1]
                    for j in 1:size(U)[2]
                        proj_list = get_projector_list(num_action, [i,j])
                        cu += U[i,j], "Proj11", control_qubits[1], "Proj11", control_qubits[2], proj_list[1], action_qubits[1],
                        proj_list[2], action_qubits[2]
                    end
                end

            elseif num_action == 3
                for i in 1:size(U)[1]
                    for j in 1:size(U)[2]
                        proj_list = get_projector_list(num_action, [i,j])
                        cu += U[i,j], "Proj11", control_qubits[1], "Proj11", control_qubits[2], proj_list[1], action_qubits[1],
                        proj_list[2], action_qubits[2], proj_list[3], action_qubits[3]
                    end
                end

            elseif num_action == 4
                for i in 1:size(U)[1]
                    for j in 1:size(U)[2]
                        proj_list = get_projector_list(num_action, [i,j])
                        cu += U[i,j], "Proj11", control_qubits[1], "Proj11", control_qubits[2], proj_list[1], action_qubits[1],
                        proj_list[2], action_qubits[2], proj_list[3], action_qubits[3], proj_list[4], action_qubits[4]
                    end
                end

            elseif num_action == 5
                for i in 1:size(U)[1]
                    for j in 1:size(U)[2]
                        proj_list = get_projector_list(num_action, [i,j])
                        cu += U[i,j], "Proj11", control_qubits[1], "Proj11", control_qubits[2], proj_list[1], action_qubits[1],
                        proj_list[2], action_qubits[2], proj_list[3], action_qubits[3], proj_list[4], action_qubits[4],
                        proj_list[5], action_qubits[5]
                    end
                end

            elseif num_action == 6
                for i in 1:size(U)[1]
                    for j in 1:size(U)[2]
                        proj_list = get_projector_list(num_action, [i,j])
                        cu += U[i,j], "Proj11", control_qubits[1], "Proj11", control_qubits[2], proj_list[1], action_qubits[1],
                        proj_list[2], action_qubits[2], proj_list[3], action_qubits[3], proj_list[4], action_qubits[4],
                        proj_list[5], action_qubits[5], proj_list[6], action_qubits[6]
                    end
                end

            elseif num_action == 7
                for i in 1:size(U)[1]
                    for j in 1:size(U)[2]
                        proj_list = get_projector_list(num_action, [i,j])
                        cu += U[i,j], "Proj11", control_qubits[1], "Proj11", control_qubits[2], proj_list[1], action_qubits[1],
                        proj_list[2], action_qubits[2], proj_list[3], action_qubits[3], proj_list[4], action_qubits[4],
                        proj_list[5], action_qubits[5], proj_list[6], action_qubits[6], proj_list[7], action_qubits[7]
                    end
                end

            elseif num_action == 8
                for i in 1:size(U)[1]
                    for j in 1:size(U)[2]
                        proj_list = get_projector_list(num_action, [i,j])
                        cu += U[i,j], "Proj11", control_qubits[1], "Proj11", control_qubits[2], proj_list[1], action_qubits[1],
                        proj_list[2], action_qubits[2], proj_list[3], action_qubits[3], proj_list[4], action_qubits[4],
                        proj_list[5], action_qubits[5], proj_list[6], action_qubits[6], proj_list[7], action_qubits[7],
                        proj_list[8], action_qubits[8]
                    end
                end

            elseif num_action > 8
                error("More than 8 action qubits currently not supported.")
            end

        #######################
        elseif num_control == 3
        #######################

            if num_action == 1
                for i in 1:size(U)[1]
                    for j in 1:size(U)[2]
                        proj_list = get_projector_list(num_action, [i,j])
                        cu += U[i,j], "Proj11", control_qubits[1], "Proj11", control_qubits[2], "Proj11", control_qubits[3],
                        proj_list[1], action_qubits[1]
                    end
                end

            elseif num_action == 2
                for i in 1:size(U)[1]
                    for j in 1:size(U)[2]
                        proj_list = get_projector_list(num_action, [i,j])
                        cu += U[i,j], "Proj11", control_qubits[1], "Proj11", control_qubits[2], "Proj11", control_qubits[3],
                        proj_list[1], action_qubits[1], proj_list[2], action_qubits[2]
                    end
                end

            elseif num_action == 3
                for i in 1:size(U)[1]
                    for j in 1:size(U)[2]
                        proj_list = get_projector_list(num_action, [i,j])
                        cu += U[i,j], "Proj11", control_qubits[1], "Proj11", control_qubits[2], "Proj11", control_qubits[3],
                        proj_list[1], action_qubits[1], proj_list[2], action_qubits[2], proj_list[3], action_qubits[3]
                    end
                end

            elseif num_action == 4
                for i in 1:size(U)[1]
                    for j in 1:size(U)[2]
                        proj_list = get_projector_list(num_action, [i,j])
                        cu += U[i,j], "Proj11", control_qubits[1], "Proj11", control_qubits[2], "Proj11", control_qubits[3],
                        proj_list[1], action_qubits[1], proj_list[2], action_qubits[2], proj_list[3], action_qubits[3],
                        proj_list[4], action_qubits[4]
                    end
                end

            elseif num_action == 5
                for i in 1:size(U)[1]
                    for j in 1:size(U)[2]
                        proj_list = get_projector_list(num_action, [i,j])
                        cu += U[i,j], "Proj11", control_qubits[1], "Proj11", control_qubits[2], "Proj11", control_qubits[3],
                        proj_list[1], action_qubits[1], proj_list[2], action_qubits[2], proj_list[3], action_qubits[3],
                        proj_list[4], action_qubits[4], proj_list[5], action_qubits[5]
                    end
                end

            elseif num_action == 6
                for i in 1:size(U)[1]
                    for j in 1:size(U)[2]
                        proj_list = get_projector_list(num_action, [i,j])
                        cu += U[i,j], "Proj11", control_qubits[1], "Proj11", control_qubits[2], "Proj11", control_qubits[3],
                        proj_list[1], action_qubits[1], proj_list[2], action_qubits[2], proj_list[3], action_qubits[3],
                        proj_list[4], action_qubits[4], proj_list[5], action_qubits[5], proj_list[6], action_qubits[6]
                    end
                end

            elseif num_action == 7
                for i in 1:size(U)[1]
                    for j in 1:size(U)[2]
                        proj_list = get_projector_list(num_action, [i,j])
                        cu += U[i,j], "Proj11", control_qubits[1], "Proj11", control_qubits[2], "Proj11", control_qubits[3],
                        proj_list[1], action_qubits[1], proj_list[2], action_qubits[2], proj_list[3], action_qubits[3],
                        proj_list[4], action_qubits[4], proj_list[5], action_qubits[5], proj_list[6], action_qubits[6],
                        proj_list[7], action_qubits[7]
                    end
                end

            elseif num_action == 8
                for i in 1:size(U)[1]
                    for j in 1:size(U)[2]
                        proj_list = get_projector_list(num_action, [i,j])
                        cu += U[i,j], "Proj11", control_qubits[1], "Proj11", control_qubits[2], "Proj11", control_qubits[3],
                        proj_list[1], action_qubits[1], proj_list[2], action_qubits[2], proj_list[3], action_qubits[3],
                        proj_list[4], action_qubits[4], proj_list[5], action_qubits[5], proj_list[6], action_qubits[6],
                        proj_list[7], action_qubits[7], proj_list[8], action_qubits[8]
                    end
                end

            elseif num_action > 8
                error("More than 8 action qubits currently not supported.")
            end

        #######################
        elseif num_control == 4
        #######################

            if num_action == 1
                for i in 1:size(U)[1]
                    for j in 1:size(U)[2]
                        proj_list = get_projector_list(num_action, [i,j])
                        cu += U[i,j], "Proj11", control_qubits[1], "Proj11", control_qubits[2], "Proj11", control_qubits[3],
                        "Proj11", control_qubits[4], proj_list[1], action_qubits[1]
                    end
                end

            elseif num_action == 2
                for i in 1:size(U)[1]
                    for j in 1:size(U)[2]
                        proj_list = get_projector_list(num_action, [i,j])
                        cu += U[i,j], "Proj11", control_qubits[1], "Proj11", control_qubits[2], "Proj11", control_qubits[3],
                        "Proj11", control_qubits[4], proj_list[1], action_qubits[1], proj_list[2], action_qubits[2]
                    end
                end

            elseif num_action == 3
                for i in 1:size(U)[1]
                    for j in 1:size(U)[2]
                        proj_list = get_projector_list(num_action, [i,j])
                        cu += U[i,j], "Proj11", control_qubits[1], "Proj11", control_qubits[2], "Proj11", control_qubits[3],
                        "Proj11", control_qubits[4], proj_list[1], action_qubits[1], proj_list[2], action_qubits[2],
                        proj_list[3], action_qubits[3]
                    end
                end

            elseif num_action == 4
                for i in 1:size(U)[1]
                    for j in 1:size(U)[2]
                        proj_list = get_projector_list(num_action, [i,j])
                        cu += U[i,j], "Proj11", control_qubits[1], "Proj11", control_qubits[2], "Proj11", control_qubits[3],
                        "Proj11", control_qubits[4], proj_list[1], action_qubits[1], proj_list[2], action_qubits[2],
                        proj_list[3], action_qubits[3], proj_list[4], action_qubits[4]
                    end
                end

            elseif num_action == 5
                for i in 1:size(U)[1]
                    for j in 1:size(U)[2]
                        proj_list = get_projector_list(num_action, [i,j])
                        cu += U[i,j], "Proj11", control_qubits[1], "Proj11", control_qubits[2], "Proj11", control_qubits[3],
                        "Proj11", control_qubits[4], proj_list[1], action_qubits[1], proj_list[2], action_qubits[2],
                        proj_list[3], action_qubits[3], proj_list[4], action_qubits[4], proj_list[5], action_qubits[5]
                    end
                end

            elseif num_action == 6
                for i in 1:size(U)[1]
                    for j in 1:size(U)[2]
                        proj_list = get_projector_list(num_action, [i,j])
                        cu += U[i,j], "Proj11", control_qubits[1], "Proj11", control_qubits[2], "Proj11", control_qubits[3],
                        "Proj11", control_qubits[4], proj_list[1], action_qubits[1], proj_list[2], action_qubits[2],
                        proj_list[3], action_qubits[3], proj_list[4], action_qubits[4], proj_list[5], action_qubits[5],
                        proj_list[6], action_qubits[6]
                    end
                end

            elseif num_action == 7
                for i in 1:size(U)[1]
                    for j in 1:size(U)[2]
                        proj_list = get_projector_list(num_action, [i,j])
                        cu += U[i,j], "Proj11", control_qubits[1], "Proj11", control_qubits[2], "Proj11", control_qubits[3],
                        "Proj11", control_qubits[4], proj_list[1], action_qubits[1], proj_list[2], action_qubits[2],
                        proj_list[3], action_qubits[3], proj_list[4], action_qubits[4], proj_list[5], action_qubits[5],
                        proj_list[6], action_qubits[6], proj_list[7], action_qubits[7]
                    end
                end

            elseif num_action == 8
                for i in 1:size(U)[1]
                    for j in 1:size(U)[2]
                        proj_list = get_projector_list(num_action, [i,j])
                        cu += U[i,j], "Proj11", control_qubits[1], "Proj11", control_qubits[2], "Proj11", control_qubits[3],
                        "Proj11", control_qubits[4], proj_list[1], action_qubits[1], proj_list[2], action_qubits[2],
                        proj_list[3], action_qubits[3], proj_list[4], action_qubits[4], proj_list[5], action_qubits[5],
                        proj_list[6], action_qubits[6], proj_list[7], action_qubits[7], proj_list[8], action_qubits[8]
                    end
                end

            elseif num_action > 8
                error("More than 8 action qubits currently not supported.")
            end

        #######################
        elseif num_control == 5
        #######################

            if num_action == 1
                for i in 1:size(U)[1]
                    for j in 1:size(U)[2]
                        proj_list = get_projector_list(num_action, [i,j])
                        cu += U[i,j], "Proj11", control_qubits[1], "Proj11", control_qubits[2], "Proj11", control_qubits[3],
                        "Proj11", control_qubits[4], "Proj11", control_qubits[5], proj_list[1], action_qubits[1]
                    end
                end

            elseif num_action == 2
                for i in 1:size(U)[1]
                    for j in 1:size(U)[2]
                        proj_list = get_projector_list(num_action, [i,j])
                        cu += U[i,j], "Proj11", control_qubits[1], "Proj11", control_qubits[2], "Proj11", control_qubits[3],
                        "Proj11", control_qubits[4], "Proj11", control_qubits[5], proj_list[1], action_qubits[1],
                        proj_list[2], action_qubits[2]
                    end
                end

            elseif num_action == 3
                for i in 1:size(U)[1]
                    for j in 1:size(U)[2]
                        proj_list = get_projector_list(num_action, [i,j])
                        cu += U[i,j], "Proj11", control_qubits[1], "Proj11", control_qubits[2], "Proj11", control_qubits[3],
                        "Proj11", control_qubits[4], "Proj11", control_qubits[5], proj_list[1], action_qubits[1],
                        proj_list[2], action_qubits[2], proj_list[3], action_qubits[3]
                    end
                end

            elseif num_action == 4
                for i in 1:size(U)[1]
                    for j in 1:size(U)[2]
                        proj_list = get_projector_list(num_action, [i,j])
                        cu += U[i,j], "Proj11", control_qubits[1], "Proj11", control_qubits[2], "Proj11", control_qubits[3],
                        "Proj11", control_qubits[4], "Proj11", control_qubits[5], proj_list[1], action_qubits[1],
                        proj_list[2], action_qubits[2], proj_list[3], action_qubits[3], proj_list[4], action_qubits[4]
                    end
                end

            elseif num_action == 5
                for i in 1:size(U)[1]
                    for j in 1:size(U)[2]
                        proj_list = get_projector_list(num_action, [i,j])
                        cu += U[i,j], "Proj11", control_qubits[1], "Proj11", control_qubits[2], "Proj11", control_qubits[3],
                        "Proj11", control_qubits[4], "Proj11", control_qubits[5], proj_list[1], action_qubits[1],
                        proj_list[2], action_qubits[2], proj_list[3], action_qubits[3], proj_list[4], action_qubits[4],
                        proj_list[5], action_qubits[5]
                    end
                end

            elseif num_action == 6
                for i in 1:size(U)[1]
                    for j in 1:size(U)[2]
                        proj_list = get_projector_list(num_action, [i,j])
                        cu += U[i,j], "Proj11", control_qubits[1], "Proj11", control_qubits[2], "Proj11", control_qubits[3],
                        "Proj11", control_qubits[4], "Proj11", control_qubits[5], proj_list[1], action_qubits[1],
                        proj_list[2], action_qubits[2], proj_list[3], action_qubits[3], proj_list[4], action_qubits[4],
                        proj_list[5], action_qubits[5], proj_list[6], action_qubits[6]
                    end
                end

            elseif num_action == 7
                for i in 1:size(U)[1]
                    for j in 1:size(U)[2]
                        proj_list = get_projector_list(num_action, [i,j])
                        cu += U[i,j], "Proj11", control_qubits[1], "Proj11", control_qubits[2], "Proj11", control_qubits[3],
                        "Proj11", control_qubits[4], "Proj11", control_qubits[5], proj_list[1], action_qubits[1],
                        proj_list[2], action_qubits[2], proj_list[3], action_qubits[3], proj_list[4], action_qubits[4],
                        proj_list[5], action_qubits[5], proj_list[6], action_qubits[6], proj_list[7], action_qubits[7]
                    end
                end

            elseif num_action == 8
                for i in 1:size(U)[1]
                    for j in 1:size(U)[2]
                        proj_list = get_projector_list(num_action, [i,j])
                        cu += U[i,j], "Proj11", control_qubits[1], "Proj11", control_qubits[2], "Proj11", control_qubits[3],
                        "Proj11", control_qubits[4], "Proj11", control_qubits[5], proj_list[1], action_qubits[1],
                        proj_list[2], action_qubits[2], proj_list[3], action_qubits[3], proj_list[4], action_qubits[4],
                        proj_list[5], action_qubits[5], proj_list[6], action_qubits[6], proj_list[7], action_qubits[7],
                        proj_list[8], action_qubits[8]
                    end
                end

            elseif num_action > 8
                error("More than 8 action qubits currently not supported.")
            end

        ######################
        elseif num_control > 5
            error("More than 5 control qubits currently not supported.")
        end
        ###

        # get sites
        sites = qc.IndexSet

        # turn into MPO
        gate = MPO(cu, sites)

        # apply to state vector, unprime
        qc.StateVector = contract(gate, qc.StateVector, maxdim=qc.MaxBondDim)#, method=qc.ContractionMethod)
        noprime!(qc.StateVector)

        # update representing matrix of quantum circuit
        if update_rep
            update_representation_arbitrary_CU!(qc, control_qubits, action_qubits, num)
        end

        # update circuit depth and bond dimension
        qc.CircuitDepth += 1
        push!(qc.BondDim, maxlinkdim(qc.StateVector))
    end
end




##########
## Test
##########

#num_qubits = 20
#control_qubits = [1, 5, 20]
#action_qubits = [7, 10, 11]
#N = length(action_qubits)
##index = [1, 3]
#for i in 1:2^N
#    for j in 1:2^N
#        #get_projector_list(N, [i,j])
#        println(i, j)
#        #println(get_projector_list_full(num_qubits, control_qubits, action_qubits, [i,j]))
#        println(get_projector_list(length(action_qubits), [i,j]))
#    end
#end
#
#c_proj = get_control_space_projectors(length(control_qubits))
#println(c_proj[1][1])


# set constants
#N = 6
#maxdim = 32
#N_meas = 100
##backend = "ED_Julia"
#backend = "MPS_ITensor"
#contmethod = "naive"
#random = false
#lintop = false
#randombond = 2
#
#qc = initialise_qcircuit(N, lintop, backend, maxdim, contmethod,
#random, randombond)
#
#PauliX!(qc, [1, 2, 6])
#
#U = [1. 0. 0. 0;
#     0. 1. 0. 0.;
#     0. 0. 0. 1.
#     0. 0. 1. 0.]
#
#X = [0. 1.;
#     1. 0.]
#
##CU_2site!(qc, U, [6, 1, 2])
#CU_general!(qc, U, [1, 6], [2, 3])

#draw(qc)

#eps = 0.01
#sample_measurement(qc, [1, 2, 3, 4, 5, 6], N_meas, eps, true, "ITensor", true)

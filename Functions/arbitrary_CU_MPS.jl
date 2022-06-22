
#########################################
## Arbitrary controlled unitary operators
#########################################


#####################
# Auxiliary functions
#####################


""" Function to update the representation of a quantum circuit qc for an
arbitrary controlled gate. The gate is generically represented by "U". """
function update_representation_arbitrary_CU!(qc::QC_IT_MPS, control_qubits, action_qubits)

    for i in 1:qc.NumQubits
        if i ∈ control_qubits
            push!(qc.Representation[i], 4) # control qubit
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
space for a given number of num_control control qubits. """
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


""" Function to implement a controlled unitary operator (on adjacent
sites). """
function CU_general!(qc::QC_IT_MPS, U, control_qubits, action_qubits, update_rep=true)

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

        #println("finished constructing MPO")

        # get sites
        sites = qc.IndexSet

        # turn into MPO
        gate = MPO(cu, sites)

        # apply to state vector, unprime
        qc.StateVector = contract(gate, qc.StateVector, maxdim=qc.MaxBondDim)#, method=qc.ContractionMethod)
        noprime!(qc.StateVector)

        # update representing matrix of quantum circuit
        if update_rep
            #update_representation_three_site!(qc, pos, 13, true, true)
            update_representation_arbitrary_CU!(qc, control_qubits, action_qubits)
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

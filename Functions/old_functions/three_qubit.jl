
#############################################################
## Three-qubit gates for Quantum Simulator - ED Julia version
#############################################################


####################
# Auxiliary functons
####################


"""  """
function get_controlled_gate_two_site(qc::QC, O::Matrix{ComplexF64}, pos::Array{Int64, 1})

    if pos[3]-pos[2] != 1
        error("Wrong ordering of positions in operator")
    end

    # get identity
    s = Index(2, "QCircuit")
    E = array(op("E", s))
    proj0 = array(op("Proj00", s)) # |0><0|
    proj1 = array(op("Proj11", s)) # |1><1|


    # control qubit above action qubits
    if pos[1] < pos[2]

        # control qubit in first line
        if pos[1] == 1
            gate0 = proj0
            gate1 = proj1
            for i in 2:pos[2]-1
                gate0 = kron(gate0, E)
                gate1 = kron(gate1, E)
            end
            gate0 = kron(gate0, E)
            gate0 = kron(gate0, E)
            gate1 = kron(gate1, O)
            for i in pos[3]+1:qc.NumQubits
                gate0 = kron(gate0, E)
                gate1 = kron(gate1, E)
            end
            gate = gate0 + gate1

        # control qubit not in first line
        else
            gate0 = E
            gate1 = E
            for i in 2:pos[1]-1
                gate0 = kron(gate0, E)
                gate1 = kron(gate1, E)
            end
            gate0 = kron(gate0, proj0)
            gate1 = kron(gate1, proj1)
            for i in pos[1]+1:pos[2]-1
                gate0 = kron(gate0, E)
                gate1 = kron(gate1, E)
            end
            gate0 = kron(gate0, E)
            gate0 = kron(gate0, E)
            gate1 = kron(gate1, O)
            for i in pos[3]+1:qc.NumQubits
                gate0 = kron(gate0, E)
                gate1 = kron(gate1, E)
            end
            gate = gate0 + gate1
        end

    # control qubit below action qubits
    else

        # action qubits on first and second line
        if pos[2] == 1
            gate0 = kron(E, E)
            gate1 = O
            for i in 3:pos[1]-1
                gate0 = kron(gate0, E)
                gate1 = kron(gate1, E)
            end
            gate0 = kron(gate0, proj0)
            gate1 = kron(gate1, proj1)
            for i in pos[1]+1:qc.NumQubits
                gate0 = kron(gate0, E)
                gate1 = kron(gate1, E)
            end
            gate = gate0 + gate1

        # action qubits not on first and second line
        else
            gate0 = E
            gate1 = E
            for i in 2:pos[2]-1
                gate0 = kron(gate0, E)
                gate1 = kron(gate1, E)
            end
            gate0 = kron(gate0, E)
            gate0 = kron(gate0, E)
            gate1 = kron(gate1, O)
            for i in pos[3]+1:pos[1]-1
                gate0 = kron(gate0, E)
                gate1 = kron(gate1, E)
            end
            gate0 = kron(gate0, proj0)
            gate1 = kron(gate1, proj1)
            for i in pos[1]+1:qc.NumQubits
                gate0 = kron(gate0, E)
                gate1 = kron(gate1, E)
            end
            gate = gate0 + gate1
        end
    end
    return gate
end


""" Auxiliary function to generate an operator O representing a three-site gate
applied on arbitrary adjacent positions given in the array pos as a Kronecker
product. """
function get_gate_three_site(qc::QC, O::Matrix{ComplexF64}, pos::Array{Int64, 1})

    # get identity
    s = Index(2, "QCircuit")
    E = array(op("E", s))

    # construct whole gate
    if pos[1] == 1
        gate = O
        for i in 4:qc.NumQubits
            gate = kron(gate, E)
        end
    else
        gate = E
        for i in 2:pos[1]-1
            gate = kron(gate, E)
        end
        gate = kron(gate, O)
        for i in (pos[3]+1):qc.NumQubits
            gate = kron(gate, E)
        end
    end

    return gate
end


""" Function to build the general matrix representation of a doubly
controlled operator O. The array pos specifies in which positions the
control qubits (pos[1], pos[2]) and the operator (pos[3]) are applied
to the qubit register. Arbitrary arrangements are supported. """
function get_doubly_controlled_gate(qc::QC, O::Matrix{ComplexF64}, pos::Array{Int64, 1})

    N = qc.NumQubits
    s = Index(2, "QCircuit")
    E = array(op("E", s))
    proj0 = array(op("Proj00", s)) # |0><0|
    proj1 = array(op("Proj11", s)) # |1><1|

    # action qubit below both control bits (e.g. pos = [1, 3, 5])
    if (pos[3] > pos[1]) && (pos[3] > pos[2])

        # first control bit in first position 1
        if pos[1] == 1

            # initialise as projectors
            gate_tmp00 = proj0
            gate_tmp01 = proj0
            gate_tmp10 = proj1
            gate_tmp11 = proj1

            # fill with identities
            for i in 2:pos[2]-1
                gate_tmp00 = kron(gate_tmp00, E)
                gate_tmp01 = kron(gate_tmp01, E)
                gate_tmp10 = kron(gate_tmp10, E)
                gate_tmp11 = kron(gate_tmp11, E)
            end

            # projection space
            gate_tmp00 = kron(gate_tmp00, proj0)
            gate_tmp01 = kron(gate_tmp01, proj1)
            gate_tmp10 = kron(gate_tmp10, proj0)
            gate_tmp11 = kron(gate_tmp11, proj1)

            # fill with identities
            for i in pos[2]+1:pos[3]-1
                gate_tmp00 = kron(gate_tmp00, E)
                gate_tmp01 = kron(gate_tmp01, E)
                gate_tmp10 = kron(gate_tmp10, E)
                gate_tmp11 = kron(gate_tmp11, E)
            end

            # apply operator in action space
            gate_tmp00 = kron(gate_tmp00, E)
            gate_tmp01 = kron(gate_tmp01, E)
            gate_tmp10 = kron(gate_tmp10, E)
            gate_tmp11 = kron(gate_tmp11, O)

            # fill with identities
            for i in pos[3]+1:N
                gate_tmp00 = kron(gate_tmp00, E)
                gate_tmp01 = kron(gate_tmp01, E)
                gate_tmp10 = kron(gate_tmp10, E)
                gate_tmp11 = kron(gate_tmp11, E)
            end

        # first control bit NOT in first position 1
        else

            # initialise as identities
            gate_tmp00 = E
            gate_tmp01 = E
            gate_tmp10 = E
            gate_tmp11 = E

            # fill with identities
            for i in 2:pos[1]-1
                gate_tmp00 = kron(gate_tmp00, E)
                gate_tmp01 = kron(gate_tmp01, E)
                gate_tmp10 = kron(gate_tmp10, E)
                gate_tmp11 = kron(gate_tmp11, E)
            end

            # projection space
            gate_tmp00 = kron(gate_tmp00, proj0)
            gate_tmp01 = kron(gate_tmp01, proj0)
            gate_tmp10 = kron(gate_tmp10, proj1)
            gate_tmp11 = kron(gate_tmp11, proj1)

            # fill with identities
            for i in pos[1]+1:pos[2]-1
                gate_tmp00 = kron(gate_tmp00, E)
                gate_tmp01 = kron(gate_tmp01, E)
                gate_tmp10 = kron(gate_tmp10, E)
                gate_tmp11 = kron(gate_tmp11, E)
            end

            # projection space
            gate_tmp00 = kron(gate_tmp00, proj0)
            gate_tmp01 = kron(gate_tmp01, proj1)
            gate_tmp10 = kron(gate_tmp10, proj0)
            gate_tmp11 = kron(gate_tmp11, proj1)

            # fill with identities
            for i in pos[2]+1:pos[3]-1
                gate_tmp00 = kron(gate_tmp00, E)
                gate_tmp01 = kron(gate_tmp01, E)
                gate_tmp10 = kron(gate_tmp10, E)
                gate_tmp11 = kron(gate_tmp11, E)
            end

            # apply operator in action space
            gate_tmp00 = kron(gate_tmp00, E)
            gate_tmp01 = kron(gate_tmp01, E)
            gate_tmp10 = kron(gate_tmp10, E)
            gate_tmp11 = kron(gate_tmp11, O)

            # fill with identities
            for i in pos[3]+1:N
                gate_tmp00 = kron(gate_tmp00, E)
                gate_tmp01 = kron(gate_tmp01, E)
                gate_tmp10 = kron(gate_tmp10, E)
                gate_tmp11 = kron(gate_tmp11, E)
            end
        end

    # action qubit between control bit (e.g. pos = [1, 5, 3])
    elseif (pos[3] > pos[1]) && (pos[3] < pos[2])

        # first control bit in first position 1
        if pos[1] == 1

            # initialise as projectors
            gate_tmp00 = proj0
            gate_tmp01 = proj0
            gate_tmp10 = proj1
            gate_tmp11 = proj1

            # fill with identities
            for i in 2:pos[3]-1
                gate_tmp00 = kron(gate_tmp00, E)
                gate_tmp01 = kron(gate_tmp01, E)
                gate_tmp10 = kron(gate_tmp10, E)
                gate_tmp11 = kron(gate_tmp11, E)
            end

            # action space
            gate_tmp00 = kron(gate_tmp00, E)
            gate_tmp01 = kron(gate_tmp01, E)
            gate_tmp10 = kron(gate_tmp10, E)
            gate_tmp11 = kron(gate_tmp11, O)


            # fill with identities
            for i in pos[3]+1:pos[2]-1
                gate_tmp00 = kron(gate_tmp00, E)
                gate_tmp01 = kron(gate_tmp01, E)
                gate_tmp10 = kron(gate_tmp10, E)
                gate_tmp11 = kron(gate_tmp11, E)
            end

            # projection space
            gate_tmp00 = kron(gate_tmp00, proj0)
            gate_tmp01 = kron(gate_tmp01, proj1)
            gate_tmp10 = kron(gate_tmp10, proj0)
            gate_tmp11 = kron(gate_tmp11, proj1)

            # fill with identities
            for i in pos[2]+1:N
                gate_tmp00 = kron(gate_tmp00, E)
                gate_tmp01 = kron(gate_tmp01, E)
                gate_tmp10 = kron(gate_tmp10, E)
                gate_tmp11 = kron(gate_tmp11, E)
            end

        # first control qubit NOT in position 1
        else

            # initialise as identities
            gate_tmp00 = E
            gate_tmp01 = E
            gate_tmp10 = E
            gate_tmp11 = E

            # fill with identities
            for i in 2:pos[1]-1
                gate_tmp00 = kron(gate_tmp00, E)
                gate_tmp01 = kron(gate_tmp01, E)
                gate_tmp10 = kron(gate_tmp10, E)
                gate_tmp11 = kron(gate_tmp11, E)
            end

            # projection space
            gate_tmp00 = kron(gate_tmp00, proj0)
            gate_tmp01 = kron(gate_tmp01, proj0)
            gate_tmp10 = kron(gate_tmp10, proj1)
            gate_tmp11 = kron(gate_tmp11, proj1)

            # fill with identities
            for i in pos[1]+1:pos[3]-1
                gate_tmp00 = kron(gate_tmp00, E)
                gate_tmp01 = kron(gate_tmp01, E)
                gate_tmp10 = kron(gate_tmp10, E)
                gate_tmp11 = kron(gate_tmp11, E)
            end

            # action space
            gate_tmp00 = kron(gate_tmp00, E)
            gate_tmp01 = kron(gate_tmp01, E)
            gate_tmp10 = kron(gate_tmp10, E)
            gate_tmp11 = kron(gate_tmp11, O)

            # fill with identities
            for i in pos[3]+1:pos[2]-1
                gate_tmp00 = kron(gate_tmp00, E)
                gate_tmp01 = kron(gate_tmp01, E)
                gate_tmp10 = kron(gate_tmp10, E)
                gate_tmp11 = kron(gate_tmp11, E)
            end

            # projection space
            gate_tmp00 = kron(gate_tmp00, proj0)
            gate_tmp01 = kron(gate_tmp01, proj1)
            gate_tmp10 = kron(gate_tmp10, proj0)
            gate_tmp11 = kron(gate_tmp11, proj1)

            # fill with identities
            for i in pos[2]+1:N
                gate_tmp00 = kron(gate_tmp00, E)
                gate_tmp01 = kron(gate_tmp01, E)
                gate_tmp10 = kron(gate_tmp10, E)
                gate_tmp11 = kron(gate_tmp11, E)
            end
        end

    # action qubit above both control bits (e.g. pos = [3, 5, 1])
    else

        # action bit in first position 1
        if pos[3] == 1

            # initialise in action space
            gate_tmp00 = E
            gate_tmp01 = E
            gate_tmp10 = E
            gate_tmp11 = O

            # fill with identities
            for i in 2:pos[1]-1
                gate_tmp00 = kron(gate_tmp00, E)
                gate_tmp01 = kron(gate_tmp01, E)
                gate_tmp10 = kron(gate_tmp10, E)
                gate_tmp11 = kron(gate_tmp11, E)
            end

            # projection space
            gate_tmp00 = kron(gate_tmp00, proj0)
            gate_tmp01 = kron(gate_tmp01, proj0)
            gate_tmp10 = kron(gate_tmp10, proj1)
            gate_tmp11 = kron(gate_tmp11, proj1)

            # fill with identities
            for i in pos[1]+1:pos[2]-1
                gate_tmp00 = kron(gate_tmp00, E)
                gate_tmp01 = kron(gate_tmp01, E)
                gate_tmp10 = kron(gate_tmp10, E)
                gate_tmp11 = kron(gate_tmp11, E)
            end

            # projection space
            gate_tmp00 = kron(gate_tmp00, proj0)
            gate_tmp01 = kron(gate_tmp01, proj1)
            gate_tmp10 = kron(gate_tmp10, proj0)
            gate_tmp11 = kron(gate_tmp11, proj1)

            # fill with identities
            for i in pos[2]+1:N
                gate_tmp00 = kron(gate_tmp00, E)
                gate_tmp01 = kron(gate_tmp01, E)
                gate_tmp10 = kron(gate_tmp10, E)
                gate_tmp11 = kron(gate_tmp11, E)
            end

        # action bit NOT in position 1
        else

            # initialise as identities
            gate_tmp00 = E
            gate_tmp01 = E
            gate_tmp10 = E
            gate_tmp11 = E

            # fill with identities
            for i in 2:pos[3]-1
                gate_tmp00 = kron(gate_tmp00, E)
                gate_tmp01 = kron(gate_tmp01, E)
                gate_tmp10 = kron(gate_tmp10, E)
                gate_tmp11 = kron(gate_tmp11, E)
            end

            # action space
            gate_tmp00 = kron(gate_tmp00, E)
            gate_tmp01 = kron(gate_tmp01, E)
            gate_tmp10 = kron(gate_tmp10, E)
            gate_tmp11 = kron(gate_tmp11, O)

            # fill with identities
            for i in pos[3]+1:pos[1]-1
                gate_tmp00 = kron(gate_tmp00, E)
                gate_tmp01 = kron(gate_tmp01, E)
                gate_tmp10 = kron(gate_tmp10, E)
                gate_tmp11 = kron(gate_tmp11, E)
            end

            # projection space
            gate_tmp00 = kron(gate_tmp00, proj0)
            gate_tmp01 = kron(gate_tmp01, proj0)
            gate_tmp10 = kron(gate_tmp10, proj1)
            gate_tmp11 = kron(gate_tmp11, proj1)

            # fill with identities
            for i in pos[1]+1:pos[2]-1
                gate_tmp00 = kron(gate_tmp00, E)
                gate_tmp01 = kron(gate_tmp01, E)
                gate_tmp10 = kron(gate_tmp10, E)
                gate_tmp11 = kron(gate_tmp11, E)
            end

            # projection space
            gate_tmp00 = kron(gate_tmp00, proj0)
            gate_tmp01 = kron(gate_tmp01, proj1)
            gate_tmp10 = kron(gate_tmp10, proj0)
            gate_tmp11 = kron(gate_tmp11, proj1)

            # fill with identities
            for i in pos[2]+1:N
                gate_tmp00 = kron(gate_tmp00, E)
                gate_tmp01 = kron(gate_tmp01, E)
                gate_tmp10 = kron(gate_tmp10, E)
                gate_tmp11 = kron(gate_tmp11, E)
            end
        end
    end

    # return complete gate
    return gate_tmp00 + gate_tmp01 + gate_tmp10 + gate_tmp11
end


""" Function to update the representation of a quantum circuit qc for a
two-site controlled gate (specified by num) applied to positions given through
the array pos. Can choose whethe to update the full or the reduced representation
(or both). """
function update_representation_three_site!(qc::QC, pos::Array{Int64, 1}, num::Int64, full::Bool, reduced::Bool)

    # adjacent controlled gate
    if abs(pos[2]-pos[1]) == 1 && abs(pos[3]-pos[3]) == 1

        if full == false && reduced == true
            for i in 1:qc.NumQubits
                if i == pos[1]
                    push!(qc.Representation[i], num) # control qubit
                elseif i == pos[2]
                    push!(qc.Representation[i], num) # control qubit
                elseif i == pos[3]
                    push!(qc.Representation[i], -num) # action qubit
                else
                    push!(qc.Representation[i], 0)
                end
            end

        elseif full == true && reduced == false
            for i in 1:qc.NumQubits
                if i == pos[1]
                    push!(qc.RepresentationFull[i], num) # control qubit
                elseif i == pos[2]
                    push!(qc.RepresentationFull[i], num) # control qubit
                elseif i == pos[3]
                    push!(qc.RepresentationFull[i], -num) # action qubit
                else
                    push!(qc.RepresentationFull[i], 0)
                end
            end

        elseif full == true && reduced == true
            for i in 1:qc.NumQubits
                if i == pos[1]
                    push!(qc.Representation[i], num) # control qubit
                    push!(qc.RepresentationFull[i], num) # control qubit
                elseif i == pos[2]
                    push!(qc.Representation[i], num) # control qubit
                    push!(qc.RepresentationFull[i], num) # control qubit
                elseif i == pos[3]
                    push!(qc.Representation[i], -num) # action qubit
                    push!(qc.RepresentationFull[i], -num) # action qubit
                else
                    push!(qc.Representation[i], 0)
                    push!(qc.RepresentationFull[i], 0)
                end
            end
        end

    # non-adjacent controlled gate
    else

        if full == false && reduced == true
            for i in 1:qc.NumQubits
                if i == pos[1]
                    push!(qc.Representation[i], num) # control qubit
                elseif i == pos[2]
                    push!(qc.Representation[i], num) # control qubit
                elseif i == pos[3]
                    push!(qc.Representation[i], -num) # action qubit
                elseif minimum(pos) < i && i < maximum(pos)
                    push!(qc.Representation[i], 15) # vertical line
                else
                    push!(qc.Representation[i], 0)
                end
            end

        elseif full == true && reduced == false
            for i in 1:qc.NumQubits
                if i == pos[1]
                    push!(qc.RepresentationFull[i], num) # control qubit
                elseif i == pos[2]
                    push!(qc.RepresentationFull[i], num) # control qubit
                elseif i == pos[3]
                    push!(qc.RepresentationFull[i], -num) # action qubit
                elseif minimum(pos) < i && i < maximum(pos)
                    push!(qc.RepresentationFull[i], 15) # vertical line
                else
                    push!(qc.RepresentationFull[i], 0)
                end
            end

        elseif full == true && reduced == true
            for i in 1:qc.NumQubits
                if i == pos[1]
                    push!(qc.Representation[i], num) # control qubit
                    push!(qc.RepresentationFull[i], num) # control qubit
                elseif i == pos[2]
                    push!(qc.Representation[i], num) # control qubit
                    push!(qc.RepresentationFull[i], num) # control qubit
                elseif i == pos[3]
                    push!(qc.Representation[i], -num) # action qubit
                    push!(qc.RepresentationFull[i], -num) # action qubit
                elseif minimum(pos) < i && i < maximum(pos)
                    push!(qc.Representation[i], 15) # vertical line
                    push!(qc.RepresentationFull[i], 15) # vertical line
                else
                    push!(qc.Representation[i], 0)
                    push!(qc.RepresentationFull[i], 0)
                end
            end
        end
    end
end


##############################
# Three-qubit gates (adjacent)
##############################


""" Function to apply the Toffoli gate to three adjacent qubits in positions
pos = [pos[1], pos[2], pos[3]] = [pos[1], pos[1]+1, pos[1]+2]. The "upper"
qubits (pos[1], pos[2]) are the control qubits, the "lower" qubit (pos[3])
the one that is acted on. This function serves mostly as an auxiliary function
in the more general toffoli! function, see below. """
function toffoli3!(qc::QC, pos::Array{Int64, 1})

   # check correct format of indices and size of circuit
   if qc.NumQubits < 3
       error("Trying to apply a 3-qubit gate to less than 3 qubits!")
   elseif length(pos) != 3
       error("Incorrect specification of indices (need 3 positions).")
   elseif pos[3]-pos[2] != 1 && pos[2]-pos[1] != 1
       error("Incorrect specification of indices (need adjecent positions).")
   end

   # get matrix representing TOFFOLI gate
   toffoli = Complex.(Matrix{Float64}(I, 8, 8))
   toffoli[7:8, 7:8] = Complex.([0. 1.; 1. 0.])

   # get gate
   gate = get_gate_three_site(qc, toffoli, pos)

   # update state vector
   qc.StateVector .= gate*qc.StateVector

   # update representing matrix of quantum circuit
   update_representation_three_site!(qc, pos, 14, true, false)

   # update circuit depth
   qc.CircuitDepth += 1
end


##################################
# three-qubit gates (not adjacent)
##################################


""" Function to implement a controlled 2-site operator (on two adjacent
sites). """
function CU_2site!(qc::QC, U, pos::Array{Int64, 1}, update_rep=true)

    # check correct format of indices and size of circuit
    if qc.NumQubits < 3
        error("Trying to apply a 3-qubit gate to less than 3 qubits!")
    elseif length(pos) != 3
        error("Incorrect specification of indices (need 3 positions).")
    elseif pos[3] - pos[2] != 1
        error("Incorrect specification of indices (last two indices must be adjacent).")
    end

    # apply general CU_2site via SWAP gates for linear topology
    if qc.LinearTopology == true

        println("Linear topology currently not implemented!")

    # master topology
    else

        # get gate
        gate = get_controlled_gate_two_site(qc, U, pos)

        # update state vector
        qc.StateVector .= gate*qc.StateVector

        # update representing matrix of quantum circuit
        update_representation_three_site!(qc, pos, 14, true, true)

        # update circuit depth
        qc.CircuitDepth += 1

    end
end


""" Function to apply the Toffoli gate in an arbitrary configuration
of the input positions. For a linear topology, the function permutes
the input positions such that they are next to each other, and the
toffoli3! function is applied subsequently. For a master topology, the
matrix representing the operation is built from scratch. """
function toffoli!(qc::QC, pos::Array{Int64, 1})

    # check correct format of indices and size of circuit
    if qc.NumQubits < 3
        error("Trying to apply a 3-qubit gate to less than 3 qubits!")
    elseif length(pos) != 3
        error("Incorrect specification of indices (need 3 positions).")
    elseif pos[1] >= pos[2]
        error("Incorrect specification of indices (check control-qubits in ascending order).")
    end

    # apply general CNOT via SWAP gates for linear topology
    if qc.LinearTopology == true

        # reference array of qubits
        qubits = [i for i in 1:qc.NumQubits]

        # get permutation of qubits from position indication
        perm = [i for i in 1:(minimum(pos)-1)]
        append!(perm, pos)
        append!(perm, setdiff(qubits, perm))

        # get non-trivial cycles that the permutation above defines
        cycs = collect(cycles(Perm(perm)))
        for cyc in cycs
            if length(cyc) == 1
                deleteat!(cycs, findall(x->x==cyc, cycs))
            end
        end

        # convert into sequence of 2-swaps on neighbouring qubits
        swaps = []
        for cyc in cycs
            swap_tmp = get_2_perms_from_cycle(cyc)
            for sw in swap_tmp
                push!(swaps, sw)
            end
        end

        # swap qubits, apply Toffoli, swap back
        for sw in swaps
            swap2!(qc::QC, sw)
        end
        toffoli3!(qc::QC, [minimum(pos), minimum(pos)+1, minimum(pos)+2])
        for sw in reverse(swaps)
            swap2!(qc::QC, sw)
        end

        # update representing matrix of quantum circuit
        update_representation_three_site!(qc, pos, 14, false, true)

    # apply general Toffoli gate via direct construction of matrix for master topology
    else

        # get X gate
        s = Index(2, "QCircuit")
        X = array(op("X", s))

        # get gate
        gate = get_doubly_controlled_gate(qc, X, pos)

        # update state vector
        qc.StateVector .= gate*qc.StateVector

        # update representing matrix of quantum circuit
        update_representation_three_site!(qc, pos, 14, true, true)

        # update circuit depth
        qc.CircuitDepth += 1
    end
end

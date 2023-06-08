
###########################################################
## Two-qubit gates for Quantum Simulator - ED Julia version
###########################################################


#####################
# Auxiliary functions
#####################


""" Auxiliary function to generate an operator O representing a two-site gate
applied on arbitrary adjacent positions given in the array pos as a Kronecker
product. E.g. specifying X on pos[2, 3] for a 4-qubit system constructs the operator
E ⊗ CNOT ⊗ E. """
function get_gate_two_site(qc::QC, O::Matrix{ComplexF64}, pos::Array{Int64, 1})

    # get identity
    s = Index(2, "QCircuit")
    E = array(op("E", s))

    # construct whole gate
    if pos[1] == 1
        gate = O

        for i in 3:qc.NumQubits
            gate = kron(gate, E)
        end
    else
        gate = E
        for i in 2:pos[1]-1
            gate = kron(gate, E)
        end
        gate = kron(gate, O)
        for i in (pos[2]+1):qc.NumQubits
            gate = kron(gate, E)
        end
    end

    return gate
end


""" Function to build the exact matrix representation applying a singly-
controlled gate O. The array pos specifies the corresponding placements
of the control- and action qubits, with pos[1] being the control qubit
and pos[2] the action qubit. The control qubit may be above or below the
action qubit. """
function get_controlled_gate(qc::QC, O::Matrix{ComplexF64}, pos::Array{Int64, 1})

    # get number of qubits, identity and projectors
    N = qc.NumQubits
    s = Index(2, "QCircuit")
    E = array(op("E", s))
    proj0 = array(op("Proj00", s)) # |0><0|
    proj1 = array(op("Proj11", s)) # |1><1|

    # control qubit "above" qubit that is acted on
    if pos[2] > pos[1]

        # qubit 1 is control bit
        if pos[1] == 1
            gate_tmp0 = proj0
            gate_tmp1 = proj1

            # fill with identities
            for i in 2:pos[2]-1
                gate_tmp0 = kron(gate_tmp0, E)
                gate_tmp1 = kron(gate_tmp1, E)
            end

            # apply operator in action space
            gate_tmp0 = kron(gate_tmp0, E)
            gate_tmp1 = kron(gate_tmp1, O)

            # fill with identities
            for i in pos[2]+1:N
                gate_tmp0 = kron(gate_tmp0, E)
                gate_tmp1 = kron(gate_tmp1, E)
            end

        # qubit 1 is NOT control bit
        else

            # initialise with identities
            gate_tmp0 = E
            gate_tmp1 = E

            # fill with identities
            for i in 2:pos[1]-1
                gate_tmp0 = kron(gate_tmp0, E)
                gate_tmp1 = kron(gate_tmp1, E)
            end

            # projection space
            gate_tmp0 = kron(gate_tmp0, proj0)
            gate_tmp1 = kron(gate_tmp1, proj1)

            # fill with identities
            for i in pos[1]+1:pos[2]-1
                gate_tmp0 = kron(gate_tmp0, E)
                gate_tmp1 = kron(gate_tmp1, E)
            end

            # apply operator in action space
            gate_tmp0 = kron(gate_tmp0, E)
            gate_tmp1 = kron(gate_tmp1, O)

            # fill with identities
            for i in pos[2]+1:N
                gate_tmp0 = kron(gate_tmp0, E)
                gate_tmp1 = kron(gate_tmp1, E)
            end
        end

    # control qubit "below" qubit that is acted on
    else

        # qubit 1 is action bit
        if pos[2] == 1
            gate_tmp0 = E
            gate_tmp1 = O

            # fill with identities
            for i in 2:pos[1]-1
                gate_tmp0 = kron(gate_tmp0, E)
                gate_tmp1 = kron(gate_tmp1, E)
            end

            # apply projectors in control space
            gate_tmp0 = kron(gate_tmp0, proj0)
            gate_tmp1 = kron(gate_tmp1, proj1)

            # fill with identities
            for i in pos[1]+1:N
                gate_tmp0 = kron(gate_tmp0, E)
                gate_tmp1 = kron(gate_tmp1, E)
            end

        # qubit 1 is NOT action bit
        else

            # initialise with identities
            gate_tmp0 = E
            gate_tmp1 = E

            # fill with identities
            for i in 2:pos[2]-1
                gate_tmp0 = kron(gate_tmp0, E)
                gate_tmp1 = kron(gate_tmp1, E)
            end

            # action space
            gate_tmp0 = kron(gate_tmp0, E)
            gate_tmp1 = kron(gate_tmp1, O)

            # fill with identities
            for i in pos[2]+1:pos[1]-1
                gate_tmp0 = kron(gate_tmp0, E)
                gate_tmp1 = kron(gate_tmp1, E)
            end

            # apply projectors in control space
            gate_tmp0 = kron(gate_tmp0, proj0)
            gate_tmp1 = kron(gate_tmp1, proj1)

            # fill with identities
            for i in pos[1]+1:N
                gate_tmp0 = kron(gate_tmp0, E)
                gate_tmp1 = kron(gate_tmp1, E)
            end
        end
    end

    # complete gate
    return gate_tmp0 + gate_tmp1
end


""" Function to update the representation of a quantum circuit qc for a
two-site controlled gate (specified by num) applied to positions given through
the array pos. Can choose whether to update the full or the reduced representation
(or both). """
function update_representation_two_site!(qc::QC, pos::Array{Int64, 1},
    num::Int64, full::Bool, reduced::Bool)

    # adjacent controlled gate
    if abs(pos[2]-pos[1]) == 1

        if full == false && reduced == true
            for i in 1:qc.NumQubits
                if i == pos[1]
                    push!(qc.Representation[i], num) # control qubit
                elseif i == pos[2]
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


""" Function to update the representation of a quantum circuit qc for a
two-site controlled gate (specified by num) applied to positions given through
the array pos. Can choose whethe to update the full or the reduced representation
(or both). """
#function update_representation_two_site!(qc::QC, pos::Array{Int64, 1}, num::Int64, full::Bool, reduced::Bool)
#
#    if full == false && reduced == true
#        for i in 1:qc.NumQubits
#            if i in pos
#                push!(qc.Representation[i], num) # control qubit
#            else
#                push!(qc.Representation[i], 0)
#            end
#        end
#
#    elseif full == true && reduced == false
#        for i in 1:qc.NumQubits
#            if i in pos
#                push!(qc.RepresentationFull[i], num) # control qubit
#            else
#                push!(qc.RepresentationFull[i], 0)
#            end
#        end
#
#    elseif full == true && reduced == true
#        for i in 1:qc.NumQubits
#            if i in pos
#                push!(qc.Representation[i], num) # control qubit
#                push!(qc.RepresentationFull[i], num) # control qubit
#            else
#                push!(qc.Representation[i], 0)
#                push!(qc.RepresentationFull[i], 0)
#            end
#        end
#    end
#end


############################
# Two-qubit gates (adjacent)
############################


""" Function to apply the CNOT-gate to two adjacent qubits in positions
pos = [pos[1], pos[2]] = [pos[1], pos[1]+1]. The "upper" qubit (pos[1])
is the control qubit, the "lower" qubit (pos[2]) the one that is acted
on. This function serves mostly as an auxiliary function in the more
general cnot! function, see below. """
function cnot2!(qc::QC, pos::Array{Int64, 1}, update_rep=true)

   if length(pos) != 2
       error("Incorrect specification of indices (need 2 positions).")
   elseif abs(pos[1]-pos[2]) != 1
       error("Incorrect specification of indices (need adjacent positions).")
   end

   # get matrix -> change to use structures from Hilbert space!
   cnot = Complex.([1. 0. 0. 0.; 0. 1. 0. 0.; 0. 0. 0. 1.; 0. 0. 1. 0.])

   # get gate
   gate = get_gate_two_site(qc, cnot, pos)

   # update state vector
   qc.StateVector .= gate*qc.StateVector

   # update representing matrix of quantum circuit
   if update_rep
       update_representation_two_site!(qc, pos, 3, true, false)
   end

   # update circuit depth
   qc.CircuitDepth += 1
end


""" Function to swap two adjacent qubits in positions pos = [pos[1], pos[2]]
= [pos[1], pos[1]+1]. """
function swap2!(qc::QC, pos::Array{Int64, 1}, update_rep=true)

   # check structure of gate indices
   if length(pos) != 2
       error("Incorrect specification of indices (need 2 positions).")
   elseif abs(pos[1]-pos[2]) != 1
       error("Incorrect specification of indices (need adjacent positions).")
   end

   # get matrices -> change to use structures from Hilbert space!
   sw = Complex.([1. 0. 0. 0.; 0. 0. 1. 0.; 0. 1. 0. 0.; 0. 0. 0. 1.])

   # get gate
   gate = get_gate_two_site(qc, sw, sort!(pos))

   # update state vector
   #println(size(gate))
   qc.StateVector .= gate*qc.StateVector

   # update representing matrix of quantum circuit
   if update_rep
      update_representation_two_site!(qc, pos, 4, true, false)
   end

   # update circuit depth
   qc.CircuitDepth += 1
end


""" Function to apply a controlled-U gate to two adjacent qubits in positions
pos = [pos[1], pos[2]] = [pos[1], pos[1]+1], where U is a unitary single-qubit
gate. The "upper" qubit (pos[1]) is the control qubit, the "lower" qubit
(pos[2]) the one that is acted on. This function serves mostly as an auxiliary
function in the more general CU! function, see below. """
function CU2!(qc::QC, U, pos::Array{Int64, 1}, update_rep=true)

   # check if always gets two indices!
   if length(pos) != 2
       error("Incorrect specification of indices (need 2 positions).")
   elseif abs(pos[1]-pos[2]) != 1
       error("Incorrect specification of indices (need adjacent positions).")
   end

   # get matrices -> change to use structures from Hilbert space!
   #E, _, _, _ = get_Pauli_matrices()
   s = Index(2, "QCircuit")
   E = array(op("E", s))
   CU = Complex.(zeros(4,4))
   CU[1:2, 1:2] = E
   CU[3:4, 3:4] = U

   # get gate
   gate = get_gate_two_site(qc, CU, pos)

   # update state vector
   qc.StateVector .= gate*qc.StateVector

   # update representing matrix of quantum circuit
   if update_rep
       update_representation_two_site!(qc, pos, 13, true, false)
   end

   # update circuit depth
   qc.CircuitDepth += 1
end


""" Function to apply a controlled-Rn gate to two adjacent qubits in positions
pos = [pos[1], pos[2]] = [pos[1], pos[1]+1], where Rn is the rotation gate as
applied in the quantum Fourier transform. The "upper" qubit (pos[1]) is the
control qubit, the "lower" qubit (pos[2]) the one that is acted on. This function
serves mostly as an auxiliary function in the more general CRn! function,
see below. """
function CRn2!(qc::QC, n::Number, pos::Array{Int64, 1}, update_rep=true)

   # check if always gets two indices!
   if length(pos) != 2
       error("Incorrect specification of indices (need 2 positions).")
   elseif abs(pos[1]-pos[2]) != 1
       error("Incorrect specification of indices (need adjacent positions).")
   end

   # get matrices -> change to use structures from Hilbert space!
   #E, _, _, _ = get_Pauli_matrices()
   s = Index(2, "QCircuit")
   E = array(op("E", s))
   CRn = Complex.(zeros(4,4))
   CRn[1:2, 1:2] = E
   CRn[3, 3] = Complex(1.)
   CRn[4, 4] = exp(sign(n)*2π*1.0im/2^abs(n))

   # get gate
   gate = get_gate_two_site(qc, CRn, pos)

   # update state vector
   qc.StateVector .= gate*qc.StateVector

   # update representing matrix of quantum circuit
   if update_rep
       update_representation_two_site!(qc, pos, 24, true, false)
   end

   # update circuit depth
   qc.CircuitDepth += 1
end


################################
# two-qubit gates (NOT adjacent)
################################


""" Function to perfom a half-swap of two qubits over an arbitrary separation,
which is given by the array pos. pos[1] is understood to be smaller than
pos[2]. The swap is obtained through an equivalent sequence of 2-swaps. Primarily
needed as an auxiliary function in the 2-qubit gates acting on non-adjacent
qubits (half-swap to bring them next to each other, then execute the gate, then
undo the half-swap). For a full swap, see the function fullSwap!. """
function swap!(qc::QC, pos::Array{Int64, 1})

    # check correct format of indices
    if length(pos) != 2
        error("Incorrect specification of indices (need 2 positions).")
    elseif pos[1] >= pos[2]
        error("Incorrect specification of indices (pos[1] should be smaller than pos[2]).")
    end

    if pos[1] == (pos[2]-1)
        swap2!(qc::QC, pos)
    else
        positions = get_2_perms_old(pos)
        for p in positions
            swap2!(qc::QC, p)
        end
    end
end


""" Function to reverse the arbitrary half-swap, performed by the swap! function.
For more information, see docstring of the swap! function. """
function unswap!(qc::QC, pos::Array{Int64, 1})

    # check correct format of indices
    if length(pos) != 2
        error("Incorrect specification of indices (need 2 positions).")
    elseif pos[1] >= pos[2]
        error("Incorrect specification of indices (pos[1] should be smaller than pos[2]).")
    end

    if pos[1] == (pos[2]-1)
        swap2!(qc::QC, pos)
    else
        positions = get_2_perms_old(pos)
        for p in reverse(positions)
            swap2!(qc::QC, p)
        end
    end
end


""" Function to apply the CNOT gate to two arbitrary qubits, defined by
positions in pos = [pos[1], pos[2]]. The first qubits in is the control
qubit, the second the qubit that is acted on. The control qubit can be
"above" or "below" the action qubit, the function implements the desired
arrangement. Depending on whether a linear or a master topology is chosen,
the CNOT gate is either constructed with auxiliary swaps or implemented directly
as a matrix. """
function cnot!(qc::QC, pos::Array{Int64, 1}, update_rep=true)

    # check correct format of indices
    if length(pos) != 2
        error("Incorrect specification of indices (need 2 positions).")
    elseif pos[1] == pos[2]
        error("Incorrect specification of indices.")
    end

    # apply general CNOT via SWAP gates for linear topology
    if qc.LinearTopology == true

        # control qubit "above" qubit that is acted on
        if pos[2] > pos[1]

            # qubits on two adjacent lines
            if pos[1] == pos[2]-1
                cnot2!(qc::QC, pos)

            # qubits not on adjacent lines
            else
                swap!(qc::QC, [pos[1], pos[2]-1])
                cnot2!(qc::QC, [pos[2]-1, pos[2]])
                unswap!(qc::QC, [pos[1], pos[2]-1])
            end

        # control qubit "below" qubit that is acted on
        else

            # qubits on two adjacent lines
            if pos[1] == (pos[2]+1)
                swap2!(qc::QC, pos)
                cnot2!(qc::QC, pos)
                swap2!(qc::QC, pos)

            # qubits on two adjacent lines
            else
                swap!(qc::QC, [pos[2], pos[1]])
                cnot2!(qc::QC, [pos[1]-1, pos[1]])
                unswap!(qc::QC, [pos[2], pos[1]])
            end
        end

        # update representing matrix of quantum circuit
        update_representation_two_site!(qc, pos, 3, false, true)

    # apply general CNOT via direct construction of matrix for master topology
    else

        # get X gate
        s = Index(2, "QCircuit")
        X = array(op("X", s))

        # get matrix representing the general CNOT gate
        gate = get_controlled_gate(qc, X, pos)

        # update state vector
        qc.StateVector .= gate*qc.StateVector

        # update representing matrix of quantum circuit
        if update_rep
            update_representation_two_site!(qc, pos, 3, true, true)
        end

        # update cicruit depth
        qc.CircuitDepth += 1

    end
end


""" Function to perform a full swap of two qubits. For a linear topology,
the swap is decomposed into the corresponding 2-swaps. In the case of a
master topology, the gate is built from composing three CNOT gates on
pos, reverse(pos) and pos. """
function fullSwap!(qc::QC, pos::Array{Int64, 1}, update_rep=true)

    # check correct format of indices
    if length(pos) != 2
        error("Incorrect specification of indices (need 2 positions).")
    elseif pos[1] >= pos[2]
        error("Incorrect specification of indices (pos[1] should be smaller than pos[2]).")
    end

    # apply full SWAP via 2-SWAP gates for linear topology
    if qc.LinearTopology == true

        # qubits on two adjacent lines
        if pos[1] == pos[2]-1
            swap2!(qc::QC, pos)

        # qubits not on adjacent lines
        else
            swap!(qc::QC, [pos[1], pos[2]])
            unswap!(qc::QC, [pos[1], pos[2]-1])
        end

        # update compressed representation
        update_representation_two_site!(qc, pos, 4, false, true)

    # for Master topology: built SWAP gate directly
    else

        # get number of qubits, identity and projectors
        N = qc.NumQubits
        s = Index(2, "QCircuit")
        E = array(op("E", s))
        proj0 = array(op("Proj00", s)) # |0><0|
        proj1 = array(op("Proj11", s)) # |1><1|
        splus = array(op("S+", s)) # |0><1|
        sminus = array(op("S-", s)) # |1><0|

        # qubit 1 is involved in SWAP
        if pos[1] == 1
            gate_tmp0 = proj0
            gate_tmp1 = splus
            gate_tmp2 = sminus
            gate_tmp3 = proj1

            # fill with identities
            for i in 2:pos[2]-1
                gate_tmp0 = kron(gate_tmp0, E)
                gate_tmp1 = kron(gate_tmp1, E)
                gate_tmp2 = kron(gate_tmp2, E)
                gate_tmp3 = kron(gate_tmp3, E)
            end

            # apply corresponding operator to second qubit involved in SWAP
            gate_tmp0 = kron(gate_tmp0, proj0)
            gate_tmp1 = kron(gate_tmp1, sminus)
            gate_tmp2 = kron(gate_tmp2, splus)
            gate_tmp3 = kron(gate_tmp3, proj1)

            # fill with identities
            for i in pos[2]+1:N
                gate_tmp0 = kron(gate_tmp0, E)
                gate_tmp1 = kron(gate_tmp1, E)
                gate_tmp2 = kron(gate_tmp2, E)
                gate_tmp3 = kron(gate_tmp3, E)
            end

        # qubit 1 is NOT involved in SWAP
        else

            # initialise with identities
            gate_tmp0 = E
            gate_tmp1 = E
            gate_tmp2 = E
            gate_tmp3 = E

            # fill with identities
            for i in 2:pos[1]-1
                gate_tmp0 = kron(gate_tmp0, E)
                gate_tmp1 = kron(gate_tmp1, E)
                gate_tmp2 = kron(gate_tmp2, E)
                gate_tmp3 = kron(gate_tmp3, E)
            end

            # first qubit involved in SWAP
            gate_tmp0 = kron(gate_tmp0, proj0)
            gate_tmp1 = kron(gate_tmp1, splus)
            gate_tmp2 = kron(gate_tmp2, sminus)
            gate_tmp3 = kron(gate_tmp3, proj1)

            # fill with identities
            for i in pos[1]+1:pos[2]-1
                gate_tmp0 = kron(gate_tmp0, E)
                gate_tmp1 = kron(gate_tmp1, E)
                gate_tmp2 = kron(gate_tmp2, E)
                gate_tmp3 = kron(gate_tmp3, E)
            end

            # second qubit involved in SWAP
            gate_tmp0 = kron(gate_tmp0, proj0)
            gate_tmp1 = kron(gate_tmp1, sminus)
            gate_tmp2 = kron(gate_tmp2, splus)
            gate_tmp3 = kron(gate_tmp3, proj1)

            # fill with identities
            for i in pos[2]+1:N
                gate_tmp0 = kron(gate_tmp0, E)
                gate_tmp1 = kron(gate_tmp1, E)
                gate_tmp2 = kron(gate_tmp2, E)
                gate_tmp3 = kron(gate_tmp3, E)
            end
        end

        # get matrix representation of gate
        gate = gate_tmp0 + gate_tmp1 + gate_tmp2 + gate_tmp3

        # update state vector
        qc.StateVector .= gate*qc.StateVector

        # update representing matrix of quantum circuit
        if update_rep
            update_representation_two_site!(qc, pos, 4, true, true)
        end

        # update cicruit depth
        qc.CircuitDepth += 1

    end
end


""" Function to apply the CRn gate to two arbitrary qubits, defined by
positions in pos = [pos[1], pos[2]]. Rn is the rotation of the QFT.
The first qubit in pos is the control qubit, the second the qubit that is
acted on. The control qubit can be "above" or "below" the action qubit,
the function implements the desired arrangement. Depending on whether a
linear or a master topology is chosen, the CRn gate is either constructed
with auxiliary swaps or implemented directly as a matrix. """
function CRn!(qc::QC, n::Number, pos::Array{Int64, 1}, update_rep=true)

    # check correct format of indices
    if length(pos) != 2
        error("Incorrect specification of indices (need 2 positions).")
    elseif pos[1] == pos[2]
        error("Incorrect specification of indices.")
    end

    # apply general CNOT via SWAP gates for linear topology
    if qc.LinearTopology == true

        # control qubit "above" qubit that is acted on
        if pos[2] > pos[1]

            # qubits on two adjacent lines
            if pos[1] == (pos[2]-1)
                CRn2!(qc::QC, n, pos)

            # qubits not on adjacent lines
            else
                swap!(qc::QC, [pos[1], pos[2]-1])
                CRn2!(qc::QC, n, [pos[2]-1, pos[2]])
                unswap!(qc::QC, [pos[1], pos[2]-1])
            end

        # control qubit "below" qubit that is acted on
        else

            # qubits on two adjacent lines
            if pos[1] == (pos[2]+1)
                swap2!(qc::QC, pos)
                CRn2!(qc::QC, n, pos)
                swap2!(qc::QC, pos)

            # qubits on two adjacent lines
            else
                swap!(qc::QC, [pos[2], pos[1]])
                CRn2!(qc::QC, n, [pos[1]-1, pos[1]])
                unswap!(qc::QC, [pos[2], pos[1]])
            end
        end

        # update representing matrix of quantum circuit
        update_representation_two_site!(qc, pos, 24, false, true)

    # apply general CRn via direct construction of matrix for master topology
    else

        # get matrix representing the general CNOT gate
        s = Index(2, "QCircuit")
        Rn = array(op("Rn", s; n=n))
        gate = get_controlled_gate(qc, Rn, pos)

        # update state vector
        qc.StateVector .= gate*qc.StateVector

        # update representing matrix of quantum circuit
        if update_rep
            update_representation_two_site!(qc, pos, 24, true, true)
        end

        # update circuit depth
        qc.CircuitDepth += 1

    end
end


""" Function to apply the CU gate to two arbitrary qubits, defined by
positions in pos = [pos[1], pos[2]]. U is a unitary single-qubit gate.
The first qubits in is the control qubit, the second the qubit that is
acted on. The control qubit can be "above" or "below" the action qubit,
the function implements the desired arrangement. Depending on whether a
linear or a master topology is chosen, the CU gate is either constructed
with auxiliary swaps or implemented directly as a matrix. """
function C_UGate!(qc::QC, U, pos::Array{Int64, 1}, update_rep=true)

    # check correct format of indices
    if length(pos) != 2
        error("Incorrect specification of indices (need 2 positions).")
    elseif pos[1] == pos[2]
        error("Incorrect specification of indices.")
    end

    # apply general CNOT via SWAP gates for linear topology
    if qc.LinearTopology == true

        # control qubit "above" qubit that is acted on
        if pos[2] > pos[1]

            # qubits on two adjacent lines
            if pos[1] == (pos[2]-1)
                CU2!(qc::QC, U, pos)

            # qubits not on adjacent lines
            else
                swap!(qc::QC, [pos[1], pos[2]-1])
                CU2!(qc::QC, U, [pos[2]-1, pos[2]])
                unswap!(qc::QC, [pos[1], pos[2]-1])
            end

        # control qubit "below" qubit that is acted on
        else

            # qubits on two adjacent lines
            if pos[1] == (pos[2]+1)
                swap2!(qc::QC, pos)
                CU2!(qc::QC, U, pos)
                swap2!(qc::QC, pos)

            # qubits on two adjacent lines
            else
                swap!(qc::QC, [pos[2], pos[1]])
                CU2!(qc::QC, U, [pos[1]-1, pos[1]])
                unswap!(qc::QC, [pos[2], pos[1]])
            end
        end

        # update representing matrix of quantum circuit
        update_representation_two_site!(qc, pos, 13, false, true)

    # apply general CU via direct construction of matrix for master topology
    else

        # get matrix representing the general CNOT gate
        gate = get_controlled_gate(qc, U, pos)

        # update state vector
        qc.StateVector .= gate*qc.StateVector

        # update representing matrix of quantum circuit
        if update_rep
            update_representation_two_site!(qc, pos, 13, true, true)
        end

        # update circuit depth
        qc.CircuitDepth += 1

    end
end


#################
# Random circuits
#################

""" Function to apply a (vertical) sequence of CNOT gates to a quantum
circuit. Use together with random unitaries to create a highly entangled
state in the fastest possible way. Can specify whether the first CNOT
gate starts on the first or on the second qubit, the others follow
automatically after that. """
function sequential_cnot!(qc::QC, start_pos::Array{Int64, 1})

    if start_pos[1] != 1 && start_pos[1] != 2
        error("Incorrect specification of starting point (pos = [1] or [2]).")
    end

    # get number of qubits and matrices
    N = qc.NumQubits
    s = Index(2, "QCircuit")
    E = array(op("E", s))
    cnot = Complex.([1. 0. 0. 0.; 0. 1. 0. 0.; 0. 0. 0. 1.; 0. 0. 1. 0.])

    # first CNOT starting on first qubit line
    if start_pos[1] == 1

        # initialise gate
        gate = cnot

        # even number of qubits
        if iseven(N)
            for i in 2:(N÷2)
                gate = kron(gate, cnot)
            end

        # odd number of qubits
        else
            for i in 2:((N-1)÷2)
                gate = kron(gate, cnot)
            end

            # need identity on last line
            gate = kron(gate, E)

        end

    # first CNOT starting on second qubit line
    elseif start_pos[1] == 2

        # initialise gate
        gate = E

        # even number of qubits
        if iseven(qc.NumQubits)
            for i in 2:(N÷2)
                gate = kron(gate, cnot)
            end

            # need identity on last line
            gate = kron(gate, E)

        # odd number of qubits
        else
            for i in 2:((N-1)÷2)+1
                gate = kron(gate, cnot)
            end
        end
    end

    # update state vector
    qc.StateVector .= gate*qc.StateVector

    # even number of qubits
    if iseven(N)

        # CNOT starting on first qubit
        if start_pos[1] == 1
            for i in 1:N
                if isodd(i)
                    push!(qc.Representation[i], 3) # control qubit
                    push!(qc.RepresentationFull[i], 3) # control qubit
                else
                    push!(qc.Representation[i], -3) # action qubit
                    push!(qc.RepresentationFull[i], -3) # action qubit
                end
            end

        # CNOT starting on second qubit
        else
            for i in 1:N
                if i == 1
                    push!(qc.Representation[i], 0)
                    push!(qc.RepresentationFull[i], 0)
                elseif iseven(i) && i != N
                    push!(qc.Representation[i], 3) # control qubit
                    push!(qc.RepresentationFull[i], 3) # control qubit
                elseif isodd(i) && i != 1
                    push!(qc.Representation[i], -3) # action qubit
                    push!(qc.RepresentationFull[i], -3) # action qubit
                elseif i == N
                    push!(qc.Representation[i], 0)
                    push!(qc.RepresentationFull[i], 0)
                end
            end
        end

    # odd number of qubits
    else

        # CNOT starting on first qubit
        if start_pos[1] == 1
            for i in 1:N
                if isodd(i) && i != N
                    push!(qc.Representation[i], 3) # control qubit
                    push!(qc.RepresentationFull[i], 3) # control qubit
                elseif iseven(i)
                    push!(qc.Representation[i], -3) # action qubit
                    push!(qc.RepresentationFull[i], -3) # action qubit
                elseif i == N
                    push!(qc.Representation[i], 0)
                    push!(qc.RepresentationFull[i], 0)
                end
            end

        # CNOT starting on second qubit
        else
            for i in 1:N
                if i == 1
                    push!(qc.Representation[i], 0)
                    push!(qc.RepresentationFull[i], 0)
                elseif iseven(i) && i != N
                    push!(qc.Representation[i], 3) # control qubit
                    push!(qc.RepresentationFull[i], 3) # control qubit
                elseif isodd(i) && i != 1
                    push!(qc.Representation[i], -3) # action qubit
                    push!(qc.RepresentationFull[i], -3) # action qubit
                elseif i == N
                    push!(qc.Representation[i], 0)
                    push!(qc.RepresentationFull[i], 0)
                end
            end
        end
    end

    # update circuit depth
    qc.CircuitDepth += 1

end


# build sequential CNOT function that can be applied to a subset
# of the whole quantum circuit



""" Function to apply a (vertical) sequence of CNOT gates to a quantum
circuit. Use together with random unitaries to create a highly entangled
state in the fastest possible way. Can specify the subregister of qubits
to which the sequential CNOT should be applied by giving the starting
position pos and the number num of qubits in the subregister. """
function sequential_cnot!(qc::QC, pos::Number, num::Number, update_rep=false)

    # get number of qubits and matrices
    N = qc.NumQubits
    s = Index(2, "QCircuit")
    E = array(op("E", s))
    cnot = Complex.([1. 0. 0. 0.; 0. 1. 0. 0.; 0. 0. 0. 1.; 0. 0. 1. 0.])

    # collect positions where CNOTs are applied
    cnot_control_pos = []
    cnot_action_pos = []

    # sequence of CNOTs starts on first qubit line
    if pos == 1

        # initialise gate
        #println("initialiased on  1, 2")
        gate = cnot
        push!(cnot_control_pos, 1)
        push!(cnot_action_pos, 2)

        # even number of qubits
        if iseven(num)

            for i in 3:2:(pos+num-1)
                #println("applying CNOT on $i, $(i+1)")
                gate = kron(gate, cnot)
                push!(cnot_control_pos, i)
                push!(cnot_action_pos, i+1)
            end

            # fill rest with identities
            for i in pos+num:N
                #println("applying E on $i")
                gate = kron(gate, E)
            end

        # odd number of qubits
        else
            for i in 3:2:(pos+num-3)
                #println("applying CNOT on $i, $(i+1)")
                gate = kron(gate, cnot)
                push!(cnot_control_pos, i)
                push!(cnot_action_pos, i+1)
            end

            # fill rest with identities
            for i in (pos+num-1):N
                #println("applying E on $i")
                gate = kron(gate, E)
            end
        end

    # sequence of CNOTs doesn't start on first qubit line
    else

        #println("initialiased on 1")

        # initialise gate
        gate = E

        # fill with identities until starting position
        for i in 2:(pos-1)
            #println("applying E on $i")
            gate = kron(gate, E)
        end

        # even number of qubits
        if iseven(num)

            for i in pos:2:(pos+num-1)
                #println("applying CNOT on $i, $(i+1)")
                gate = kron(gate, cnot)
                push!(cnot_control_pos, i)
                push!(cnot_action_pos, i+1)
            end

            # fill rest with identities
            for i in (pos+num):N
                #println("applying E on $i")
                gate = kron(gate, E)
            end

        # odd number of qubits
        else
            for i in pos:2:(pos+num-3)
                #println("applying CNOT on $i, $(i+1)")
                gate = kron(gate, cnot)
                push!(cnot_control_pos, i)
                push!(cnot_action_pos, i+1)
            end

            # fill rest with identities
            for i in (pos+num-1):N
                #println("applying E on $i")
                gate = kron(gate, E)
            end
        end
    end

    # update state vector
    qc.StateVector .= gate*qc.StateVector

    # update representation, if desired
    if update_rep
        for i in 1:N
            if i ∈ cnot_control_pos
                push!(qc.Representation[i], 3) # control qubit
                push!(qc.RepresentationFull[i], 3) # control qubit
            elseif i ∈ cnot_action_pos
                push!(qc.Representation[i], -3) # action qubit
                push!(qc.RepresentationFull[i], -3) # action qubit
            else
                push!(qc.Representation[i], 0)
                push!(qc.RepresentationFull[i], 0)
            end
        end
    end

    # update circuit depth
    qc.CircuitDepth += 1

end

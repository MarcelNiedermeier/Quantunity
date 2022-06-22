
############################################################################
## Two-qubit gates for Quantum Simulator - MPS ITensor version with full MPO
############################################################################


#####################
# Auxiliary functions
#####################


""" Auxiliary function to apply a given two-site gate "gateName" (as
specified in the Hilbert space) to a quantum circuit in the positions
given in the array pos. """
function apply_two_site_gate!(qc::QC_IT_MPS, pos::Array{Int64, 1}, gateName::String)

    if gateName ∉ ["CNOT", "SWAP", "CU"]
        error("Wrong gate name: not available or mispelled (choose CNOT, SWAP, CU). ")
    end

    # get sites
    sites = qc.IndexSet

    # construct two-site gate
    gate = op(gateName, sites[pos[1]], sites[pos[2]])

    orthogonalize!(qc.StateVector, pos[1])
    wf = (qc.StateVector[pos[1]]*qc.StateVector[pos[2]])*gate
    noprime!(wf)
    ind = uniqueinds(qc.StateVector[pos[1]], qc.StateVector[pos[2]])
    U, S, V = svd(wf, ind, cutoff=1E-8)
    qc.StateVector[pos[1]] = U
    qc.StateVector[pos[2]] = S*V

end


""" Auxiliary function to apply a given two-site gate "gateName" (as
specified in the Hilbert space) to a quantum circuit in the positions
given in the array pos. To be used with a gate U that is applied depending
on a control qubit. """
function apply_two_site_gate!(qc::QC_IT_MPS, pos::Array{Int64, 1}, gateName::String, U::Matrix{ComplexF64})

    if gateName ∉ ["CNOT", "SWAP", "CU"]
        error("Wrong gate name: not available or mispelled (choose CNOT, SWAP, CU). ")
    end

    # get sites
    sites = qc.IndexSet

    # get elements of operator U to be applied
    U11 = U[1,1]
    U12 = U[1,2]
    U21 = U[2,1]
    U22 = U[2,2]

    # construct two-site gate
    #gate = op(gateName, sites[pos[1]], sites[pos[2]]; U[1,1], U[1,2], U[2,1], U[2,2])
    gate = op(gateName, sites[pos[1]], sites[pos[2]]; U11, U12, U21, U22)

    orthogonalize!(qc.StateVector, pos[1])
    wf = (qc.StateVector[pos[1]]*qc.StateVector[pos[2]])*gate
    noprime!(wf)
    ind = uniqueinds(qc.StateVector[pos[1]], qc.StateVector[pos[2]])
    U, S, V = svd(wf, ind, cutoff=1E-8)
    qc.StateVector[pos[1]] = U
    qc.StateVector[pos[2]] = S*V

end


""" Function to update the representation of a quantum circuit qc for a
two-site controlled gate (specified by num) applied to positions given through
the array pos. Can choose whethe to update the full or the reduced representation
(or both). """
function update_representation_two_site!(qc::QC_IT_MPS, pos::Array{Int64, 1}, num::Int64, full::Bool, reduced::Bool)

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


###########################
# Two-site gates (adjacent)
###########################


""" Function to apply the CNOT-gate to two adjacent qubits in positions
pos = [pos[1], pos[2]] = [pos[1], pos[1]+1]. The "upper" qubit (pos[1])
is the control qubit, the "lower" qubit (pos[2]) the one that is acted
on. This function serves mostly as an auxiliary function in the more
general cnot! function, see below. """
function cnot2!(qc::QC_IT_MPS, pos::Array{Int64, 1}, cutoff=1E-8)

   if length(pos) != 2
       error("Incorrect specification of indices (need 2 positions).")
   elseif abs(pos[1]-pos[2]) != 1
       error("Incorrect specification of indices (need adjacent positions).")
   end

   # update state vector with two-site CNOT gate
   apply_two_site_gate!(qc, pos, "CNOT")

   # update representing matrix of quantum circuit
   update_representation_two_site!(qc, pos, 3, true, false)

   # update circuit depth and bond dimension
   qc.CircuitDepth += 1
   push!(qc.BondDim, maxlinkdim(qc.StateVector))

end


""" Function to swap two adjacent qubits in positions pos = [pos[1], pos[2]]
= [pos[1], pos[1]+1]. """
function swap2!(qc::QC_IT_MPS, pos::Array{Int64, 1})

   # check structure of gate indices
   if length(pos) != 2
       error("Incorrect specification of indices (need 2 positions).")
   elseif abs(pos[1]-pos[2]) != 1
       error("Incorrect specification of indices (need adjacent positions).")
   end

   # update state vector with two-site CNOT gate
   apply_two_site_gate!(qc, pos, "SWAP")

   # update representing matrix of quantum circuit
   update_representation_two_site!(qc, pos, 4, true, false)

   # update circuit depth and bond dimension
   qc.CircuitDepth += 1
   push!(qc.BondDim, maxlinkdim(qc.StateVector))
end


""" Function to apply a controlled-U gate to two adjacent qubits in positions
pos = [pos[1], pos[2]] = [pos[1], pos[1]+1], where U is a unitary single-qubit
gate. The "upper" qubit (pos[1]) is the control qubit, the "lower" qubit
(pos[2]) the one that is acted on. This function serves mostly as an auxiliary
function in the more general CU! function, see below. """
function CU2!(qc::QC_IT_MPS, U::Matrix{ComplexF64}, pos::Array{Int64, 1})

   # check if always gets two indices!
   if length(pos) != 2
       error("Incorrect specification of indices (need 2 positions).")
   elseif abs(pos[1]-pos[2]) != 1
       error("Incorrect specification of indices (need adjacent positions).")
   end

   # update state vector with two-site CNOT gate
   apply_two_site_gate!(qc, pos, "CU", U)

   # update representing matrix of quantum circuit
   update_representation_two_site!(qc, pos, 13, true, false)

   # update circuit depth and bond dimension
   qc.CircuitDepth += 1
   push!(qc.BondDim, maxlinkdim(qc.StateVector))

end


################################
# Two-qubit gates (NOT adjacent)
################################


""" Function to perfom a half-swap of two qubits over an arbitrary separation,
which is given by the array pos. pos[1] is understood to be smaller than
pos[2]. The swap is obtained through an equivalent sequence of 2-swaps. Primarily
needed as an auxiliary function in the 2-qubit gates acting on non-adjacent
qubits (half-swap to bring them next to each other, then execute the gate, then
undo the half-swap). For a full swap, see the function fullSwap!. """
function swap!(qc::QC_IT_MPS, pos::Array{Int64, 1})

    # check correct format of indices
    if length(pos) != 2
        error("Incorrect specification of indices (need 2 positions).")
    elseif pos[1] >= pos[2]
        error("Incorrect specification of indices (pos[1] should be smaller than pos[2]).")
    end


    """
    sites_tmp1 = qc.IndexSet[pos[1]]
    sites_tmp2 = qc.IndexSet[pos[2]]

    qc.IndexSet[pos[1]] = sites_tmp2
    qc.IndexSet[pos[2]] = sites_tmp1

    qc.StateVector = permute(qc.StateVector, qc.IndexSet)
    """

    if pos[1] == (pos[2]-1)
        swap2!(qc::QC_IT_MPS, pos)
    else
        positions = get_2_perms_old(pos)
        for p in positions
            println("swap pos p ", p)
            swap2!(qc::QC_IT_MPS, p)
        end
    end
end


""" Function to reverse the arbitrary half-swap, performed by the swap! function.
For more information, see docstring of the swap! function. """
function unswap!(qc::QC_IT_MPS, pos::Array{Int64, 1})

    # check correct format of indices
    if length(pos) != 2
        error("Incorrect specification of indices (need 2 positions).")
    elseif pos[1] >= pos[2]
        error("Incorrect specification of indices (pos[1] should be smaller than pos[2]).")
    end

    if pos[1] == (pos[2]-1)
        println("doing this unswap")
        swap2!(qc::QC_IT_MPS, pos)
    else
        positions = get_2_perms_old(pos)
        for p in reverse(positions)
            println("unswap pos p ", p)
            swap2!(qc::QC_IT_MPS, p)
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
function cnot!(qc::QC_IT_MPS, pos::Array{Int64, 1})

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
                cnot2!(qc::QC_IT_MPS, pos)

            # qubits not on adjacent lines
            else
                swap!(qc::QC_IT_MPS, [pos[1], pos[2]-1])
                cnot2!(qc::QC_IT_MPS, [pos[2]-1, pos[2]])
                unswap!(qc::QC_IT_MPS, [pos[1], pos[2]-1])
            end

        # control qubit "below" qubit that is acted on
        else

            # qubits on two adjacent lines
            if pos[1] == (pos[2]+1)
                swap2!(qc::QC_IT_MPS, pos)
                cnot2!(qc::QC_IT_MPS, pos)
                swap2!(qc::QC_IT_MPS, pos)

            # qubits on two adjacent lines
            else
                swap!(qc::QC_IT_MPS, [pos[2], pos[1]])
                cnot2!(qc::QC_IT_MPS, [pos[1]-1, pos[1]])
                unswap!(qc::QC_IT_MPS, [pos[2], pos[1]])
            end
        end

        # update representing matrix of quantum circuit
        update_representation_two_site!(qc, pos, 3, false, true)

    # for Master topology: construct MPO representing CNOT
    else

        # initialise CNOT
        cnot = OpSum()
        cnot += "Proj00", pos[1]
        cnot += "Proj11", pos[1], "X", pos[2]

        # get sites
        sites = qc.IndexSet

        # turn into MPO
        gate = MPO(cnot, sites)

        # apply to state vector
        qc.StateVector = contract(gate, qc.StateVector, maxdim=qc.MaxBondDim)#, method="naive")

        # update representing matrix of quantum circuit
        update_representation_two_site!(qc, pos, 3, true, true)

        # update cicruit depth and bond dimension
        qc.CircuitDepth += 1
        push!(qc.BondDim, maxlinkdim(qc.StateVector))

    end
end


""" Function to perform a full swap of two qubits. For a linear topology,
the swap is decomposed into the corresponding 2-swaps. In the case of a
master topology, the gate is built from composing three CNOT gates on
pos, reverse(pos) and pos. """
function fullSwap!(qc::QC_IT_MPS, pos::Array{Int64, 1})

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
            swap2!(qc::QC_IT_MPS, pos)

        # qubits not on adjacent lines
        else
            swap!(qc::QC_IT_MPS, [pos[1], pos[2]])
            unswap!(qc::QC_IT_MPS, [pos[1], pos[2]-1])
        end

        # update compressed representation
        update_representation_two_site!(qc, pos, 4, false, true)

    # for Master topology: construct MPO representing full SWAP
    else

        # construct SWAP
        swap = OpSum()
        swap += "Proj00", pos[1], "Proj00", pos[2]
        swap += "S+", pos[1], "S-", pos[2]
        swap += "S-", pos[1], "S+", pos[2]
        swap += "Proj11", pos[1], "Proj11", pos[2]

        # get sites
        sites = qc.IndexSet

        # turn into MPO
        gate = MPO(swap, sites)

        # apply to state vector
        qc.StateVector = contract(gate, qc.StateVector, maxdim=qc.MaxBondDim, method="naive")

        # update representing matrix of quantum circuit
        update_representation_two_site!(qc, pos, 4, true, true)

        # update cicruit depth and bond dimension
        qc.CircuitDepth += 1
        push!(qc.BondDim, maxlinkdim(qc.StateVector))

    end
end


""" Function to apply the CU gate to two arbitrary qubits, defined by
positions in pos = [pos[1], pos[2]]. U is a unitary single-qubit gate.
The first qubits in is the control qubit, the second the qubit that is
acted on. The control qubit can be "above" or "below" the action qubit,
the function implements the desired arrangement. Depending on whether a
linear or a master topology is chosen, the CU gate is either constructed
with auxiliary swaps or implemented directly as a matrix. """
function CU!(qc::QC_IT_MPS, U, pos::Array{Int64, 1})

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
                CU2!(qc::QC_IT_MPS, U, pos)

            # qubits not on adjacent lines
            else
                swap!(qc::QC_IT_MPS, [pos[1], pos[2]-1])
                CU2!(qc::QC_IT_MPS, U, [pos[2]-1, pos[2]])
                unswap!(qc::QC_IT_MPS, [pos[1], pos[2]-1])
            end

        # control qubit "below" qubit that is acted on
        else

            # qubits on two adjacent lines
            if pos[1] == (pos[2]+1)
                swap2!(qc::QC_IT_MPS, pos)
                CU2!(qc::QC_IT_MPS, U, pos)
                swap2!(qc::QC_IT_MPS, pos)

            # qubits on two adjacent lines
            else
                swap!(qc::QC_IT_MPS, [pos[2], pos[1]])
                CU2!(qc::QC_IT_MPS, U, [pos[1]-1, pos[1]])
                unswap!(qc::QC_IT_MPS, [pos[2], pos[1]])
            end
        end

        # update representing matrix of quantum circuit
        update_representation_two_site!(qc, pos, 13, false, true)

    # for Master topology: construct MPO representing CU
    else

        # initialise CU
        cu = OpSum()
        cu += "Proj00", pos[1]
        cu += U[1,1], "Proj11", pos[1], "Proj00", pos[2]
        cu += U[1,2], "Proj11", pos[1], "S+", pos[2]
        cu += U[2,1], "Proj11", pos[1], "S-", pos[2]
        cu += U[2,2], "Proj11", pos[1], "Proj11", pos[2]

        # get sites
        sites = qc.IndexSet

        # turn into MPO
        gate = MPO(cu, sites)

        # apply to state vector
        qc.StateVector = contract(gate, qc.StateVector, maxdim=qc.MaxBondDim, method="naive")

        # update representing matrix of quantum circuit
        update_representation_two_site!(qc, pos, 13, true, true)

        # update circuit depth and bond dimension
        qc.CircuitDepth += 1
        push!(qc.BondDim, maxlinkdim(qc.StateVector))

    end
end

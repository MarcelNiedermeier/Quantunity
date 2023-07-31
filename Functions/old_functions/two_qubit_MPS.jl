
############################################################################
## Two-qubit gates for Quantum Simulator - MPS ITensor version with full MPO
############################################################################


#####################
# Auxiliary functions
#####################


""" Auxiliary function to apply a given two-site gate "gateName" (as
specified in the Hilbert space) to a quantum circuit in the positions
given in the array pos. """
function apply_two_site_gate!(qc::QC_IT_MPS, pos::Array{Int64, 1},
    gateName::String, recordFidelity=false)

    if gateName ∉ ["CNOT", "SWAP"]
        error("Wrong gate name: not available or mispelled (choose CNOT, SWAP). ")
    end

    # get sites, construct two-site gate
    sites = qc.IndexSet
    #println(sites)
    #println("should be here")
    gate = op(gateName, sites[pos[1]], sites[pos[2]])

    # get untruncated MPS
    orthogonalize!(qc.StateVector, pos[1])
    wf = (qc.StateVector[pos[1]]*qc.StateVector[pos[2]])*gate
    noprime!(wf)
    ind = uniqueinds(qc.StateVector[pos[1]], qc.StateVector[pos[2]])
    #U_full, S_full, V_full = svd(wf, ind) # full SVD
    U, S, V = svd(wf, ind, maxdim=qc.MaxBondDim)
    #println("U: ", size(U_full))
    #println("S: ", size(S_full))
    #println("V: ", size(V_full))
    #psi = deepcopy(qc.StateVector)
    #psi[pos[1]] = U_full
    #psi[pos[2]] = S_full*V_full

    # apply, do SVD (optimised for truncated tensor)
    #orthogonalize!(qc.StateVector, pos[1])
    #wf = (qc.StateVector[pos[1]]*qc.StateVector[pos[2]])*gate
    #noprime!(wf)
    #ind = uniqueinds(qc.StateVector[pos[1]], qc.StateVector[pos[2]])
    #U, S, V = svd(wf, ind, maxdim=qc.MaxBondDim)

    #println("U trunc: ", size(U))
    #println("S trunc: ", size(S))
    #println("V trunc: ", size(V))

    # calculate fidelity: find two parts of sum of singular values
    sumS_1 = sum(diag(array(S)*array(S)))
    sumS_2 = norm(wf - U*S*V)^2
    fid = (sumS_1)/(sumS_1 + sumS_2)

    # replace updated tensors in MPS
    qc.StateVector[pos[1]] = U
    qc.StateVector[pos[2]] = S*V

    # check fidelity calculation by finding exact overlap
    #sum_tot = sum(diag(array(S_full)*array(S_full)))
    #fid_ex = (sumS_1)/(sum_tot)

    #if recordFidelity == true
    #    push!(qc.TwoQubitFidelityExactSelected, fid_ex)
    #    NfidSel_ex = 1.0
    #    for f in qc.TwoQubitFidelityExactSelected
    #        NfidSel_ex *= f
    #    end
    #    push!(qc.NQubitFidelityExactSelected, NfidSel_ex)
    #    push!(qc.AverageTwoQubitFidelityExactSelected, NfidSel_ex^(1/length(qc.TwoQubitFidelityExactSelected)))
    #end

    # calculate N-qubit fidelity and av. 2-qubit fidelity, update quantum circuit
    push!(qc.TwoQubitFidelity, fid)
    Nfid = 1.0
    for f in qc.TwoQubitFidelity
        Nfid *= f
    end
    push!(qc.NQubitFidelity, Nfid)
    push!(qc.AverageTwoQubitFidelity, Nfid^(1/length(qc.TwoQubitFidelity)))

    # keep second list of selected fidelities
    if recordFidelity == true
        push!(qc.TwoQubitFidelitySelected, fid)
        NfidSel = 1.0
        for f in qc.TwoQubitFidelitySelected
            NfidSel *= f
        end
        push!(qc.NQubitFidelitySelected, NfidSel)
        push!(qc.AverageTwoQubitFidelitySelected, NfidSel^(1/length(qc.TwoQubitFidelitySelected)))
    end

end


""" Auxiliary function to apply a given two-site gate "gateName" (as
specified in the Hilbert space) to a quantum circuit in the positions
given in the array pos. To be used with a gate U that is applied depending
on a control qubit. """
function apply_two_site_gate_U!(qc::QC_IT_MPS, pos::Array{Int64, 1},
    gateName::String, U::Matrix{ComplexF64})

    if gateName ∉ ["CU"]
        error("Wrong gate name: not available or mispelled (choose CU). ")
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
    gate = op(gateName, sites[pos[1]], sites[pos[2]]; U11=U11, U12=U12, U21=U21, U22=U22)

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
given in the array pos. To be used with an additional paramter n. """
function apply_two_site_gate_n!(qc::QC_IT_MPS, pos::Array{Int64, 1},
    gateName::String, n::Number)

    #if gateName ∉ ["CNOT", "SWAP", "CRn"]
    if gateName ∉ ["CRn"]
        error("Wrong gate name: not available or mispelled (choose CRn). ")
    end

    # get sites
    sites = qc.IndexSet
    #println("here")
    #println(sites[pos[1]], sites[pos[2]])
    #println(gateName)

    # get elements of operator U to be applied
    #U11 = U[1,1]
    #U12 = U[1,2]
    #U21 = U[2,1]
    #U22 = U[2,2]

    # construct two-site gate
    #gate = op(gateName, sites[pos[1]], sites[pos[2]]; U[1,1], U[1,2], U[2,1], U[2,2])
    gate = op(gateName, sites[pos[1]], sites[pos[2]]; n=n)

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
   apply_two_site_gate!(qc, sort!(pos), "SWAP")

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
   apply_two_site_gate_U!(qc, pos, "CU", U)

   # update representing matrix of quantum circuit
   update_representation_two_site!(qc, pos, 13, true, false)

   # update circuit depth and bond dimension
   qc.CircuitDepth += 1
   push!(qc.BondDim, maxlinkdim(qc.StateVector))

end


""" Function to apply a controlled-U gate to two adjacent qubits in positions
pos = [pos[1], pos[2]] = [pos[1], pos[1]+1], where U is a unitary single-qubit
gate. The "upper" qubit (pos[1]) is the control qubit, the "lower" qubit
(pos[2]) the one that is acted on. This function serves mostly as an auxiliary
function in the more general CU! function, see below. """
function CRn2!(qc::QC_IT_MPS, n::Number, pos::Array{Int64, 1})

   # check if always gets two indices!
   if length(pos) != 2
       error("Incorrect specification of indices (need 2 positions).")
   elseif abs(pos[1]-pos[2]) != 1
       error("Incorrect specification of indices (need adjacent positions).")
   end

   # update state vector with two-site CNOT gate
   apply_two_site_gate_n!(qc, pos, "CRn", n)

   # update representing matrix of quantum circuit
   update_representation_two_site!(qc, pos, 24, true, false)

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
            #println("swap pos p ", p)
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
        #println("doing this unswap")
        swap2!(qc::QC_IT_MPS, pos)
    else
        positions = get_2_perms_old(pos)
        for p in reverse(positions)
            #println("unswap pos p ", p)
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
function U_2site!(qc::QC_IT_MPS, U, pos::Array{Int64, 1}, update_rep=true)

    # check correct format of indices
    if length(pos) != 2
        error("Incorrect specification of indices (need 2 positions).")
    elseif pos[2] - pos[1] != 1
        error("Incorrect specification of indices.")
    end

    # apply general CNOT via SWAP gates for linear topology
    if qc.LinearTopology == true

        println("Linear topology currently not implemented")

        ## control qubit "above" qubit that is acted on
        #if pos[2] > pos[1]

        #    # qubits on two adjacent lines
        #    if pos[1] == (pos[2]-1)
        #        cnot2!(qc::QC_IT_MPS, pos)

        #    # qubits not on adjacent lines
        #    else
        #        swap!(qc::QC_IT_MPS, [pos[1], pos[2]-1])
        #        cnot2!(qc::QC_IT_MPS, [pos[2]-1, pos[2]])
        #        unswap!(qc::QC_IT_MPS, [pos[1], pos[2]-1])
        #    end

        ## control qubit "below" qubit that is acted on
        #else

        #    # qubits on two adjacent lines
        #    if pos[1] == (pos[2]+1)
        #        swap2!(qc::QC_IT_MPS, pos)
        #        cnot2!(qc::QC_IT_MPS, pos)
        #        swap2!(qc::QC_IT_MPS, pos)

        #    # qubits on two adjacent lines
        #    else
        #        swap!(qc::QC_IT_MPS, [pos[2], pos[1]])
        #        cnot2!(qc::QC_IT_MPS, [pos[1]-1, pos[1]])
        #        unswap!(qc::QC_IT_MPS, [pos[2], pos[1]])
        #    end
        #end

        # update representing matrix of quantum circuit
        #if update_rep
        #    update_representation_two_site!(qc, pos, 3, false, true)
        #end

    # for Master topology: construct MPO representing CNOT
    else

        # initialise U_2site
        u_2site = OpSum()
        u_2site += U[1,1], "Proj00", pos[1], "Proj00", pos[2]
        u_2site += U[1,2], "Proj00", pos[1], "Proj01", pos[2]
        u_2site += U[1,3], "Proj01", pos[1], "Proj00", pos[2]
        u_2site += U[1,4], "Proj01", pos[1], "Proj01", pos[2]

        u_2site += U[2,1], "Proj00", pos[1], "Proj10", pos[2]
        u_2site += U[2,2], "Proj00", pos[1], "Proj11", pos[2]
        u_2site += U[2,3], "Proj01", pos[1], "Proj10", pos[2]
        u_2site += U[2,4], "Proj01", pos[1], "Proj11", pos[2]

        u_2site += U[3,1], "Proj10", pos[1], "Proj00", pos[2]
        u_2site += U[3,2], "Proj10", pos[1], "Proj01", pos[2]
        u_2site += U[3,3], "Proj11", pos[1], "Proj00", pos[2]
        u_2site += U[3,4], "Proj11", pos[1], "Proj01", pos[2]

        u_2site += U[4,1], "Proj10", pos[1], "Proj10", pos[2]
        u_2site += U[4,2], "Proj10", pos[1], "Proj11", pos[2]
        u_2site += U[4,3], "Proj11", pos[1], "Proj10", pos[2]
        u_2site += U[4,4], "Proj11", pos[1], "Proj11", pos[2]

        # get sites
        sites = qc.IndexSet

        # turn into MPO
        gate = MPO(u_2site, sites)

        # apply to state vector, unprime
        qc.StateVector = contract(gate, qc.StateVector, maxdim=qc.MaxBondDim)#, method="naive")
        noprime!(qc.StateVector)

        # update representing matrix of quantum circuit
        if update_rep
            update_representation_two_site!(qc, pos, 25, true, true)
        end

        # update cicruit depth and bond dimension
        qc.CircuitDepth += 1
        push!(qc.BondDim, maxlinkdim(qc.StateVector))

    end
end



""" Function to apply the CNOT gate to two arbitrary qubits, defined by
positions in pos = [pos[1], pos[2]]. The first qubits in is the control
qubit, the second the qubit that is acted on. The control qubit can be
"above" or "below" the action qubit, the function implements the desired
arrangement. Depending on whether a linear or a master topology is chosen,
the CNOT gate is either constructed with auxiliary swaps or implemented directly
as a matrix. """
function cnot!(qc::QC_IT_MPS, pos::Array{Int64, 1}, update_rep=true)

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
        if update_rep
            update_representation_two_site!(qc, pos, 3, false, true)
        end

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

        # apply to state vector, unprime
        qc.StateVector = contract(gate, qc.StateVector, maxdim=qc.MaxBondDim)#, method="naive")
        noprime!(qc.StateVector)

        # update representing matrix of quantum circuit
        if update_rep
            update_representation_two_site!(qc, pos, 3, true, true)
        end

        # update cicruit depth and bond dimension
        qc.CircuitDepth += 1
        push!(qc.BondDim, maxlinkdim(qc.StateVector))

    end
end


""" Function to perform a full swap of two qubits. For a linear topology,
the swap is decomposed into the corresponding 2-swaps. In the case of a
master topology, the gate is built from composing three CNOT gates on
pos, reverse(pos) and pos. """
function fullSwap!(qc::QC_IT_MPS, pos::Array{Int64, 1}, update_rep=true)

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

        # apply to state vector, unprime
        qc.StateVector = contract(gate, qc.StateVector, maxdim=qc.MaxBondDim)#, method="naive")
        noprime!(qc.StateVector)

        # update representing matrix of quantum circuit
        if update_rep
            update_representation_two_site!(qc, pos, 4, true, true)
        end

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
function CU!(qc::QC_IT_MPS, U, pos::Array{Int64, 1}, update_rep=true)

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

        # apply to state vector, unprime
        qc.StateVector = contract(gate, qc.StateVector, maxdim=qc.MaxBondDim, method="naive")
        noprime!(qc.StateVector)

        # update representing matrix of quantum circuit
        if update_rep
            update_representation_two_site!(qc, pos, 13, true, true)
        end

        # update circuit depth and bond dimension
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
function CRn!(qc::QC_IT_MPS, n::Number, pos::Array{Int64, 1}, update_rep=true)

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
                CRn2!(qc::QC_IT_MPS, n, pos)

            # qubits not on adjacent lines
            else
                swap!(qc::QC_IT_MPS, [pos[1], pos[2]-1])
                CRn2!(qc::QC_IT_MPS, n, [pos[2]-1, pos[2]])
                unswap!(qc::QC_IT_MPS, [pos[1], pos[2]-1])
            end

        # control qubit "below" qubit that is acted on
        else

            # qubits on two adjacent lines
            if pos[1] == (pos[2]+1)
                swap2!(qc::QC_IT_MPS, pos)
                CRn2!(qc::QC_IT_MPS, n, pos)
                swap2!(qc::QC_IT_MPS, pos)

            # qubits on two adjacent lines
            else
                swap!(qc::QC_IT_MPS, [pos[2], pos[1]])
                CRn2!(qc::QC_IT_MPS, n, [pos[1]-1, pos[1]])
                unswap!(qc::QC_IT_MPS, [pos[2], pos[1]])
            end
        end

        # update representing matrix of quantum circuit
        update_representation_two_site!(qc, pos, 24, false, true)

    # for Master topology: construct MPO representing CU
    else

        # initialise CRn
        crn = OpSum()
        crn += "Proj00", pos[1]
        crn += Complex(1.), "Proj11", pos[1], "Proj00", pos[2]
        crn += Complex(0.), "Proj11", pos[1], "S+", pos[2]
        crn += Complex(0.), "Proj11", pos[1], "S-", pos[2]
        crn += exp(sign(n)*2π*1.0im/2^abs(n)), "Proj11", pos[1], "Proj11", pos[2]

        # get sites
        sites = qc.IndexSet

        # turn into MPO
        gate = MPO(crn, sites)
        #println(gate)

        # apply to state vector, unprime
        qc.StateVector = contract(gate, qc.StateVector, maxdim=qc.MaxBondDim)#, method="naive")
        noprime!(qc.StateVector)

        # update representing matrix of quantum circuit
        if update_rep
            update_representation_two_site!(qc, pos, 24, true, true)
        end

        # update circuit depth and bond dimension
        qc.CircuitDepth += 1
        push!(qc.BondDim, maxlinkdim(qc.StateVector))

    end
end



#################
# Random Circuits
#################


""" Function to apply a (vertical) sequence of CNOT gates to a quantum
circuit. Use together with random unitaries to create a highly entangled
state in the fastest possible way. Can specify whether the first CNOT
gate starts on the first or on the second qubit, the others follow
automatically after that. """
function sequential_cnot!(qc::QC_IT_MPS, start_pos::Array{Int64, 1})

    if start_pos[1] != 1 && start_pos[1] != 2
        error("Incorrect specification of starting point (pos = [1] or [2]).")
    end

    # get number of qubits and matrices
    N = qc.NumQubits
    s = Index(2, "QCircuit")
    E = array(op("E", s))

    # first CNOT starting on first qubit line
    if start_pos[1] == 1

        # even number of qubits, record 2-qubit fidelity in middle
        if iseven(N)

            # get "middle" positions in bulk
            middle_pos = [Int(N/2)-1, Int(N/2)]
            for i in 1:2:N
                if i in middle_pos
                    apply_two_site_gate!(qc, [i, i+1], "CNOT", true)
                else
                    apply_two_site_gate!(qc, [i, i+1], "CNOT", false)
                end
            end

        # odd number of qubits, record 2-qubit fidelity in middle
        else

            # get "middle" positions in bulk
            middle_pos = [Int((N-1)/2)-1, Int((N-1)/2)]
            for i in 1:2:(N-1)
                if i in middle_pos
                    apply_two_site_gate!(qc, [i, i+1], "CNOT", true)
                else
                    apply_two_site_gate!(qc, [i, i+1], "CNOT", false)
                end
            end
        end

    # first CNOT starting on second qubit line
    elseif start_pos[1] == 2

        # even number of qubits, record 2-qubit fidelity in middle
        if iseven(qc.NumQubits)

            # get "middle" positions in bulk
            middle_pos = [Int(N/2)-1, Int(N/2)]
            for i in 2:2:(N-1)
                if i in middle_pos
                    apply_two_site_gate!(qc, [i, i+1], "CNOT", true)
                else
                    apply_two_site_gate!(qc, [i, i+1], "CNOT", false)
                end
            end

        # odd number of qubits, record 2-qubit fidelity in middle
        else

            # get "middle" positions in bulk
            middle_pos = [Int((N-1)/2)-1, Int((N-1)/2)]
            for i in 2:2:N
                if i in middle_pos
                    apply_two_site_gate!(qc, [i, i+1], "CNOT", true)
                else
                    apply_two_site_gate!(qc, [i, i+1], "CNOT", false)
                end
            end
        end
    end

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

    # update circuit depth and bond dimension
    qc.CircuitDepth += 1
    push!(qc.BondDim, maxlinkdim(qc.StateVector))

end


## second sequential CNOT with specified start and end positions


""" Function to apply a (vertical) sequence of CNOT gates to a quantum
circuit. Use together with random unitaries to create a highly entangled
state in the fastest possible way. Can specify whether the first CNOT
gate starts on the first or on the second qubit, the others follow
automatically after that. """
function sequential_cnot!(qc::QC_IT_MPS, pos::Number, num::Number,
     update_rep=false, record_fidelity=false)


    # get number of qubits and matrices
    N = qc.NumQubits
    s = Index(2, "QCircuit")
    E = array(op("E", s))

    # collect positions where CNOTs are applied
    cnot_control_pos = []
    cnot_action_pos = []

    # approximated middle position if fidelity should be recorded
    middle_pos = [round((pos+num-1)/2), round((pos+num-1)/2)+1]

    # first CNOT starting on first qubit line
    if pos == 1

        # even number of qubits, record 2-qubit fidelity in middle
        if iseven(num)

            # get "middle" positions in bulk
            #middle_pos = [Int((pos+num)/2)-1, Int((pos+num)/2)]
            # apply sequetial CNOTs, record fidelity in bulk if desired
            for i in 1:2:(pos+num-1)
                if i in middle_pos
                    #println("would be recording fidelity in pos [$i, $(i+1)]")
                    apply_two_site_gate!(qc, [i, i+1], "CNOT", record_fidelity)
                    push!(cnot_control_pos, i)
                    push!(cnot_action_pos, i+1)
                else
                    apply_two_site_gate!(qc, [i, i+1], "CNOT", false)
                    push!(cnot_control_pos, i)
                    push!(cnot_action_pos, i+1)
                end
            end

        # odd number of qubits, record 2-qubit fidelity in middle
        else

            # get "middle" positions in bulk
            #middle_pos = [Int((pos+num-1)/2)-1, Int((pos+num-1)/2)]
            for i in 1:2:(pos+num-3)
                if i in middle_pos
                    #println("would be recording fidelity in pos [$i, $(i+1)]")
                    apply_two_site_gate!(qc, [i, i+1], "CNOT", record_fidelity)
                    push!(cnot_control_pos, i)
                    push!(cnot_action_pos, i+1)
                else
                    apply_two_site_gate!(qc, [i, i+1], "CNOT", false)
                    push!(cnot_control_pos, i)
                    push!(cnot_action_pos, i+1)
                end
            end
        end

    # sequence of CNOTs doesn't start on first qubit line
    else

        # even number of qubits, record 2-qubit fidelity in middle
        if iseven(num)

            # get "middle" positions in bulk
            #middle_pos = [Int((pos+num-1)/2)-1, Int((pos+num-1)/2)]
            for i in pos:2:(pos+num-1)
                if i in middle_pos
                    #println("would be recording fidelity in pos [$i, $(i+1)]")
                    apply_two_site_gate!(qc, [i, i+1], "CNOT", record_fidelity)
                    push!(cnot_control_pos, i)
                    push!(cnot_action_pos, i+1)
                else
                    apply_two_site_gate!(qc, [i, i+1], "CNOT", false)
                    push!(cnot_control_pos, i)
                    push!(cnot_action_pos, i+1)
                end
            end

        # odd number of qubits, record 2-qubit fidelity in middle
        else

            # get "middle" positions in bulk
            #middle_pos = [Int((pos+num-1-1)/2)-1, Int((pos+num-1-1)/2)]
            for i in pos:2:(pos+num-3)
                if i in middle_pos
                    #println("would be recording fidelity in pos [$i, $(i+1)]")
                    apply_two_site_gate!(qc, [i, i+1], "CNOT", record_fidelity)
                    push!(cnot_control_pos, i)
                    push!(cnot_action_pos, i+1)
                else
                    apply_two_site_gate!(qc, [i, i+1], "CNOT", false)
                    push!(cnot_control_pos, i)
                    push!(cnot_action_pos, i+1)
                end
            end
        end
    end

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

    # update circuit depth and bond dimension
    qc.CircuitDepth += 1
    push!(qc.BondDim, maxlinkdim(qc.StateVector))

end


############################################################################
## Two-qubit gates for Quantum Simulator - MPS ITensor version with full MPO
############################################################################


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

   """
   # get matrices
   #E, _, _, _ = get_Pauli_matrices()
   cnot = Complex.([1. 0. 0. 0.; 0. 1. 0. 0.; 0. 0. 0. 1.; 0. 0. 1. 0.])
   cnot_reverse = Complex.([1. 0. 0. 0.; 0. 0. 0. 1.; 0. 0. 1. 0.; 0. 1. 0. 0.])

   # outgoing indices
   #sites_out = siteinds(2, N)

   # get MPS and corresponding indices
   #sites = qc.IndexSet
   #psi = qc.StateVector

   # control qubit above action qubit
   if pos[2]-1 == pos[1]
       gate = ITensor(cnot, qc.IndexSet[pos[1]], ITensors.prime(qc.IndexSet[pos[1]]), qc.IndexSet[pos[2]], ITensors.prime(qc.IndexSet[pos[2]]))
       #orthogonalize!(qc.StateVector, pos[1])
   # control qubit below action qubit
   else
       gate = ITensor(cnot_reverse, qc.IndexSet[pos[1]], ITensors.prime(qc.IndexSet[pos[1]]), qc.IndexSet[pos[2]], ITensors.prime(qc.IndexSet[pos[2]]))
       #orthogonalize!(qc.StateVector, pos[2])
   end

   orthogonalize!(qc.StateVector, pos[1])
   updatedTensor = (qc.StateVector[pos[1]]*qc.StateVector[pos[2]])*gate
   noprime!(updatedTensor)
   inds = uniqueinds(qc.StateVector[pos[1]], qc.StateVector[pos[2]])
   U, S, V = svd(updatedTensor, inds, cutoff=cutoff)
   qc.StateVector[pos[1]] = U
   qc.StateVector[pos[2]] = S*V
   """

   """
   # apply identities
   for i in 1:pos[1]-1
       gate = ITensor(E, qc.IndexSet[i], sites_out[i])
       qc.StateVector = gate*qc.StateVector
       qc.IndexSet[i] = sites_out[i]
   end

   # apply two-site gate
   gate = ITensor(cnot, qc.IndexSet[pos[1]], sites_out[pos[1]], qc.IndexSet[pos[2]], sites_out[pos[2]])
   qc.StateVector = gate*qc.StateVector
   qc.IndexSet[pos[1]] = sites_out[pos[1]]
   qc.IndexSet[pos[2]] = sites_out[pos[2]]

   # apply identities
   for i in pos[2]+1:qc.NumQubits
       gate = ITensor(E, qc.IndexSet[i], sites_out[i])
       qc.StateVector = gate*qc.StateVector
       qc.IndexSet[i] = sites_out[i]
   end
   """

   """
   if pos[1] == 1
       gate = cnot

       for i in 3:qc.NumQubits
           gate = kron(gate, E)
       end
   else
       gate = E
       for i in 2:pos[1]-1
           gate = kron(gate, E)
       end
       gate = kron(gate, cnot)
       for i in (pos[2]+1):qc.NumQubits
           gate = kron(gate, E)
       end
   end

   # create a new set of outgoing indices and define ITensor from array
   sites_out = siteinds(2, N)
   gate = ITensor(gate, qc.IndexSet, sites_out)

   # update state vector and outgoing indices
   qc.StateVector = gate*qc.StateVector
   qc.IndexSet = sites_out
   """

   # initialise CNOT
   cnot = OpSum()
   cnot += "Proj00", pos[1]
   cnot += "Proj11", pos[1], "X", pos[2]

   # get sites
   sites = qc.IndexSet

   # turn into MPO
   gate = MPO(cnot, sites)

   # apply to state vector
   qc.StateVector = contract(gate, qc.StateVector, maxdim=qc.MaxBondDim)#, method=qc.ContractionMethod)

   # update representing matrix of quantum circuit
   for i in 1:qc.NumQubits
       if i == pos[1]
           push!(qc.RepresentationFull[i], 3) # control qubit
       elseif i == pos[2]
           push!(qc.RepresentationFull[i], -3) # action qubit
       else
           push!(qc.RepresentationFull[i], 0)
       end
   end

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

   """
   # get matrices
   #E, _, _, _ = get_Pauli_matrices()
   #swap = Complex.([1. 0. 0. 0.; 0. 0. 1. 0.; 0. 1. 0. 0.; 0. 0. 0. 1.])

   sites_tmp1 = qc.IndexSet[pos[1]]
   sites_tmp2 = qc.IndexSet[pos[2]]

   qc.IndexSet[pos[1]] = sites_tmp2
   qc.IndexSet[pos[2]] = sites_tmp1

   qc.StateVector = permute(qc.StateVector, qc.IndexSet)
   """

   """
   # outgoing indices
   sites_out = siteinds(2, N)

   println("here1")

   # apply identities
   for i in 1:pos[1]-1
       println("here2")
       gate = ITensor(E, qc.IndexSet[i], sites_out[i])
       qc.StateVector = gate*qc.StateVector
       qc.IndexSet[i] = sites_out[i]
   end

   # apply two-site gate
   println("here3")
   gate = ITensor(swap, qc.IndexSet[pos[1]], sites_out[pos[1]], qc.IndexSet[pos[2]], sites_out[pos[2]])
   qc.StateVector = gate*qc.StateVector
   qc.IndexSet[pos[1]] = sites_out[pos[1]]
   qc.IndexSet[pos[2]] = sites_out[pos[2]]

   # apply identities
   for i in pos[2]+1:qc.NumQubits
       println("here4")
       gate = ITensor(E, qc.IndexSet[i], sites_out[i])
       qc.StateVector = gate*qc.StateVector
       qc.IndexSet[i] = sites_out[i]
   end
   """

   """
   if pos[1] == 1
       gate = swap

       for i in 3:qc.NumQubits
           gate = kron(gate, E)
       end
   else
       gate = E
       for i in 2:pos[1]-1
           gate = kron(gate, E)
       end
       gate = kron(gate, swap)
       for i in (pos[2]+1):qc.NumQubits
           gate = kron(gate, E)
       end
   end

   # create a new set of outgoing indices and define ITensor from array
   sites_out = siteinds(2, N)
   gate = ITensor(gate, qc.IndexSet, sites_out)

   # update state vector and outgoing indices
   qc.StateVector = gate*qc.StateVector
   qc.IndexSet = sites_out
   """

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
   qc.StateVector = contract(gate, qc.StateVector, maxdim=qc.MaxBondDim)#, method=qc.ContractionMethod)

   # update representing matrix of quantum circuit
   for i in 1:qc.NumQubits
       if i in pos
           push!(qc.RepresentationFull[i], 4)
       else
           push!(qc.RepresentationFull[i], 0)
       end
   end

   # update circuit depth and bond dimension
   qc.CircuitDepth += 1
   push!(qc.BondDim, maxlinkdim(qc.StateVector))
end


""" Function to apply a controlled-U gate to two adjacent qubits in positions
pos = [pos[1], pos[2]] = [pos[1], pos[1]+1], where U is a unitary single-qubit
gate. The "upper" qubit (pos[1]) is the control qubit, the "lower" qubit
(pos[2]) the one that is acted on. This function serves mostly as an auxiliary
function in the more general CU! function, see below. """
function CU2!(qc::QC_IT_MPS, U, pos::Array{Int64, 1})

   # check if always gets two indices!
   if length(pos) != 2
       error("Incorrect specification of indices (need 2 positions).")
   elseif abs(pos[1]-pos[2]) != 1
       error("Incorrect specification of indices (need adjacent positions).")
   end

   """
   # get matrices
   E, _, _, _ = get_Pauli_matrices()
   CU = Complex.(zeros(4,4))
   CU[1:2, 1:2] = E
   CU[3:4, 3:4] = U

   # outgoing indices
   sites_out = siteinds(2, N)

   # apply identities
   for i in 1:pos[1]-1
       gate = ITensor(E, qc.IndexSet[i], sites_out[i])
       qc.StateVector = gate*qc.StateVector
       qc.IndexSet[i] = sites_out[i]
   end

   # apply two-site gate
   gate = ITensor(CU, qc.IndexSet[pos[1]], sites_out[pos[1]], qc.IndexSet[pos[2]], sites_out[pos[2]])
   qc.StateVector = gate*qc.StateVector
   qc.IndexSet[pos[1]] = sites_out[pos[1]]
   qc.IndexSet[pos[2]] = sites_out[pos[2]]

   # apply identities
   for i in pos[2]+1:qc.NumQubits
       gate = ITensor(E, qc.IndexSet[i], sites_out[i])
       qc.StateVector = gate*qc.StateVector
       qc.IndexSet[i] = sites_out[i]
   end
   """

   """
   if pos[1] == 1
       gate = CU

       for i in 3:qc.NumQubits
           gate = kron(gate, E)
       end
   else
       gate = E
       for i in 2:pos[1]-1
           gate = kron(gate, E)
       end
       gate = kron(gate, CU)
       for i in (pos[2]+1):qc.NumQubits
           gate = kron(gate, E)
       end
   end

   # create a new set of outgoing indices and define ITensor from array
   sites_out = siteinds(2, N)
   gate = ITensor(gate, qc.IndexSet, sites_out)

   # update state vector and outgoing indices
   qc.StateVector = gate*qc.StateVector
   qc.IndexSet = sites_out
   """


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
   qc.StateVector = contract(gate, qc.StateVector, maxdim=qc.MaxBondDim)#, method=qc.ContractionMethod)

   # update representing matrix of quantum circuit
   for i in 1:qc.NumQubits
       if i == pos[1]
           push!(qc.RepresentationFull[i], 13) # control qubit
       elseif i == pos[2]
           push!(qc.RepresentationFull[i], -13) # action qubit
       else
           push!(qc.RepresentationFull[i], 0)
       end
   end

   # update circuit depth and bond dimension
   qc.CircuitDepth += 1
   push!(qc.BondDim, maxlinkdim(qc.StateVector))

end



################################
# two-qubit gates (not adjacent)
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

    """
    # update representing matrix of quantum circuit
    for i in 1:qc.NumQubits
        if i in pos
            push!(qc.Representation[i], 4)
        else
            push!(qc.Representation[i], 0)
        end
    end

    # update circuit depth
    qc.CircuitDepth += 1
    """

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
        #println("Feature currently not implemented.")

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
        for i in 1:qc.NumQubits
            if i == pos[1]
                push!(qc.Representation[i], 3) # control qubit
            elseif i == pos[2]
                push!(qc.Representation[i], -3) # action qubit
            elseif minimum(pos) < i && i < maximum(pos)
                push!(qc.Representation[i], 15) # vertical line
            else
                push!(qc.Representation[i], 0)
            end
        end

    # apply general CNOT via direct construction of matrix for master topology
    else

        # initialise CNOT
        #cnot = AutoMPO()
        cnot = OpSum()
        cnot += "Proj00", pos[1]
        cnot += "Proj11", pos[1], "X", pos[2]

        # get sites
        sites = qc.IndexSet

        # turn into MPO
        gate = MPO(cnot, sites)

        # apply to state vector
        qc.StateVector = contract(gate, qc.StateVector, maxdim=qc.MaxBondDim, method="naive")

        # update representing matrix of quantum circuit
        for i in 1:qc.NumQubits
            if i == pos[1]
                push!(qc.Representation[i], 3) # control qubit
                push!(qc.RepresentationFull[i], 3) # control qubit
            elseif i == pos[2]
                push!(qc.Representation[i], -3) # action qubit
                push!(qc.RepresentationFull[i], -3) # action qubit
            elseif minimum(pos) < i && i < maximum(pos)
                push!(qc.Representation[i], 15) # vertical line
                push!(qc.RepresentationFull[i], 15) # vertical line
            else
                push!(qc.Representation[i], 0)
                push!(qc.RepresentationFull[i], 0)
            end
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
        for i in 1:qc.NumQubits
            if i == pos[1]
                push!(qc.Representation[i], 4) # swap qubit
            elseif i == pos[2]
                push!(qc.Representation[i], 4) # swap qubit
            elseif minimum(pos) < i && i < maximum(pos)
                push!(qc.Representation[i], 15) # vertical line
            else
                push!(qc.Representation[i], 0)
            end
        end

    # for Master topology: compose full SWAP with CNOT gates
    else

        """
        # elementary decomposition of SWAP gate
        cnot!(qc::QC_IT_ED, pos)
        cnot!(qc::QC_IT_ED, reverse(pos))
        cnot!(qc::QC_IT_ED, pos)

        # delete CNOTs from representation and replace with SWAP
        for i in 1:qc.NumQubits
            if i == pos[1]
                #deleteat!(qc.Representation[i], [end-2, end-1, end])
                pop!(qc.Representation[i])
                pop!(qc.Representation[i])
                pop!(qc.Representation[i])
                push!(qc.Representation[i], 4) # swap qubit
            elseif i == pos[2]
                #deleteat!(qc.Representation[i], [end-2, end-1, end])
                pop!(qc.Representation[i])
                pop!(qc.Representation[i])
                pop!(qc.Representation[i])
                push!(qc.Representation[i], 4) # swap qubit
            elseif minimum(pos) < i && i < maximum(pos)
                #deleteat!(qc.Representation[i], [end-2, end-1, end])
                pop!(qc.Representation[i])
                pop!(qc.Representation[i])
                pop!(qc.Representation[i])
                push!(qc.Representation[i], 15) # vertical line
            else
                #deleteat!(qc.Representation[i], [end-1, end])
                pop!(qc.Representation[i])
                pop!(qc.Representation[i])
            end
        end
        """

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
        for i in 1:qc.NumQubits
            if i == pos[1]
                push!(qc.Representation[i], 4) # swap qubit
                push!(qc.RepresentationFull[i], 4) # swap qubit
            elseif i == pos[2]
                push!(qc.Representation[i], 4) # swap qubit
                push!(qc.RepresentationFull[i], 4) # swap qubit
            elseif minimum(pos) < i && i < maximum(pos)
                push!(qc.Representation[i], 15) # vertical line
                push!(qc.RepresentationFull[i], 15) # vertical line
            else
                push!(qc.Representation[i], 0)
                push!(qc.RepresentationFull[i], 0)
            end
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
function CU!(qc::QC_IT_MPS, U, pos::Array{Int64, 1})

    # check correct format of indices
    if length(pos) != 2
        error("Incorrect specification of indices (need 2 positions).")
    elseif pos[1] == pos[2]
        error("Incorrect specification of indices.")
    end

    # apply general CNOT via SWAP gates for linear topology
    if qc.LinearTopology == true

        #println("Currently not implemented.")

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
        for i in 1:qc.NumQubits
            if i == pos[1]
                push!(qc.Representation[i], 13) # control qubit
            elseif i == pos[2]
                push!(qc.Representation[i], -13) # action qubit
            elseif minimum(pos) < i && i < maximum(pos)
                push!(qc.Representation[i], 15) # vertical line
            else
                push!(qc.Representation[i], 0)
            end
        end

    # apply general CU via direct construction of matrix for master topology
    else

        """
        N = qc.NumQubits
        E, _, _, _  = get_Pauli_matrices()

        # projectors |0><0| and |1><1|
        proj0 = Complex.([1. 0.; 0. 0.])
        proj1 = Complex.([0. 0.; 0. 1.])

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

                # apply X gate in action space
                gate_tmp0 = kron(gate_tmp0, E)
                gate_tmp1 = kron(gate_tmp1, U)

                # fill with identities
                for i in pos[2]+1:N
                    gate_tmp0 = kron(gate_tmp0, E)
                    gate_tmp1 = kron(gate_tmp1, E)
                end

                # complete CU gate
                gate = gate_tmp0 + gate_tmp1

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

                # apply X gate in action space
                gate_tmp0 = kron(gate_tmp0, E)
                gate_tmp1 = kron(gate_tmp1, U)

                # fill with identities
                for i in pos[2]+1:N
                    gate_tmp0 = kron(gate_tmp0, E)
                    gate_tmp1 = kron(gate_tmp1, E)
                end

                # complete CU gate
                gate = gate_tmp0 + gate_tmp1
            end

        # control qubit "below" qubit that is acted on
        else

            # qubit 1 (in position 2) is action bit
            if pos[2] == 1
                gate_tmp0 = E
                gate_tmp1 = U

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

                # complete CU gate
                gate = gate_tmp0 + gate_tmp1

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

                gate_tmp0 = kron(gate_tmp0, E)
                gate_tmp1 = kron(gate_tmp1, U)

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

                # complete CU gate
                gate = gate_tmp0 + gate_tmp1
            end
        end
        """

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
        for i in 1:qc.NumQubits
            if i == pos[1]
                push!(qc.Representation[i], 13) # control qubit
                push!(qc.RepresentationFull[i], 13) # control qubit
            elseif i == pos[2]
                push!(qc.Representation[i], -13) # action qubit
                push!(qc.RepresentationFull[i], -13) # action qubit
            elseif minimum(pos) < i && i < maximum(pos)
                push!(qc.Representation[i], 15) # vertical line
                push!(qc.RepresentationFull[i], 15) # vertical line
            else
                push!(qc.Representation[i], 0)
                push!(qc.RepresentationFull[i], 0)
            end
        end

        # update circuit depth and bond dimension
        qc.CircuitDepth += 1
        push!(qc.BondDim, maxlinkdim(qc.StateVector))

    end
end

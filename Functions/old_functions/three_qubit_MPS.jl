
##############################################################################
## Three-qubit gates for Quantum Simulator - MPS ITensor version with full MPO
##############################################################################


####################
# Auxiliary functons
####################





""" Function to update the representation of a quantum circuit qc for a
two-site controlled gate (specified by num) applied to positions given through
the array pos. Can choose whether to update the full or the reduced representation
(or both). """
function update_representation_three_site!(qc::QC_IT_MPS, pos::Array{Int64, 1}, num::Int64, full::Bool, reduced::Bool)

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
function toffoli3!(qc::QC_IT_MPS, pos::Array{Int64, 1})

   # check correct format of indices and size of circuit
   if qc.NumQubits < 3
       error("Trying to apply a 3-qubit gate to less than 3 qubits!")
   elseif length(pos) != 3
       error("Incorrect specification of indices (need 3 positions).")
   elseif pos[3]-pos[2] != 1 && pos[2]-pos[1] != 1
       error("Incorrect specification of indices (need adjecent positions).")
   end

   # initialise TOFFOLI
   toffoli = OpSum()
   toffoli += "Proj00", pos[1], "Proj00", pos[2]
   toffoli += "Proj00", pos[1], "Proj11", pos[2]
   toffoli += "Proj11", pos[1], "Proj00", pos[2]
   toffoli += "Proj11", pos[1], "Proj11", pos[2], "X", pos[3]

   # get sites
   sites = qc.IndexSet

   # turn into MPO
   gate = MPO(toffoli, sites)

   # apply to state vector
   qc.StateVector = contract(gate, qc.StateVector, maxdim=qc.MaxBondDim)#, method=qc.ContractionMethod)

   # update representing matrix of quantum circuit
   update_representation_three_site!(qc, pos, 14, true, false)

   # update circuit depth and bond dimension
   qc.CircuitDepth += 1
   push!(qc.BondDim, maxlinkdim(qc.StateVector))
end


##################################
# Three-qubit gates (NOT adjacent)
##################################


""" Function to implement a controlled 2-site operator (on two adjacent
sites). """
function CU_2site!(qc::QC_IT_MPS, U, pos::Array{Int64, 1}, update_rep=true)

    # check correct format of indices and size of circuit
    if qc.NumQubits < 3
        error("Trying to apply a 3-qubit gate to less than 3 qubits!")
    elseif length(pos) != 3
        error("Incorrect specification of indices (need 3 positions).")
    elseif pos[3] - pos[2] != 1
        error("Incorrect specification of indices (last two indices must be adjacent).")
    end

    # apply general CNOT via SWAP gates for linear topology
    if qc.LinearTopology == true

        println("Linear topology currently not implemented!")

        # reference array of qubits
        #qubits = [i for i in 1:qc.NumQubits]

        # get permutation of qubits from position indication
        #perm = [i for i in 1:(minimum(pos)-1)]
        #append!(perm, pos)
        #append!(perm, setdiff(qubits, perm))

        # get non-trivial cycles that the permutation above defines
        #cycs = collect(cycles(Perm(perm)))
        #for cyc in cycs
        #    if length(cyc) == 1
        #        deleteat!(cycs, findall(x->x==cyc, cycs))
        #    end
        #end

        # convert into sequence of 2-swaps on neighbouring qubits
        #swaps = []
        #for cyc in cycs
        #    swap_tmp = get_2_perms_from_cycle(cyc)
        #    for sw in swap_tmp
        #        push!(swaps, sw)
        #    end
        #end

        # swap qubits, apply Toffoli, swap back
        #for sw in swaps
        #    swap2!(qc::QC_IT_MPS, sw)
        #end
        #toffoli3!(qc::QC_IT_MPS, [minimum(pos), minimum(pos)+1, minimum(pos)+2])
        #for sw in reverse(swaps)
        #    swap2!(qc::QC_IT_MPS, sw)
        #end

        # update representing matrix of quantum circuit
        #update_representation_three_site!(qc, pos, 14, false, true)
        #for i in 1:qc.NumQubits
        #    if i == pos[1]
        #        push!(qc.Representation[i], 14) # control qubit
        #    elseif i == pos[2]
        #        push!(qc.Representation[i], 14) # control qubit
        #    elseif i == pos[3]
        #        push!(qc.Representation[i], -14) # action qubit
        #    elseif minimum(pos) < i && i < maximum(pos)
        #        push!(qc.Representation[i], 15) # vertical line
        #    else
        #        push!(qc.Representation[i], 0)
        #    end
        #end

    # apply general Toffoli gate via direct construction of matrix for master topology
    else

        # initialise CU_2site
        cu_2site = OpSum()
        cu_2site += "Proj00", pos[1] # do nothing
        cu_2site += U[1,1], "Proj11", pos[1], "Proj00", pos[2], "Proj00", pos[3]
        cu_2site += U[1,2], "Proj11", pos[1], "Proj00", pos[2], "Proj01", pos[3]
        cu_2site += U[1,3], "Proj11", pos[1], "Proj01", pos[2], "Proj00", pos[3]
        cu_2site += U[1,4], "Proj11", pos[1], "Proj01", pos[2], "Proj01", pos[3]

        cu_2site += U[2,1], "Proj11", pos[1], "Proj00", pos[2], "Proj10", pos[3]
        cu_2site += U[2,2], "Proj11", pos[1], "Proj00", pos[2], "Proj11", pos[3]
        cu_2site += U[2,3], "Proj11", pos[1], "Proj01", pos[2], "Proj10", pos[3]
        cu_2site += U[2,4], "Proj11", pos[1], "Proj01", pos[2], "Proj11", pos[3]

        cu_2site += U[3,1], "Proj11", pos[1], "Proj10", pos[2], "Proj00", pos[3]
        cu_2site += U[3,2], "Proj11", pos[1], "Proj10", pos[2], "Proj01", pos[3]
        cu_2site += U[3,3], "Proj11", pos[1], "Proj11", pos[2], "Proj00", pos[3]
        cu_2site += U[3,4], "Proj11", pos[1], "Proj11", pos[2], "Proj01", pos[3]

        cu_2site += U[4,1], "Proj11", pos[1], "Proj10", pos[2], "Proj10", pos[3]
        cu_2site += U[4,2], "Proj11", pos[1], "Proj10", pos[2], "Proj11", pos[3]
        cu_2site += U[4,3], "Proj11", pos[1], "Proj11", pos[2], "Proj10", pos[3]
        cu_2site += U[4,4], "Proj11", pos[1], "Proj11", pos[2], "Proj11", pos[3]

        # get sites
        sites = qc.IndexSet

        # turn into MPO
        gate = MPO(cu_2site, sites)

        # apply to state vector, unprime
        qc.StateVector = contract(gate, qc.StateVector, maxdim=qc.MaxBondDim)#, method=qc.ContractionMethod)
        noprime!(qc.StateVector)

        # update representing matrix of quantum circuit
        if update_rep
            update_representation_three_site!(qc, pos, 13, true, true)
        end

        # update circuit depth and bond dimension
        qc.CircuitDepth += 1
        push!(qc.BondDim, maxlinkdim(qc.StateVector))
    end
end


""" Function to apply the Toffoli gate in an arbitrary configuration
of the input positions. For a linear topology, the function permutes
the input positions such that they are next to each other, and the
toffoli3! function is applied subsequently. For a master topology, the
matrix representing the operation is built from scratch. """
function toffoli!(qc::QC_IT_MPS, pos::Array{Int64, 1}, update_rep=true)

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
            swap2!(qc::QC_IT_MPS, sw)
        end
        toffoli3!(qc::QC_IT_MPS, [minimum(pos), minimum(pos)+1, minimum(pos)+2])
        for sw in reverse(swaps)
            swap2!(qc::QC_IT_MPS, sw)
        end

        # update representing matrix of quantum circuit
        update_representation_three_site!(qc, pos, 14, false, true)
        #for i in 1:qc.NumQubits
        #    if i == pos[1]
        #        push!(qc.Representation[i], 14) # control qubit
        #    elseif i == pos[2]
        #        push!(qc.Representation[i], 14) # control qubit
        #    elseif i == pos[3]
        #        push!(qc.Representation[i], -14) # action qubit
        #    elseif minimum(pos) < i && i < maximum(pos)
        #        push!(qc.Representation[i], 15) # vertical line
        #    else
        #        push!(qc.Representation[i], 0)
        #    end
        #end

    # apply general Toffoli gate via direct construction of matrix for master topology
    else

        # initialise TOFFOLI
        toffoli = OpSum()
        toffoli += "Proj00", pos[1], "Proj00", pos[2]
        toffoli += "Proj00", pos[1], "Proj11", pos[2]
        toffoli += "Proj11", pos[1], "Proj00", pos[2]
        toffoli += "Proj11", pos[1], "Proj11", pos[2], "X", pos[3]

        # get sites
        sites = qc.IndexSet

        # turn into MPO
        gate = MPO(toffoli, sites)

        # apply to state vector, unprime
        qc.StateVector = contract(gate, qc.StateVector, maxdim=qc.MaxBondDim)#, method=qc.ContractionMethod)
        noprime!(qc.StateVector)

        # update representing matrix of quantum circuit
        if update_rep
            update_representation_three_site!(qc, pos, 14, true, true)
        end

        # update circuit depth and bond dimension
        qc.CircuitDepth += 1
        push!(qc.BondDim, maxlinkdim(qc.StateVector))
    end
end

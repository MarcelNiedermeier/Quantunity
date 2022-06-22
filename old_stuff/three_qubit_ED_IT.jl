
###############################################################
## Three-qubit gates for Quantum Simulator - ED ITensor version
###############################################################


""" Function to apply the Toffoli gate to three adjacent qubits in positions
pos = [pos[1], pos[2], pos[3]] = [pos[1], pos[1]+1, pos[1]+2]. The "upper"
qubits (pos[1], pos[2]) are the control qubits, the "lower" qubit (pos[3])
the one that is acted on. This function serves mostly as an auxiliary function
in the more general toffoli! function, see below. """
function toffoli3!(qc::QC_IT_ED, pos::Array{Int64, 1})

   # check correct format of indices and size of circuit
   if qc.NumQubits < 3
       error("Trying to apply a 3-qubit gate to less than 3 qubits!")
   elseif length(pos) != 3
       error("Incorrect specification of indices (need 3 positions).")
   elseif pos[3]-pos[2] != 1 && pos[2]-pos[1] != 1
       error("Incorrect specification of indices (need adjecent positions).")
   end

   # get matrices
   E, _, _, _ = get_Pauli_matrices()
   toffoli = Complex.(Matrix{Float64}(I, 8, 8))
   toffoli[7:8, 7:8] = Complex.([0. 1.; 1. 0.])

   if pos[1] == 1
       gate = toffoli
       for i in 4:qc.NumQubits
           gate = kron(gate, E)
       end
   else
       gate = E
       for i in 2:pos[1]-1
           gate = kron(gate, E)
       end
       gate = kron(gate, toffoli)
       for i in (pos[3]+1):qc.NumQubits
           gate = kron(gate, E)
       end
   end

   # create a new set of outgoing indices and define ITensor from Hadamard array
   #sites_out = siteinds(d, N)
   #gate = ITensor(gate, sites, sites_out)
   #qc.StateVector .= gate*qc.StateVector

   # create a new set of outgoing indices and define ITensor from array
   sites_out = siteinds(2, N)
   gate = ITensor(gate, qc.IndexSet, sites_out)

   # update state vector and outgoing indices
   qc.StateVector = gate*qc.StateVector
   qc.IndexSet = sites_out

   # update representing matrix of quantum circuit
   for i in 1:qc.NumQubits
       if i == pos[1]
           #push!(qc.Representation[i], 12)
           push!(qc.RepresentationFull[i], 14) # control qubit
       elseif i == pos[2]
           push!(qc.RepresentationFull[i], -14) # action qubit
       else
           #push!(qc.Representation[i], 0)
           push!(qc.RepresentationFull[i], 0)
       end
   end

   # update circuit depth
   qc.CircuitDepth += 1
end


##################################
# three-qubit gates (not adjacent)
##################################


""" Function to apply the Toffoli gate in an arbitrary configuration
of the input positions. For a linear topology, the function permutes
the input positions such that they are next to each other, and the
toffoli3! function is applied subsequently. For a master topology, the
matrix representing the operation is built from scratch. """
function toffoli!(qc::QC_IT_ED, pos::Array{Int64, 1})

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
            swap2!(qc::QC_IT_ED, sw)
        end
        toffoli3!(qc::QC_IT_ED, [minimum(pos), minimum(pos)+1, minimum(pos)+2])
        for sw in reverse(swaps)
            swap2!(qc::QC_IT_ED, sw)
        end

        # update representing matrix of quantum circuit
        for i in 1:qc.NumQubits
            if i == pos[1]
                push!(qc.Representation[i], 14) # control qubit
            elseif i == pos[2]
                    push!(qc.Representation[i], 14) # control qubit
            elseif i == pos[3]
                push!(qc.Representation[i], -14) # action qubit
            elseif minimum(pos) < i && i < maximum(pos)
                push!(qc.Representation[i], 15) # vertical line
            else
                push!(qc.Representation[i], 0)
            end
        end

    # apply general Toffoli gate via direct construction of matrix for master topology
    else

        N = qc.NumQubits
        E, X, _, _  = get_Pauli_matrices()

        # projectors |0><0| and |1><1|
        proj0 = Complex.([1. 0.; 0. 0.])
        proj1 = Complex.([0. 0.; 0. 1.])

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

                # apply X gate in action space
                gate_tmp00 = kron(gate_tmp00, E)
                gate_tmp01 = kron(gate_tmp01, E)
                gate_tmp10 = kron(gate_tmp10, E)
                gate_tmp11 = kron(gate_tmp11, X)

                # fill with identities
                for i in pos[3]+1:N
                    gate_tmp00 = kron(gate_tmp00, E)
                    gate_tmp01 = kron(gate_tmp01, E)
                    gate_tmp10 = kron(gate_tmp10, E)
                    gate_tmp11 = kron(gate_tmp11, E)
                end

                # complete Toffoli gate
                gate = gate_tmp00 + gate_tmp01 + gate_tmp10 + gate_tmp11

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

                # apply X gate in action space
                gate_tmp00 = kron(gate_tmp00, E)
                gate_tmp01 = kron(gate_tmp01, E)
                gate_tmp10 = kron(gate_tmp10, E)
                gate_tmp11 = kron(gate_tmp11, X)

                # fill with identities
                for i in pos[3]+1:N
                    gate_tmp00 = kron(gate_tmp00, E)
                    gate_tmp01 = kron(gate_tmp01, E)
                    gate_tmp10 = kron(gate_tmp10, E)
                    gate_tmp11 = kron(gate_tmp11, E)
                end

                # complete Toffoli gate
                gate = gate_tmp00 + gate_tmp01 + gate_tmp10 + gate_tmp11
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
                gate_tmp11 = kron(gate_tmp11, X)


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

                # complete Toffoli gate
                gate = gate_tmp00 + gate_tmp01 + gate_tmp10 + gate_tmp11

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
                gate_tmp11 = kron(gate_tmp11, X)

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

                # complete Toffoli gate
                gate = gate_tmp00 + gate_tmp01 + gate_tmp10 + gate_tmp11
            end

        # action qubit above both control bits (e.g. pos = [3, 5, 1])
        else

            # action bit in first position 1
            if pos[3] == 1

                # initialise in action space
                gate_tmp00 = E
                gate_tmp01 = E
                gate_tmp10 = E
                gate_tmp11 = X

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

                # complete Toffoli gate
                gate = gate_tmp00 + gate_tmp01 + gate_tmp10 + gate_tmp11

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
                gate_tmp11 = kron(gate_tmp11, X)

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

                # complete Toffoli gate
                gate = gate_tmp00 + gate_tmp01 + gate_tmp10 + gate_tmp11
            end
        end

        # update state vector
        #qc.StateVector .= gate*qc.StateVector

        # create a new set of outgoing indices and define ITensor from array
        sites_out = siteinds(2, N)
        gate = ITensor(gate, qc.IndexSet, sites_out)

        # update state vector and outgoing indices
        qc.StateVector = gate*qc.StateVector
        qc.IndexSet = sites_out

        # update representing matrix of quantum circuit
        for i in 1:qc.NumQubits
            if i == pos[1]
                push!(qc.Representation[i], 14) # control qubit
                push!(qc.RepresentationFull[i], 14) # control qubit
            elseif i == pos[2]
                push!(qc.Representation[i], -14) # action qubit
                push!(qc.RepresentationFull[i], -14) # action qubit
            elseif minimum(pos) < i && i < maximum(pos)
                push!(qc.Representation[i], 15) # vertical line
                push!(qc.RepresentationFull[i], 15) # vertical line
            else
                push!(qc.Representation[i], 0)
                push!(qc.RepresentationFull[i], 0)
            end
        end

        # update circuit depth
        qc.CircuitDepth += 1
    end
end

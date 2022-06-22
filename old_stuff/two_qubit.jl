
########################################
## Two-qubit gates for Quantum Simulator
########################################

""" Function to apply the CNOT-gate to two adjacent qubits in positions
pos = [pos[1], pos[2]] = [pos[1], pos[1]+1]. The "upper" qubit (pos[1])
is the control qubit, the "lower" qubit (pos[2]) the one that is acted
on. This function serves mostly as an auxiliary function in the more
general cnot! function, see below. """
function cnot2!(qc::QC, pos::Array{Int64, 1})

   if length(pos) != 2
       error("Incorrect specification of indices (need 2 positions).")
   elseif abs(pos[1]-pos[2]) != 1
       error("Incorrect specification of indices (need adjacent positions).")
   end

   # get matrices
   E, _, _, _ = get_Pauli_matrices()
   cnot = Complex.([1. 0. 0. 0.; 0. 1. 0. 0.; 0. 0. 0. 1.; 0. 0. 1. 0.])

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

   # create a new set of outgoing indices and define ITensor from Hadamard array
   #sites_out = siteinds(d, N)
   #gate = ITensor(gate, sites, sites_out)
   #return gate*psi
   qc.StateVector .= gate*qc.StateVector

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

   # update circuit depth
   qc.CircuitDepth += 1
end

""" Function to swap two adjacent qubits in positions pos = [pos[1], pos[2]]
= [pos[1], pos[1]+1]. """
function swap2!(qc::QC, pos::Array{Int64, 1})

   # check structure of gate indices
   if length(pos) != 2
       error("Incorrect specification of indices (need 2 positions).")
   elseif abs(pos[1]-pos[2]) != 1
       error("Incorrect specification of indices (need adjacent positions).")
   end

   # get matrices
   E, _, _, _ = get_Pauli_matrices()
   swap = Complex.([1. 0. 0. 0.; 0. 0. 1. 0.; 0. 1. 0. 0.; 0. 0. 0. 1.])

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

   # create a new set of outgoing indices and define ITensor from Hadamard array
   #sites_out = siteinds(d, N)
   #gate = ITensor(gate, sites, sites_out)
   #return gate*psi
   #psi = gate*psi
   qc.StateVector .= gate*qc.StateVector

   # update representing matrix of quantum circuit
   for i in 1:qc.NumQubits
       if i in pos
           push!(qc.RepresentationFull[i], 4)
       else
           push!(qc.RepresentationFull[i], 0)
       end
   end

   # update circuit depth
   qc.CircuitDepth += 1
end

""" Function to apply a controlled-U gate to two adjacent qubits in positions
pos = [pos[1], pos[2]] = [pos[1], pos[1]+1], where U is a unitary single-qubit
gate. The "upper" qubit (pos[1]) is the control qubit, the "lower" qubit
(pos[2]) the one that is acted on. This function serves mostly as an auxiliary
function in the more general CU! function, see below. """
function CU2!(qc::QC, pos::Array{Int64, 1})

   # check if always gets two indices!
   if length(pos) != 2
       error("Incorrect specification of indices (need 2 positions).")
   elseif abs(pos[1]-pos[2]) != 1
       error("Incorrect specification of indices (need adjacent positions).")
   end

   # get matrices
   E, _, _, _ = get_Pauli_matrices()
   CU = Complex.(zeros(4,4))
   CU[1:2, 1:2] = E
   CU[3:4, 3:4] = U

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

   # create a new set of outgoing indices and define ITensor from Hadamard array
   #sites_out = siteinds(d, N)
   #gate = ITensor(gate, sites, sites_out)
   #return gate*psi
   qc.StateVector .= gate*qc.StateVector

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

   # update circuit depth
   qc.CircuitDepth += 1
end


################################
# two-qubit gates (not adjacent)
################################

""" Function to perfom a swap of two qubits over an arbitrary separation,
which is given by the array pos. pos[1] is understood to be smaller than
pos[2]. The swap is obtained through an equivalent sequence of 2-swaps. """
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

""" Function to reverse the arbitrary swap, performed by the swap! function."""
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


# the two functions above perform "half-swaps", still need function
# for full arbitrary swap!




""" Function to apply the CNOT gate to two arbitrary qubits, defined by
positions in pos = [pos[1], pos[2]]. The first qubits in is the control
qubit, the second the qubit that is acted on. The control qubit can be
"above" or "below" the action qubit, the function implements the desired
arrangement. """
function cnot!(qc::QC, pos::Array{Int64, 1})

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

        println("am here")
        println("position: ", pos)

        N = qc.NumQubits
        E, X, _, - = get_Pauli_matrices()

        # projectors |0><0| and |1><1|
        proj0 = Complex.([1. 0.; 0. 0.])
        proj1 = Complex.([0. 0.; 0. 1.])
        println(proj0)
        println(proj1)

        # control qubit "above" qubit that is acted on
        if pos[2] > pos[1]

            # qubit 1 is control bit
            if pos[1] == 1
                gate_tmp0 = proj0
                gate_tmp1 = proj1

                # fill with identities
                println("now here")
                println("position: ", pos)
                println(typeof(pos))
                for i in 2:pos[2]-1
                    println(i)
                end
                println("now here")
                println("position: ", pos)
                for i in 2:pos[2]-1
                    println("now here")
                    gate_tmp0 = kron(gate_tmp0, E)
                    gate_tmp1 = kron(gate_tmp1, E)
                end

                # apply X gate in action space
                gate_tmp0 = kron(gate_tmp0, E)
                gate_tmp1 = kron(gate_tmp1, X)

                # fill with identities
                for i in pos[2]+1:N
                    gate_tmp0 = kron(gate_tmp0, E)
                    gate_tmp1 = kron(gate_tmp1, E)
                end

                # complete CNOT gate
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

                gate_tmp0 = kron(gate_tmp0, proj0)
                gate_tmp1 = kron(gate_tmp1, proj1)

                # fill with identities
                for i in pos[1]:pos[2]-1
                    gate_tmp0 = kron(gate_tmp0, E)
                    gate_tmp1 = kron(gate_tmp1, E)
                end

                # apply X gate in action space
                gate_tmp0 = kron(gate_tmp0, E)
                gate_tmp1 = kron(gate_tmp1, X)

                # fill with identities
                for i in pos[2]+1:N
                    gate_tmp0 = kron(gate_tmp0, E)
                    gate_tmp1 = kron(gate_tmp1, E)
                end

                # complete CNOT gate
                gate = gate_tmp0 + gate_tmp1
            end

        # control qubit "below" qubit that is acted on
        else

            # qubit 1 is action bit
            if pos[1] == 1
                gate_tmp0 = E
                gate_tmp1 = X

                # fill with identities
                for i in 2:pos[2]-1
                    gate_tmp0 = kron(gate_tmp0, E)
                    gate_tmp1 = kron(gate_tmp1, E)
                end

                # apply projectors in control space
                gate_tmp0 = kron(gate_tmp0, proj0)
                gate_tmp1 = kron(gate_tmp1, proj1)

                # fill with identities
                for i in pos[2]+1:N
                    gate_tmp0 = kron(gate_tmp0, E)
                    gate_tmp1 = kron(gate_tmp1, E)
                end

                # complete CNOT gate
                gate = gate_tmp0 + gate_tmp1

            # qubit 1 is NOT action bit
            else

                # initialise with identities
                gate_tmp0 = E
                gate_tmp1 = E

                # fill with identities
                for i in 2:pos[1]-1
                    gate_tmp0 = kron(gate_tmp0, E)
                    gate_tmp1 = kron(gate_tmp1, E)
                end

                gate_tmp0 = kron(gate_tmp0, E)
                gate_tmp1 = kron(gate_tmp1, X)

                # fill with identities
                for i in pos[1]:pos[2]-1
                    gate_tmp0 = kron(gate_tmp0, E)
                    gate_tmp1 = kron(gate_tmp1, E)
                end

                # apply projectors in control space
                gate_tmp0 = kron(gate_tmp0, proj0)
                gate_tmp1 = kron(gate_tmp1, proj1)

                # fill with identities
                for i in pos[2]+1:N
                    gate_tmp0 = kron(gate_tmp0, E)
                    gate_tmp1 = kron(gate_tmp1, E)
                end

                # complete CNOT gate
                gate = gate_tmp0 + gate_tmp1
            end
        end

        # update statevector
        qc.StateVector .= gate*qc.StateVector

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

        # update cicruit depth
        qc.CircuitDepth += 1

    end
end

""" Function to apply the CU gate to two arbitrary qubits, defined by
positions in pos = [pos[1], pos[2]]. U is a unitary single-qubit gate.
The first qubits in is the control qubit, the second the qubit that is
acted on. The control qubit can be "above" or "below" the action qubit,
the function implements the desired arrangement. """
function CU!(qc::QC, pos::Array{Int64, 1}, U)

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
                CU2!(qc::QC, pos, U)

            # qubits not on adjacent lines
            else
                swap!(qc::QC, [pos[1], pos[2]-1])
                CU2!(qc::QC, [pos[2]-1, pos[2]], U)
                unswap!(qc::QC, [pos[1], pos[2]-1])
            end

        # control qubit "below" qubit that is acted on
        else

            # qubits on two adjacent lines
            if pos[1] == (pos[2]+1)
                swap2!(qc::QC, pos)
                CU2!(qc::QC, pos, U)
                swap2!(qc::QC, pos)

            # qubits on two adjacent lines
            else
                swap!(qc::QC, [pos[2], pos[1]])
                CU2!(qc::QC, [pos[1]-1, pos[1]], U)
                unswap!(qc::QC, [pos[2], pos[1]])
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

                gate_tmp0 = kron(gate_tmp0, proj0)
                gate_tmp1 = kron(gate_tmp1, proj1)

                # fill with identities
                for i in pos[1]:pos[2]-1
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

            # qubit 1 is action bit
            if pos[1] == 1
                gate_tmp0 = E
                gate_tmp1 = U

                # fill with identities
                for i in 2:pos[2]-1
                    gate_tmp0 = kron(gate_tmp0, E)
                    gate_tmp1 = kron(gate_tmp1, E)
                end

                # apply projectors in control space
                gate_tmp0 = kron(gate_tmp0, proj0)
                gate_tmp1 = kron(gate_tmp1, proj1)

                # fill with identities
                for i in pos[2]+1:N
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
                for i in 2:pos[1]-1
                    gate_tmp0 = kron(gate_tmp0, E)
                    gate_tmp1 = kron(gate_tmp1, E)
                end

                gate_tmp0 = kron(gate_tmp0, E)
                gate_tmp1 = kron(gate_tmp1, U)

                # fill with identities
                for i in pos[1]:pos[2]-1
                    gate_tmp0 = kron(gate_tmp0, E)
                    gate_tmp1 = kron(gate_tmp1, E)
                end

                # apply projectors in control space
                gate_tmp0 = kron(gate_tmp0, proj0)
                gate_tmp1 = kron(gate_tmp1, proj1)

                # fill with identities
                for i in pos[2]+1:N
                    gate_tmp0 = kron(gate_tmp0, E)
                    gate_tmp1 = kron(gate_tmp1, E)
                end

                # complete CU gate
                gate = gate_tmp0 + gate_tmp1
            end
        end

        # update statevector
        qc.StateVector .= gate*qc.StateVector

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

        # update circuit depth
        qc.CircuitDepth += 1

    end
end

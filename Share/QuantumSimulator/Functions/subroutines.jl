
###################################
## Subroutines for quantum circuits
###################################

# collect a few subroutines which are important building blocks for
# bigger quantum circuits


#####################
# Auxiliary functions
#####################


""" Function to add a "block" in the representation of a quantum circuit qc,
representing a given subroutine. """
function update_block_representation!(qc, pos, num, gatestring::String,
    rep_number::Number)

    # add information about the gate block in the conversion dictionary
    qc.ConversionTable[rep_number] = gatestring[1:4]
    qc.ConversionTable[rep_number+1] = gatestring[5:8]

    # update representation of quantum circuit
    for j in 1:2
        if j == 1
            for i in 1:qc.NumQubits
                if i == pos
                    push!(qc.Representation[i], rep_number)
                    push!(qc.RepresentationFull[i], rep_number)
                elseif i == pos+num-1
                    push!(qc.Representation[i], 22)
                    push!(qc.RepresentationFull[i], 22)
                elseif pos < i && i < pos+num-1
                    push!(qc.Representation[i], 20) # vertical line
                    push!(qc.RepresentationFull[i], 20) # vertical line
                else
                    push!(qc.Representation[i], 0)
                    push!(qc.RepresentationFull[i], 0)
                end
            end
        else
            for i in 1:qc.NumQubits
                if i == pos
                    push!(qc.Representation[i], rep_number+1)
                    push!(qc.RepresentationFull[i], rep_number+1)
                elseif i == pos+num-1
                    push!(qc.Representation[i], 23)
                    push!(qc.RepresentationFull[i], 23)
                elseif pos < i && i < pos+num-1
                    push!(qc.Representation[i], 21) # vertical line
                    push!(qc.RepresentationFull[i], 21) # vertical line
                else
                    push!(qc.Representation[i], 0)
                    push!(qc.RepresentationFull[i], 0)
                end
            end
        end
    end
end


########################
# Random Quantum Circuit
########################


""" Function to apply a sequence of gates generating a random quantum
state to the circuit qc, starting at position pos and covering a number
num of qubits. The sequence is of depth D. Every four layers of gates,
the state vector is renormalised. By default, the representation is changed
to a "block". """
function randomCircuit!(qc, pos, num, D, compact_rep=true)

    # check if size of subregister is compatible with circuit size
    if (pos+num-1) > qc.NumQubits
        error("You are trying to apply the transformation between qubits
        [$(pos), $(pos+num-1)], this however exceeds the number of qubits
        in the quantum circuit.")
    end

    # set representation updates to false if compact rep desired
    if compact_rep
        update_rep = false
    else
        update_rep = true
    end

    # build random circuit
    for j in 1:(D÷2)

        # apply random unitaries to each site
        #println("random single sites")
        random_single_site_gates!(qc, [i for i in pos:pos+num-1], update_rep)

        # apply vertically stacked CNOT gates
        if isodd(j)
            #println("applying sequential CNOTs between qubits $pos and $(pos+num-1)")
            sequential_cnot!(qc, pos, num, update_rep)
        else
            #println("applying sequential CNOTs between qubits $(pos+1) and $(pos+num-1)")
            sequential_cnot!(qc, pos+1, num-1, update_rep)
        end

        # renormalise every four layers
        if j%4 == 0
            #println("renormalising")
            qc.StateVector = (1/real(norm(qc.StateVector)))*qc.StateVector
        end
    end

    if compact_rep
        update_block_representation!(qc, pos, num, "-RANDOM-", 200)
    end

end


###########################
# Quantum Fourier Transform
###########################


""" Apply the QFT to num qubits of the quantum circuit qc, starting at
position pos. """
function QFT!(qc, pos, num, compact_rep=true)

    # check if size of subregister is compatible with circuit size
    if (pos+num-1) > qc.NumQubits
        error("You are trying to apply the transformation between qubits
        [$(pos), $(pos+num-1)], this however exceeds the number of qubits
        in the quantum circuit.")
    end

    # if compact representation desired: suppress representations of gate functions
    if compact_rep
        update_rep = false
    else
        update_rep = true
    end

    # Hadamard gates and controlled rotations
    for i in 1:num
        hadamard!(qc, [pos+i-1], update_rep)
        k = 0
        for j in pos+num-i:-1:pos+1
            #println("rotation by $(j-pos+1) for CRn from $(pos+num-1-k) to $(pos+i-1)")
            CRn!(qc, (j-pos+1), [pos+num-1-k, pos+i-1], update_rep)
            k = k+1
        end
    end

    # Swaps in the end
    for i in 0:(num÷2-1)
        fullSwap!(qc, [pos+i, pos+num-1-i], update_rep)
    end

    # update compact representation
    if compact_rep
        update_block_representation!(qc, pos, num, "-+QFT+--", 300)
    end

end


""" Apply the inverse QFT to num qubits of the quantum circuit qc, starting
at position pos. """
function invQFT!(qc, pos, num, compact_rep=true)

    # check if size of subregister is compatible with circuit size
    if (pos+num-1) > qc.NumQubits
        error("You are trying to apply the transformation between qubits
        [$(pos), $(pos+num-1)], this however exceeds the number of qubits
        in the quantum circuit.")
    end

    # if compact representation desired: suppress representations of gate functions
    if compact_rep
        update_rep = false
    else
        update_rep = true
    end

    # Swaps in the beginning
    for i in (num÷2-1):-1:0
        fullSwap!(qc, [pos+i, pos+num-1-i], update_rep)
    end

    # Hadamard gates and hermitian conjugate of controlled rotations
    for i in num:-1:1

        k = 2
        for j in pos+i:pos+num-1
            #println("rotation by -$k for CRn from $j to $(pos+i-1)")
            CRn!(qc, -k, [j, pos+i-1], update_rep)
            k = k+1
        end

        # initial code
        #for j in i+1:pos+num-1
        #    CRn!(qc, -(j-i+1), [j, i], update_rep)
        #    #CRn!(qc, -(j-i+1), [j, pos+i-1], update_rep)
        #end

        hadamard!(qc, [pos+i-1], update_rep)
    end

    # update compact representation
    if compact_rep
        update_block_representation!(qc, pos, num, "-invQFT-", 400)
    end

end


##########################
# Quantum Phase Estimation
##########################

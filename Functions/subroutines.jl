
###################################
## Subroutines for quantum circuits
###################################

# collect a few subroutines which are important building blocks for
# bigger quantum circuits


#####################
# Auxiliary functions
#####################


""" Function to add a "block" in the representation of a quantum circuit qc,
representing a given subroutine. pos """
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


######################
# Entanglement Entropy
######################


""" Function to compute the bipartite entanglement entropy (for every
bond) of the MPS representing the state of a quantum circuit. By
convention, the von Neumann entropy is evaluated. If α is set, the
corresponding Rényi entropy is calculated instead. A cutoff for the
SVD yielding the singular values at the bipartition can be selected. """
function entanglement_entropy(qc::QC_IT_MPS; α=1, cutoff=1E-3)

    # check input
    if α < 0
        error("Choose α as a positive real number.")
    end

    # get data of quantum circuit
    state = qc.StateVector
    N = qc.NumQubits

    # save entanglement entropies
    EEs = zeros(N-1)

    # loop through orthocenters of MPS, calculate entanglement entropy
    for b in 1:N-1
        orthogonalize!(state, b)
        if b == 1 # left boundary / upper qubit
            U,S,V = svd(state[b], (siteind(state,b)), cutoff=cutoff)
        else # bulk
            U,S,V = svd(state[b], (linkind(state, b-1), siteind(state,b)),
                cutoff=cutoff)
        end

        # von Neumann entropy
        if α == 1
            SvN = 0.0
            for n=1:ITensors.dim(S, 1)
                p = S[n,n]^2
                #println("p = $p")
                SvN -= p * log(p)
            end
            EEs[b] = SvN

        # Rényi entropy
        else # α ≠ 1
            SR = 0.0
            for n=1:ITensors.dim(S, 1)
                p = S[n,n]^2
                SR += p^α
            end
            EEs[b] = 1/(1-α) * log(SR)
        end
    end

    return EEs
end


#############
# Bell States
#############


""" Function to create the Bell state on two qubits in pos[1], pos[2]. The
states are defined by the integer num, as:
1 -> |00⟩ + |11⟩;
2 -> |00⟩ - |11⟩
3 -> |01⟩ + |10⟩
4 -> |01⟩ - |10⟩. """
function BellState!(qc, num, pos; compact_rep=true)

    if num ∉ [1, 2, 3, 4]
        error("Incorrect specification of Bell state.")
    elseif pos[2] < pos[1]
        error("Wrong order of positions.")
    elseif pos[2] > qc.NumQubits
        error("Given position exceeds number of qubits in circuit.")
    end

    # set representation updates to false if compact rep desired
    if compact_rep
        update_rep = false
    else
        update_rep = true
    end

    # first Bell state
    if num == 1
        Hadamard!(qc, [pos[1]], update_rep=update_rep)
        Cnot!(qc, [pos[1], pos[2]], update_rep=update_rep)

    # second Bell state
    elseif num == 2
        PauliX!(qc2, [pos[1]], update_rep=update_rep)
        Hadamard!(qc2, [pos[1]], update_rep=update_rep)
        Cnot!(qc2, [pos[1], pos[2]], update_rep=update_rep)

    # third Bell state
    elseif num == 3
        Hadamard!(qc2, [pos[1]], update_rep=update_rep)
        PauliX!(qc2, [pos[2]], update_rep=update_rep)
        Cnot!(qc2, [pos[1], pos[2]], update_rep=update_rep)

    # fourth Bell state
    else # num == 4
        Hadamard!(qc2, [pos[1]], update_rep=update_rep)
        PauliX!(qc2, [pos[2]], update_rep=update_rep)
        PauliZ!(qc2, [pos[1], pos[2]], update_rep=update_rep)
        Cnot!(qc2, [pos[1], pos[2]], update_rep=update_rep)

    end

    if compact_rep
        if num == 1
            update_block_representation!(qc, pos[1], pos[2], "-BELL-1-", 251)
        elseif num == 2
            update_block_representation!(qc, pos[1], pos[2], "-BELL-2-", 262)
        elseif num == 3
            update_block_representation!(qc, pos[1], pos[2], "-BELL-3-", 273)
        else # num == 4
            update_block_representation!(qc, pos[1], pos[2], "-BELL-4-", 284)
        end
    end
end


########################
# Random Quantum Circuit
########################


""" Function to apply a (different) random unitary single-site gate to a
quantum circuit in the positions given in the array pos. """
function random_single_site_gates!(qc::QC_IT_MPS, pos::Array; update_rep=true)

    # loop through given positions, generate random unitary and apply
    for i in pos

        # get random angles
        α = 2π * rand()
        β = 2π * rand()
        γ = 2π * rand()
        δ = 2π * rand()

        # sample random unitary 2x2 matrix, apply
        gate = op("RandU", qc.IndexSet[i]; α=α, β=β, γ=γ, δ=δ)
        updatedTensor = gate*qc.StateVector[i]
        noprime!(updatedTensor)
        qc.StateVector[i] = updatedTensor
    end

    # update representing matrix of quantum circuit
    if update_rep
      update_representation_single_site!(qc, pos, 19)
   end

    # update circuit depth and bond dimension
    qc.CircuitDepth += 1
    push!(qc.BondDim, maxlinkdim(qc.StateVector))

end


""" Function to apply a layer of CNOT gates on neighbouring qubit lines
to a quantum circuit. Can indicate the starting position. """
function cnot_layer!(qc::QC_IT_MPS, start_pos, end_pos; update_rep=false, recordEE=true,
    α=α, cutoff=cutoff)
    for i in start_pos:2:(end_pos-1)
        Cnot!(qc, [i, i+1], update_rep=update_rep, recordEE=recordEE,
            α=α, cutoff=cutoff)
    end
end


""" Function to apply a layer of CNOT gates on neighbouring qubit lines
to a quantum circuit. Can indicate the starting position. """
function Cnot_layer!(qc::Union{QC,QC_DM}, start_pos, end_pos; update_rep=false)
    for i in start_pos:2:(end_pos-1)
        Cnot!(qc, [i, i+1], update_rep=update_rep)
    end
end


""" Function to apply a sequence of gates generating a random quantum
state to the circuit qc, starting at position pos and covering a number
num of qubits. The sequence is of depth D, counting the number of 2-qubit
gates. Every four layers of gates, the state vector is renormalised. By
default, the representation is changed to a "block". """
function randomCircuit!(qc::QC_IT_MPS, pos, num, D; compact_rep=true, recordEE=false,
    α=1, cutoff=1E-3)

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
    for j in 1:D

        # apply random unitaries to each site
        random_single_site_gates!(qc, [i for i in pos:pos+num-1], update_rep=update_rep)

        # apply vertically stacked CNOT gates
        if isodd(j)
            cnot_layer!(qc, pos, pos+num-1, update_rep=update_rep,
                recordEE=recordEE, α=α, cutoff=cutoff)
        else
            cnot_layer!(qc, pos+1, pos+num-1, update_rep=update_rep,
                recordEE=recordEE, α=α, cutoff=cutoff)
        end
    end

    if compact_rep
        update_block_representation!(qc, pos, num, "-RANDOM-", 200)
    end
end


""" Function to apply a sequence of gates generating a random quantum
state to the circuit qc, starting at position pos and covering a number
num of qubits. The sequence is of depth D, counting the number of 2-qubit
gates. By default, the representation is changed to a "block". """
function randomCircuit!(qc::Union{QC,QC_DM}, pos, num, D; compact_rep=true)

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
    for j in 1:D

        # apply random unitaries to each site
        #random_single_site_gates!(qc, [i for i in pos:pos+num-1], update_rep=update_rep)
        for i in [i for i in pos:pos+num-1]
            RandomU!(qc, [i], update_rep=update_rep)
        end

        # apply vertically stacked CNOT gates
        if isodd(j)
            Cnot_layer!(qc, pos, pos+num-1, update_rep=update_rep)
        else
            Cnot_layer!(qc, pos+1, pos+num-1, update_rep=update_rep)
        end
    end

    if compact_rep
        update_block_representation!(qc, pos, num, "-RANDOM-", 200)
    end
end


""" Function to apply a sequence of gates generating a random quantum
state to the circuit qc, starting at position pos and covering a number
num of qubits. The sequence is of depth D. Every four layers of gates,
the state vector is renormalised. By default, the representation is changed
to a "block". """
function randomCircuit2!(qc, pos, num, D; compact_rep=true, recordEE=false,
    α=1, cutoff=1E-3)

    println("WARNING: Deprecated function.")

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
        random_single_site_gates!(qc, [i for i in pos:pos+num-1], update_rep=update_rep)
        #if recordEE
        #    push!(qc.EntanglementEntropy, entanglement_entropy(qc))
        #end

        # apply vertically stacked CNOT gates
        if isodd(j)
            #println("applying sequential CNOTs between qubits $pos and $(pos+num-1)")
            sequential_cnot!(qc, pos, num, update_rep)
        else
            #println("applying sequential CNOTs between qubits $(pos+1) and $(pos+num-1)")
            sequential_cnot!(qc, pos+1, num-1, update_rep)
        end
        if recordEE
            push!(qc.EntanglementEntropy, entanglement_entropy(qc, α=α, cutoff=cutoff))
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
function QFT!(qc::Union{QC,QC_DM}, pos, num; compact_rep=true, no_rep=false)

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
        Hadamard!(qc, [pos+i-1], update_rep=update_rep)

        k = 0
        for j in pos+num-i:-1:pos+1
            #println("rotation by $(j-pos+1) for CRn from $(pos+num-1-k) to $(pos+i-1)")
            CRn!(qc, (j-pos+1), [pos+num-1-k, pos+i-1], update_rep=update_rep)
            k = k+1
        end
    end

    # Swaps in the end
    for i in 0:(num÷2-1)
        #fullSwap!(qc, [pos+i, pos+num-1-i], update_rep)
        Swap!(qc, [pos+i, pos+num-1-i], update_rep=update_rep)
    end

    # update compact representation
    if compact_rep && no_rep==false
        update_block_representation!(qc, pos, num, "-+QFT+--", 300)
    end

end


""" Apply the QFT to num qubits of the quantum circuit qc, starting at
position pos. """
function QFT!(qc::QC_IT_MPS, pos, num; compact_rep=true, no_rep=false, recordEE=false,
    α=1, cutoff=1E-3)

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
        Hadamard!(qc, [pos+i-1], update_rep=update_rep)

        k = 0
        for j in pos+num-i:-1:pos+1
            #println("rotation by $(j-pos+1) for CRn from $(pos+num-1-k) to $(pos+i-1)")
            CRn!(qc, (j-pos+1), [pos+num-1-k, pos+i-1], update_rep=update_rep,
                recordEE=recordEE, α=α, cutoff=cutoff)
            k = k+1
        end
    end

    # Swaps in the end
    for i in 0:(num÷2-1)
        #fullSwap!(qc, [pos+i, pos+num-1-i], update_rep)
        Swap!(qc, [pos+i, pos+num-1-i], update_rep=update_rep,
            recordEE=recordEE, α=α, cutoff=cutoff)
    end

    # update compact representation
    if compact_rep && no_rep==false
        update_block_representation!(qc, pos, num, "-+QFT+--", 300)
    end

end


""" Apply the AQFT to num qubits of the quantum circuit qc, starting at
position pos. m is the cutoff parameter: only controlled rotations up
to R_m will be applied (therefore needs to be smaller than subregister
size). """
function AQFT!(qc::Union{QC,QC_DM}, m, pos, num; compact_rep=true, no_rep=false)

    # check if size of subregister is compatible with circuit size
    if (pos+num-1) > qc.NumQubits
        error("You are trying to apply the transformation between qubits
        [$(pos), $(pos+num-1)], this however exceeds the number of qubits
        in the quantum circuit.")
    elseif m > num
        error("Incorrect choice for parameter m (needs to be smaller
        than subregister size).)")
    end

    # if compact representation desired: suppress representations of gate functions
    if compact_rep
        update_rep = false
    else
        update_rep = true
    end

    # Hadamard gates and controlled rotations
    for i in 1:num
        Hadamard!(qc, [pos+i-1], update_rep=update_rep)
        k = 0
        for j in pos+num-i:-1:pos+1
            #println("rotation by $(j-pos+1) for CRn from $(pos+num-1-k) to $(pos+i-1)")
            if (j-pos+1)  <= m
                #println("rotation by $(j-pos+1) for CRn from $(pos+num-1-k) to $(pos+i-1)")
                CRn!(qc, (j-pos+1), [pos+num-1-k, pos+i-1], update_rep=update_rep)
            end
            k = k+1
        end
    end

    # Swaps in the end
    for i in 0:(num÷2-1)
        #fullSwap!(qc, [pos+i, pos+num-1-i], update_rep)
        Swap!(qc, [pos+i, pos+num-1-i], update_rep=update_rep)
    end

    # update compact representation
    if compact_rep && no_rep==false
        update_block_representation!(qc, pos, num, "-AQFT+--", 300)
    end
end


""" Apply the AQFT to num qubits of the quantum circuit qc, starting at
position pos. m is the cutoff parameter: only controlled rotations up
to R_m will be applied (therefore needs to be smaller than subregister
size). """
function AQFT!(qc::QC_IT_MPS, m, pos, num; compact_rep=true, no_rep=false, recordEE=false,
    α=1, cutoff=1E-3)

    # check if size of subregister is compatible with circuit size
    if (pos+num-1) > qc.NumQubits
        error("You are trying to apply the transformation between qubits
        [$(pos), $(pos+num-1)], this however exceeds the number of qubits
        in the quantum circuit.")
    elseif m > num
        error("Incorrect choice for parameter m (needs to be smaller
        than subregister size).)")
    end

    # if compact representation desired: suppress representations of gate functions
    if compact_rep
        update_rep = false
    else
        update_rep = true
    end

    # Hadamard gates and controlled rotations
    for i in 1:num
        Hadamard!(qc, [pos+i-1], update_rep=update_rep)
        k = 0
        for j in pos+num-i:-1:pos+1
            #println("rotation by $(j-pos+1) for CRn from $(pos+num-1-k) to $(pos+i-1)")
            if (j-pos+1)  <= m
                #println("rotation by $(j-pos+1) for CRn from $(pos+num-1-k) to $(pos+i-1)")
                CRn!(qc, (j-pos+1), [pos+num-1-k, pos+i-1], update_rep=update_rep,
                    recordEE=recordEE, α=α, cutoff=cutoff)
            end
            k = k+1
        end
    end

    # Swaps in the end
    for i in 0:(num÷2-1)
        #fullSwap!(qc, [pos+i, pos+num-1-i], update_rep)
        Swap!(qc, [pos+i, pos+num-1-i], update_rep=update_rep,
            recordEE=recordEE, α=α, cutoff=cutoff)
    end

    # update compact representation
    if compact_rep && no_rep==false
        update_block_representation!(qc, pos, num, "-AQFT+--", 300)
    end
end


""" Apply the inverse QFT to num qubits of the quantum circuit qc, starting
at position pos. """
function invQFT!(qc::QC_IT_MPS, pos, num; compact_rep=true, no_rep=false, recordEE=false,
    α=1, cutoff=1E-3)

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
        #fullSwap!(qc, [pos+i, pos+num-1-i], update_rep)
        Swap!(qc, [pos+i, pos+num-1-i], update_rep=update_rep,
            recordEE=recordEE, α=α, cutoff=cutoff)
    end

    # Hadamard gates and hermitian conjugate of controlled rotations
    for i in num:-1:1

        k = 2
        for j in pos+i:pos+num-1
            #println("rotation by -$k for CRn from $j to $(pos+i-1)")
            CRn!(qc, -k, [j, pos+i-1], update_rep=update_rep,
                recordEE=recordEE, α=α, cutoff=cutoff)
            k = k+1
        end

        # initial code
        #for j in i+1:pos+num-1
        #    CRn!(qc, -(j-i+1), [j, i], update_rep)
        #    #CRn!(qc, -(j-i+1), [j, pos+i-1], update_rep)
        #end

        Hadamard!(qc, [pos+i-1], update_rep=update_rep)
    end

    # update compact representation
    if compact_rep && no_rep==false
        update_block_representation!(qc, pos, num, "-invQFT-", 400)
    end
end


""" Apply the inverse QFT to num qubits of the quantum circuit qc, starting
at position pos. """
function invQFT!(qc::Union{QC,QC_DM}, pos, num; compact_rep=true, no_rep=false)

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
        #fullSwap!(qc, [pos+i, pos+num-1-i], update_rep)
        Swap!(qc, [pos+i, pos+num-1-i], update_rep=update_rep)
    end

    # Hadamard gates and hermitian conjugate of controlled rotations
    for i in num:-1:1

        k = 2
        for j in pos+i:pos+num-1
            #println("rotation by -$k for CRn from $j to $(pos+i-1)")
            CRn!(qc, -k, [j, pos+i-1], update_rep=update_rep)
            k = k+1
        end

        # initial code
        #for j in i+1:pos+num-1
        #    CRn!(qc, -(j-i+1), [j, i], update_rep)
        #    #CRn!(qc, -(j-i+1), [j, pos+i-1], update_rep)
        #end

        Hadamard!(qc, [pos+i-1], update_rep=update_rep)
    end

    # update compact representation
    if compact_rep && no_rep==false
        update_block_representation!(qc, pos, num, "-invQFT-", 400)
    end
end


""" Apply the inverse AQFT to num qubits of the quantum circuit qc, starting
at position pos. m is the cutoff parameter: only controlled rotations up
to R_m will be applied (therefore needs to be smaller than subregister
size)."""
function invAQFT!(qc::Union{QC,QC_DM}, m, pos, num; compact_rep=true, no_rep=false)

    # check if size of subregister is compatible with circuit size
    if (pos+num-1) > qc.NumQubits
        error("You are trying to apply the transformation between qubits
        [$(pos), $(pos+num-1)], this however exceeds the number of qubits
        in the quantum circuit.")
    elseif m > num
        error("Incorrect choice for parameter m (needs to be smaller
        than subregister size).)")
    end

    # if compact representation desired: suppress representations of gate functions
    if compact_rep
        update_rep = false
    else
        update_rep = true
    end

    # Swaps in the beginning
    for i in (num÷2-1):-1:0
        #fullSwap!(qc, [pos+i, pos+num-1-i], update_rep)
        Swap!(qc, [pos+i, pos+num-1-i], update_rep=update_rep)
    end

    # Hadamard gates and hermitian conjugate of controlled rotations
    for i in num:-1:1

        k = 2
        for j in pos+i:pos+num-1
            #println("rotation by -$k for CRn from $j to $(pos+i-1)")
            if k  <= m
                CRn!(qc, -k, [j, pos+i-1], update_rep=update_rep)
            end
            k = k+1
        end

        Hadamard!(qc, [pos+i-1], update_rep=update_rep)
    end

    # update compact representation
    if compact_rep && no_rep==false
        update_block_representation!(qc, pos, num, "-inAQFT-", 400)
    end
end


""" Apply the inverse AQFT to num qubits of the quantum circuit qc, starting
at position pos. m is the cutoff parameter: only controlled rotations up
to R_m will be applied (therefore needs to be smaller than subregister
size)."""
function invAQFT!(qc::QC_IT_MPS, m, pos, num; compact_rep=true, no_rep=false,
    recordEE=false, α=1, cutoff=1E-3)

    # check if size of subregister is compatible with circuit size
    if (pos+num-1) > qc.NumQubits
        error("You are trying to apply the transformation between qubits
        [$(pos), $(pos+num-1)], this however exceeds the number of qubits
        in the quantum circuit.")
    elseif m > num
        error("Incorrect choice for parameter m (needs to be smaller
        than subregister size).)")
    end

    # if compact representation desired: suppress representations of gate functions
    if compact_rep
        update_rep = false
    else
        update_rep = true
    end

    # Swaps in the beginning
    for i in (num÷2-1):-1:0
        #fullSwap!(qc, [pos+i, pos+num-1-i], update_rep)
        Swap!(qc, [pos+i, pos+num-1-i], update_rep=update_rep,
            recordEE=recordEE, α=α, cutoff=cutoff)
    end

    # Hadamard gates and hermitian conjugate of controlled rotations
    for i in num:-1:1

        k = 2
        for j in pos+i:pos+num-1
            #println("rotation by -$k for CRn from $j to $(pos+i-1)")
            if k  <= m
                CRn!(qc, -k, [j, pos+i-1], update_rep=update_rep,
                    recordEE=recordEE, α=α, cutoff=cutoff)
            end
            k = k+1
        end

        Hadamard!(qc, [pos+i-1], update_rep=update_rep)
    end

    # update compact representation
    if compact_rep && no_rep==false
        update_block_representation!(qc, pos, num, "-inAQFT-", 400)
    end
end


##########################
# Quantum Phase Estimation
##########################


""" Calculate U^n in computational basis with dummy unitary map that has
eigenvalues θ1 and θ2 for the eigenvectors |+⟩, |-⟩. """
function U_n(θ1, θ2, n)
    H = 1/√2 * Complex.([1. 1.; 1. -1.])
    U_tmp = Complex.([exp(1.0im*2π*n*θ1) 0.; 0. exp(1.0im*2π*n*θ2)])
    return H*U_tmp*H
end


""" Calculate U^n in computational basis with dummy unitary 2-site map that has
eigenvalues θ1, θ2, θ3 and θ4 for the Bell states as eigenvectors. """
function U_n_2site(θ1, θ2, θ3, θ4, n)
    T = 1/√2 * Complex.([1. 1. 0. 0.; 0. 0. 1. 1.; 0. 0. 1. -1.; 1. -1. 0. 0.])
    T_inv = 1/√2 * Complex.([1. 0. 0. 1.; 1. 0. 0. -1.; 0. 1. 1. 0.; 0. 1. -1. 0.])
    U_tmp = Complex.([exp(1.0im*2π*n*θ1) 0. 0. 0.; 0. exp(1.0im*2π*n*θ2) 0. 0.; 0. 0. exp(1.0im*2π*n*θ3) 0.; 0. 0. 0. exp(1.0im*2π*n*θ4)])
    return T*U_tmp*T_inv
end


""" Calculate estimated phase in QPE algorithm, given a measurement result
of qubits. """
function recover_phase_estimate(phase_array::Array)
    N = length(phase_array)
    phase = 0.
    for j in 1:N
        #phase += phase_array[end-j+1] * 2.0^(-N-1+j)
        phase += (1/2)^(j) * phase_array[j]
        #println("$(phase_array[j]), $((1/2)^(j))")
    end
    return phase
end


""" Function to find the measurement result with the highest proportion
in the whole sample. Returns the bitstring and its corresponding
probability. """
function get_highest_prob_measurement(qc)

    # find index of highest probability
    p_max_ind = argmax(collect(values(qc.ClassicalBitsProportion)))

    # find highest probability and corresponding bitstring
    state_max = collect(keys(qc.ClassicalBitsProportion))[p_max_ind]
    p_max = collect(values(qc.ClassicalBitsProportion))[p_max_ind]

    return state_max, p_max
end


""" Function to read out the k measurements with the highest probabilities. """
function get_highest_prob_measurement(qc, k::Int)

    # get probabilities and obtain list of k highest values
    probas_tmp = collect(values(qc.ClassicalBitsProportion))
    sort!(probas_tmp, rev=true)
    if length(probas_tmp) >= k
        max_probas = probas_tmp[1:k]
    else
        max_probas = probas_tmp[:]
    end

    # find indices of k highest probabilities and find corresponding bitstrings
    ind = findall(x -> x ∈ max_probas, collect(values(qc.ClassicalBitsProportion)))
    measurements_max = collect(keys(qc.ClassicalBitsProportion))[ind]

    return measurements_max, max_probas
end


""" Function to get the phase from a QPE circuit after measurement. """
function QPE_get_phase(qc, k=1)

    # get highest prob measured bitstrings
    states_max, probs_max = get_highest_prob_measurement(qc, k)

    # safeguard in case less than k phases are actually measured
    if length(states_max) >= k
        num_phases = k
    else
        num_phases = length(states_max)
    end

    # convert bitstrings to phases, rescale to [0, 2π]
    phases = zeros(num_phases)
    for i in 1:num_phases
        phases[i] = 2π * recover_phase_estimate(states_max[i])
    end

    return phases, probs_max
end


""" Function to apply a quantum phase estimation of operator U to a
specified part of a quantum circuit. Can define precision of estimated
phase. Preparation of input state must be taken care of outside of this
subroutine."""
function QPE!(qc, U, pos, num_bin_digits; compact_rep=true, α=1,
    cutoff=1E-3, recordEE=recordEE)

    println("WARNING: this function needs to be updated to account for an
        operator in sub-circuit structure.")

    # get number of qubits for operator
    num_qubits_U = Int(log2(size(U)[1]))

    # check if size of subregister is compatible with circuit size
    if (pos+num_bin_digits+num_qubits_U-1) > qc.NumQubits
        error("You are trying to apply the transformation between qubits
        [$(pos), $(pos+num_bin_digits+num_qubits_U-1)], this
        however exceeds the number of qubits in the quantum circuit.")
    end

    # if compact representation desired: suppress representations of gate functions
    if compact_rep
        update_rep = false
    else
        update_rep = true
    end

    # initialise phase estimation register
    Hadamard!(qc, [i for i in pos:(pos+num_bin_digits-1)], update_rep=update_rep)

    # do controlled rotations
    local n = 1
    #n = 1
    for i in (pos+num_bin_digits-1):-1:(pos)
        println("QPE, n = $n")
        #println("i = $i")
        #println("action qubits: $([j for j in (pos+num_bin_digits):
        #    (pos+num_bin_digits+num_qubits_U-1)])")
        # construct controlled operator
        CU_general!(qc, U^n, [i], [j for j in (pos+num_bin_digits):
            (pos+num_bin_digits+num_qubits_U-1)], update_rep)
        n = 2*n
    end

    # do inverse QFT (switch off representation)
    no_rep = true
    invQFT!(qc, pos, num_bin_digits, compact_rep=compact_rep, no_rep=no_rep,
        α=α, cutoff=cutoff, recordEE=recordEE)

    # update compact representation
    if compact_rep
        update_block_representation!(qc, pos, num_bin_digits+num_qubits_U, "-+QPE+--", 500)
    end
end


""" Function to apply an inverse quantum phase estimation of operator U to a
specified part of a quantum circuit. Can define precision of estimated
phase. Need for example in HHL algorithm. """
function invQPE!(qc, U, pos, num_bin_digits; compact_rep=true, α=1,
    cutoff=1E-3, recordEE=recordEE)

    println("WARNING: this function needs to be updated to account for an
        operator in sub-circuit structure.")

    # get number of qubits for operator
    num_qubits_U = Int(log2(size(U)[1]))

    # check if size of subregister is compatible with circuit size
    if (pos+num_bin_digits+num_qubits_U-1) > qc.NumQubits
        error("You are trying to apply the transformation between qubits
        [$(pos), $(pos+num_bin_digits+num_qubits_U-1)], this
        however exceeds the number of qubits in the quantum circuit.")
    end

    # if compact representation desired: suppress representations of gate functions
    if compact_rep
        update_rep = false
    else
        update_rep = true
    end

    # start with QFT (switch off representation)
    no_rep = true
    QFT!(qc, pos, num_bin_digits, compact_rep=compact_rep, no_rep=no_rep,
        α=α, cutoff=cutoff, recordEE=recordEE)

    # do controlled rotations
    local n = 2^(num_bin_digits-1)
    #n = 2^(num_bin_digits-1)
    for i in pos:(pos+num_bin_digits-1)
        #println("i = $i")
        #println("inv QPE, n = $n")
        #println("action qubits: $([j for j in (pos+num_bin_digits):
        #    (pos+num_bin_digits+num_qubits_U-1)])")
        # construct controlled operator
        CU_general!(qc, U^n, [i], [j for j in (pos+num_bin_digits):
            (pos+num_bin_digits+num_qubits_U-1)], update_rep)
        n = n÷2
    end

    # finalise with Hadamards ("reverse initialisation")
    Hadamard!(qc, [i for i in pos:(pos+num_bin_digits-1)], update_rep=update_rep)

    # update compact representation
    if compact_rep
        update_block_representation!(qc, pos, num_bin_digits+num_qubits_U, "-invQPE-", 600)
    end
end


###############
# Grover Search
###############


""" Function implementing the Grover diffusion operator H(2|s⟩⟨s| - I)H.
Operator starts at qubit line given by start_pos and extends over
num_qubits qubits. """
function Grover_diffusor!(qc, start_pos, num_qubits; recordEE=false,
    compact_rep=false, α=1, cutoff=1E-3)

    # if compact representation desired: suppress representations of gate functions
    if compact_rep
        update_rep = false
    else
        update_rep = true
    end

    # register over which Grover diffusor is applied
    register = [i for i in start_pos:(start_pos+num_qubits-1)]

    # wrap in Hadamards
    Hadamard!(qc, register, update_rep=update_rep)

    # 2|0⟩⟨0| - I
    PauliX!(qc, register, update_rep=update_rep)

    action_qubits = [start_pos+num_qubits-1]
    control_register = [i for i in start_pos:(start_pos+num_qubits-2)]
    C_PauliZ!(qc, control_register, action_qubits, update_rep=update_rep,
        recordEE=recordEE, α=α, cutoff=cutoff)

    PauliX!(qc, register, update_rep=update_rep)

    # wrap in Hadamards
    Hadamard!(qc, register, update_rep=update_rep)

    # update compact representation
    if compact_rep
        update_block_representation!(qc, start_pos, num_qubits, "-+Diff--", 900)
    end

end


""" Function implementing the Grover diffusion operator H(2|s⟩⟨s| - I)H.
Operator starts at qubit line given by start_pos and extends over
num_qubits qubits. It is controlled by an additional qubit. """
function C_Grover_diffusor!(qc, control_qubit, start_pos,
    num_qubits; compact_rep=false, recordEE=false, α=1, cutoff=1E-3)

    if length(control_qubit) > 1
        error("Too many control qubits given (only one).")
    end

    # if compact representation desired: suppress representations of gate functions
    if compact_rep
        update_rep = false
    else
        update_rep = true
    end

    # register over which Grover diffusor is applied
    register = [i for i in start_pos:(start_pos+num_qubits-1)]

    # wrap in controlled Hadamards
    for i in 1:length(register)
        C_Hadamard!(qc, control_qubit, [register[i]], update_rep=update_rep,
            recordEE=recordEE, α=α, cutoff=cutoff)
    end

    # 2|0⟩⟨0| - I
    for i in 1:length(register)
        C_PauliX!(qc, control_qubit, [register[i]], update_rep=update_rep,
            recordEE=recordEE, α=α, cutoff=cutoff)
    end

    action_qubits = [start_pos+num_qubits-1]
    control_register = copy(control_qubit)
    for i in start_pos:(start_pos+num_qubits-2)
        push!(control_register, i)
    end
    C_PauliZ!(qc, control_register, action_qubits, update_rep=update_rep,
        recordEE=recordEE, α=α, cutoff=cutoff)

    for i in 1:length(register)
        C_PauliX!(qc, control_qubit, [register[i]], update_rep=update_rep,
            recordEE=recordEE, α=α, cutoff=cutoff)
    end

    # wrap in Hadamards
    for i in 1:length(register)
        C_Hadamard!(qc, control_qubit, [register[i]], update_rep=update_rep,
            recordEE=recordEE, α=α, cutoff=cutoff)
    end

end


""" Function to mark a single computational basis element in a
Grover-type circuit. E.g. if |001⟩ is supposed to be marked,
its sign in some superposition state is inverted whereas the
operator acts as an identity on all other computational basis
states. """
function phase_oracle_single!(qc, marked_element, start_pos; recordEE=false,
    compact_rep=false, α=1, cutoff=1E-3)

    # if compact representation desired: suppress representations of gate functions
    if compact_rep
        update_rep = false
    else
        update_rep = true
    end

    num_qubits = length(marked_element)

    # invert list representing marked element to map to qubits
    #marked_register = reverse(marked_element)
    marked_register = marked_element # reversing doesn't seem to be necessary

    # construct phase oracle
    Pauli_register = []
    for i in 1:length(marked_register)
        if marked_register[i] == 0
            push!(Pauli_register, start_pos+i-1)
        end
    end

    PauliX!(qc, Pauli_register, update_rep=update_rep)

    control_register = [i for i in start_pos:(start_pos+num_qubits-2)]
    action_qubits = [start_pos+num_qubits-1]
    C_PauliZ!(qc, control_register, action_qubits, update_rep=update_rep,
        recordEE=recordEE, α=α, cutoff=cutoff)

    PauliX!(qc, Pauli_register, update_rep=update_rep)

    # update compact representation
    if compact_rep
        update_block_representation!(qc, start_pos, num_qubits, "-Phase--", 800)
    end

end


""" Function to mark a single computational basis element in a
Grover-type circuit. E.g. if |001⟩ is supposed to be marked,
its sign in some superposition state is inverted whereas the
operator acts as an identity on all other computational basis
states. It is controlled by an additional qubit. """
function C_phase_oracle_single!(qc, control_qubit, marked_element,
    start_pos; compact_rep=false, recordEE=false, α=1, cutoff=1E-3)

    # if compact representation desired: suppress representations of gate functions
    if compact_rep
        update_rep = false
    else
        update_rep = true
    end

    num_qubits = length(marked_element)

    # invert list representing marked element to map to qubits
    #marked_register = reverse(marked_element)
    marked_register = marked_element # reversing doesn't seem to be necessary

    # construct phase oracle
    Pauli_register = []
    for i in 1:length(marked_register)
        if marked_register[i] == 0
            push!(Pauli_register, start_pos+i-1)
        end
    end

    for i in 1:length(Pauli_register)
        C_PauliX!(qc, control_qubit, [Pauli_register[i]], update_rep=update_rep,
            recordEE=recordEE, α=α, cutoff=cutoff)
    end

    action_qubits = [start_pos+num_qubits-1]
    control_register = copy(control_qubit)
    for i in start_pos:(start_pos+num_qubits-2)
        push!(control_register, i)
    end
    C_PauliZ!(qc, control_register, action_qubits, update_rep=update_rep,
        recordEE=recordEE, α=α, cutoff=cutoff)

    #PauliX!(qc, Pauli_register)
    for i in 1:length(Pauli_register)
        C_PauliX!(qc, control_qubit, [Pauli_register[i]], update_rep=update_rep,
            recordEE=recordEE, α=α, cutoff=cutoff)
    end
end


""" Function to mark multiple computational basis elements in a
Grover-type circuit. Should be given a list of marked elements. """
function phase_oracle_multiple!(qc, marked_elements, start_pos; recordEE=false,
    compact_rep=false, α=1, cutoff=1E-3)

    # apply separate phase oracle for each marked element
    for marked_element in marked_elements
        phase_oracle_single!(qc, marked_element, start_pos, recordEE=recordEE,
            compact_rep=compact_rep, α=α, cutoff=cutoff)
    end
end


""" Function to mark multiple computational basis elements in a
Grover-type circuit. Should be given a list of marked elements.
It is controlled by an additional qubit. """
function C_phase_oracle_multiple!(qc, control_qubit, marked_elements,
    start_pos; recordEE=false, compact_rep=false, α=1, cutoff=1E-3)

    # apply separate phase oracle for each marked element
    for marked_element in marked_elements
        C_phase_oracle_single!(qc, control_qubit, marked_element, start_pos,
            recordEE=recordEE, compact_rep=compact_rep, α=α, cutoff=cutoff)
    end
end


""" Function to apply the whole Grover operator to a quantum circuit. """
function Grover_operator!(qc, marked_elements, start_pos; recordEE=false,
    compact_rep=false, α=1, cutoff=1E-3)

    # controlled phase oracle for one or more marked elements
    if length(marked_elements) == 1
        phase_oracle_single!(qc, marked_elements[1], start_pos,
            recordEE=recordEE, compact_rep=compact_rep, α=α, cutoff=cutoff)
    else
        phase_oracle_multiple!(qc, marked_elements, start_pos,
            recordEE=recordEE, compact_rep=compact_rep, α=α, cutoff=cutoff)
    end

    # controlled Grover diffusor
    num_qubits = length(marked_elements[1])
    Grover_diffusor!(qc, start_pos, num_qubits, recordEE=recordEE,
        compact_rep=compact_rep, α=α, cutoff=cutoff)

    # update compact representation
    #if compact_rep
    #    update_block_representation!(qc, pos, num_qubits, "-Grover-", 800)
    #end
end


""" Function to apply the whole Grover operator to a quantum circuit,
with one control qubit (useful for quantum counting!). """
function C_Grover_operator!(qc, control_qubit, marked_elements, start_pos;
    recordEE=false, compact_rep=false, α=1, cutoff=1E-3)

    # controlled phase oracle for one or more marked elements
    if length(marked_elements) == 1
        C_phase_oracle_single!(qc, control_qubit, marked_elements[1], start_pos,
            recordEE=recordEE, compact_rep=compact_rep, α=α, cutoff=cutoff)
    else
        C_phase_oracle_multiple!(qc, control_qubit, marked_elements, start_pos,
            recordEE=recordEE, compact_rep=compact_rep, α=α, cutoff=cutoff)
    end

    # controlled Grover diffusor
    num_qubits = length(marked_elements[1])
    C_Grover_diffusor!(qc, control_qubit, start_pos, num_qubits,
        recordEE=recordEE, compact_rep=compact_rep, α=α, cutoff=cutoff)
end


##################
# Shor's Algorithm
##################


""" Function to find how many qubits are required to factorise the number
N with Short algorithm (in the upper register, twice the number of binary
digits). """
function qubits_needed(N)
    q = 0
    for i in 1:100
        if N^2 <= 2^i && 2^i <= 2*N^2
            #println("you need $i qubits")
            q = i
        end
    end
    println("You need $q qubits: N² = $(N^2), 2^q = $(2^q), 2N² = $(2*N^2), number states ≈ $(2^q/N)")
    return q
end


""" Function to find the order r of x mod N (x^r mod N = 1). """
function find_order(x, N)

    # numbers not co-prime
    if gcd(BigInt(x), N) != 1
        return 0
    end
    for i in 1:N
        if mod(BigInt(x)^i, N) == 1
            return i
        end
    end
end


""" Function to find the largest order for some x < N (returns all
values of x which have the largest order). """
function find_largest_order(N)
    orders = zeros(N)
    xs = [i for i in 1:N]
    for x in xs
        orders[x] = find_order(x, N)
    end
    max_order = maximum(orders)
    x_max = findall(x->x==max_order, orders)
    return x_max, max_order
end


""" Function to find all the states (as decimal numbers) which occur
in a superposition in the upper register in Shor's algorithm, for a
given N, x, order and value of a to collapse the lower register. """
function find_upper_register(q, x, N, a, max_order)

    #chosen_a = rand(0:2^(q)-1)
    chosen_a = a
    #println("chosen a ", chosen_a)
    z = mod(BigInt(x)^chosen_a, BigInt(N))
    println("measured z = $z")
    a_vals_temp = []
    for a in 0:2^(2*q)-1
        #println("a = $a")
        if mod(BigInt(x)^a, BigInt(N)) == z
            push!(a_vals_temp, a)
        end
        if length(a_vals_temp) == 3
            break
        end
    end
    #println(a_vals_temp)
    l = a_vals_temp[1]
    r = a_vals_temp[2] - a_vals_temp[1]
    #println("l = $l")
    #println("r = $r")

    a_vals = []
    A = ceil(2^q/max_order)+1
    #println("A = $A")
    for j in 0:A
        if l+j*r <= (2^q)-1
            push!(a_vals, Int64(l+j*r))
        end
    end
    return a_vals
end


""" Function to prepare the state of the upper register in Shor's algorithm
(after the controlled operations) as an MPS, given that the lower register has
collapsed it into a specific superposition (set by giving a value of a). It
is automatically assumed that for a given N, a value of x is chosen such that
the order of x is maximal (and therefore the number of states in the resulting
superposition minimal). """
function prepare_Shor_superposition(N, a, backend)

    # find the largest order for the given integer and number of qubits
    x_max, max_order = find_largest_order(N)
    q = qubits_needed(N)
    println("period: ", max_order)
    x = x_max[1]
    println("x with max order selected: $x")
    #println("smallest number of states in Fourier register: $((2^q)/max_order)")

    # find the states in the upper register after measurement
    println("Measuring upper register")
    a_states = find_upper_register(q, x, N, a, max_order)
    println("States: ", a_states)

    # convert a_states into initial state for upper register
    println("Constructing superposition state ($((2^q)/max_order) elements)")
    psi = initial_state_arbitrary_superposition(a_states, q, backend)

    return psi, q
end


########################
# Hamiltonian Simulation
########################


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


#############
# Bell States
#############


""" Function to create the Bell state on two qubits in pos[1], pos[2]. The
states are defined by the integer num, as:
1 -> |00⟩ + |11⟩;
2 -> |00⟩ - |11⟩
3 -> |01⟩ + |10⟩
4 -> |01⟩ - |10⟩. """
function BellState!(qc, num::Int, pos, compact_rep=true)

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
        hadamard!(qc, [pos[1]], update_rep)
        cnot!(qc, [pos[1], pos[2]], update_rep)

    # second Bell state
    elseif num == 2
        PauliX!(qc2, [pos[1]], update_rep)
        hadamard!(qc2, [pos[1]], update_rep)
        cnot!(qc2, [pos[1], pos[2]], update_rep)

    # third Bell state
    elseif num == 3
        hadamard!(qc2, [pos[1]], update_rep)
        PauliX!(qc2, [pos[2]], update_rep)
        cnot!(qc2, [pos[1], pos[2]], update_rep)

    # fourth Bell state
    else # num == 4
        hadamard!(qc2, [pos[1]], update_rep)
        PauliX!(qc2, [pos[2]], update_rep)
        PauliZ!(qc2, [pos[1], pos[2]], update_rep)
        cnot!(qc2, [pos[1], pos[2]], update_rep)

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
function QFT!(qc, pos, num, compact_rep=true, no_rep=false)

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
    if compact_rep && no_rep==false
        update_block_representation!(qc, pos, num, "-+QFT+--", 300)
    end

end


""" Apply the inverse QFT to num qubits of the quantum circuit qc, starting
at position pos. """
function invQFT!(qc, pos, num, compact_rep=true, no_rep=false)

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
    if compact_rep && no_rep==false
        update_block_representation!(qc, pos, num, "-invQFT-", 400)
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
    p_max = maximum(keys(qc.ClassicalBitsProportion))
    return qc.ClassicalBitsProportion[p_max], p_max
end


function get_highest_prob_measurement(qc, k::Int)

    println(keys(qc.ClassicalBitsProportion)[1])
    keys_tmp = Array(keys(qc.ClassicalBitsProportion))
    sort!(keys_tmp, rev=true)

    max_probas = keys_temp[1:k]
    measurements_max = []

    for i in 1:k
        append!(measurements_max, qc.ClassicalBitsProportion[max_probas[i]])
    end

    return measurements_max, max_probas

end

""" Function to apply a quantum phase estimation of operator U to a
specified part of a quantum circuit. Can define precision of estimated
phase. Preparation of input state must be taken care of outside of this
subroutine."""
function QPE!(qc, U, pos, num_bin_digits, compact_rep=true)

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
    hadamard!(qc, [i for i in pos:(pos+num_bin_digits-1)], update_rep)

    # do controlled rotations
    local n = 1
    for i in (pos+num_bin_digits-1):-1:(pos)
        println("i = $i")
        println("action qubits: $([j for j in (pos+num_bin_digits):
            (pos+num_bin_digits+num_qubits_U-1)])")
        # construct controlled operator
        CU_general!(qc, U^n, [i], [j for j in (pos+num_bin_digits):
            (pos+num_bin_digits+num_qubits_U-1)], update_rep)
        n = 2*n
    end

    # do inverse QFT (switch off representation)
    no_rep = true
    invQFT!(qc, pos, num_bin_digits, compact_rep, no_rep)

    # update compact representation
    if compact_rep
        update_block_representation!(qc, pos, num_bin_digits+num_qubits_U, "-+QPE+--", 500)
    end
end

""" Function to apply an inverse quantum phase estimation of operator U to a
specified part of a quantum circuit. Can define precision of estimated
phase. Need for example in HHL algorithm. """
function invQPE!(qc, U, pos, num_bin_digits, compact_rep=true)

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
    QFT!(qc, pos, num_bin_digits, compact_rep, no_rep)

    # do controlled rotations
    local n = 2^(num_bin_digits-1)
    for i in pos:(pos+num_bin_digits-1)
        println("i = $i")
        println("n = $n")
        println("action qubits: $([j for j in (pos+num_bin_digits):
            (pos+num_bin_digits+num_qubits_U-1)])")
        # construct controlled operator
        CU_general!(qc, U^n, [i], [j for j in (pos+num_bin_digits):
            (pos+num_bin_digits+num_qubits_U-1)], update_rep)
        n = n÷2
    end

    # finalise with Hadamards ("reverse initialisation")
    hadamard!(qc, [i for i in pos:(pos+num_bin_digits-1)], update_rep)

    # update compact representation
    if compact_rep
        update_block_representation!(qc, pos, num_bin_digits+num_qubits_U, "-invQPE-", 600)
    end
end


########################
# Hamiltonian Simulation
########################

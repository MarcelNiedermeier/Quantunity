
########################
## Custom Gate Functions
########################


#####################
# Auxiliary functions
#####################


""" Function to apply a custom gate at the desired position. pos is
the index of the first qubit of the total qubit register that the
custom gate acts on. """
function apply_custom_gate(qc, CGate::custom_gate, pos::Int64)

    for qubit in collect(keys(CGate.Representation))
        #println("qubit: $(qubit)")
        for gate in CGate.Representation[qubit]
            #println("gate: $(gate)")

            # gates with 4 parameters
            if length(gate) == 5
                UGate!(qc, [qubit+pos-1], gate[2], gate[3], gate[4], gate[5], false)

            # gates with 1 parameter
            elseif length(gate) == 2
                if gate[1] == 9.
                    RXGate!(qc, [qubit+pos-1], gate[2], false)
                elseif gate[1] == 10.
                    RYGate!(qc, [qubit+pos-1], gate[2], false)
                elseif gate[1] == 11.
                    RZGate!(qc, [qubit+pos-1], gate[2], false)
                elseif gate[1] == 12.
                    PhaseShift!(qc, [qubit+pos-1], gate[2], false)
                end

            # fixed gates
            elseif length(gate) == 1
                if gate[1] == 1.
                    hadamard!(qc, [qubit+pos-1], false)
                elseif gate[1] == 2.
                    PauliX!(qc, [qubit+pos-1], false)
                elseif gate[1] == 6.
                    PauliY!(qc, [qubit+pos-1], false)
                elseif gate[1] == 5.
                    PauliZ!(qc, [qubit+pos-1], false)
                elseif gate[1] == 18.
                    SqrtX!(qc, [qubit+pos-1], false)
                elseif gate[1] == 7.
                    SGate!(qc, [qubit+pos-1], false)
                elseif gate[1] == 8.
                    TGate!(qc, [qubit+pos-1], false)
                end

            # two-site gates
            elseif length(gate) == 3
                if gate[3] == 1.
                    #println("here", [Int64(gate[1]), Int64(gate[2])])
                    cnot!(qc, [Int64(gate[1]), Int64(gate[2])], false)
                end
            end
        end
    end

    # update other circuit quantities
    qc.HasCustomGate = true

    # generate name string
    if length(CGate.GateName) == 1
        gatestring = "-++"*CGate.GateName[1]*"++--"
    elseif length(CGate.GateName) == 2
        gatestring = "-+"*CGate.GateName[1:2]*"++--"
    elseif length(CGate.GateName) == 3
        gatestring = "-+"*CGate.GateName[1:3]*"+--"
    elseif length(CGate.GateName) == 4
        gatestring = "-"*CGate.GateName[1:4]*"+--"
    else
        gatestring = "-"*CGate.GateName*"--"
    end

    qc.ConversionTable[100] = gatestring[1:4]
    qc.ConversionTable[101] = gatestring[5:8]

    # update representation of quantum circuit (general "block")
    for j in 1:2
        if j == 1
            for i in 1:qc.NumQubits
                if i == pos
                    push!(qc.Representation[i], 100)
                    push!(qc.RepresentationFull[i], 100)
                elseif i == pos+CGate.NumQubitsGate-1
                    push!(qc.Representation[i], 22)
                    push!(qc.RepresentationFull[i], 22)
                elseif pos < i && i < pos+CGate.NumQubitsGate-1
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
                    push!(qc.Representation[i], 101)
                    push!(qc.RepresentationFull[i], 101)
                elseif i == pos+CGate.NumQubitsGate-1
                    push!(qc.Representation[i], 23)
                    push!(qc.RepresentationFull[i], 23)
                elseif pos < i && i < pos+CGate.NumQubitsGate-1
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


""" Function to update the representation of a custom gate CGate for a
gate (specified by num) applied  in positions given through the array pos.
Needs to be supplied with data about the gate, i.e. the gate type and possible
parameters, which are saved in the representation dictionary as a list. """
function update_representation_custom_gate!(CGate::custom_gate, pos::Array{Int64, 1},
    gate_info::Array, two_site=false)

    if two_site == false
       for i in 1:CGate.NumQubitsGate
           if i in pos
               push!(CGate.Representation[i], gate_info)
           else
               push!(CGate.Representation[i], [0.])
           end
       end
    else
       for i in 1:CGate.NumQubitsGate
           if i == pos[1]
                push!(CGate.Representation[i], gate_info)
           else
                push!(CGate.Representation[i], [0.])
           end
       end
    end
end


####################
# Single-qubit gates
####################


""" Function to apply the (single-qubit) Hadamard gate to a custom gate CGate
in every position specified in the array pos. """
function hadamard!(CGate::custom_gate, pos::Array{Int64, 1})
   update_representation_custom_gate!(CGate, pos, [1.])
end


""" Function to apply the (single-qubit) Pauli X gate to a custom gate CGate
in every position specified in the array pos. """
function PauliX!(CGate::custom_gate, pos::Array{Int64, 1})
   update_representation_custom_gate!(CGate, pos, [2.])
end


""" Function to apply the (single-qubit) Pauli Y gate to a custom gate CGate
in every position specified in the array pos. """
function PauliY!(CGate::custom_gate, pos::Array{Int64, 1})
   update_representation_custom_gate!(CGate, pos, [6.])
end


""" Function to apply the (single-qubit) Pauli Z gate to acustom gate CGate
in every position specified in the array pos. """
function PauliZ!(CGate::custom_gate, pos::Array{Int64, 1})
   update_representation_custom_gate!(CGate, pos, [5.])
end


""" Function to apply the (single-qubit) √X gate to a custom gate CGate
in every position specified in the array pos. """
function SqrtX!(CGate::custom_gate, pos::Array{Int64, 1})
   update_representation_custom_gate!(CGate, pos, [18.])
end


""" Function to apply the (single-qubit) S gate to a custom gate CGate
in every position specified in the array pos. """
function SGate!(CGate::custom_gate, pos::Array{Int64, 1})
   update_representation_custom_gate!(CGate, pos, [7.])
end


""" Function to apply the (single-qubit) T gate to a custom gate CGate
in every position specified in the array pos. """
function TGate!(CGate::custom_gate, pos::Array{Int64, 1})
   update_representation_custom_gate!(CGate, pos, [8.])
end


""" Function to apply the (single-qubit) Rx gate to a custom gate CGate
in every position specified in the array pos. Performs a rotation around
the x-axis specified by theta."""
function RXGate!(CGate::custom_gate, pos::Array{Int64, 1}, theta::Number)
   update_representation_custom_gate!(CGate, pos, [9., theta])
end


""" Function to apply the (single-qubit) Ry gate to a custom gate CGate
in every position specified in the array pos. Performs a rotation around
the x-axis specified by theta."""
function RYGate!(CGate::custom_gate, pos::Array{Int64, 1}, theta::Number)
   update_representation_custom_gate!(CGate, pos, [10., theta])
end


""" Function to apply the (single-qubit) Rz gate to a custom gate CGate
in every position specified in the array pos. Performs a rotation around
the x-axis specified by theta."""
function RZGate!(CGate::custom_gate, pos::Array{Int64, 1}, theta::Number)
   update_representation_custom_gate!(CGate, pos, [11., theta])
end


""" Function to apply a phase shift by the angle theta to a custom gate CGate
in every position specified in the array pos. """
function PhaseShift!(CGate::custom_gate, pos::Array{Int64, 1}, theta::Number)
   update_representation_custom_gate!(CGate, pos, [12., theta])
end


""" Function to apply custom-defined unitary operator U to a custom gate CGate
in every position specified in the array pos. The matrix is defined as
U = [α β; γ δ] with the corresponding parameters. """
function UGate!(CGate::custom_gate, pos::Array{Int64, 1}, α::Number,
   β::Number, γ::Number, δ::Number)
   update_representation_custom_gate!(CGate, pos, [19., α, β, γ, δ])
end


#################
# Two-qubit gates
#################


""" Function to apply the CNOT gate to a custom gate CGate in positions
specified in the array pos. """
function cnot!(CGate::custom_gate, pos::Array{Int64, 1})
   update_representation_custom_gate!(CGate, pos, [Float64(pos[1]),
      Float64(pos[2]), 1.], true)
end


###################
# multi-qubit gates
###################

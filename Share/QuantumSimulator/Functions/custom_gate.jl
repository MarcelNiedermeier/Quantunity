
########################
## Custom Gate Functions
########################


#####################
# Auxiliary functions
#####################


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

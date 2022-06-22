
################################################################
## Single qubit gates for Quantum Simulator - ED ITensor version
################################################################


""" Function to apple the (single-qubit) Hadamard gate to a quantum circuit qc
in every position specified in the array pos. """
function hadamard!(qc::QC_IT_ED, pos::Array{Int64, 1})

   # get matrices
   h = 1/sqrt(2) * Complex.([1. 1.; 1. -1.])
   E, _, _, _ = get_Pauli_matrices()

   # apply single qubit gates at indicated positions, otherwise apply identity
   sites_out = siteinds(2, N)
   for i in 1:qc.NumQubits
       if i in pos
           gate = ITensor(h, qc.IndexSet[i], sites_out[i])
       else
           gate = ITensor(E, qc.IndexSet[i], sites_out[i])
       end
       qc.StateVector = gate*qc.StateVector
       qc.IndexSet[i] = sites_out[i]
   end

   # create a new set of outgoing indices and define ITensor from array
   #sites_out = siteinds(2, N)
   #gate = ITensor(gate, qc.IndexSet, sites_out)

   # update state vector and outgoing indices
   #qc.StateVector = gate*qc.StateVector
   #qc.IndexSet = sites_out

   # update representing matrix of quantum circuit
   for i in 1:qc.NumQubits
       if i in pos
           push!(qc.Representation[i], 1)
           push!(qc.RepresentationFull[i], 1)
       else
           push!(qc.Representation[i], 0)
           push!(qc.RepresentationFull[i], 0)
       end
   end

   # update circuit depth
   qc.CircuitDepth += 1
end


""" Function to apple the (single-qubit) Pauli X gate to a quantum circuit qc
in every position specified in the array pos. """
function PauliX!(qc::QC_IT_ED, pos::Array{Int64, 1})

   # get matrices
   E, X, _, _ = get_Pauli_matrices()

   # apply single qubit gates at indicated positions, otherwise apply identity
   sites_out = siteinds(2, N)
   for i in 1:qc.NumQubits
       if i in pos
           gate = ITensor(X, qc.IndexSet[i], sites_out[i])
       else
           gate = ITensor(E, qc.IndexSet[i], sites_out[i])
       end
       qc.StateVector = gate*qc.StateVector
       qc.IndexSet[i] = sites_out[i]
   end

   # create a new set of outgoing indices and define ITensor from array
   #sites_out = siteinds(2, N)
   #gate = ITensor(gate, qc.IndexSet, sites_out)

   # update state vector and outgoing indices
   #qc.StateVector = gate*qc.StateVector
   #qc.IndexSet = sites_out


   # update representing matrix of quantum circuit
   for i in 1:qc.NumQubits
       if i in pos
           push!(qc.Representation[i], 2)
           push!(qc.RepresentationFull[i], 2)
       else
           push!(qc.Representation[i], 0)
           push!(qc.RepresentationFull[i], 0)
       end
   end

   # update circuit depth
   qc.CircuitDepth += 1
end


""" Function to apple the (single-qubit) Pauli Y gate to a quantum circuit qc
in every position specified in the array pos. """
function PauliY!(qc::QC_IT_ED, pos::Array{Int64, 1})

   # get matrices
   E, _, Y, _ = get_Pauli_matrices()

   # apply single qubit gates at indicated positions, otherwise apply identity
   sites_out = siteinds(2, N)
   for i in 1:qc.NumQubits
       if i in pos
           gate = ITensor(Y, qc.IndexSet[i], sites_out[i])
       else
           gate = ITensor(E, qc.IndexSet[i], sites_out[i])
       end
       qc.StateVector = gate*qc.StateVector
       qc.IndexSet[i] = sites_out[i]
   end

   # create a new set of outgoing indices and define ITensor from array
   #sites_out = siteinds(2, N)
   #gate = ITensor(gate, qc.IndexSet, sites_out)

   # update state vector and outgoing indices
   #qc.StateVector = gate*qc.StateVector
   #qc.IndexSet = sites_out

   # update representing matrix of quantum circuit
   for i in 1:qc.NumQubits
       if i in pos
           push!(qc.Representation[i], 6)
           push!(qc.RepresentationFull[i], 6)
       else
           push!(qc.Representation[i], 0)
           push!(qc.RepresentationFull[i], 0)
       end
   end

   # update circuit depth
   qc.CircuitDepth += 1
end


""" Function to apple the (single-qubit) Pauli Z gate to a quantum circuit qc
in every position specified in the array pos. """
function PauliZ!(qc::QC_IT_ED, pos::Array{Int64, 1})

   # get matrices
   E, _, _, Z = get_Pauli_matrices()

   # apply single qubit gates at indicated positions, otherwise apply identity
   sites_out = siteinds(2, N)
   for i in 1:qc.NumQubits
       if i in pos
           gate = ITensor(Z, qc.IndexSet[i], sites_out[i])
       else
           gate = ITensor(E, qc.IndexSet[i], sites_out[i])
       end
       qc.StateVector = gate*qc.StateVector
       qc.IndexSet[i] = sites_out[i]
   end

   # create a new set of outgoing indices and define ITensor from array
   #sites_out = siteinds(2, N)
   #gate = ITensor(gate, qc.IndexSet, sites_out)

   # update state vector and outgoing indices
   #qc.StateVector = gate*qc.StateVector
   #qc.IndexSet = sites_out

   # update representing matrix of quantum circuit
   for i in 1:qc.NumQubits
       if i in pos
           push!(qc.Representation[i], 5)
           push!(qc.RepresentationFull[i], 5)
       else
           push!(qc.Representation[i], 0)
           push!(qc.RepresentationFull[i], 0)
       end
   end

   # update circuit depth
   qc.CircuitDepth += 1
end


""" Function to apple the (single-qubit) S gate to a quantum circuit qc
in every position specified in the array pos. """
function SGate!(qc::QC_IT_ED, pos::Array{Int64, 1})

   # get matrices
   E, _, _, _ = get_Pauli_matrices()
   S = Complex.([1. 0.; 0. 1.0im])

   # apply single qubit gates at indicated positions, otherwise apply identity
   sites_out = siteinds(2, N)
   for i in 1:qc.NumQubits
       if i in pos
           gate = ITensor(S, qc.IndexSet[i], sites_out[i])
       else
           gate = ITensor(E, qc.IndexSet[i], sites_out[i])
       end
       qc.StateVector = gate*qc.StateVector
       qc.IndexSet[i] = sites_out[i]
   end

   # create a new set of outgoing indices and define ITensor from array
   #sites_out = siteinds(2, N)
   #gate = ITensor(gate, qc.IndexSet, sites_out)

   # update state vector and outgoing indices
   #qc.StateVector = gate*qc.StateVector
   #qc.IndexSet = sites_out

   # update representing matrix of quantum circuit
   for i in 1:qc.NumQubits
       if i in pos
           push!(qc.Representation[i], 7)
           push!(qc.RepresentationFull[i], 7)
       else
           push!(qc.Representation[i], 0)
           push!(qc.RepresentationFull[i], 0)
       end
   end

   # update circuit depth
   qc.CircuitDepth += 1
end


""" Function to apple the (single-qubit) T gate to a quantum circuit qc
in every position specified in the array pos. """
function TGate!(qc::QC_IT_ED, pos::Array{Int64, 1})

   # get matrices
   E, _, _, _ = get_Pauli_matrices()
   T = Complex.([1. 0.; 0. exp(1.0im*pi/4)])

   # apply single qubit gates at indicated positions, otherwise apply identity
   sites_out = siteinds(2, N)
   for i in 1:qc.NumQubits
       if i in pos
           gate = ITensor(T, qc.IndexSet[i], sites_out[i])
       else
           gate = ITensor(E, qc.IndexSet[i], sites_out[i])
       end
       qc.StateVector = gate*qc.StateVector
       qc.IndexSet[i] = sites_out[i]
   end

   # create a new set of outgoing indices and define ITensor from array
   #sites_out = siteinds(2, N)
   #gate = ITensor(gate, qc.IndexSet, sites_out)

   # update state vector and outgoing indices
   #qc.StateVector = gate*qc.StateVector
   #qc.IndexSet = sites_out

   # update representing matrix of quantum circuit
   for i in 1:qc.NumQubits
       if i in pos
           push!(qc.Representation[i], 8)
           push!(qc.RepresentationFull[i], 8)
       else
           push!(qc.Representation[i], 0)
           push!(qc.RepresentationFull[i], 0)
       end
   end

   # update circuit depth
   qc.CircuitDepth += 1
end


""" Function to apple the (single-qubit) Rx gate to a quantum circuit qc
in every position specified in the array pos. Performs a rotation around
the x-axis specified by theta."""
function RXGate!(qc::QC_IT_ED, pos::Array{Int64, 1}, theta)

   # get matrices
   E, X, _, _ = get_Pauli_matrices()
   RX = exp(-1.0im * X * theta/2)

   # apply single qubit gates at indicated positions, otherwise apply identity
   sites_out = siteinds(2, N)
   for i in 1:qc.NumQubits
       if i in pos
           gate = ITensor(RX, qc.IndexSet[i], sites_out[i])
       else
           gate = ITensor(E, qc.IndexSet[i], sites_out[i])
       end
       qc.StateVector = gate*qc.StateVector
       qc.IndexSet[i] = sites_out[i]
   end

   # create a new set of outgoing indices and define ITensor from array
   #sites_out = siteinds(2, N)
   #gate = ITensor(gate, qc.IndexSet, sites_out)

   # update state vector and outgoing indices
   #qc.StateVector = gate*qc.StateVector
   #qc.IndexSet = sites_out

   # update representing matrix of quantum circuit
   for i in 1:qc.NumQubits
       if i in pos
           push!(qc.Representation[i], 9)
           push!(qc.RepresentationFull[i], 9)
       else
           push!(qc.Representation[i], 0)
           push!(qc.RepresentationFull[i], 0)
       end
   end

   # update circuit depth
   qc.CircuitDepth += 1
end


""" Function to apple the (single-qubit) Ry gate to a quantum circuit qc
in every position specified in the array pos. Performs a rotation around
the x-axis specified by theta."""
function RYGate!(qc::QC_IT_ED, pos::Array{Int64, 1}, theta)

   # get matrices
   E, _, Y, _ = get_Pauli_matrices()
   RY = exp(-1.0im * Y * theta/2)

   # apply single qubit gates at indicated positions, otherwise apply identity
   sites_out = siteinds(2, N)
   for i in 1:qc.NumQubits
       if i in pos
           gate = ITensor(RY, qc.IndexSet[i], sites_out[i])
       else
           gate = ITensor(E, qc.IndexSet[i], sites_out[i])
       end
       qc.StateVector = gate*qc.StateVector
       qc.IndexSet[i] = sites_out[i]
   end

   # create a new set of outgoing indices and define ITensor from array
   #sites_out = siteinds(2, N)
   #gate = ITensor(gate, qc.IndexSet, sites_out)

   # update state vector and outgoing indices
   #qc.StateVector = gate*qc.StateVector
   #qc.IndexSet = sites_out

   # update representing matrix of quantum circuit
   for i in 1:qc.NumQubits
       if i in pos
           push!(qc.Representation[i], 10)
           push!(qc.RepresentationFull[i], 10)
       else
           push!(qc.Representation[i], 0)
           push!(qc.RepresentationFull[i], 0)
       end
   end

   # update circuit depth
   qc.CircuitDepth += 1
end


""" Function to apple the (single-qubit) Rz gate to a quantum circuit qc
in every position specified in the array pos. Performs a rotation around
the x-axis specified by theta."""
function RZGate!(qc::QC_IT_ED, pos::Array{Int64, 1}, theta)

   # get matrices
   E, _, _, Z = get_Pauli_matrices()
   RZ = exp(-1.0im * Z * theta/2)

   # apply single qubit gates at indicated positions, otherwise apply identity
   sites_out = siteinds(2, N)
   for i in 1:qc.NumQubits
       if i in pos
           gate = ITensor(RZ, qc.IndexSet[i], sites_out[i])
       else
           gate = ITensor(E, qc.IndexSet[i], sites_out[i])
       end
       qc.StateVector = gate*qc.StateVector
       qc.IndexSet[i] = sites_out[i]
   end

   # create a new set of outgoing indices and define ITensor from array
   #sites_out = siteinds(2, N)
   #gate = ITensor(gate, qc.IndexSet, sites_out)

   # update state vector and outgoing indices
   #qc.StateVector = gate*qc.StateVector
   #qc.IndexSet = sites_out

   # update representing matrix of quantum circuit
   for i in 1:qc.NumQubits
       if i in pos
           push!(qc.Representation[i], 11)
           push!(qc.RepresentationFull[i], 11)
       else
           push!(qc.Representation[i], 0)
           push!(qc.RepresentationFull[i], 0)
       end
   end

   # update circuit depth
   qc.CircuitDepth += 1
end


""" Function to apply a phase shift by the angle theta to a quantum circuit
qc in every position specified in the array pos. """
function PhaseShift!(qc::QC_IT_ED, pos::Array{Int64, 1}, theta)

   # get matrices
   E, _, _, _ = get_Pauli_matrices()
   P = E * exp(1.0im * theta)

   # apply single qubit gates at indicated positions, otherwise apply identity
   sites_out = siteinds(2, N)
   for i in 1:qc.NumQubits
       if i in pos
           gate = ITensor(P, qc.IndexSet[i], sites_out[i])
       else
           gate = ITensor(E, qc.IndexSet[i], sites_out[i])
       end
       qc.StateVector = gate*qc.StateVector
       qc.IndexSet[i] = sites_out[i]
   end

   # create a new set of outgoing indices and define ITensor from array
   #sites_out = siteinds(2, N)
   #gate = ITensor(gate, qc.IndexSet, sites_out)

   # update state vector and outgoing indices
   #qc.StateVector = gate*qc.StateVector
   #qc.IndexSet = sites_out

   # update representing matrix of quantum circuit
   for i in 1:qc.NumQubits
       if i in pos
           push!(qc.Representation[i], 12)
           push!(qc.RepresentationFull[i], 12)
       else
           push!(qc.Representation[i], 0)
           push!(qc.RepresentationFull[i], 0)
       end
   end

   # update circuit depth
   qc.CircuitDepth += 1
end


###########################################
## Single qubit gates for Quantum Simulator
###########################################


""" Function to apple the (single-qubit) Hadamard gate to a state psi
(with N qubits) in every position specified in the array pos. """
function hadamard!(qc::QC, pos::Array{Int64, 1})

   # get matrices
   h = 1/sqrt(2) * Complex.([1. 1.; 1. -1.])
   E, _, _, _ = get_Pauli_matrices()

   # initialise first matrix
   if pos[1] == 1
       gate = h
   else
       gate = E
   end

   # construct whole operator
   for i in 2:qc.NumQubits
       if i in pos
           gate = kron(gate, h)
       else
           gate = kron(gate, E)
       end
   end

   # create a new set of outgoing indices and define ITensor from Hadamard array
   #sites_out = siteinds(d, N)
   #gate = ITensor(gate, sites, sites_out)

   #return gate*psi
   qc.StateVector .= gate*qc.StateVector
   #psi .= gate*psi

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

""" Function to apple the (single-qubit) Pauli-X gate to a state psi
(with N qubits) in every position specified in the array pos. """
function PauliX!(qc::QC, pos::Array{Int64, 1})

   # get matrices
   E, X, _, _ = get_Pauli_matrices()

   # initialise first matrix
   if pos[1] == 1
       gate = X
   else
       gate = E
   end

   # construct whole operator
   for i in 2:qc.NumQubits
       if i in pos
           gate = kron(gate, X)
       else
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

""" Function to apple the (single-qubit) Pauli-Y gate to a state psi
(with N qubits) in every position specified in the array pos. """
function PauliY!(qc::QC, pos::Array{Int64, 1})

   # get matrices
   E, _, Y, _ = get_Pauli_matrices()

   # initialise first matrix
   if pos[1] == 1
       gate = Y
   else
       gate = E
   end

   # construct whole operator
   for i in 2:qc.NumQubits
       if i in pos
           gate = kron(gate, Y)
       else
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

""" Function to apple the (single-qubit) Pauli-Z gate to a state psi
(with N qubits) in every position specified in the array pos. """
function PauliZ!(qc::QC, pos::Array{Int64, 1})

   # get matrices
   E, _, _, Z = get_Pauli_matrices()

   # initialise first matrix
   if pos[1] == 1
       gate = Z
   else
       gate = E
   end

   # construct whole operator
   for i in 2:qc.NumQubits
       if i in pos
           gate = kron(gate, Z)
       else
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

""" Function to apple the (single-qubit) S-gate to a state psi
(with N qubits) in every position specified in the array pos. """
function SGate!(qc::QC, pos::Array{Int64, 1})

   # get matrices
   E, _, _, _ = get_Pauli_matrices()
   S = Complex.([1. 0.; 0. 1.0im])

   # initialise first matrix
   if pos[1] == 1
       gate = S
   else
       gate = E
   end

   # construct whole operator
   for i in 2:qc.NumQubits
       if i in pos
           gate = kron(gate, S)
       else
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

""" Function to apple the (single-qubit) T-gate to a state psi
(with N qubits) in every position specified in the array pos. """
function TGate!(qc::QC, pos::Array{Int64, 1})

   # get matrices
   E, _, _, _ = get_Pauli_matrices()
   T = Complex.([1. 0.; 0. exp(1.0im*pi/4)])

   # initialise first matrix
   if pos[1] == 1
       gate = T
   else
       gate = E
   end

   # construct whole operator
   for i in 2:qc.NumQubits
       if i in pos
           gate = kron(gate, T)
       else
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

""" Function to apple the (single-qubit) Rx-gate to a state psi
(with N qubits) in every position specified in the array pos. Performs
a rotation around the x-axis specified by theta."""
function RXGate!(qc::QC, pos::Array{Int64, 1}, theta)

   # get matrices
   E, X, _, _ = get_Pauli_matrices()
   RX = exp(-1.0im * X * theta/2)

   # initialise first matrix
   if pos[1] == 1
       gate = RX
   else
       gate = E
   end

   # construct whole operator
   for i in 2:qc.NumQubits
       if i in pos
           gate = kron(gate, RX)
       else
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

""" Function to apple the (single-qubit) Ry-gate to a state psi
(with N qubits) in every position specified in the array pos. Performs
a rotation around the y-axis specified by theta."""
function RYGate!(qc::QC, pos::Array{Int64, 1}, theta)

   # get matrices
   E, _, Y, _ = get_Pauli_matrices()
   RY = exp(-1.0im * Y * theta/2)

   # initialise first matrix
   if pos[1] == 1
       gate = RY
   else
       gate = E
   end

   # construct whole operator
   for i in 2:qc.NumQubits
       if i in pos
           gate = kron(gate, RY)
       else
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

""" Function to apple the (single-qubit) Rz-gate to a state psi
(with N qubits) in every position specified in the array pos. Performs
a rotation around the z-axis specified by theta."""
function RZGate!(qc::QC, pos::Array{Int64, 1}, theta)

   # get matrices
   E, _, _, Z = get_Pauli_matrices()
   RZ = exp(-1.0im * Z * theta/2)

   # initialise first matrix
   if pos[1] == 1
       gate = RZ
   else
       gate = E
   end

   # construct whole operator
   for i in 2:qc.NumQubits
       if i in pos
           gate = kron(gate, RZ)
       else
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

""" Function to apply a phase shift by the angle theta to a state psi
(with N qubits) in every position specified in the array pos. """
function PhaseShift!(qc::QC, pos::Array{Int64, 1}, theta)

   # get matrices
   E, _, _, _ = get_Pauli_matrices()

   # initialise first matrix
   if pos[1] == 1
       gate = E * exp(1.0im * theta)
   else
       gate = E
   end

   # construct whole operator
   for i in 2:qc.NumQubits
       if i in pos
           gate = kron(gate, E * exp(1.0im * theta))
       else
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

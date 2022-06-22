
#########################################################
# File to test the different quantum simulator components
#########################################################


#include("quantum_simulator.jl")
include("QSim.jl")

# set properties
N = 5
N_meas = 100
register = [1, 3]
theta = π/6
U = Complex.([1. 0.; 0. 1.0im])
backend = "ED_Julia"
#backend = "ED_ITensor"

# initialise circuit (set true for linear topology, else master topology)
#qcircuit_master = initialise_qcircuit(N)#, true)
#qcircuit_linear = initialise_qcircuit(N, true)

qcircuit_master = initialise_qcircuit(N, false, backend)#, true)
qcircuit_linear = initialise_qcircuit(N, true, backend)

# check attributes
println(qcircuit_linear.StateVector)
println(qcircuit_linear.NumQubits)
println(qcircuit_linear.CircuitDepth)
println(qcircuit_linear.ClassicalBits)
println("lin top: ", qcircuit_linear.LinearTopology)
println("Rep ", qcircuit_linear.Representation)
println("Rep full ", qcircuit_linear.RepresentationFull)
println("Conversion table ", qcircuit_linear.ConversionTable)

# check single-qubit gates
hadamard!(qcircuit_linear, [1, 2, 3, 4, 5])
PauliX!(qcircuit_linear, [1, 2, 3, 5])
PauliY!(qcircuit_linear, [1, 4, 5])
PauliZ!(qcircuit_linear, [1, 2, 3, 4, 5])
SGate!(qcircuit_linear, [4, 5])
TGate!(qcircuit_linear, [2, 3, 4, 5])
RXGate!(qcircuit_linear, [1, 5], theta)
RYGate!(qcircuit_linear, [1, 2, 4,], theta)
RZGate!(qcircuit_linear, [1, 3, 5], theta)
PhaseShift!(qcircuit_linear, [4, 5], theta)

hadamard!(qcircuit_master , [1, 2, 3, 4, 5])
PauliX!(qcircuit_master , [1, 2, 3, 5])
PauliY!(qcircuit_master , [1, 4, 5])
PauliZ!(qcircuit_master , [1, 2, 3, 4, 5])
SGate!(qcircuit_master , [4, 5])
TGate!(qcircuit_master , [2, 3, 4, 5])
RXGate!(qcircuit_master , [1, 5], theta)
RYGate!(qcircuit_master , [1, 2, 4,], theta)
RZGate!(qcircuit_master , [1, 3, 5], theta)
PhaseShift!(qcircuit_master , [4, 5], theta)

# check measurements
#sample_measurement(qcircuit, [1, 2], N_meas)
#sample_measurement(qcircuit, [3, 4, 5], N_meas)
#sample_measurement(qcircuit_linear, [1, 2, 3, 4, 5], N_meas)


# check drawings of circuit
draw_circuit(qcircuit_linear)
draw_circuit(qcircuit_master)


# study 2-site gates
cnot2!(qcircuit_master, [1, 2])
swap2!(qcircuit_master, [3, 4])
CU2!(qcircuit_master, U, [2, 3])

#swap!(qcircuit_master, [1, 5])
#unswap!(qcircuit_master, [1, 5])
#draw_circuit(qcircuit_master, true)


cnot2!(qcircuit_linear, [1, 2])
swap2!(qcircuit_linear, [3, 4])
CU2!(qcircuit_linear, U, [2, 3])

#swap!(qcircuit_linear, [1, 5])
#unswap!(qcircuit_linear, [1, 5])
#draw_circuit(qcircuit_linear, true)


# check CU gates
CU!(qcircuit_linear, U, [2, 5])
CU!(qcircuit_linear, U, [4, 1])
draw_circuit(qcircuit_linear, true)


CU!(qcircuit_master, U, [2, 4])
CU!(qcircuit_master, U, [5, 2])
draw_circuit(qcircuit_master, true)

# check CNOT gates
cnot!(qcircuit_linear, [3, 5])
cnot!(qcircuit_linear, [4, 2])
draw_circuit(qcircuit_linear, true)


cnot!(qcircuit_master, [2, 4])
cnot!(qcircuit_master, [5, 2])
draw_circuit(qcircuit_master, true)

# test Toffoli
toffoli3!(qcircuit_linear, [3, 4, 5])
toffoli3!(qcircuit_master, [3, 4, 5])

toffoli!(qcircuit_linear, [3, 4, 5])
toffoli!(qcircuit_linear, [1, 3, 5])
toffoli!(qcircuit_linear, [2, 5, 3])
toffoli!(qcircuit_linear, [2, 4, 1])
draw_circuit(qcircuit_linear)

toffoli!(qcircuit_master, [1, 3, 5])
toffoli!(qcircuit_master, [2, 5, 3])
toffoli!(qcircuit_master, [2, 4, 1])
draw_circuit(qcircuit_master, true)

# test full swap
fullSwap!(qcircuit_linear, [1, 5])
fullSwap!(qcircuit_master, [1, 5])

draw_circuit(qcircuit_linear)
draw_circuit(qcircuit_master)


#println(qcircuit.StateVector)

#measure!(qcircuit, [1, 2])
#println(qcircuit.ClassicalBits)

#draw_circuit(qcircuit)

#name = "test_circuit.csv"
#save_circuit(qcircuit, name)
#load_circuit(name)


N = 4

backend1 = "ED_Julia"
backend2 = "ED_ITensor"

qc1 = initialise_qcircuit(N, true, backend1)
qc2 = initialise_qcircuit(N, true, backend2)

PauliX!(qc1, [1, 4])
PauliX!(qc2, [1, 4])
#hadamard!(qc2, [1, 2, 4, 5])
#cnot!(qc1, [1, 4])
#cnot!(qc2, [1, 4])
cnot2!(qc1, [1, 2])
cnot2!(qc2, [1, 2])
cnot2!(qc1, [2, 3])
cnot2!(qc2, [2, 3])
cnot2!(qc1, [3, 4])
cnot2!(qc2, [3, 4])
#fullSwap!(qc1, [1, 3])
println("now")
#fullSwap!(qc2, [1, 3])
#swap!(qc2, [1, 3])
#unswap!(qc2, [1, 2])
#swap2!(qc2, [1, 2])

"""
println(qc1.StateVector)
swap2!(qc1, [1, 2])
println(qc1.StateVector)
swap2!(qc1, [2, 3])
println(qc1.StateVector)
swap2!(qc1, [1, 2])
println(qc1.StateVector)

println("\n")

println(reshape(array(qc2.StateVector), (1, 2^4)))
swap2!(qc2, [1, 2])
println(reshape(array(qc2.StateVector), (1, 2^4)))
swap2!(qc2, [2, 3])
println(reshape(array(qc2.StateVector), (1, 2^4)))
#swap2!(qc2, [2, 3])
#println(reshape(array(qc2.StateVector), (1, 2^4)))
swap2!(qc2, [1, 2])
println(reshape(array(qc2.StateVector), (1, 2^4)))
"""



sample_measurement(qc1, [1, 2, 3, 4], N_meas)
sample_measurement(qc2, [1, 2, 3, 4], N_meas)

draw_circuit(qc1, true)
draw_circuit(qc2, true)

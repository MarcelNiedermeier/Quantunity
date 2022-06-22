


include("QSim.jl")

#backend = "MPS_ITensor"
backend = "ED_Julia"

N = 5
N_meas = 5
lintop = false
maxdim = 40
theta = π/8

qc = initialise_qcircuit(N, lintop, backend, maxdim)

println(qc.StateVector)

hadamard!(qc, [i for i in 1:qc.NumQubits])
#PauliX!(qc, [1, 2, 7])
#PauliY!(qc, [4, 7, 9])
#PauliZ!(qc, [1, 10, 16])
#SGate!(qc, [3, 12, 19])
#TGate!(qc, [6, 15])
#SqrtX!(qc, [7, 13, 20])
#PhaseShift!(qc, [6, 9, 11, 15], theta)
#RXGate!(qc, [4, 14, 18], theta)
#RYGate!(qc, [2, 6, 17], theta)
#RZGate!(qc, [1, 8, 13], theta)

PauliX!(qc, [1, 2, 3, 4, 5])
PauliY!(qc, [1, 2, 3, 4, 5])
PauliZ!(qc, [1, 2, 3, 4, 5])
SGate!(qc, [1, 2, 3, 4, 5])
TGate!(qc, [1, 2, 3, 4, 5])
SqrtX!(qc, [1, 2, 3, 4, 5])
PhaseShift!(qc, [1, 2, 3, 4, 5], theta)
RXGate!(qc, [1, 2, 3, 4, 5], theta)
RYGate!(qc, [1, 2, 3, 4, 5], theta)
RZGate!(qc, [1, 2, 3, 4, 5], theta)

draw_circuit(qc)

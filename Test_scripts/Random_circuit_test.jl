
######################
## Random Circuit Test
######################


using DelimitedFiles
using Statistics
include("../QSim.jl")


# set constants
N = 10
maxdim = 20
maxdims = [4, 8, 12, 16, 20, 24, 28, 32, 36, 40, 44, 48, 52, 56, 60, 64]# 64, 128]
N_sample = 4
#N_meas = 10
#backend = "ED_Julia"
backend = "MPS_ITensor"
contmethod = "naive"
random = false
lintop = false
randombond = 2
#randombonds = [1, 2, 4]#, 8, 16]


# intialise quantum circuit
qc = initialise_qcircuit(N, lintop, backend, maxdim, contmethod, random, randombond)
#initial_state = deepcopy(qc.StateVector)

# construct random circuit
pos = 1
num = 10
D = 20
randomCircuit!(qc, pos, num, D, false)

initial_state = deepcopy(qc.StateVector)

#sequential_cnot!(qc, pos, num, true)
#QFT!(qc, pos, num, true)
#invQFT!(qc, pos, num, true)
#
#fid = abs(ITensors.dot(initial_state, qc.StateVector))^2
#println("fidelity: $fid")

draw(qc, true)

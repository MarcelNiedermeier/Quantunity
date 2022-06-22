
############################################
## Tutorial for Bernstein-Vazirani algorithm
############################################

# the Bernstein-Vazirani algorithm can find a secret bitstring
# which is encoded in the oracle (which consists of CNOT gates).

# load functions for quantum circuit calculations
include("../QSim.jl")

# set constants
N = 4
N_meas = 1000
backend = "ED_Julia" # feel free to test it for exact state vectors as well!
#backend = "MPS_ITensor"
lintop = false

# set up circuit
qc = initialise_qcircuit(N, lintop, backend)

# prepare initial superposition
hadamard!(qc, [1, 2, 3, 4])
PauliZ!(qc, [4])

# construct/apply oracle
# this specific oracle encodes the bitstring [1, 0, 1]
cnot!(qc, [1, 4])
cnot!(qc, [3, 4])

# you could also construct the oracle as a custom gate
#U_BV = initialise_custom_gate(qc, N_reg, "BV")
#
#cnot!(U_BV, [1, 4])
#cnot!(U_BV, [3, 4])
#
#apply_custom_gate(qc, U_BV, 1)

# wrap up circuit
hadamard!(qc, [1, 2, 3])

# final measurement
println("Investigating result of Bernstein-Vazirani circuit: ")
register = [1, 2, 3]
sample_measurement(qc, register, N_meas)

# the measument yields the "secret" bitstring with 100% certainty,
# as desired :)

# let's have a look at the algorithm!
draw(qc)


############################
# Bernstein-Vazirani circuit
############################

include("../QSim.jl")


# set constants
N = 4
N_meas = 100
backend_ED = "ED_Julia"
backend = "MPS_ITensor"
lintop = false

# test custom gate
N_reg = 4
theta = 0.3843453
α = 1/sqrt(2)
β = 1/sqrt(2)
γ = 1/sqrt(2)
ϵ = -1/sqrt(2)



# set up circuit
qc = initialise_qcircuit(N, lintop, backend)

# prepare initial superposition
hadamard!(qc, [1, 2, 3, 4])
PauliZ!(qc, [4])

# Bernstein-Vazirani custom gates
U_BV = initialise_custom_gate(qc, N_reg, "BV")

cnot!(U_BV, [1, 4])
cnot!(U_BV, [3, 4])

apply_custom_gate(qc, U_BV, 1)


# construct/apply oracle
#cnot!(qc, [1, 4])
#cnot!(qc, [3, 4])

# wrap up circuit
hadamard!(qc, [1, 2, 3])
fullSwap!(qc, [1, 3])

# final measurement
println("Investigating result of Bernstein-Vazirani circuit: ")
register = [1, 2, 3]
sample_measurement(qc, register, N_meas)

draw(qc)

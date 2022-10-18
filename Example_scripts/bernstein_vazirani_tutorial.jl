
############################################
## Tutorial for Bernstein-Vazirani algorithm
############################################

using JLD2
using Plots

# the Bernstein-Vazirani algorithm can find a secret bitstring
# which is encoded in the oracle (which consists of CNOT gates).

# load functions for quantum circuit calculations
include("../QSim.jl")

# set constants
N = 4
N_meas = 1000
#backend = "ED_Julia" # feel free to test it for exact state vectors as well!
backend = "MPS_ITensor"
maxbond = 10

# set up circuit
qc = initialise_qcircuit(N, backend, maxdim=maxbond)

# prepare initial superposition
hadamard!(qc, [1, 2, 3, 4])
PauliZ!(qc, [4])

# construct/apply oracle
# this specific oracle encodes the bitstring [1, 0, 1]
#cnot!(qc, [1, 4])
#cnot!(qc, [3, 4])

# you could also construct the oracle as a custom gate
N_reg = 4
U_BV = initialise_custom_gate(qc, N_reg, "BV")

cnot!(U_BV, [1, 4])
cnot!(U_BV, [3, 4])

println(U_BV.Representation)

apply_custom_gate(qc, U_BV, 1)

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

#println(qc.ClassicalBitsProportion)
#states = collect(keys(qc.ClassicalBitsProportion))
#probs = collect(values(qc.ClassicalBitsProportion))

#dict = Dict(
#    "N" => N,
#    "N_meas" => N_meas,
#    "backend" => backend,
#    "states" => states,
#    "probs" => probs
#    )

#dict["new data"] = [1. , 2.0]

#println(d)

#save_data(qc, "Test_data.jld2", [])

#data = load_object("Test_data.jld2", d)
#data = load("Test_data.jld2")#["single_stored_object"]
#data = load_data("Test_data.jld2")

#println(data)

#println(qc.MeasurementResult)


# get measurement results and corresponding frequencies
#states, freqs = get_measurement_histogram(qc, t)
#states = collect(keys(qc.MeasurementResult))
#freqs = collect(values(qc.MeasurementResult))

#states_dec = bit_array_to_int.(states)

#println(states_dec)
#println(freqs)

# make (rudimentary) plot; save
#p = bar(states_dec, freqs, label="measurement")
#xlabel!("states (decimal expansion)")
#ylabel!("frequency")
#title!(title)

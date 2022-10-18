
##############################################
## Tutorial 2 for Bernstein-Vazirani algorithm
##############################################


# Here, let's go to somewhat bigger circuits and discover the
# MPS approach in the Bernstein-Vazirani algorithm!


# load functions for quantum circuit calculations
include("../QSim.jl")

# set constants
N = 20
N_meas = 1000
#backend = "ED_Julia" # feel free to test it for exact state vectors as well!
backend = "MPS_ITensor"
lintop = false
maxdim = 10

""" Function to generate a secret bitstring which will be
determined by the Bernstein-Vazirani algorithm. """
function generate_secret_bitstring(N)

    bitstring = []
    for i in 1:N
        push!(bitstring, rand([0, 1]))
    end

    return bitstring
end


""" Function to apply the oracle encoding a given secret
bitstring to a Bernstein-Vazirani circuit. """
function BV_oracle!(qc, bitstring)

    N = qc.NumQubits

    if length(bitstring) > N - 1
        error("Length of bitstring exceeds what can be encoded in
        quantum circuit!")
    end

    for i in 1:length(bitstring)
        if bitstring[i] == 1
            cnot!(qc, [i, N])
        end
    end
end



# set up circuit
qc = initialise_qcircuit(N, lintop, backend, maxdim)

# prepare initial superposition
hadamard!(qc, [i for i in 1:N])
PauliZ!(qc, [N])

# let's generate a secret bitstring
bitstring = generate_secret_bitstring(N-1)

# construct oracle according to this bitstring
BV_oracle!(qc, bitstring)

# wrap up circuit
hadamard!(qc, [i for i in 1:N])

# final measurement
register = [i for i in 1:N]
sample_measurement(qc, register, N_meas, save_measurement=false)

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
#
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

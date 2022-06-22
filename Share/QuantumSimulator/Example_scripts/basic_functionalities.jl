
##############################################
## Overview over basic functionalities of QSim
##############################################

# "mother file" for all functions
include("../QSim.jl")

# let's set some constants
N = 8 # total number of qubits
maxdim = 20 # maximum allowed bond dimension for MPS
contmethod = "naive" # contraction method for MPS (not so important right now)
random = false # no random intial state for circuit
lintop = false # possibility to use linear qubit topology, also don't touch for now :)
randombond = 10 # bond dimension of random initial MPS, if needed

# choose the method of calculating the quantum circuit
# every function described below exists for both backends and has the
# exact same structure! So you decide if you want to work with exact
# vectors or with MPS by simply setting the backend option, and the
# remaining code will look exactly the same!
backend_ED = "ED_Julia" # calculate with exact state vectors (up to 14ish qubits for now)
backend_MPS = "MPS_ITensor" # calculate with MPS (any higher number of qubits)

# now let's build some quantum circuits!
# need to start by initialising a quantum circuit object. Check both ED and MPS
qc_ED = initialise_qcircuit(N, lintop, backend_ED, maxdim,
contmethod, random, randombond)

qc_MPS = initialise_qcircuit(N, lintop, backend_MPS, maxdim,
contmethod, random, randombond)


#########################################################
# now we can apply gates. Let's do some single-site gates
#########################################################


# apply a Hadamard gate to every qubit in qc_ED
#"!" is important in the function call
hadamard!(qc_ED, [i for i in 1:N])

# then look at what has happened to the circuit (output on terminal)
draw(qc_ED)

# in the same fashion can apply every single-site gate. Basic gates:
PauliX!(qc_MPS, [1, 3, 5, 7])
PauliY!(qc_MPS, [1])
PauliZ!(qc_MPS, [4, 5, 6])
TGate!(qc_MPS, [4, 7, 8])
SGate!(qc_MPS, [1,2, 5])
SqrtX!(qc_MPS, [3, 4, 8])

# some gates depend on a parameter, e.g. rotations around x, y, z
θ = π/6 # random angle
RXGate!(qc_MPS, [2], θ)
RYGate!(qc_MPS, [5, 7], θ)
RZGate!(qc_MPS, [3, 8], θ)
PhaseShift!(qc_MPS, [4], θ)

# can also define a custom 2x2 unitary gate, U = [α β; γ δ]
# Note: the function will automatically check if the matrix is unitary!
α = exp(1.0im*2π*1/√3)
β = Complex(0.0)
γ = Complex(0.0)
δ = exp(1.0im*2π*1/√5)
UGate!(qc_MPS, [i for i in 1:N], α, β, γ, δ)

# look how the representation has changed!
draw(qc_MPS)

# Explanation of numbers: at the bottom there is the MAXIMUM bond
# dimension after each gate layer. As we apply only single-site gates
# here, the state remains a product state, i.e. the bond dimension
# is simply 1. On the right, there are the (possibly different) bond
# dimensions of the final state. The first number is the bond dim. between
# qubit 1 and qubit 2, and so on until the last number which is the bond
# dim. between qubit N-1 and qubit N.

# Note also: when the quantum circuit gets too big, the representation
# in the terminal is automatically cut off. So the draw feauture is
# most useful to check if small quantum circuits are constructed
# correctly.


#################################
# Next apply some two-site gates!
#################################

# Most importantly: CNOT gate. Always needs two positions, first
# is the "control qubit", second the "action qubit"
cnot!(qc_ED, [1, 4])

# look at circuit representation to see what has happened
draw(qc_ED)

# other two-site gates include a SWAP (sorry for the naming convetion below)
fullSwap!(qc_MPS, [3, 6])

# There are also controlled rotations and controlled arbitrary unitaries,
# but not relevant for us now


####################
# Custom gate blocks
####################

# if you need to re-apply a whole collection of gates multiple (many)
# times, you can group those gates into a custom gate block and
# then simply apply the whole block

# declare that you want a custom gate, here with 5 qubits
N_reg = 5
customgate = initialise_custom_gate(qc_ED, N_reg, "cust")

# let's put some gates into our custom gate block
# Note: exactly same syntax as when applying gates to quantum circuit!
cnot!(customgate, [1, 4])
cnot!(customgate, [3, 4])
PauliX!(customgate,[1, 2, 3])

# finally, apply custom gate to quantum circuit, starting at qubit 1
apply_custom_gate(qc_ED, customgate, 1)

# now look at what happenend to the quantum circuit!
draw(qc_ED)

# now we could apply the custom many times via a loop
# can change position where custom gate is appied!
for i in 1:4
    apply_custom_gate(qc_ED, customgate, i)
end
draw(qc_ED)


##############
# measurements
##############

# there are two types of measurements that we can perform in this framework:
# we can measure the circuit once, which leads to a collapse of the
# wave function. Afterwards, we could process the state further by applying
# more gates to it. The other possibility is to sample a large number
# of measurements to immediately obtain a statistical distribution of the
# outcomes. This is the kind of measurement that would only be performed
# at the end of the circuit. Those outcomes could then be further used
# to perform some classical post-processing on them.


# Let's try it. First we create a new circuit and construct the initial
# state |010101⟩

N = 6 # new number of qubits
N_meas = 1000 # number of measurements performed on quantum circuit

qc = initialise_qcircuit(N, lintop, backend_MPS, maxdim,
contmethod, random, randombond)

# create state by applying X-gates
PauliX!(qc, [2, 4, 6])

# now we sample a measurement of this state
register = [i for i in 1:N] # define subregister of qubits to be measured
meas = sample_measurement(qc, register, N_meas)

# meas is a dictionary containing the different outcomes with their
# corresponding frequency of occurrence
println("Measurement outcomes: ", meas)

# as expected, we find the state the system is in with 100% certainty,
# since it is just a computational basis state!

# let's draw it again. The measument is represented by M:
draw(qc)


# now try a different setup, with a quntum state in a perfect superposition:
qc2 = initialise_qcircuit(N, lintop, backend_MPS, maxdim,
contmethod, random, randombond)
hadamard!(qc2, [i for i in 1:N])

# measure only first four qubits, just to see what happens
register = [1, 2, 3, 4]
meas2 = sample_measurement(qc2, register, N_meas)

# get all 16 possible configurations, with roughly equal proportions!

# finally, let's see what a projective measurement does. We prepare
# the same state again, then do a projective measurement followed
# by a large sample of measurements. Those latter measurements should
# then all yield the same state, namely the one the circuit collapsed
# into

qc3 = initialise_qcircuit(N, lintop, backend_MPS, maxdim,
contmethod, random, randombond)
hadamard!(qc3, [i for i in 1:N])
register = [i for i in 1:N]
projective_measurement!(qc3, register)
meas3 = sample_measurement(qc3, register, N_meas)

# indeed, that is the outcome which is observed :)

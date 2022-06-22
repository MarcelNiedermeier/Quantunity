
#####################################
## Quantum Fourier Transform Tutorial
#####################################

# let's load again the functions
include("../QSim.jl")

# set constants
N = 5
maxdim = 30
N_meas = 100
#backend = "ED_Julia"
backend = "MPS_ITensor"
lintop = false
contmethod = "naive"
random = true # now, let's use a random initial state!
randombond = 3

# set up circuit
qc = initialise_qcircuit(N, lintop, backend, maxdim,
contmethod, random, randombond)

# save the initial state for later comparison
# the state is encoded as the "StateVector" subtype in the qc object
# use deepcopy to make sure to hold a true copy of the state,
# and not just a new pointer to the same memory address
initial_state = deepcopy(qc.StateVector)

# We also need the intial state as a vector, such that we can compute
# the Fourier coefficients directly. For this, we have to contract the MPS to
# recover a "normal" state vector:

initial_wf = deepcopy(get_wavefunction(qc))


#####################
# build QFT algorithm
#####################

# CRn implements a controlled rotation, where the rotation matrix is
# [1. 1.; 1. exp(sign(n)*2π*1.0im/2^abs(n))]. Specifying a negative
# n therefore leads to a controlled rotation in the other direction.
# The first number in the CRn! function is the value of n, afterwards
# the positions of the gate are specified as for CNOT, by giving the
# array ["control qubit", "action qubit"].

# first qubit
hadamard!(qc, [1])
CRn!(qc, 5, [5, 1])
CRn!(qc, 4, [4, 1])
CRn!(qc, 3, [3, 1])
CRn!(qc, 2, [2, 1])

# second qubit
hadamard!(qc, [2])
CRn!(qc, 4, [5, 2])
CRn!(qc, 3, [4, 2])
CRn!(qc, 2, [3, 2])

# third qubit
hadamard!(qc, [3])
CRn!(qc, 3, [5, 3])
CRn!(qc, 2, [4, 3])

# fourth qubit
hadamard!(qc, [4])
CRn!(qc, 2, [5, 4])

# fifth qubit
hadamard!(qc, [5])

# final swaps
fullSwap!(qc, [1, 5])
fullSwap!(qc, [2, 4])


# now let's have a look at the quantum circuit. Observe how the bond
# dimension changes throughout. In the representation, "Rn" is a rotation
# around the z-axis by a number "n" (unfortunately the n doesn't change,
# so it's always "Rn", even for different degrees of rotation)
draw(qc)


#####################
# Check result of QFT
#####################

# Now we should find out whether our circuit has implemented the quantum
# Fourier transform correctly. For this we compute the components exactly
# and compare!

""" Function to calculate the coefficients of a Fourier-transformed
state directly. """
function QFT_coeff(j, N, initial_state)
    ω = exp(1.0im*2*π/(2.0^N))
    b_j = Complex(0.0)
    for k in 0:(2^N-1)
        b_j += initial_state[k+1] * ω^(j*k)
    end
    return b_j/√(2^N)
end

# calculate Fourier coefficients exactly by hand
coeffs = []
for j in 0:(2^N-1)
    push!(coeffs, round(QFT_coeff(j, N, initial_wf), digits=16))
end

# next we compare the coefficients explicitly for the calculation "by hand"
# and the MPS
for i in 1:2^N
    println("QFT: $(get_wavefunction(qc)[i]),        direct calc: $(coeffs[i])")
end

# the terminal will be a bit full now, but you can easily see that they
# are both the same, up to numerical precision!

# to be completely sure, let's calculate also the norm of the difference
# QFT(|ψ_{exact}⟩) - QFT(|ψ_{MPS}⟩)
println("Check norm difference: ", norm(coeffs - get_wavefunction(qc)))

# it is indeed zero (again up to numerical precision) :)


#####################################
# Subroutines for QFT and inverse QFT
#####################################

# it would be annoying if we had to type in all the gates needed for
# a QFT every time. Therefore, there are subroutines to perform the
# QFT and the inverse QFT on an arbitrary subregister of qubits. All
# we have to do is specify the starting point of the QFT and the length
# of the subregister.

# To see it in action, let's apply the inverse QFT subroutine to the
# quantum circuit above: we state that we start at the first qubit
# and inverse Fourier transform a register of length N (i.e. the whole
# quantum circuit)
invQFT!(qc, 1, N)

# Let's see what we get!
draw(qc)

# the final state now should be (almost) identical to the initial state
# that we have saved in the very beginning. Almost, because we are restricting
# the bond dimension, i.e. we loose some information along the way, but
# hopefully not too much. Let's compare, by calculating
# || |ψ_{initial}⟩ - QFT^{-1}(QFT(|ψ_{initial}⟩)) ||:

# note that here we are taking the difference of two MPS, not two normal
# vectors. ITensor can handle this difference, and computes the norm
# in MPS-fasion
print("Check norm difference: ", norm(initial_state - qc.StateVector))


#########
# Summary
#########

# using those subroutines, Fourier transforming any subregister of a
# potentially large quantum circuit becomes a no-brainer! For concreteness,
# let's take a larger circuit and do a QFT and an inverse QFT

# set constants
N = 30
maxdim = 20
N_meas = 100
random = false

# set up circuit
qc2 = initialise_qcircuit(N, lintop, backend, maxdim,
contmethod, random, randombond)

# set quantum circuit to initial state |11....1⟩
PauliX!(qc2, [i for i in 1:N])

# do QFT and inverse QFT, e.g. on a 20-qubit subregister starting at qubit 5
# Note: this will take a some time to run, probably a few seconds up to a minute
QFT!(qc2, 5, 20)
invQFT!(qc2, 5, 20)

# perform a measurement to confirm that we indeed recover the initial
# state, as we should!
meas = sample_measurement(qc2, [i for i in 1:N], N_meas)

# check again what it looks like
draw(qc2)

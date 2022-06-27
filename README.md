# QuantumSimulator

This package is a general-purpose, high-level quantum circuit simulator. It combines both an exact state vector and a matrix product state (MPS)-based approach (add REF for MPS). The MPS functionalities are built on top of the existing Julia package `ITensors` (REF ITensors). The user can construct any quantum circuit consisting of the usual quantum gates, with the possibility to define custom gates acting on up to eight qubits. By adjusting the maximum bond dimension if MPS are chosen as a backend, the user can place an upper bound on the possible entanglement in the quantum circuit and study different quantum algorithms as a function of the size of the accessible Hilbert space. It has been shown that this approach is comparable to running quantum circuits on NISQ quantum computing devices (REF to PRX paper); thus MPS-based quantum computation gives a means to analyse to computational power of currently available quantum computers.

## Overview 

The main building blocks of the QuantumSimulator are contained in the file `QSim.jl`. Each quantum circuit is initiliased as a struct, which keeps track of the ciruit's attributes. The file `Hilbert_Space.jl` contains predefined operators for a (multi-) qubit system (see the documentation of ITensors for more details). Different measurement functions are contained in the `measurement.jl` file. In particular, one may perform either an "end-of-circuit" measurement with large statistics (many samples of the quantum state are evaluated), or a single "mid-circuit" measurement, leading to a subsequent collapse of the wave function. 
All the quantum gates are contained in the directory `functions`. Each gate function is named intuitively after the gate it is modeling, and acts on a quantum circuit object by updating its state vector. By the virtue of Julia's multiple dispatch funtionalities, it is automatically recognised whether the state vector is a "true" vector or an MPS; therefore all MPS-based quantum circuits can be written with the exact same high-level commands as the exact quantum circuits. In addition, each quantum circuit object holds information about which gates have already been applied to it. The `draw()` function is able to read out this information and print it on the terminal, such that the user can easily verify if a quantum circuit under construction is implemented correctly. In order to get started with using the package, we recommend to have a look at the scripts in the `Example_scripts` directory. They contain various basic quantum circuits with a many explanations, allowing the user to quickly learn how to use the main functions of this package.

## Installation

The QuantumSimulator is entirely built on top of the existing `ITensors` (add link) Julia package. Note that `ITensors` currently requires at least the Julia 1.6 distribution.

-- more details to come --

## Usage and Examples

Below we will give an overview of the most important functions and how to use them. More details can be found in the files in the `Example_scripts` directory. 

### A Minimal Working Example: The Bernstein-Vazirani Circuit

To give you a first flavour of how to use the QuantumSimulator, we will implement the Bernstein-Vazirani algorithm. The purpose of this section is not to make use of the underlying MPS-backend, but simply to get an idea of how the high-level functions and different building block of this package are implemented. 

We don't provide a full explanation of the Bernstein-Vazirani algorithm below, for more details see (REF). However, the basic idea is to determine a secret bitstring. This bitsring needs to be encoded into a quantum oracle, and the algorithm will be able to find the bitstring with only a single call to the oracle.

First we need to import the package:
```
include("YOUR_PATH_TO_PACKAGE/QSim.jl")
```

Then, we need to declare the main parameters defining the quantum circuit:
```
N = 4
N_meas = 1000
backend = "ED_Julia" # or "MPS_ITensor"
lintop = false # remnant from previous version, for now simply keep as "false"
```
The `backend` parameter determines whether the computation is done on exact state vectors or MPS. After that, the code is *exactly* the same for both cases. Neat!

Quantum circuit objects are declared as follows:
```
qc = initialise_qcircuit(N, lintop, backend)
```

Here, you can see already that the main parameters are the number of qubits and the backend (full state vectors *or* MPS) needed for the simulation. Now, let's build the quantum circuit:
```
# prepare initial superposition
hadamard!(qc, [1, 2, 3, 4])
PauliZ!(qc, [4])

# construct/apply oracle: this specific oracle encodes the bitstring [1, 0, 1]
cnot!(qc, [1, 4])
cnot!(qc, [3, 4])

# wrap up circuit
hadamard!(qc, [1, 2, 3])
```
The above code already tells you everything about how to apply (multiple) single-site gates and singly controlled gates: for any single site gate, you can specify an abitrary list of qubit positions where the gate(s) should be applied. As for the CNOT gate (and later more general controlled gates), the placement is always described by a list of two positions, the first one being the control qubit and the second one the "action qubit".

In order to recover the secret bitstring, we will of course need to measure the final state. A statistical measurement can be performed by calling:
```
register = [1, 2, 3]
sample_measurement(qc, register, N_meas)
```
Similarly to applying single-site gates, the measurement can be done on an arbitrary sub-register of the quantum circuit by choosing the desired qubits in the `register` argument ("arbitrary register" also means non-adjacent qubits).

Finally, we may draw the quantum circuit to the terminal by executing:
```
draw(qc)
```
This yields the representation
```
1   |0⟩ -H-------○-------H-----M
2   |0⟩ -H-------|-------H-----M
3   |0⟩ -H-------|---○---H-----M
4   |0⟩ -H---Z---+---+----------
```
The symbols used in the circuit above are hopefully self-explanatory. Due to the finite width of the terminal, the `draw()` function automatically truncated the graphical representation if the circuit becomes too deep to fit into a single line. Therefore, this feature should be seen rather as a cross-check for developing (smaller parts of) a quantum circuit than a fully-fledged graphics output. 

Another feature we would like to introduce at this stage is the possibility to implement custom sub-circuits. This will be especially useful in the context of quantum algorithms where a subroutine is applied/looped over many times throughout the circuit, e.g. for a trotterised time evolution. Here, one could declare the oracle as a stand-alone subroutine by replacing the CNOTs in the above circuit by the following commands:
```
#cnot!(qc, [1, 4])
#cnot!(qc, [3, 4])
N_reg = 4
U_BV = initialise_custom_gate(qc, N_reg, "BV")
cnot!(U_BV, [1, 4])
cnot!(U_BV, [3, 4])
apply_custom_gate(qc, U_BV, 1)
```
The gates are applied to the subcircuit in the exact same fashion as they would be applied to the "full" circuit. Note that upon applying the subciruit, one has to specify the position to which the first qubit in the subcircuit corresponds (as the subcircuit may have less qubits than the full circuit). Correspondingly, the graphical representation of the quantum circuit now contains the custom block:
```
1   |0⟩ -H-------+BV++---H-----M
2   |0⟩ -H-------|   |---H-----M
3   |0⟩ -H-------|   |---H-----M
4   |0⟩ -H---Z---+++++----------
```
In the top line of the custom gate block, the name string is included (please choose a maximum length of five characters).


### Quantum Fourier Transform

Here, we will show how one can implement a quantum Fourier transform, which an important building block of many more advanced quantum algorithms. To check the accuracy of our results, we will compare it to an exact evaluation of the Fourier coefficients. Let's start again by importing the simulator and setting the initial parameters:
```
include("YOUR_PATH_TO_PACKAGE/QSim.jl")
N = 5
maxdim = 30 # maximum allowed bond dimension
N_meas = 100
backend = "MPS_ITensor" 
lintop = false
contmethod = "naive" # also no need to worry about
random = true # now, let's use a random initial state!
randombond = 3 # bond dimension of random state
```

Then, we can initialise the quantum circuit with a random initial state and copy the wave functions for later comparisons:
```
qc = initialise_qcircuit(N, lintop, backend, maxdim,
contmethod, random, randombond)

# the state is encoded as the "StateVector" subtype in the qc object
# use deepcopy to make sure to hold a true copy of the state,
# and not just a new pointer to the same memory address
initial_state = deepcopy(qc.StateVector)

# we also need the intial state as a vector, such that we can compute
# the Fourier coefficients directly. For this, we have to contract the MPS to
# recover a "normal" state vector:
initial_wf = deepcopy(get_wavefunction(qc))
```

Now we will build the quantum circuit of the QFT for five qubits. For further details on the structure of the algorithm, feel free to consult the usual texts on quantum computation.
```
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
```
Here, the `CRn!`'s implement controlled rotations, where the rotation matrix is \[1. 1.; 1. exp(sign(n)*2π*1.0im/2^abs(n))]. Specifying a negative
n therefore leads to a controlled rotation in the other direction. The first number in the `CRn!` function is the value of n, afterwards the positions of the gate are specified as for `CNOT`.

Let's check again the graphical representation of our quantum circuit:
```
1   |0⟩ -H---Rn--Rn--Rn--Rn------------------------------------------⊗------  10
2   |0⟩ -----|---|---|---○---H---Rn--Rn--Rn--------------------------|---⊗--  7
3   |0⟩ -----|---|---○-----------|---|---○---H---Rn--Rn--------------|---|--  4
4   |0⟩ -----|---○---------------|---○-----------|---○---H---Rn------|---⊗--  2
5   |0⟩ -----○-------------------○---------------○-----------○---H---⊗------
Chi:   3   3   8   14  8   16  16  16  15  14  14  14  16  16  16  16  9   10 
```
Some new feautures have appeared here. Firstly, note that for the moment the initial state is still called "|0>", even though we have declared it as a random state. More interesting, however, are the numbers appearing below the circuit. They give you some information about the evolution of the bond dimension. The numbers below the circuit represent the maximum bond dimension after each gate layer (recall that the initial state was chosen with a bond dimension of 3, that's why we start at 3 here). The numbers at the end of the circuit describe the bond dimensions of the final state (i.e. 10 between qubit 1 and qubit 2, 7 between qubit 2 and qubit 3 and so on).


So far so good, now we hold the Fourier transform of our initial quantum. But crucially, we would like to know *how good* the transformation was, given that we have restricted the accuracy by putting an upper bound on the bond dimension (i.e. the accessible entanglement) in the algorithm. Well, let's calculate the QFT analytically and see what we get. 
```
""" Function to calculate the coefficients of a Fourier-transformed
state analytically. """
function QFT_coeff(j, N, initial_state)
    ω = exp(1.0im*2*π/(2.0^N))
    b_j = Complex(0.0)
    for k in 0:(2^N-1)
        b_j += initial_state[k+1] * ω^(j*k)
    end
    return b_j/√(2^N)
end

coeffs = []
for j in 0:(2^N-1)
    push!(coeffs, round(QFT_coeff(j, N, initial_wf), digits=16))
end
```

In order to quantify the difference between the analytical solution and the MPS-based quantum circuit, we can evaluate the norm difference:
```
# QFT(|ψ_{exact}⟩) - QFT(|ψ_{MPS}⟩)
println("Check norm difference: ", norm(coeffs - get_wavefunction(qc)))
```
Play around with some different numbers of qubits in the algorithms and different bond dimensions to see how the errors depend on both!

Since QFTs and inverse QFTs are essential building blocks of many quantum algorithms, it would be annoying if we had to type them by hand every time and for different qubit numbers. Therefore, we have implemented subroutines which create the corresponding quantum circuits with the following commands:
```
QFT!(qc, 3, 6)
invQFT!(qc, 3, 6)
```
The positional arguments should be read as follows: we create an (inverse) QFT ranging over a subregister of length 6, starting at qubit 3. This allows us to place the (inverse) QFT at arbitrary positions in the quantum circuit (of course you should make sure that your quantum circuit has been declared big enough to host a QFT with the given positional arguments).




### Quantum Phase Estimation

-- fill in --

### Overview: Available Quantum Gates

For completeness, below we list all gates which are currently implemented with dedicated function calls. The names and definitions follow the usual conventions for quantum gates.
```
# Basic single-site gates
hadamard!(qc, [...])
PauliX!(qc, [...])
PauliY!(qc, [...])
PauliZ!(qc, [...])
TGate!(qc, [...])
SGate!(qc, [...])
SqrtX!((qc, [...])

# some gates depend on a parameter, e.g. rotations around x, y, z
θ = π/6 # random angle
RXGate!(qc, [...], θ)
RYGate!(qc, [...], θ)
RZGate!(qc, [...], θ)
PhaseShift!(qc, [...], θ)

# can also define a custom 2x2 unitary gate, U = [α β; γ δ]
# Note: this function will automatically check if the matrix U is unitary!
α = exp(1.0im*2π*1/√3)
β = Complex(0.0)
γ = Complex(0.0)
δ = exp(1.0im*2π*1/√5)
UGate!(qc, [...], α, β, γ, δ)

# two-site gates
cnot!(qc, [...])
fullSwap!(qc, [...])
#
#

# three-site gates
#
#
```

## Disclaimer

This Julia package is still under development as a part of the author's doctoral research in the groups for *Correlated Quantum Materials* and *Quantum Transport* at Aalto University in Espoo, Finland. If you discover any bugs, inconsitencies or other issues while using it, or if you have suggestions for features which should be implemented but aren't, please contact me at `marcel.niedermeier[at]aalto.fi`. I would be very glad for any feedback!

## License, Citation

-- yet to come --

## References

-- yet to come --

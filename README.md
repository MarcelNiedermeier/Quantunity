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
#cnot!(qc, [3, 4
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

-- fill in --

### Quantum Phase Estimation

-- fill in --

### Overview: Available Quantum Gates

-- fill in --

## Disclaimer

This Julia package is still under development as a part of the author's doctoral research in the groups for *Correlated Quantum Materials* and *Quantum Transport* at Aalto University in Espoo, Finland. If you discover any bugs, inconsitencies or other issues while using it, or if you have suggestions for features which should be implemented but aren't, please contact me at `marcel.niedermeier[at]aalto.fi`. I would be very glad for any feedback!

## License, Citation

-- yet to come --

## References

-- yet to come --

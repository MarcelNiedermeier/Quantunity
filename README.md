# QuantumSimulator

This package is a general-purpose, high-level quantum circuit simulator. It combines both an exact state vector and a matrix product state (MPS)-based approach (add REF for MPS). The MPS functionalities are built on top of the existing Julia package `ITensors` (REF ITensors). The user can construct any quantum circuit consisting of the usual quantum gates, with the possibility to define custom gates acting on up to eight qubits. By adjusting the maximum bond dimension if MPS are chosen as a backend, the user can place an upper bound on the possible entanglement in the quantum circuit and study different quantum algorithms as a function of the size of the accessible Hilbert space. It has been shown that this approach is comparable to running quantum circuits on NISQ quantum computing devices (REF to PRX paper); thus MPS-based quantum computation gives a means to analyse to computational power of currently available quantum computers.

## Overview 

The main building blocks of the QuantumSimulator are contained in the file `QSim.jl`. Each quantum circuit is initiliased as a struct, which keeps track of the ciruit's attributes. The file `Hilbert_Space.jl` contains predefined operators for a (multi-) qubit system (see the documentation of ITensors for more details). Different measurement functions are contained in the `measurement.jl` file. In particular, one may perform either an "end-of-circuit" measurement with large statistics (many samples of the quantum state are evaluated), or a single "mid-circuit" measurement, leading to a subsequent collapse of the wave function. 
All the quantum gates are contained in the directory `functions`. Each gate function is named intuitively after the gate it is modeling, and acts on a quantum circuit object by updating its state vector. By the virtue of Julia's multiple dispatch funtionalities, it is automatically recognised whether the state vector is a "true" vector or an MPS; therefore all MPS-based quantum circuits can be written with the exact same high-level commands as the exact quantum circuits. In addition, each quantum circuit object holds information about which gates have already been applied to it. The `draw()` function is able to read out this information and print it on the terminal, such that the user can easily verify if a quantum circuit under construction is implemented correctly. In order to get started with using the package, we recommend to have a look at the scripts in the `Example_scripts` directory. They contain various basic quantum circuits with a many explanations, allowing the user to quickly learn how to use the main functions of this package.

## Installation

The QuantumSimulator is entirely built on top of the existing `ITensors` (add link) Julia package. Note that `ITensors` currently requires at least the Julia 1.6 distribution.

-- more details to come --

## Usage and Examples

Below we will give an overview of the most important functions and how to use them. More details can be found in the `Example_scripts` directory. 

### Example 1

-- fill in --

### Example 2

-- fill in --

### Example 3

-- fill in --

## Disclaimer

This Julia package is still under development as a part of the author's doctoral research in the groups for *Correlated Quantum Materials* and *Quantum Transport* at Aalto University in Espoo, Finland. If you discover any bugs, inconsitencies or other issues while using it, or if you have suggestions for features which should be implemented but aren't, please contact me at `marcel.niedermeier[at]aalto.fi`. I would be very glad for any feedback!

## License, Citation

-- yet to come --

## References

-- yet to come --

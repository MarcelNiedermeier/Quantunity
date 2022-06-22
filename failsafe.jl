
##########################
## Old and dummy functions
##########################

# test

##############
# custom gates
##############

""" Function to initialise a custom gate object.
Only gets meta-information about gate object, doesn't do construction!
MPS VERSION! """
function initialise_custom_gate(qc::QC_IT_MPS, N_reg::Int64, gate_name::String)#, backend="MPS_ITensor", maxdim=100, contmethod="naive")

    if length(gate_name) > 5
        error("Chosen gate name is too long (max. 5 characters)!")
    end

    # get information about quantum circuit
    N = qc.NumQubits
    #lintop = qc.LinearTopology
    #maxdim = qc.MaxBondDim
    #contmethod = qc.ContractionMethod

    # initialise custom gate object
    #CGate = custom_gate_MPS(N, N_reg, 0, maxdim, contmethod,
    #SortedDict(), SortedDict(), Dict(), lintop)
    CGate = custom_gate_MPS(N, N_reg, 0, gate_name, SortedDict())

    # set up representation of custom gate
    for i in 1:N_reg
        CGate.Representation[i] = []
    end
    return CGate
end


""" Function to apply a custom gate at the desired position. pos is
the index of the first qubit of the total qubit register that the
custom gate acts on.
MPS VERSION! """
function apply_custom_gate(qc::QC_IT_MPS, CGate::custom_gate_MPS, pos::Int64)

    for qubit in collect(keys(CGate.Representation))
        for gate in CGate.Representation[qubit]

            # gates with 4 parameters
            if length(gate) == 5
                UGate!(qc, [qubit+pos-1], gate[2], gate[3], gate[4], gate[5], false)

            # gates with 1 parameter
            elseif length(gate) == 2
                if gate[1] == 9.
                    RXGate!(qc, [qubit+pos-1], gate[2], false)
                elseif gate[1] == 10.
                    RYGate!(qc, [qubit+pos-1], gate[2], false)
                elseif gate[1] == 11.
                    RZGate!(qc, [qubit+pos-1], gate[2], false)
                elseif gate[1] == 12.
                    PhaseShift!(qc, [qubit+pos-1], gate[2], false)
                end

            # fixed gates
            elseif length(gate) == 1
                if gate[1] == 1.
                    hadamard!(qc, [qubit+pos-1], false)
                elseif gate[1] == 2.
                    PauliX!(qc, [qubit+pos-1], false)
                elseif gate[1] == 6.
                    PauliY!(qc, [qubit+pos-1], false)
                elseif gate[1] == 5.
                    PauliZ!(qc, [qubit+pos-1], false)
                elseif gate[1] == 18.
                    SqrtX!(qc, [qubit+pos-1], false)
                elseif gate[1] == 7.
                    SGate!(qc, [qubit+pos-1], false)
                elseif gate[1] == 8.
                    TGate!(qc, [qubit+pos-1], false)
                end

            # two-site gates
            elseif length(gate) == 3
                if gate[3] == 1.
                    cnot!(qc, [Int64(gate[1]), Int64(gate[2])], false)
                end
            end
        end
    end

    # update other circuit quantities

    # generate name string
    if length(CGate.GateName) == 1
        gatestring = "-++"*CGate.GateName[1]*"++--"
    elseif length(CGate.GateName) == 2
        gatestring = "-+"*CGate.GateName[1:2]*"++--"
    elseif length(CGate.GateName) == 3
        gatestring = "-+"*CGate.GateName[1:3]*"+--"
    elseif length(CGate.GateName) == 4
        gatestring = "-"*CGate.GateName[1:4]*"+--"
    else
        gatestring = "-"*CGate.GateName*"--"
    end

    qc.ConversionTable[100] = gatestring[1:4]
    qc.ConversionTable[101] = gatestring[5:8]

    # update representation of quantum circuit (general "block")
    for j in 1:2
        if j == 1
            for i in 1:qc.NumQubits
                if i == pos
                    push!(qc.Representation[i], 100)
                    push!(qc.RepresentationFull[i], 100)
                elseif i == pos+CGate.NumQubitsGate-1
                    push!(qc.Representation[i], 22)
                    push!(qc.RepresentationFull[i], 22)
                elseif pos < i && i < pos+CGate.NumQubitsGate-1
                    push!(qc.Representation[i], 20) # vertical line
                    push!(qc.RepresentationFull[i], 20) # vertical line
                else
                    push!(qc.Representation[i], 0)
                    push!(qc.RepresentationFull[i], 0)
                end
            end
        else
            for i in 1:qc.NumQubits
                if i == pos
                    push!(qc.Representation[i], 101)
                    push!(qc.RepresentationFull[i], 101)
                elseif i == pos+CGate.NumQubitsGate-1
                    push!(qc.Representation[i], 23)
                    push!(qc.RepresentationFull[i], 23)
                elseif pos < i && i < pos+CGate.NumQubitsGate-1
                    push!(qc.Representation[i], 21) # vertical line
                    push!(qc.RepresentationFull[i], 21) # vertical line
                else
                    push!(qc.Representation[i], 0)
                    push!(qc.RepresentationFull[i], 0)
                end
            end
        end
    end
end


#################
## in development
#################


""" Function to save a quantum circuit into a .csv file, such that a complex
circuit could be used within another file without the need for including it
explicitly in the new script. Need to be provided with the path/filename
to store the .csv file. """
function save_circuit(qc, name)

    # check if path is string
    if typeof(name) != String
        error("Wrong data type for file path (need string).")
    end

    # get data from circuit
    N = qc.NumQubits
    rep = qc.Representation
    maxdim = qc.MaxBondDim
    repFull = qc.RepresentationFull
    len = length(collect(values(rep))[1])
    lenFull = length(collect(values(repFull))[1])

    matQC = zeros(Int, N, len)
    matQCFull = zeros(Int, N, lenFull)

    # write to matrix
    for i in 1:N
        for j in 1:len
            matQC[i, j] = Int(rep[i][j])
        end
        for k in 1:lenFull
            matQCFull[i, k] = Int(repFull[i][k])
        end
    end

    # save in CSV file
    # header line with information
    # matrices with representations
    open("Data/Random_circuits/random_circuit_rand_MPS_fid_D_$(depth)_bond_$(maxdim)_samples_$(N_samples).csv", "w") do io
        writedlm(io, matQCFull, ", ")
    end

    #open(name, "w") do io
    #    writedlm(io, matQCFull)
    #end

end


function load_circuit(name)

    matQC = readdlm(name)
    println(matQC)

    N = size(matQC)[1]

end


####################
## - end development
####################


#############
# measurement
#############


""" Function to sample the measurements of the quantum circuit qc N_meas times.
Register specifies the sequence of qubits to be measured; can be an arbitrary
subregister of the full state. Prints out the different outcomes with their
corresponding frequency of occurrence ("empirical probability"). To be used at
the end of a quantum circuit to evaluate the result; doesn't perform collapse
of the wave function (as multiple samples are usually desired). Only shows the
measured states which have a frequency of occurrence bigger than eps (in order
to keep the output more readable and to suppress low-probability states).
ED VERSION FOR JULIA ARRAYS!"""
function sample_measurement(qc::QC, register::Array{Int64, 1}, N_meas=100,
    eps=0.005, verbose=true, save_measurement=true)

    # check correct ordering of values in register
    for i in 1:length(register)-1
        if register[i] >= register[i+1]
            error("Wrong ordering of qubit indices in register.")
        end
    end

    # obtain full probability distribution from wavefunction
    N_qubits = qc.NumQubits
    states = [reverse(digits(i, base=2, pad=N_qubits)) for i in 0:(2^N_qubits-1)]
    probabilities = abs.(qc.StateVector).^2

    # create dictionary with different states and corresponding probabilities
    probs = Dict()
    for i in 1:2^N_qubits
        probs[states[i]] = probabilities[i]
    end

    # if only measuring subregister: caculate marginalised probability distribution
    N_register = 2^length(register)
    states_marg = [reverse(digits(i, base=2, pad=length(register))) for i in 0:N_register-1]
    prob_list = []
    register_marg = Dict()

    # create a lists of equivalent states with same qubits values in selected registers
    for st in states_marg
        equiv_states = filter_binary_numbers(states, st, register)

        # calculate probability of corresponding subspaces
        prob_marg = 0
        for e_st in equiv_states
            prob_marg += probs[e_st]
        end
        push!(prob_list, prob_marg)
        register_marg[st] = (prob_marg, equiv_states)
    end

    # calculate weights corresponding to marginalised probs, sample number
    # of measurements
    weights = Weights(Vector{Float64}(collect(values(prob_list))))
    samp = StatsBase.sample(Random.GLOBAL_RNG, states_marg, weights, N_meas)
    if save_measurement
        for i in 1:N_meas
            qc.ClassicalBits[i] = samp[i]
        end
    end

    samp = proportionmap(samp)
    for key in sort!(collect(keys(samp)))
        #qc.ClassicalBitsProportion[key] = samp[key]
        qc.ClassicalBitsProportion[samp[key]] = key
    end

    # summarise result of measurements
    if verbose
        println("\n")
        println("Full qubit register: $([i for i in 1:N_qubits])")
        println("Doing $N_meas measurements of the register $register: ")
        println("Obtain the following states with their corresponding frequencies: ")
        println("Ignore states which have occured with a frequency of less than $eps: ")
        for key in sort!(collect(keys(samp)))
            item = samp[key]
            if item > eps
                println("State: ", key, ", p = ", item)
            end
        end
        println("\n")
    end

    # update representing matrix of quantum circuit
    for i in 1:qc.NumQubits
        if i in register
            push!(qc.Representation[i], 17)
            push!(qc.RepresentationFull[i], 17)
        else
            push!(qc.Representation[i], 0)
            push!(qc.RepresentationFull[i], 0)
        end
    end
end



""" Function to sample the measurements of the quantum circuit qc N_meas times.
Register specifies the sequence of qubits to be measured; can be an arbitrary
subregister of the full state. Prints out the different outcomes with their
corresponding frequency of occurrence ("empirical probability"). To be used at
the end of a quantum circuit to evaluate the result; doesn't perform collapse
of the wave function (as multiple samples are usually desired). Only shows the
measured states which have a frequency of occurrence bigger than eps (in order
to keep the output more readable and to suppress low-probability states). Can
choose between two algorithms to draw the sample of the MPS. Unless drawing
very few samples of a large register of a state with very high bond dimension,
the SVD-based algorithm should be preferred (especially if large statistics
are desired).
MPS VERSION FOR ITENSOR! """
function sample_measurement(qc::QC_IT_MPS, register::Array{Int64, 1}, N_meas=100,
    eps=0.005, verbose=true, algorithm="SVDbased", save_measurement=true)

    if algorithm ∉ ["SVDbased", "DirectSampling"]
        error("Invalid choice of algorithm (either SVDbased or DirectSampling).")
    end

    # check correct ordering of values in register
    for i in 1:length(register)-1
        if register[i] >= register[i+1]
            error("Wrong ordering of qubit indices in register.")
        end
    end

    # sample measurements of chosen register and get corresponding frequencies
    if algorithm == "SVDbased"
        measurements = sampleMPS!(qc.StateVector, register, N_meas)
    else # direct sampling
        measurements = sampleMPS2!(qc.StateVector, register, N_meas)
    end
    res = []
    for i in 1:N_meas
        push!(res, measurements[i, :])
        if save_measurement
            qc.ClassicalBits[i] = measurements[i, :]
        end
    end

    freq = proportionmap(res)
    for key in sort!(collect(keys(freq)))
        #qc.ClassicalBitsProportion[key] = freq[key]
        qc.ClassicalBitsProportion[freq[key]] = key
    end

    # ignore measurement results if probability smaller than eps
    for (key, value) in freq
        if value < eps
            delete!(freq, key)
        end
    end

    # summarise result of measurements
    if verbose
        println("\n")
        println("Full qubit register: $([i for i in 1:qc.NumQubits])")
        println("Doing $N_meas measurements of the register $register: ")
        println("Obtain the following states with their corresponding frequencies: ")
        println("Ignore states which have occured with a frequency of less than $eps: ")
        for key in sort!(collect(keys(freq)))
            println("State: ", key, ", p = ", freq[key])
        end
        println("\n")
    end

    # update representing matrix of quantum circuit
    for i in 1:qc.NumQubits
        if i in register
            push!(qc.Representation[i], 17)
            push!(qc.RepresentationFull[i], 17)
        else
            push!(qc.Representation[i], 0)
            push!(qc.RepresentationFull[i], 0)
        end
    end

    return freq
end


""" Function to perform a measurement of the quantum circuit qc.
Register specifies the sequence of qubits to be measured; can be an arbitrary
subregister of the full state. In contrast to the function sample_measurement(),
this function only evaluates a single sample and collapses the wavefunction
accordingly. To be used for a measurement within the circuit. The results of the
measurement are saved as classical bits in the ClassicalBits dictionary of the
quantum circuit object. Note that the explicit contruction of the projector
requires projecting the state onto a 2^(N-length(register))-dimensional subspace.
As the projected state is constructed sequentially (and compressed after each)
step, this doesn't lead to a breakdown of the code. However it is discouraged
to project a large qubit register onto a small subregister. Can choose between
two algorithms to draw the sample of the MPS. In case a large subregister of an
MPS with high bond dimension is sampled, the direct sampling should be preferred,
otherwise the SVD-based sampling.
MPS VERSION FOR ITENSOR! """
function projective_measurement!(qc::QC_IT_MPS, register::Array{Int64, 1},
    verbose=true, algorithm="SVDbased")

    if algorithm ∉ ["SVDbased", "DirectSampling"]
        error("Invalid choice of algorithm (either SVDbased or DirectSampling).")
    end

    # check correct ordering of values in register
    for i in 1:length(register)-1
        if register[i] >= register[i+1]
            error("Wrong ordering of qubit indices in register.")
        end
    end

    # get number of qubits and states
    N_qubits = qc.NumQubits
    states = [reverse(digits(i, base=2, pad=N_qubits)) for i in 0:(2^N_qubits-1)]

    # sample measurement of chosen register and get basis of subspace
    if algorithm == "SVDbased"
        measurement = [sampleMPS!(qc.StateVector, register, 1)[1, :]]
    else # direct sampling algorithm
        measurement = [sampleMPS2!(qc.StateVector, register, 1)[1, :]]
    end
    subspace_basis = filter_binary_numbers(states, measurement[1], register)

    # get coefficients for weighted sum
    coeff_list = zeros(Complex, length(subspace_basis))
    for i in 1:length(subspace_basis)
        coeff_list[i] = get_coefficient(qc, subspace_basis[i])
    end

    # construct state projected onto subspace
    out_state = coeff_list[1]*MPS_computationalBasis(qc.IndexSet, subspace_basis[1])
    for i in 2:length(subspace_basis)
        if verbose
            if i%100 == 0
                println("updating out state $i")
            end
        end
        out_state = out_state + coeff_list[i]*MPS_computationalBasis(qc.IndexSet,
                                                subspace_basis[i])
    end

    # save and renormalise state projected onto subspace
    qc.StateVector = 1/(norm(out_state)) * out_state

    # summarise result of measurements
    if verbose
        println("\n")
        println("Full qubit register: $([i for i in 1:N_qubits])")
        println("Measured register $register: ")
        println("The qubit register has collapsed to $(measurement[1]). ")
        println("\n")
    end

    # update representing matrix of quantum circuit
    for i in 1:qc.NumQubits
        if i in register
            push!(qc.Representation[i], 17)
            push!(qc.RepresentationFull[i], 17)
        else
            push!(qc.Representation[i], 0)
            push!(qc.RepresentationFull[i], 0)
        end
    end

    # save classical bits obtained from measurement
    for i in 1:length(register)
        qc.ClassicalBits[register[i]] = measurement[1][i]
    end

    # update cicruit depth and bond dimension
    qc.CircuitDepth += 1
    push!(qc.BondDim, maxlinkdim(qc.StateVector))
end


""" Function to perform a measurement of the quantum circuit qc.
Register specifies the sequence of qubits to be measured; can be an arbitrary
subregister of the full state. In contrast to the function sample_measurement(),
this function only evaluates a single sample and collapses the wavefunction
accordingly. To be used for a measurement within the circuit. The results of the
measurement are saved as classical bits in the ClassicalBits dictionary of the
quantum circuit object. Note that the explicit contruction of the projector
requires projecting the state onto a 2^(N-length(register))-dimensional subspace!
ED VERSION FOR JULIA ARRAYS! """
function projective_measurement!(qc::QC, register::Array{Int64, 1}, verbose=true)

    # check correct ordering of values in register
    for i in 1:length(register)-1
        if register[i] >= register[i+1]
            error("Wrong ordering of qubit indices in register.")
        end
    end

    # obtain full probability distribution from wavefunction
    N_qubits = qc.NumQubits
    states = [reverse(digits(i, base=2, pad=N_qubits)) for i in 0:(2^N_qubits-1)]
    probabilities = abs.(qc.StateVector).^2

    # create dictionary with different states and corresponding probabilities
    probs = Dict()
    for i in 1:2^N_qubits
        probs[states[i]] = probabilities[i]
    end

    # if only measuring subregister: caculate marginalised probability distribution
    N_register = 2^length(register)
    states_marg = [reverse(digits(i, base=2, pad=length(register))) for i in 0:N_register-1]
    prob_list = []
    register_marg = Dict()

    # create a lists of equivalent states with same qubits values in selected registers
    for st in states_marg
        equiv_states = filter_binary_numbers(states, st, register)

        # calculate probability of corresponding subspaces
        prob_marg = 0
        for e_st in equiv_states
            prob_marg += probs[e_st]
        end
        push!(prob_list, prob_marg)
        register_marg[st] = (prob_marg, equiv_states)
    end

    # calculate weights corresponding to marginalised probs, sample state
    weights = Weights(Vector{Float64}(prob_list))
    samp = StatsBase.sample(Random.GLOBAL_RNG, states_marg, weights, 1)

    # contruct projector onto measured subspace
    proj = Complex.(zeros(2^N_qubits, 2^N_qubits))
    for eq_st in register_marg[samp[1]][2]
        ind = bit_array_to_int(eq_st) + 1
        proj[ind, ind] = 1.0
    end

    # apply projector to state, renormalise
    qc.StateVector .= proj*qc.StateVector
    qc.StateVector .= qc.StateVector/sqrt(sum(abs.(qc.StateVector).^2))

    # save classical bits obtained from measurement
    for i in 1:length(register)
        qc.ClassicalBits[register[i]] = samp[1][i]
    end

    #sort!(collect(keys(qc.ClassicalBits)))

    # summarise result of measurements
    if verbose
        println("\n")
        println("Full qubit register: $([i for i in 1:N_qubits])")
        println("Measured register $register: ")
        println("The qubit register has collapsed to $(samp[1]). ")
        println("\n")
    end

    # update representing matrix of quantum circuit
    for i in 1:qc.NumQubits
        if i in register
            push!(qc.Representation[i], 16)
            push!(qc.RepresentationFull[i], 16)
        else
            push!(qc.Representation[i], 0)
            push!(qc.RepresentationFull[i], 0)
        end
    end
end


""" Function to compute the wavefunction coefficient of a given bitstring
in a quantum state represented by an MPS (by direct sequential contraction).
The bitstring is given as a list representing the configuration of a
computational basis vector, e.g. [0, 1, 0, 1]. """
function get_coefficient(qc::QC_IT_MPS, tensor_components::Array)

    if length(tensor_components) != qc.NumQubits
        error("Number of given configuration ($(length(tensor_components))
        elements) doesn't match number of qubits ($(qc.NumQubits))!")
    end

    # initialise coefficient, contract MPS
    V = ITensor(1.)
    for j in 1:qc.NumQubits

        #println(qc.StateVector[j])
        #println(state(qc.IndexSet[j], tensor_components[j]+1))

        # state(Index, 1) = [1, 0] and state(Index, 2) = [0, 1]
        V_tmp = (qc.StateVector[j]*state(qc.IndexSet[j], tensor_components[j]+1))
        #println(V_tmp)
        #V *= (qc.StateVector[j]*state(qc.IndexSet[j], tensor_components[j]+1))
        V *= V_tmp
    end

    # save result
    return scalar(V)

end


""" Function to compute the wavefunction coefficient of a given bitstring
in a quantum state represented by an MPS (by direct sequential contraction).
Note that the bitstring must be given in decimal form, i.e. for a quantum
circuit with N qubits the states of the computational basis are numbered
by |0>, ... |2^N-1>. """
function get_bitstring_coefficient(qc::QC_IT_MPS, bitstring::Int64)

    # get number of qubits, check validity of given bitstring
    N = qc.NumQubits
    if bitstring < 0 || bitstring > 2^N-1
        error("Given bitstring is not a possible measurement result for $N qubits.")
    end

    # convert to list specifying which basis states the MPS has to be cotracted with
    tensor_components = reverse(digits(bitstring, base=2, pad=N))

    # initialise coefficient, contract MPS
    V = ITensor(1.)
    for j in 1:N

        # state(Index, 1) = [1, 0] and state(Index, 2) = [0, 1]
        V *= (qc.StateVector[j]*state(qc.IndexSet[j], tensor_components[j]+1))
    end

    # save result
    return scalar(V)
end


""" Function to compute the wavefunction coefficient of a given bitstring
in a quantum state represented by an MPS (by direct sequential contraction).
Note that the bitstring must be given in decimal form, i.e. for a quantum
circuit with N qubits the states of the computational basis are numbered
by |0>, ... |2^N-1>. """
function get_bitstring_coefficient(qc::QC, bitstring::Int64)

    # get number of qubits, check validity of given bitstring
    N = qc.NumQubits
    if bitstring < 0 || bitstring > 2^N-1
        error("Given bitstring is not a possible measurement result for $N qubits.")
    end
    return qc.StateVector[bitstring]
end


""" Function to calculate the wavefunction corresponding to a quantum
circuit built from an MPS by explicit contraction of the MPS. Very
inefficient scaling of O(2^N)*O(N*χ²) (need to perform the contraction
exponentially many times), therefore only suitable for check with small
systems which can also be tackled with the exact state vector approach. """
function get_wavefunction(qc::QC_IT_MPS)

    # get number of qubits and initialise container
    N = qc.NumQubits
    wavefunction = Complex.(zeros(2^N))

    # loop through all possible components/bitstrings
    for i in 0:2^N-1
        wavefunction[i+1] = get_bitstring_coefficient(qc, i)
    end
    return wavefunction
end


""" Function to get the exact list of probabilities corresponding to
a wave function in the MPS representation of a quantum circuit. Can choose
between the direct or the cumulative probability distribution. Use only
for small quantum circuits. """
function get_probabilities(qc::QC_IT_MPS, cumulative=false)
    probs = abs.(get_wavefunction(qc)).^2
    if cumulative==false
        return probs
    else
        return cumsum(probs)
    end
end


""" Function to return the probalities corresponding to a wave function
in the exact representation of a quantum circuit. Can choose between the
direct or the cumulative probability distribution. The calculation is trivial,
this function merely serves as a structural counterpart for the more involved
MPS-based calculation. """
function get_probabilities(qc::QC, cumulative=false)
    probs = abs.(qc.StateVector).^2
    if cumulative==false
        return probs
    else
        return cumsum(probs)
    end
end


""" Sample given number N_meas of measurements from the probability
distribution defined by an MPS. Particularly effective for a large
number of measurements, as only a single sweep through the MPS is
required. Specify a register from the sample is evaluated directly. """
function sampleMPS!(M::MPS, register::Array{Int64, 1}, N_meas::Int64)

    N = length(M)
    N_reg = length(register)
    result = zeros(Int, N_meas, N_reg)

    for i in 1:N_reg

        # move orthocenter to qubit in question
        orthogonalize!(M, register[i])
        s = siteind(M, register[i])
        A = M[register[i]]

        # get list of probalities marginalised for qubit i
        probs = zeros(2)
        for n in 1:2
            An = ITensor()
            projn = ITensor(s)
            projn[s => n] = 1.0
            An = A * dag(projn)
            probs[n] = real(scalar(dag(An) * An))
        end

        # draw sample
        w = Weights((probs))
        result[:, i] = StatsBase.sample(Random.GLOBAL_RNG, [0, 1], w, N_meas)
    end
    return result
end


""" Sample given number N_meas of measurements from the probability
distribution defined by an MPS. Particularly effective for a large
number of measurements, as only a single sweep through the MPS is
required. """
function sampleMPS!(M::MPS, N_meas::Int64)

    N = length(M)
    result = zeros(Int, N_meas, N)

    for i in 1:N

        # move orthocenter to qubit in question
        orthogonalize!(M, i)
        s = siteind(M, i)
        A = M[i]

        # get list of probalities marginalised for qubit i
        probs = zeros(2)
        for n in 1:2
            An = ITensor()
            projn = ITensor(s)
            projn[s => n] = 1.0
            An = A * dag(projn)
            probs[n] = real(scalar(dag(An) * An))
        end

        # draw sample
        w = Weights((probs))
        result[:, i] = StatsBase.sample(Random.GLOBAL_RNG, [0, 1], w, N_meas)

    end
    return result
end


""" Adapted version of ITensor.sample(), which draws N_meas samples from
the probability distribution defined through an MPS representing a wave
function. Can specify a subregister from which to draw the samples. """
function sampleMPS2!(M::MPS, register::Array{Int64, 1}, N_meas::Int)

    #if abs(1.0 - norm(M)) > 1E-8
    #    error("sample: MPS is not normalized")
    #end

    orthogonalize!(M, 1)
    N = length(M)
    N_reg = length(register)
    result = zeros(Int, N_meas, N_reg)

    for l in 1:N_meas

        result_tmp = []
        A = M[1]

        for j in 1:N

            s = siteind(M, j)

            # Compute the probability of each state
            # one-by-one and stop when the random
            # number r is below the total prob so far

            # Will need n,An, and pn below
            k = 0
            proj1 = ITensor(s)
            proj1[s => 1] = 1.0
            An = A * dag(proj1)
            pn = real(scalar(dag(An) * An))

            if rand() > pn
                k = 1
                proj2 = ITensor(s)
                proj2[s => 2] = 1.0
                An = A * dag(proj2)
                pn += real(scalar(dag(An) * An))
            end

            # record measurement result only if in subregister
            if j in register
                push!(result_tmp, k)
            end

            if j < N
              A = M[j + 1] * An
              A *= (1.0 / sqrt(pn))
            end

        end
        result[l, :] = result_tmp[:]
    end

    return result
end


###############
# Old functions
###############


""" Function to extract the 2-permutations from a general array of poisitions,
pos = [pos[1], pos[2]], such that the permutation represented by pos is
equivalent to the sequence of 2-permutations. """
function get_2_perms(pos)
    perms_forward = []
    perms_backward = []
    for i in (pos[1]+1):pos[2]
        push!(perms, [i-1, i])
    end
    return perms
end


""" Function to sample the measurements of the quantum circuit qc N_meas times.
Register specifies the sequence of qubits to be measured; can be an arbitrary
subregister of the full state. Prints out the different outcomes with their
corresponding frequency of occurrence ("empirical probability"). To be used at
the end of a quantum circuit to evaluate the result; doesn't perform collapse
of the wave function (as multiple samples are usually desired). Only shows the
measured states which have a frequency of occurrence bigger than eps (in order
to keep the output more readable and to suppress low-probability states).
MPS VERSION FOR ITENSOR! """
function sample_measurement_old(qc::QC_IT_MPS, register::Array{Int64, 1}, N_meas=100, eps=0.001, verbose=true)

    # check correct ordering of values in register
    for i in 1:length(register)-1
        if register[i] >= register[i+1]
            error("Wrong ordering of qubit indices in register.")
        end
    end

    # get number of qubits and states
    N_qubits = qc.NumQubits
    states = [reverse(digits(i, base=2, pad=N_qubits)) for i in 0:(2^N_qubits-1)]

    # start by sampling from full wave function
    measurements = []
    for i in 1:N_meas
        s = ITensors.sample(orthogonalize!(qc.StateVector, 1))
        s_new = s - ones(Int32, length(s))
        push!(measurements, s_new)
    end

    # dictionary of states with their relative frequency
    freq = proportionmap(measurements)

    # if only measuring subregister: caculate marginalised probability distribution
    N_register = 2^length(register)
    states_marg = [reverse(digits(i, base=2, pad=length(register))) for i in 0:N_register-1]
    prob_list = []
    register_marg = Dict()

    # create a lists of equivalent states with same qubits values in selected registers
    for st in states_marg
        equiv_states = filter_binary_numbers(states, st, register)

        # calculate probability of corresponding subspaces
        prob_marg = 0
        for e_st in equiv_states
            if e_st in keys(freq)
                prob_marg += freq[e_st]
            end
        end
        push!(prob_list, prob_marg)
        register_marg[st] = (prob_marg, equiv_states)
    end

    # summarise result of measurements
    if verbose
        println("\n")
        println("Full qubit register: $([i for i in 1:N_qubits])")
        println("Doing $N_meas measurements of the register $register: ")
        println("Obtain the following states with their corresponding frequencies: ")
        println("Ignore states which have occured with a frequency of less than $eps: ")
        for key in sort!(collect(keys(register_marg)))
            item = register_marg[key]
            if item[1] > eps
                println("State: ", key, ", p = ", item[1])
            end
        end
        println("\n")
    end

    # update representing matrix of quantum circuit
    for i in 1:qc.NumQubits
        if i in register
            push!(qc.Representation[i], 17)
            push!(qc.RepresentationFull[i], 17)
        else
            push!(qc.Representation[i], 0)
            push!(qc.RepresentationFull[i], 0)
        end
    end

    return register_marg
end


""" test to benchmark performance (probably worse than other function)"""
function projective_measurement2!(qc::QC_IT_MPS, register::Array{Int64, 1}, verbose=true)

    # check correct ordering of values in register
    for i in 1:length(register)-1
        if register[i] >= register[i+1]
            error("Wrong ordering of qubit indices in register.")
        end
    end

    # get number of qubits and states
    N_qubits = qc.NumQubits
    states = [reverse(digits(i, base=2, pad=N_qubits)) for i in 0:(2^N_qubits-1)]

    # sample measurement of chosen register and get basis of subspace
    measurement = [sampleMPS!(qc.StateVector, register, 1)[1, :]]
    subspace_basis = filter_binary_numbers(states, measurement[1], register)

    # try to sum weighted states directly
    coeff_list = zeros(Complex, length(subspace_basis))
    #state_list = []
    for i in 1:length(subspace_basis)
        if verbose
            if i%100 == 0
                println("subspace basis $i")
            end
        end
        coeff_list[i] = get_coefficient(qc, subspace_basis[i])
        #push!(state_list, get_coefficient(qc, subspace_basis[i])*MPS_computationalBasis(qc.IndexSet, subspace_basis[i]))
    end
    #out_state = add(state_list...; cutoff=1e-4)
    out_state = coeff_list[1]*MPS_computationalBasis(qc.IndexSet, subspace_basis[1])
    for i in 2:length(subspace_basis)
        if verbose
            if i%100 == 0
                println("updating out state $i")
            end
        end

        # try compressing only after every fourth step
        # gives no effective difference in performance
        if i%4 == 0
            out_state = add(out_state, coeff_list[i]*MPS_computationalBasis(qc.IndexSet,
                            subspace_basis[i]))
        else
            out_state = add(out_state, coeff_list[i]*MPS_computationalBasis(qc.IndexSet,
                            subspace_basis[i]), method="naive")
        end
    end
    qc.StateVector = 1/(norm(out_state)) * out_state

    # construct projector
    #projector_list = []
    #for i in 1:length(subspace_basis)
    #    #println("subspace basis $i")
    #    basis_state = MPS_computationalBasis(qc.IndexSet, subspace_basis[i])
    #    proj = projector(basis_state, normalize=true)
    #    push!(projector_list, proj)
    #end
    ##proj_tot = add(projector_list...; cutoff=1e-4)
    #proj_tot = add(projector_list[1], projector_list[2]; cutoff=1e-8)
    #for i in 3:length(projector_list)
    #    #println("constructing projector $i")
    #    proj_tot = add(proj_tot, projector_list[i]; cutoff=1e-8)
    #end

    ## apply projector to state
    #qc.StateVector = contract(proj_tot, qc.StateVector, maxdim=qc.MaxBondDim)#, method="naive")
    #noprime!(qc.StateVector)

    # summarise result of measurements
    if verbose
        println("\n")
        println("Full qubit register: $([i for i in 1:N_qubits])")
        println("Measured register $register: ")
        println("The qubit register has collapsed to $(measurement[1]). ")
        println("\n")
    end

    # update representing matrix of quantum circuit
    for i in 1:qc.NumQubits
        if i in register
            push!(qc.Representation[i], 17)
            push!(qc.RepresentationFull[i], 17)
        else
            push!(qc.Representation[i], 0)
            push!(qc.RepresentationFull[i], 0)
        end
    end

    # save classical bits obtained from measurement
    for i in 1:length(register)
        qc.ClassicalBits[register[i]] = measurement[1][i]
    end

    # update cicruit depth and bond dimension
    qc.CircuitDepth += 1
    push!(qc.BondDim, maxlinkdim(qc.StateVector))
end


""" Adapted version of ITensor.sample(), which draws N_meas samples from
the probability distribution defined through an MPS representing a wave
function. FUNCTION TO PLAY AROUND WITH, DON'T USE IN CODE!"""
function sampleMPS_old(m::MPS, N_meas::Int)
  N = length(m)

  #if orthocenter(m) != 1
    #error("sample: MPS m must have orthocenter(m)==1")
  #end
  if abs(1.0 - norm(m[1])) > 1E-8
    error("sample: MPS is not normalized")
  end

  result_matrix = zeros(Int, N_meas, N)

  for l in 1:N_meas

      result = zeros(Int, N)
      #result_matrix = zeros(Int, N_meas, N)
      A = m[1]
      d = 2

      for j in 1:N


        s = siteind(m, j)
        #d = ITensors.dim(s)

        # Compute the probability of each state
        # one-by-one and stop when the random
        # number r is below the total prob so far
        pdisc = 0.0
        r = rand()
        #r_vec = rand(N_meas)

        # Will need n,An, and pn below
        #k = 1
        An = ITensor()
        pn = 0.0
        #pn_vec = zeros(Float, N_meas)



        #k = n-1 # correct automatically for state 0, 1 results
        k = 0
        #println("State ", 1)
        proj1 = ITensor(s)
        proj1[s => 1] = 1.0
        An = A * dag(proj1)
        pn = real(scalar(dag(An) * An))
        pdisc += pn
        #println("p1: ", pn)
        #println("pdisc: ", pdisc)
        #(r < pdisc) && break
        #n += 1

        if r > pdisc
            k = 1
            #println("State ", 2)
            proj2 = ITensor(s)
            proj2[s => 2] = 1.0
            An = A * dag(proj2)
            pn = real(scalar(dag(An) * An))
            pdisc += pn
            #println("p2: ", pn)
            #println("pdisc: ", pdisc)

        end


        #while n <= d
        #for n = 1:d

        #  k = n-1 # correct automatically for state 0, 1 results
        #  println("State ", n)
        #  projn = ITensor(s)
        #  projn[s => n] = 1.0
        #  An = A * dag(projn)
        #  pn = real(scalar(dag(An) * An))
        #  pdisc += pn
        #  println("pn: ", pn)
        #  println("pdisc: ", pdisc)
        #  (r < pdisc) && break
        #  #n += 1
        #
        #end

        #result[j] = n
        result[j] = k

        if j < N
          A = m[j + 1] * An
          A *= (1.0 / sqrt(pn))
        end

      end
    result_matrix[l, :] = result[:]
  end

  return result_matrix
end

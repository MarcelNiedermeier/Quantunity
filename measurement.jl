
#################################################
## Functions related to sampling and measurements
#################################################


##########################
# Statistical Measurements
##########################

""" Function to sample the measurements of the quantum circuit qc N_meas times.
Register specifies the sequence of qubits to be measured; can be an arbitrary
subregister of the full state. Prints out the different outcomes with their
corresponding frequency of occurrence ("empirical probability"). To be used at
the end of a quantum circuit to evaluate the result; doesn't perform collapse
of the wave function (as multiple samples are usually desired). Only shows the
measured states which have a frequency of occurrence bigger than eps (in order
to keep the output more readable and to suppress low-probability states).
ED VERSION FOR JULIA ARRAYS!"""
function sample_measurement(qc::QC, register::Array{Int64, 1}, N_meas=100;
    eps=0.005, verbose=true, save_measurement=true, plot=true)

    # check correct ordering of values in register
    for i in 1:length(register)-1
        if register[i] >= register[i+1]
            error("Wrong ordering of qubit indices in register.")
        end
    end

    # save number of measurements in quantum circuit
    qc.NumMeasurements = N_meas

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
    #println("samp ", samp)
    for key in sort!(collect(keys(samp)))
        qc.ClassicalBitsProportion[key] = samp[key]
        #qc.ClassicalBitsProportion[samp[key]] = key
    end
    #println("qc.ClassicalBitsProportion ", qc.ClassicalBitsProportion)

    # save full measurement results
    found_states = collect(keys(samp))
    for st in states_marg
        if st in found_states
            qc.MeasurementResult[st] = samp[st]
        else
            qc.MeasurementResult[st] = 0.
        end
    end
    #println("full measurement results ", qc.MeasurementResult)

    # make (rudimentary) plot
    if plot
        freqs = collect(values(qc.MeasurementResult))
        states_dec = bit_array_to_int.(collect(keys(qc.MeasurementResult)))
        p = bar(states_dec, freqs, label="measurement")
        xlabel!("states (decimal expansion)")
        ylabel!("frequency")
        display(p)
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
function sample_measurement(qc::QC_IT_MPS, register::Array{Int64, 1}, N_meas=100;
    eps=0.005, verbose=true, algorithm="ITensor", save_measurement=true, plot=true)

    if algorithm ∉ ["ITensor", "SVDbased", "DirectSampling"]
        error("Invalid choice of algorithm (either SVDbased or DirectSampling).")
    end

    # check correct ordering of values in register
    for i in 1:length(register)-1
        if register[i] >= register[i+1]
            error("Wrong ordering of qubit indices in register.")
        end
    end

    # save number of measurements in quantum circuit
    qc.NumMeasurements = N_meas

    # sample array of measurements
    if algorithm == "SVDbased"
        measurements = sampleMPS!(qc.StateVector, register, N_meas)
    elseif algorithm == "DirectSampling" # direct sampling
        measurements = sampleMPS2!(qc.StateVector, register, N_meas)
    else # ITensor
        measurements = sampleMPS3!(qc.StateVector, register, N_meas)
    end

    # save measurements in qc if desired
    res = []
    for i in 1:N_meas
        push!(res, measurements[i, :])
        if save_measurement
            qc.ClassicalBits[i] = measurements[i, :]
        end
    end

    # obtain frequency of occurence of each distinct measurement result
    freq = proportionmap(res)
    #println("freq ", freq)
    for key in sort!(collect(keys(freq)))
        #println("key (state): ", key)
        #println("frequency: ", freq[key])
        qc.ClassicalBitsProportion[key] = freq[key]
        #qc.ClassicalBitsProportion[freq[key]] = key
    end

    if save_measurement
        N_register = 2^length(register)
        states_marg = [reverse(digits(i, base=2, pad=length(register))) for i in 0:N_register-1]

        # save full measurement results
        found_states = collect(keys(freq))
        for st in states_marg
            if st in found_states
                qc.MeasurementResult[st] = freq[st]
            else
                qc.MeasurementResult[st] = 0.
            end
        end

        # make (rudimentary) plot
        if plot
            freqs = collect(values(qc.MeasurementResult))
            states_dec = bit_array_to_int.(collect(keys(qc.MeasurementResult)))
            p = bar(states_dec, freqs, label="measurement")
            xlabel!("states (decimal expansion)")
            ylabel!("frequency")
            annotate!((100, 0.2, "χ = $(maxdim)"))
            display(p)
        end
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


#########################
# Projective Measurements
#########################

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
function projective_measurement!(qc::QC_IT_MPS, register::Array{Int64, 1};
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
function projective_measurement!(qc::QC, register::Array{Int64, 1}; verbose=true)

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


############################
# Samling, MPS Manipulations
############################

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
    println("observer norm before: ", norm(M))

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
    println("observer norm after: ", norm(M))
    return result
end


""" Sample given number N_meas of measurements from the probability
distribution defined by an MPS. Particularly effective for a large
number of measurements, as only a single sweep through the MPS is
required. """
function sampleMPS!(M::MPS, N_meas::Int64)

    N = length(M)
    result = zeros(Int, N_meas, N)

    println("observer norm before: ", norm(M))

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
    println("observer norm after: ", norm(M))
    return result
end


""" Adapted version of ITensor.sample(), which draws N_meas samples from
the probability distribution defined through an MPS representing a wave
function. Can specify a subregister from which to draw the samples. """
function sampleMPS2!(M::MPS, register::Array{Int64, 1}, N_meas::Int)

    #if abs(1.0 - norm(M)) > 1E-8
    #    error("sample: MPS is not normalized")
    #end

    println("observer norm before: ", norm(M))

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
    println("observer norm after: ", norm(M))

    return result
end


""" Adapted version of ITensor.sample(), which draws N_meas samples from
the probability distribution defined through an MPS representing a wave
function. Can specify a subregister from which to draw the samples. """
function sampleMPS3!(M::MPS, register::Array{Int64, 1}, N_meas::Int)

    # renormalise is needed
    if abs(1.0 - norm(M)) > 1E-8
        M = 1/norm(M) * M
    end

    # container for results
    N_reg = length(register)
    result = zeros(Int, N_meas, N_reg)

    # draw desired number of samples from full MPS
    samples = []
    for i in 1:N_meas
        sample_full = ITensors.sample!(M)
        result[i, :] = getindex(sample_full, register)[:]
    end

    return result - ones(Int, N_meas, N_reg)

end

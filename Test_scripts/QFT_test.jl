
###########
## QFT test
###########

using LinearAlgebra
include("../QSim.jl")


# set constants
N = 13
maxdim = 20
#maxdims = [4, 8, 12, 16, 20, 24, 28, 32, 36, 40, 44, 48, 52, 56, 60, 64]# 64, 128]
N_sample = 4
N_meas = 1000
#backend = "ED_Julia"
backend = "MPS_ITensor"
contmethod = "naive"
random = false
lintop = false
randombond = 10
#randombonds = [1, 2, 4]#, 8, 16]


"""
function QFT_coeff_had(j, N)

    ω = exp(1.0im*2*π/(2.0^N))
    b_j = 0.0

    for k in 0:(2^N-1)
        b_j += ω^(j*k)
    end

    return b_j/2^N
end
"""


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


""" Function to calculate the coefficients of a inverse Fourier-transformed
state directly. """
function invQFT_coeff(j, N, initial_state)
    ω = exp(-1.0im*2*π/(2.0^N))
    b_j = Complex(0.0)
    for k in 0:(2^N-1)
        b_j += initial_state[k+1] * ω^(j*k)
    end
    return b_j/√(2^N)
end


# intialise quantum circuit
qc = initialise_qcircuit(N, lintop, backend, maxdim, contmethod,
random, randombond)



# set up initial random state
D = 20
randomCircuit!(qc, 1, N, D, true)
println(qc.StateVector)

# save intial state
#initial_state = deepcopy(qc.StateVector)
initial_state = deepcopy(get_wavefunction(qc))

# do QFT
QFT!(qc, 1, N)

# measure
#sample_measurement(qc, [i for i in 1:N])

# draw
draw(qc)

#println(qc.StateVector)

# calculate Fourier coefficients by hand
coeffs = []
j = 0
for j in 0:(2^N-1)
    #println("coeff for j=$j: $(round(QFT_coeff_had(j, N), digits=5))")
    #println("coeff for j=$j: $(round(QFT_coeff(j, N, initial_state), digits=5))")
    push!(coeffs, round(QFT_coeff(j, N, initial_state), digits=16))
end


#for i in 1:2^N
    #println("QFT: $(qc.StateVector[i]),        direct calc: $(coeffs[i])")
    #println("QFT: $(get_wavefunction(qc)[i]),        direct calc: $(coeffs[i])")
#end
#println(coeffs)


println("Check norm: ", sum(coeffs .* conj.(coeffs)))
#println("Check difference: ", norm(coeffs - qc.StateVector))
println("Check difference: ", norm(coeffs - get_wavefunction(qc)))

# check fidelity between quantum circuit QFT and mathematical QFT
fid = abs(dot(get_wavefunction(qc), coeffs))^2
println("fidelity: $fid")

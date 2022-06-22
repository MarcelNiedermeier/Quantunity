
#########################################################
# File to test the different quantum simulator components
#########################################################


include("quantum_simulator.jl")

"""
N = 4
N_meas = 100
theta = pi/6

psi = initialise_state(N)
println("initial ", psi)
do_measurement(psi, N_meas)

pos = [4]
pos2 = [4, 2]
#pos3 = [1, 2]
#pos4 = [1,3]
PauliX!(psi, pos, N)
#SGate!(psi, pos, N)
#TGate!(psi, pos, N)
#RXGate!(psi, [2], N, theta)
#RYGate!(psi, [3], N, theta)
#RZGate!(psi, [4], N, theta)
println("here")
println("after X ", psi)
do_measurement(psi, N_meas)
#psi = cnot(psi, pos2, N)
cnot!(psi, pos2, N)
println("after cnot ", psi)

#cnot!(psi, pos3, N)
#cnot!(psi, pos4, N)

#pos = [1, 2, 3]
#pos = [1]
#psi = hadamard(psi, pos, N)
#println(psi)

#psi = PauliY(psi, pos, N)
#println(psi)

do_measurement(psi, N_meas)

#println(reverse(get_2_perms(reverse([3,1]))))


# check if sampling works
#println(samplepair(Random.GLOBAL_RNG, 4))
#println(StatsBase.sample(Random.GLOBAL_RNG, [1, 2, 3], Weights([1/3, 1/3, 1/3])))

#weigths = Weights(probs)
#items = [1, 2, 3, 4, 5, 6, 7, 8]
#samp = sample(Random.GLOBAL_RNG, items, weights, N_meas)
#println(samp)
"""



N = 3
N_meas = 100
register = [1, 3]

psi = initialise_state(N)
#PauliX!(psi, [1, 4], N)
#toffoli!(psi, [1, 5, 3], N)
hadamard!(psi, [1, 2, 3], N)
#do_measurement(psi, N_meas)
sample_measurement(psi, [1, 2, 3], N_meas)
measure!(psi, register)
sample_measurement(psi, [1, 2, 3], N_meas)

psi2 = initialise_state(N)
hadamard!(psi2, [2], N)
PauliX!(psi2, [3], N)
sample_measurement(psi2, [1, 2, 3], N_meas)

psi3 = initialise_state(N)
#hadamard!(psi2, [2], N)
PauliX!(psi3, [2, 3], N)
sample_measurement(psi3, [1, 2, 3], N_meas)


#pos = [3, 5, 2]
#toffoli!(psi, pos, N)
#println("hi")

#bin_nums = [reverse(digits(i, base=2, pad=4)) for i in 0:15]

#out = filter_binary_numbers(bin_nums, [0, 0], [1, 3])
#println(out)

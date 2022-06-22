
###################################
# test measurement function for MPS
###################################

include("../QSim.jl")
using BenchmarkTools


# set constants
N = 10
N_meas = 50
backend_ED = "ED_Julia"
backend = "MPS_ITensor"
lintop = false
random = true
bond = 12
verbose = false
eps = 0.005

# set up circuit
qc = initialise_qcircuit(N, lintop, backend, 250, "naive", random, 25)
#qc_ED = initialise_qcircuit(N, lintop, backend_ED)

# apply some dummy gates
PauliX!(qc, [1, 3, 5])
#PauliX!(qc_ED, [1, 3, 5, 7, 9])
#hadamard!(qc, [i for i in 1:N])
#hadamard!(qc_ED, [i for i in 1:N])

# draw cicruits
draw(qc)
#draw(qc_ED)

# measure
register = [i for i in 1:N]
#register = [1, 3, 5]
#sample_measurement(qc, register, N_meas)
#sample_measurement(qc_ED, register, N_meas)

# test sampleMPS
#res1 = sampleMPS!(qc.StateVector, register, N_meas)
#println("after measurement: ", qc.StateVector)
#res2 = sampleMPS2!(qc.StateVector, N_meas)
#println("after measurement: ", qc.StateVector)
#res3 = sampleMPS3!(qc.StateVector, N_meas, bond)
#println("after measurement: ", qc.StateVector)
#println(res3)
#res4 = sampleMPS4!(qc.StateVector, N_meas, bond)
#println("after measurement: ", qc.StateVector)
#res5 = sampleMPS3!(qc.StateVector, N_meas, bond)

#println(res1==res2)
#println(res2==res3)

#println(res1)

#res1 = sampleMPS!(qc.StateVector, register, N_meas)
#res2 = sampleMPS2!(qc.StateVector, register, N_meas)
#println(res1==res2)
#println("res1: ", res1)
#println("res2: ", res2)

#timing1 = @btime sampleMPS!(qc.StateVector, register, N_meas)
#timing2 = @btime sampleMPS2!(qc.StateVector, register, N_meas)

# store excution time in variable
#timing1 = @belapsed sampleMPS!(qc.StateVector, register, N_meas)
#timing2 = @belapsed sampleMPS2!(qc.StateVector, register, N_meas)

#timing1 = @benchmark sampleMPS!(qc.StateVector, register, N_meas)
#timing2 = @benchmark sampleMPS2!(qc.StateVector, register, N_meas)

#println(timing1)
#println(timing2)


freq = sample_measurement(qc, register, N_meas, verbose)
projective_measurement!(qc, register, verbose)
freq = sample_measurement(qc, register, N_meas, verbose)
println(qc.StateVector)

#@btime projective_measurement!(qc, register, verbose)
#@btime projective_measurement2!(qc, register, verbose)



#@btime sample_measurement(qc, register, N_meas, eps, verbose)
#@btime sample_measurement2(qc, register, N_meas, eps, verbose)
#@btime sampleMPS!(qc.StateVector, N_meas)
#@btime sampleMPS2!(qc.StateVector, N_meas)
#@btime sampleMPS4!(qc.StateVector, N_meas, bond)

#freq_ED = sample_measurement(qc_ED, register, N_meas)
#projective_measurement!(qc_ED, register)
#freq_ED = sample_measurement(qc_ED, register, N_meas)

#println(qc.ClassicalBits)
#println(qc_ED.ClassicalBits)

"""
N = 3
sites = siteinds("QCircuit", N)
configuration = [0, 0, 1]
basis_state = MPS_computationalBasis(sites, configuration)
for i in 1:N
    println(basis_state[i])
end

basis_proj = projector(basis_state)
for i in 1:N
    println(basis_proj[i])
end
"""

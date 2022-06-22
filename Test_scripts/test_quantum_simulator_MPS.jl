
include("QSim.jl")

backend = "MPS_ITensor"
#backend = "ED_Julia"

N = 20
N_meas = 5
lintop = false
maxdim = 40

qc = initialise_qcircuit(N, lintop, backend, maxdim)

println(qc.StateVector)
#println(measure(qc))

hadamard!(qc, [i for i in 1:qc.NumQubits])
PauliX!(qc, [1, 2, 7])

#println(measure(qc))
#sample_measurement(qc, [i for i in 1:N], N_meas)

#PauliX!(qc, [1, 2, 3, 4])

#println(measure(qc))
#draw_circuit(qc)

#hadamard!(qc, [1, 2, 3, 4])

#println(measure(qc))
#draw_circuit(qc)

#println(maxlinkdim(qc.StateVector))

#cnot2!(qc, [1, 2])
#cnot!(qc, [1, 8])
#cnot!(qc, [4, 3])

#println(measure(qc))
#sample_measurement(qc, [1, 2], N_meas)
draw_circuit(qc)

fullSwap!(qc, [1, 2])

#println(measure(qc))
#sample_measurement(qc, [1, 2], N_meas)
draw_circuit(qc)
println(qc.BondDim)

toffoli!(qc, [1, 7, 3])
toffoli!(qc, [3, 4, 8])
toffoli!(qc, [2, 6, 1])
toffoli!(qc, [4, 7, 2])
toffoli!(qc, [3, 6, 2])
toffoli!(qc, [4, 7, 5])

fullSwap!(qc, [3, 8])
fullSwap!(qc, [2, 4])
fullSwap!(qc, [1, 5])

#for i in 1:qc.NumQubits-1
#    cnot!(qc, [i, i+1])
#end


#println(measure(qc))
#sample_measurement(qc, [1, 2, 3, 4], N_meas)
draw_circuit(qc, false)
draw_circuit(qc, true)

println(qc.BondDim)
println(qc.MaxBondDim)


N = 4
s = siteinds(2,N)
chi = 4
psi = randomMPS(s;linkdims=chi)

#orthogonalize!(qc.StateVector, 1)
#println(maxlinkdim(psi))
#println(ITensors.sample!(psi))

result = zeros(Int, N)
A = psi[1]
println(A)

for j in 1:2

    si = siteind(psi, j)
    d = ITensors.dim(si)
    println(si, d)

    pdisc = 0.0
    r = rand()

    n = 1
    An = ITensor()
    pn = 0.0

    while n <= d
        projn = ITensor(si)
        println(projn)

        projn[si=>n] = 1.0
        An = A*dag(projn)
        println(An)

        n += 1
    end

end






N = 10

# initialise CNOT
cnot = AutoMPO()
cnot += "ProjUp", 1
cnot += "ProjDn", 1, "X", 2

# get sites
sites = siteinds("Qubit",N)
#sites = qc.IndexSet

# turn into MPO
gate = MPO(cnot, sites)

println(gate)


# check single site gates as MPO
xtest = OpSum()
xtest += "X", 5

# get sites
sites = siteinds("Qubit",N)
#sites = qc.IndexSet

# turn into MPO
gateX = MPO(xtest, sites)

println("checking X gate: ", gateX)
println(gateX[5])

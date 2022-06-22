include("../QSim.jl")


N = 16+1 # total number of qubits
maxdim = 20 # maximum allowed bond dimension for MPS
contmethod = "naive" # contraction method for MPS (not so important right now)
random = false # no random intial state for circuit
lintop = false # possibility to use linear qubit topology, also don't touch for now :)
randombond = 10 # bond dimension of random initial MPS, if needed'
dt = 0.1
last_qubit = N

backend_MPS = "MPS_ITensor" # calculate with MPS (any higher number of qubits)


qc_MPS = initialise_qcircuit(N, lintop, backend_MPS, maxdim,
contmethod, random, randombond)

qc_MPS2 = initialise_qcircuit(N, lintop, backend_MPS, maxdim,
contmethod, random, randombond)


for j in 0:3	
	hadamard!(qc_MPS, [4 + 4j, 2 + mod(4+4j, 16)])

	for i in 3 : 5
		cnot!(qc_MPS, [last_qubit, 1 + mod(i + 4j, 16)])
	end


	#some random gate
	N_reg = 1
	customgate_exp1 = initialise_custom_gate(qc_MPS, N_reg, "cust")
	α = exp(1.0im*2π*1/4*dt)
	β = Complex(0.0)
	γ = Complex(0.0)
	δ = exp(1.0im*2π*1/4*dt)
	UGate!(customgate_exp1, [last_qubit], α, β, γ, δ)
	apply_custom_gate(qc_MPS, customgate_exp1, 2)

end

draw(qc_MPS)


for j in 0:3	
	hadamard!(qc_MPS2, [2+4j, 3 + mod(4+4j, 16)])

	for i in 1 : 6
		cnot!(qc_MPS2, [last_qubit, 1 + mod(i + 4j, 16)])
	end

	#some random gate
	N_reg = 1
	customgate_exp1 = initialise_custom_gate(qc_MPS, N_reg, "cust")
	α = exp(1.0im*2π*1/4*dt)
	β = Complex(0.0)
	γ = Complex(0.0)
	δ = exp(1.0im*2π*1/4*dt)
	UGate!(customgate_exp1, [last_qubit], α, β, γ, δ)
	apply_custom_gate(qc_MPS2, customgate_exp1, 2)

end

draw(qc_MPS2)
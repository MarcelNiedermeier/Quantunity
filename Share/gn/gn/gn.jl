include("../QSim.jl")


N = 32+1 # total number of qubits
maxdim = 20 # maximum allowed bond dimension for MPS
contmethod = "naive" # contraction method for MPS (not so important right now)
random = false # no random intial state for circuit
lintop = false # possibility to use linear qubit topology, also don't touch for now :)
randombond = 10 # bond dimension of random initial MPS, if needed'
dt = 0.1 #timestep 
last_qubit = N

#some constants
r = 1
k = 0.5
m = 1
g = 1.01

backend_MPS = "MPS_ITensor" # calculate with MPS (any higher number of qubits)


qc_MPS = initialise_qcircuit(N, lintop, backend_MPS, maxdim,
contmethod, random, randombond)



for t in 1:1  #timesteps
	for j in 0:6 
		
		######1
		hadamard!(qc_MPS, [2+4j, 4 + mod(4+4j, N)])

		for i in 1:7
			cnot!(qc_MPS, [last_qubit, 1 + mod(i + 4j, N)])
		end


		N_reg = 1
		customgate_exp1 = initialise_custom_gate(qc_MPS, N_reg, "cust")
		α = exp(1.0im*1/4*t*dt)
		β = Complex(0.0)
		γ = Complex(0.0)
		δ = exp(1.0im*1/4*t*dt)
		UGate!(customgate_exp1, [last_qubit], α, β, γ, δ)
		apply_custom_gate(qc_MPS, customgate_exp1, last_qubit)


		for i in reverse(1:7)
			cnot!(qc_MPS, [last_qubit, i+4j])
		end
		hadamard!(qc_MPS, [2+4j, 8+4j])



		#########2
		hadamard!(qc_MPS, [1+4j, 3 + mod(4+4j, N)])

		for i in 0:6
			cnot!(qc_MPS, [last_qubit, 1 + mod(i + 4j, N)])
		end

		apply_custom_gate(qc_MPS, customgate_exp1, last_qubit)

		for i in reverse(0:6)
			cnot!(qc_MPS, [last_qubit, 1 + mod(i + 4j, N)])
		end
		hadamard!(qc_MPS, [1+4j, 3 + mod(4+4j, N)])



		#########3

		N_reg = 1
		customgateR = initialise_custom_gate(qc_MPS, N_reg, "cust")
		α = Complex(1/sqrt(2))
		β = -im/sqrt(2)
		γ = im/sqrt(2)
		δ = Complex(-1/sqrt(2))
		UGate!(customgateR, [last_qubit], α, β, γ, δ)
		apply_custom_gate(qc_MPS, customgateR, 2+4j)
		apply_custom_gate(qc_MPS, customgateR, 4 + mod(4+4j, N))


		for i in 1:7
			cnot!(qc_MPS, [last_qubit, 1 + mod(i + 4j, N)])
		end

		apply_custom_gate(qc_MPS, customgate_exp1, last_qubit)


		for i in reverse(1:7)
			cnot!(qc_MPS, [last_qubit, 1 + mod(i + 4j, N)])
		end

		apply_custom_gate(qc_MPS, customgateR, 2+4j)
		apply_custom_gate(qc_MPS, customgateR, 4 + mod(4+4j, N))

		

		#########4
		apply_custom_gate(qc_MPS, customgateR, 1+4j)
		apply_custom_gate(qc_MPS, customgateR, 3 + mod(4+4j, N))

		for i in 0:6
			cnot!(qc_MPS, [last_qubit, 1 + mod(i + 4j, N)])
		end

		apply_custom_gate(qc_MPS, customgate_exp1, last_qubit)

		for i in reverse(0:6)
			cnot!(qc_MPS, [last_qubit, 1 + mod(i + 4j, N)])
		end
		apply_custom_gate(qc_MPS, customgateR, 1+4j)
		apply_custom_gate(qc_MPS, customgateR, 3 + mod(4+4j, N))


		
		#######5
		hadamard!(qc_MPS, [4+4j, 2 + mod(4+4j, N)])

		for i in 3:5
			cnot!(qc_MPS, [last_qubit, 1 + mod(i + 4j, N)])
		end

		N_reg = 1
		customgate_exp2 = initialise_custom_gate(qc_MPS, N_reg, "cust")
		α = exp(-1.0im*r/4*t*dt)
		β = Complex(0.0)
		γ = Complex(0.0)
		δ = exp(-1.0im*r/4*t*dt)
		UGate!(customgate_exp2, [last_qubit], α, β, γ, δ)
		apply_custom_gate(qc_MPS, customgate_exp2, last_qubit)

		for i in reverse(3:5)
			cnot!(qc_MPS, [last_qubit, 1 + mod(i + 4j, N)])
		end
		
		hadamard!(qc_MPS, [4+4j, 2 + mod(4+4j, N)])


		######6
		hadamard!(qc_MPS, [3+4j, 1 + mod(4+4j, N)])

		for i in 2:4
			cnot!(qc_MPS, [last_qubit, 1 + mod(i + 4j, N)])
		end


		apply_custom_gate(qc_MPS, customgate_exp2, last_qubit)


		for i in reverse(2:4)
			cnot!(qc_MPS, [last_qubit, 1 + mod(i + 4j, N)])
		end

		hadamard!(qc_MPS, [3+4j, 1 + mod(4+4j, N)])


		########7 
		apply_custom_gate(qc_MPS, customgateR, 4+4j)
		apply_custom_gate(qc_MPS, customgateR, 2 + mod(4+4j, N))

		for i in 3:5
			cnot!(qc_MPS, [last_qubit, 1 + mod(i + 4j, N)])
		end

		apply_custom_gate(qc_MPS, customgate_exp2, last_qubit)


		for i in reverse(3:5)
			cnot!(qc_MPS, [last_qubit, 1 + mod(i + 4j, N)])
		end

		apply_custom_gate(qc_MPS, customgateR, 4+4j)
		apply_custom_gate(qc_MPS, customgateR, 2 + mod(4+4j, N))


		#######8

		apply_custom_gate(qc_MPS, customgateR, 3+4j)
		apply_custom_gate(qc_MPS, customgateR, 1 + mod(4+4j, N))

		for i in 2:4
			cnot!(qc_MPS, [last_qubit, 1 + mod(i + 4j, N)])
		end

		apply_custom_gate(qc_MPS, customgate_exp2, last_qubit)


		for i in reverse(2:4)
			cnot!(qc_MPS, [last_qubit, 1 + mod(i + 4j, N)])
		end

		apply_custom_gate(qc_MPS, customgateR, 3+4j)
		apply_custom_gate(qc_MPS, customgateR, 1 + mod(4+4j, N))


		######9
		hadamard!(qc_MPS, [4 + 4j, 4 + mod(4+4j, N)])

		for i in 3:7
			cnot!(qc_MPS, [last_qubit, 1 + mod(i + 4j, N)])
		end


		N_reg = 1
		customgate_exp3 = initialise_custom_gate(qc_MPS, N_reg, "cust")
		α = -exp(1.0im*r/4*t*dt)
		β = Complex(0.0)
		γ = Complex(0.0)
		δ = -exp(1.0im*r/4*t*dt)
		UGate!(customgate_exp3, [last_qubit], α, β, γ, δ)
		apply_custom_gate(qc_MPS, customgate_exp3, last_qubit)


		for i in reverse(3:7)
			cnot!(qc_MPS, [last_qubit, 1 + mod(i + 4j, N)])
		end

		hadamard!(qc_MPS, [4 + 4j, 4 + mod(4+4j, N)])


		######10
		hadamard!(qc_MPS, [3 + 4j, 3 + mod(4+4j, N)])

		for i in 2:6
			cnot!(qc_MPS, [last_qubit, 1 + mod(i + 4j, N)])
		end


		apply_custom_gate(qc_MPS, customgate_exp3, last_qubit)

		for i in reverse(2:6)
			cnot!(qc_MPS, [last_qubit, 1 + mod(i + 4j, N)])
		end

		hadamard!(qc_MPS, [3 + 4j, 3 + mod(4+4j, N)])



		########11
		apply_custom_gate(qc_MPS, customgateR, 4+4j)
		apply_custom_gate(qc_MPS, customgateR, 4 + mod(4+4j, N))

		for i in 3:7
			cnot!(qc_MPS, [last_qubit, 1 + mod(i + 4j, N)])
		end

		apply_custom_gate(qc_MPS, customgate_exp3, last_qubit)


		for i in reverse(3:7)
			cnot!(qc_MPS, [last_qubit, 1 + mod(i + 4j, N)])
		end
		apply_custom_gate(qc_MPS, customgateR, 4+4j)
		apply_custom_gate(qc_MPS, customgateR, 4 + mod(4+4j, N))



		########12
		apply_custom_gate(qc_MPS, customgateR, 3+4j)
		apply_custom_gate(qc_MPS, customgateR, 3 + mod(4+4j, N))

		for i in 2:6
			cnot!(qc_MPS, [last_qubit, 1 + mod(i + 4j, N)])
		end

		apply_custom_gate(qc_MPS, customgate_exp3, last_qubit)


		for i in reverse(2:6)
			cnot!(qc_MPS, [last_qubit, 1 + mod(i + 4j, N)])
		end
		apply_custom_gate(qc_MPS, customgateR, 3+4j)
		apply_custom_gate(qc_MPS, customgateR, 3 + mod(4+4j, N))




		#######13
		hadamard!(qc_MPS, [2+4j, 2 + mod(4+4j, N)])

		for i in 1:5
			cnot!(qc_MPS, [last_qubit, 1 + mod(i + 4j, N)])
		end


		N_reg = 1
		customgate_exp4 = initialise_custom_gate(qc_MPS, N_reg, "cust")
		α = exp(1.0im*r/4*t*dt)
		β = Complex(0.0)
		γ = Complex(0.0)
		δ = exp(1.0im*r/4*t*dt)
		UGate!(customgate_exp4, [last_qubit], α, β, γ, δ)
		apply_custom_gate(qc_MPS, customgate_exp4, last_qubit)

		for i in reverse(1:5)
			cnot!(qc_MPS, [last_qubit, 1 + mod(i + 4j, N)])
		end

		hadamard!(qc_MPS, [2+4j, 2 + mod(4+4j, N)])



		#######14
		hadamard!(qc_MPS, [1+4j, 1 + mod(4+4j, N)])

		for i in 0:4
			cnot!(qc_MPS, [last_qubit, 1 + mod(i + 4j, N)])
		end


		apply_custom_gate(qc_MPS, customgate_exp4, last_qubit)

		for i in reverse(0:4)
			cnot!(qc_MPS, [last_qubit, 1 + mod(i + 4j, N)])
		end
		hadamard!(qc_MPS, [1+4j, 1 + mod(4+4j, N)])



		########15
		apply_custom_gate(qc_MPS, customgateR, 2+4j)
		apply_custom_gate(qc_MPS, customgateR, 2 + mod(4+4j, N))

		for i in 1:5
			cnot!(qc_MPS, [last_qubit, 1 + mod(i + 4j, N)])
		end

		apply_custom_gate(qc_MPS, customgate_exp4, last_qubit)


		for i in reverse(1:5)
			cnot!(qc_MPS, [last_qubit, 1 + mod(i + 4j, N)])
		end

		apply_custom_gate(qc_MPS, customgateR, 2+4j)
		apply_custom_gate(qc_MPS, customgateR, 2 + mod(4+4j, N))



		########16
		apply_custom_gate(qc_MPS, customgateR, 1+4j)
		apply_custom_gate(qc_MPS, customgateR, 1 + mod(4+4j, N))

		for i in 0:4
			cnot!(qc_MPS, [last_qubit, 1 + mod(i + 4j, N)])
		end

		apply_custom_gate(qc_MPS, customgate_exp4, last_qubit)


		for i in reverse(0:4)
			cnot!(qc_MPS, [last_qubit, 1 + mod(i + 4j, N)])
		end

		apply_custom_gate(qc_MPS, customgateR, 1+4j)
		apply_custom_gate(qc_MPS, customgateR, 1 + mod(4+4j, N))


		#######17

		cnot!(qc_MPS, [last_qubit, 8+4j])
		cnot!(qc_MPS, [last_qubit, 7+4j])

		N_reg = 1
		customgate_exp5 = initialise_custom_gate(qc_MPS, N_reg, "cust")
		α = exp(1.0im*k*g^2/4*t*dt)
		β = Complex(0.0)
		γ = Complex(0.0)
		δ = exp(1.0im*k*g^2/4*t*dt)
		UGate!(customgate_exp5, [last_qubit], α, β, γ, δ)
		apply_custom_gate(qc_MPS, customgate_exp5, last_qubit)

		cnot!(qc_MPS, [last_qubit, 7+4j])
		cnot!(qc_MPS, [last_qubit, 8+4j])


		cnot!(qc_MPS, [last_qubit, 8+4j])
		cnot!(qc_MPS, [last_qubit, 6+4j])

		N_reg = 1
		customgate_exp6 = initialise_custom_gate(qc_MPS, N_reg, "cust")
		α = exp(-1.0im*k*g^2/4*t*dt)
		β = Complex(0.0)
		γ = Complex(0.0)
		δ = exp(-1.0im*k*g^2/4*t*dt)
		UGate!(customgate_exp6, [last_qubit], α, β, γ, δ)
		apply_custom_gate(qc_MPS, customgate_exp6, last_qubit)

		cnot!(qc_MPS, [last_qubit, 6+4j])
		cnot!(qc_MPS, [last_qubit, 8+4j])


		cnot!(qc_MPS, [last_qubit, 8+4j])
		cnot!(qc_MPS, [last_qubit, 5+4j])
		apply_custom_gate(qc_MPS, customgate_exp6, last_qubit)
		cnot!(qc_MPS, [last_qubit, 5+4j])
		cnot!(qc_MPS, [last_qubit, 8+4j])


		cnot!(qc_MPS, [last_qubit, 7+4j])
		cnot!(qc_MPS, [last_qubit, 6+4j])
		apply_custom_gate(qc_MPS, customgate_exp6, last_qubit)
		cnot!(qc_MPS, [last_qubit, 6+4j])
		cnot!(qc_MPS, [last_qubit, 7+4j])


		cnot!(qc_MPS, [last_qubit, 7+4j])
		cnot!(qc_MPS, [last_qubit, 5+4j])
		apply_custom_gate(qc_MPS, customgate_exp6, last_qubit)
		cnot!(qc_MPS, [last_qubit, 5+4j])
		cnot!(qc_MPS, [last_qubit, 7+4j])


		cnot!(qc_MPS, [last_qubit, 6+4j])
		cnot!(qc_MPS, [last_qubit, 5+4j])
		apply_custom_gate(qc_MPS, customgate_exp5, last_qubit)
		cnot!(qc_MPS, [last_qubit, 5+4j])
		cnot!(qc_MPS, [last_qubit, 6+4j])

		cnot!(qc_MPS, [last_qubit, 8+4j])
		N_reg = 1
		customgate_exp7 = initialise_custom_gate(qc_MPS, N_reg, "cust")
		α = exp(-1.0im*(m*k+r)/2*t*dt)
		β = Complex(0.0)
		γ = Complex(0.0)
		δ = exp(-1.0im*(m*k+r)/2*t*dt)
		UGate!(customgate_exp7, [last_qubit], α, β, γ, δ)
		apply_custom_gate(qc_MPS, customgate_exp7, last_qubit)
		cnot!(qc_MPS, [last_qubit, 8+4j])

		cnot!(qc_MPS, [last_qubit, 7+4j])
		apply_custom_gate(qc_MPS, customgate_exp7, last_qubit)
		cnot!(qc_MPS, [last_qubit, 7+4j])

		cnot!(qc_MPS, [last_qubit, 6+4j])
		customgate_exp8 = initialise_custom_gate(qc_MPS, N_reg, "cust")
		α = exp(1.0im*(m*k+r)/2*t*dt)
		β = Complex(0.0)
		γ = Complex(0.0)
		δ = exp(1.0im*(m*k+r)/2*t*dt)
		UGate!(customgate_exp8, [last_qubit], α, β, γ, δ)
		apply_custom_gate(qc_MPS, customgate_exp8, last_qubit)
		cnot!(qc_MPS, [last_qubit, 6+4j])

		cnot!(qc_MPS, [last_qubit, 5+4j])
		apply_custom_gate(qc_MPS, customgate_exp8, last_qubit)
		cnot!(qc_MPS, [last_qubit, 5+4j])
	end
end
draw(qc_MPS)




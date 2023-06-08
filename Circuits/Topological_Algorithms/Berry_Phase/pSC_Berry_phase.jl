
########################
## Berry Phase p-wave SC
########################

include("../../QuantumSimulator_DEV/QSim.jl")
include("topological_quantities.jl")


# parameters
Nk = 1000 # 1000
Nk_ex = 40
dt = 0.02
N_meas = 10000
N_qubits = 4
n_prec = 3
maxdim = 10


# parameter space for p-wave SC
mus = [-3, -1, 1, 3]
#mus = [1]
kys = collect(range(-π, π, length=20))
Berry_phases_pSC = zeros(length(mus), length(kys))
Berry_phases_pSC_ex = zeros(length(mus), length(kys))
Berry_phases_pSC_QPE = zeros(length(mus), length(kys))

# loop through parameter space and calculate Berry phase
for i in 1:length(mus)
   for j in 1:length(kys)
      println("doing μ = $(mus[i]), ky = $(kys[j])")
      params = [kys[j], 1, 1, mus[i]]
      Berry_phases_pSC_ex[i, j] = Berry_phase(H_p_wave, Nk_ex, params)
      Berry_phases_pSC[i, j] = measure_Berry_phase_1D(H_p_wave, params, N_meas, Nk, dt)

      #phases, probs_max = Berry_phase_QPE(H_p_wave, params, Nk, dt, N_qubits, n_prec, N_meas, backend="MPS_ITensor", maxbond=maxdim)
      #println("phases, probs_max: ", phases, probs_max)

      # map to right range
      #if phases[1] < π
      #   Berry_phases_pSC_QPE[i, j] = phases[1]
      #else
      #   Berry_phases_pSC_QPE[i, j] = phases[1] - 2π
      #end

      #println("Berry phase from QPE: ", Berry_phases_pSC_QPE[i, j])

   end
end


# save data
datafile[1, 1] = Nk
datafile[1, 2] = Nk_ex
datafile[1, 3] = dt
datafile[1, 4] = N_meas
datafile[1, 5] = N_qubits
datafile[1, 6] = n_prec
datafile[1, 7] = maxdim
datafile[2:2+length(mus)-1, :] = Berry_phases_pSC[:, :]
datafile[2+length(mus):2+2*length(mus)-1, :] = Berry_phases_pSC_QPE[:, :]
datafile[2+2*length(mus):2+3*length(mus)-1, :] = Berry_phases_pSC_ex[:, :]

open("Data/data_pSC_Berry_Nk_$(Nk)_dt_$(dt)_N_qubits_$(N_qubits)_maxdim_$(maxdim).csv", "w") do io
    writedlm(io, datafile, ", ")
end



p1 = plot!(kys, Berry_phases_pSC_ex[1, :], label="μ = $(mus[1]), exact", legend=true)
plot!(kys, Berry_phases_pSC[1, :], label="μ = $(mus[1]), H-test", marker=:star, seriestype=:scatter)
#plot!(kys, Berry_phases_pSC_QPE[1, :], label="μ = $(mus[1])", marker=:star, seriestype=:scatter)
for i in 2:length(mus)
   plot!(kys, Berry_phases_pSC[i, :], label="μ = $(mus[i]), H-test", marker=:star, seriestype=:scatter)
   #plot!(kys, Berry_phases_pSC_QPE[i, :], label="μ = $(mus[i])", marker=:star, seriestype=:scatter)
   plot!(kys, Berry_phases_pSC_ex[i, :], label="μ = $(mus[i]), exact")
end
xlabel!("ky")
ylabel!("θ")
savefig(p1, "Plots/berry_pSC_Berry_Nk_$(Nk)_dt_$(dt)_N_qubits_$(N_qubits)_maxdim_$(maxdim).png")


########################
## Berry Phase QWZ model
########################

"""
Script to calculate the Berry phase of the QWZ model using a Quantum
Phase Estimation algorithm and a Hadamard Test algorithm. The algorithm
is looped over the parameters [u, k_y] of the computation (with an adiabatic
loop executed over k_x). A backend of matrix product states is used.
"""

using ArgParse
include("../../QuantumSimulator_DEV/QSim.jl")


# parameters
#Nk = 2000 # 1000
#Nk_ex = 20
#dt = 0.02
#N_meas = 20000
#N_qubits = 7
#n_prec = 5
#maxdim = 20


""" Function to parse arguments from command line. """
function parse_commandline()
   s = ArgParseSettings()
   @add_arg_table s begin
      "--N_qubits"
         help = "number of qubits"
         default = 10
      "--n_prec"
         help = "number of qubits measured in QPE"
         default = 7
      "--Nk"
         help = "number of steps in adiabatic loop"
         default = 1000
      "--dt"
         help = "increment in time evolution"
         default = 0.02
      "--N_meas"
         help = "number of measurements taken"
         default = 20000
      "--maxdim"
         help = "maximum bond dimension"
         default = 20
      "--u"
         help = "parameter u of QWZ model"
         default = -1
      "--ky"
         help = "parameter ky of QWZ model"
         default = 0
   end
   return parse_args(s)
end

# get arguments from command line input
@show parsed_args = parse_commandline()

# set constants and input arguments
Nk = parsed_args["Nk"] # 1000
dt = parsed_args["dt"]
N_meas = parsed_args["N_meas"]
N_qubits = parsed_args["N_qubits"]
n_prec = parsed_args["n_prec"]
maxdim = parsed_args["maxdim"]
u = parsed_args["u"]
ky = parsed_args["ky"]


# parameter space for QWZ model
#us = [-3, -1, 1, 3]
#us = [-1]
#kys = collect(range(-π, π, length=Nk_ex))
#kys = [-2]
#Berry_phases_QWZ = zeros(length(us), length(kys))
#Berry_phases_QWZ_ex = zeros(length(us), length(kys))
#Berry_phases_QWZ_QPE = zeros(length(us), length(kys))
#Berry_phases_QWZ_QPE_av = zeros(length(us), length(kys))

# save results
Berry_phase_QWZ = 0.
Berry_phase_QWZ_ex = 0.
Berry_phase_QWZ_QPE = 0.
Berry_phases_QWZ_QPE = []
Berry_probs_QWZ_QPE = []


###############
function main()
###############

   # set parameters of Hamiltonian
   params = [ky, u]

   # exact calculation, Hadamard test, QPE
   Berry_phase_QWZ_ex = Berry_phase(QWZ_model, Nk_ex, params)
   Berry_phase_QWZ = measure_Berry_phase_1D(QWZ_model, params, N_meas, Nk, dt)
   phases, probs_max = Berry_phase_QPE(QWZ_model, params, Nk, dt, N_qubits, n_prec, N_meas, backend="MPS_ITensor", maxbond=maxdim)

   # recover phase estimate from QPE
   phase_est = phases[1]
   #phase_est_av = phases[1]*probs_max[1] + phases[2]*probs_max[2]
   #Berry_phase_QWZ_QPE[i, j] = phase_est-π
   Berry_phase_QWZ_QPE = π-phase_est # where does this come from??
   #Berry_phase_QWZ_QPE_av[i, j] = π-phase_est_av # where does this come from??

   # save also remaining data from QPE
   for i in 1:length(phases)
      push!(Berry_phases_QWZ_QPE, π-phases[i])
      push!(Berry_probs_QWZ_QPE, probs_max[i])
   end
end


# do main calculation
main()


# save output in datafile
path = "/scratch/work/niederm1/calculations/QuantumSimulations/Data/Berry_phase/QWZ_Berry_phase_Nqubits_$(N_qubits)_Nk_$(Nk)_maxdim_$(maxdim)_u_$(u)_ky_$(ky).jld2"
other_data = [:N_qubits, :n_prec, :Nk, :dt, :N_meas, :maxdim, :u, :ky, :Berry_phase_QWZ_ex, :Berry_phase_QWZ, :Berry_phase_QWZ_QPE, :Berry_phases_QWZ_QPE, :Berry_probs_QWZ_QPE]
save_data(path, other_data)


#println("Berry phase exact: ", Berry_phases_QWZ_ex)
#println("Berry phase from QPE: ", Berry_phases_QWZ_QPE)
#println("Berry phase from QPE av: ", Berry_phases_QWZ_QPE_av)

# datafile for output
#datafile = zeros(3*length(us)+1, length(kys))

# save data
#datafile[1, 1] = Nk
#datafile[1, 2] = Nk_ex
#datafile[1, 3] = dt
#datafile[1, 4] = N_meas
#datafile[1, 5] = N_qubits
#datafile[1, 6] = n_prec
#datafile[1, 7] = maxdim
#datafile[2:2+length(us)-1, :] = Berry_phases_QWZ[:, :]
#datafile[2+length(us):2+2*length(us)-1, :] = Berry_phases_QWZ_QPE[:, :]
#datafile[2+2*length(us):2+3*length(us)-1, :] = Berry_phases_QWZ_ex[:, :]

#open("Data/data_QWZ_Berry_Nk_$(Nk)_dt_$(dt)_N_qubits_$(N_qubits)_maxdim_$(maxdim).csv", "w") do io
#    writedlm(io, datafile, ", ")
#end

# plotting
#p1 = plot!(kys, Berry_phases_QWZ_ex[1, :], label="u = $(us[1]), exact", legend=:bottomleft)
#plot!(kys, Berry_phases_QWZ[1, :], label="u = $(us[1]), H-test", marker=:circle, seriestype=:scatter)
#plot!(kys, Berry_phases_QWZ_QPE[1, :], label="u = $(us[1]), QPE", marker=:star, seriestype=:scatter)
##plot!(kys, Berry_phases_QWZ_QPE_av[1, :], label="u = $(us[1]), QPE av", marker=:circle, seriestype=:scatter)
#for i in 2:length(us)
#   plot!(kys, Berry_phases_QWZ_ex[i, :], label="u = $(us[i]), exact")
#   plot!(kys, Berry_phases_QWZ[i, :], label="u = $(us[i]), H-test", marker=:circle, seriestype=:scatter)
#   plot!(kys, Berry_phases_QWZ_QPE[i, :], label="u = $(us[i]), QPE", marker=:star, seriestype=:scatter)
#   #plot!(kys, Berry_phases_QWZ_QPE_av[i, :], label="u = $(us[i]), QPE av", marker=:circle, seriestype=:scatter)
#end
#xlabel!("ky")
#ylabel!("θ")
#savefig(p1, "Plots/berry_QWZ_Berry_Nk_$(Nk)_dt_$(dt)_N_qubits_$(N_qubits)_maxdim_$(maxdim).png")

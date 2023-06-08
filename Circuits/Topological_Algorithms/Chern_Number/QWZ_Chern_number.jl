
#########################
## Chern Number QWZ model
#########################

include("../../QuantumSimulator_DEV/QSim.jl")
include("topological_quantities.jl")



#######################################
# Quantum Circuit adiabatic double loop
#######################################


# Brillouin zone
Nk = 20
δk = 2π/Nk
N_link = 200 # steps for adiabatic loop per link, check for different numbers!
dt = 0.02 # 0.1

# quantum circuit
N_meas = 100000 # very small Berry phases, need lot of measurement precision

# parameter space
Nmu = 4 # 1, 20
mu_min = -3
mu_max = 3
#mus = LinRange(mu_min, mu_max, Nmu)
mus = [-3, -1, 1, 3]
#mus = [-1, 1]
Nmu_ex = 100
mus_ex = LinRange(mu_min, mu_max, Nmu_ex)

# save Chern numbers
Clist_ex = zeros(Nmu_ex)
Clist_adiabatic_loop = zeros(Nmu)

# exact
for a in 1:Nmu_ex
    params = [mus_ex[a]]
    _, Clist_ex[a] = Chern_number(QWZ_model, params, Nk)
end

# quantum
for a in 1:Nmu
    println("Doing μ = $(mus[a])")
    params = [mus[a]]
    #Berry_flux_ex, Clist_ex[a] = Chern_number(QWZ_model, params, Nk)
    _, Clist_adiabatic_loop[a] = measure_Chern_number_adiabatic_loop(QWZ_model, params, N_link, Nk, dt, N_meas)
    #Berry_flux, Berry_flux_ex, Clist_adiabatic_loop[a], Clist_ex[a] = measure_Chern_number_plaquette_comparison(QWZ_model, params, N_link, Nk, dt, N_meas)
end

println("Chern number from adiabatic loop: ", Clist_adiabatic_loop)
println("Chern number exact: ", Clist_ex)


#################
# Data processing
#################

# plotting
p1 = plot!(mus_ex, Clist_ex, label="exact", legend=:bottomright)
#plot!(mus, Clist_link, label="quantum link", markershape=:circle)
plot!(mus, Clist_adiabatic_loop, label="quantum adiabatic", markershape=:cross)
xlabel!("μ")
ylabel!("C")
savefig(p1, "Plots/Chern_number_QWZ_Nk_$(Nk)_Nlink_$(N_link)_dt_$(dt)_N_meas_$(N_meas).png")


# save data
datafile = zeros(5, Nmu_ex)
datafile[1, 1] = Nk
datafile[1, 1] = δk
datafile[1, 1] = N_link
datafile[1, 1] = dt
datafile[1, 1] = N_meas

datafile[2, 1:Nmu] = mus[:]
datafile[3, :] = mus_ex[:]
datafile[4, 1:Nmu] = Clist_adiabatic_loop[:]
datafile[5, :] = Clist_ex[:]

#open("Chern_number_p_wave.csv", "w") do io
open("Data/Chern_number_QWZ_Nk_$(Nk)_Nlink_$(N_link)_dt_$(dt)_N_meas_$(N_meas).csv", "w") do io
    writedlm(io, datafile, ", ")
end



#################################
# Quantum Circuit link estimation
#################################

"""
# Brillouin zone
Nk = 16

# quantum circuit
N_meas = 5000

# parameter space
Nmu = 4 # 1, 20
mu_min = -3
mu_max = 3
#mus = LinRange(mu_min, mu_max, Nmu)
mus = [-3, -1, 1, 3]


# save Chern numbers
Clist_link = zeros(Nmu)


for a in 1:Nmu
    params = [mus[a]]
    Clist_link[a] = measure_Chern_number(QWZ_model, params, Nk, N_meas)
end

println("Chern numbers for QWZ model (quantum link estimation): ", Clist_link)
"""

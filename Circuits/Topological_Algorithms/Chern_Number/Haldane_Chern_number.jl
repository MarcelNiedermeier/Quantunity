
#############################
## Chern Number Haldane model
#############################

include("../../QuantumSimulator_DEV/QSim.jl")
include("topological_quantities.jl")


#######################################
# Quantum Circuit adiabatic double loop
#######################################


# Brillouin zone
Nk = 20
δk = 2π/Nk
N_link = 100 # steps for adiabatic loop per link
dt = 0.02 # 0.1

# quantum circuit
N_meas = 100000 # very small Berry phases, need lot of measurement precision

# parameter space
Nm = 10 # 1, 20
Nϕ = 10
Nm_ex = 200
Nϕ_ex = 200
m_min = -3*sqrt(3)
m_max = 3*sqrt(3)
ϕ_min = -π
ϕ_max = π
ms = LinRange(m_min, m_max, Nm)
ϕs = LinRange(ϕ_min, ϕ_max, Nϕ)
ms_ex = LinRange(m_min, m_max, Nm_ex)
ϕs_ex = LinRange(ϕ_min, ϕ_max, Nϕ_ex)

# Haldane model parameters
t1 = 1
t2 = 1
#m = params[3]
#ϕ = params[4]
a = 1

# save Chern numbers
C_ex = zeros(Nm_ex, Nϕ_ex)
C_adiabatic_loop = zeros(Nm, Nϕ)

# exact
#for i in 1:Nm_ex
#    for j in 1:Nϕ_ex
#        params = [t1, t2, ms_ex[i], ϕs_ex[j], a]
#        C_ex[i, j] = Chern_number(Haldane_model, params, Nk)
#    end
#end

# quantum
for i in 1:Nm
    for j in 1:Nϕ
        println("Doing m = $(ms[i]) and ϕ = $(ϕs[j])")
        params = [t1, t2, ms[i], ϕs[j], a]
        _, C_adiabatic_loop[i, j] = measure_Chern_number_adiabatic_loop(Haldane_model, params, N_link, Nk, dt, N_meas)
    end
end



#################
# Data processing
#################

# plotting
#p1 = heatmap!(ϕs_ex, ms_ex, C_ex')
p1 = heatmap!(ϕs, ms, C_adiabatic_loop')
#plot!(mus, Clist_link, label="quantum link", markershape=:circle)
#plot!(mus, Clist_adiabatic_loop, label="quantum adiabatic", markershape=:cross)
xlabel!("ϕ")
ylabel!("m")
savefig(p1, "Plots/Chern_number_Haldane_ad_loop_Nk_$(Nk)_Nlink_$(N_link)_dt_$(dt)_N_meas_$(N_meas).png")


# save data
datafile = zeros(Nm_ex+3, Nϕ_ex)
datafile[1, 1] = Nk
datafile[1, 1] = δk
datafile[1, 1] = N_link
datafile[1, 1] = dt
datafile[1, 1] = N_meas
datafile[2, :] = ms_ex[:]
datafile[3, :] = ϕs_ex[:]
datafile[4:, :] = C_ex[:]

open("Data/Chern_number_QWZ_Nk_$(Nk)_Nlink_$(N_link)_dt_$(dt)_N_meas_$(N_meas).csv", "w") do io
    writedlm(io, datafile, ", ")
end



# quantum
#for a in 1:Nmu
#    println("Doing μ = $(mus[a])")
#    params = [mus[a]]
#    Clist_adiabatic_loop[a] = measure_Chern_number_adiabatic_loop(QWZ_model, params, N_link, Nk, dt, N_meas)
#end

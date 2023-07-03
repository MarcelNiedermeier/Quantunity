
#########################
## Chern Number QWZ model
#########################

include("/scratch/work/niederm1/apps/QuantumSimulator/QSim.jl")


#######################################
# Quantum Circuit adiabatic double loop
#######################################


function main_Chern_QWZ(ind)

    # hard parameters
    Nk = 20
    δk = 2π/Nk
    #dt = 0.02 # 0.1
    Nu = 40 # 1, 20
    u_min = -3
    u_max = 3
    #N_meas = 100000 # very small Berry phases, need lot of measurement precision

    # get parameter space
    us = LinRange(u_min, u_max, Nu)
    N_measurements = [10000, 20000, 50000, 100000, 200000, 500000]
    N_links = [100, 200, 500, 1000, 2000, 5000]
    dts = [0.001, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1]
    parameter_space = collect(Iterators.product(us, N_measurements, N_links, dts))

    ps = parameter_space[ind]
    u = ps[1]
    N_meas = ps[2]
    N_link = ps[3]
    dt = ps[4]

    # parameters of Hamiltonian
    params = [u]

    # exact calculation, HT
    _, C_ad_loop_ex = Chern_number(QWZ_model, params, Nk)
    _, C_ad_loop_HT = measure_Chern_number_adiabatic_loop_HT(QWZ_model, params, N_link, Nk, dt, N_meas)

    return (Nk, Nu, u, N_meas, N_link, dt, C_ad_loop_ex, C_ad_loop_HT)
end
